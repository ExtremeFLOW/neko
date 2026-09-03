! Copyright (c) 2025-2026 The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> ALE Manager: Handles Mesh Motion
module ale_manager
  use num_types, only : rp, dp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use field, only : field_t
  use coefs, only : coef_t
  use space, only : space_t
  use ax_product, only : ax_t, ax_helm_allocator
  use krylov, only : ksp_t, ksp_monitor_t, krylov_solver_factory
  use precon, only : pc_t, precon_allocator, precon_destroy
  use bc_list, only : bc_list_t
  use checkpoint, only : chkp_t
  use zero_dirichlet, only : zero_dirichlet_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use dofmap, only : dofmap_t
  use jacobi, only : jacobi_t
  use hsmg, only : hsmg_t
  use phmg, only : phmg_t
  use device_jacobi, only : device_jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use profiler, only : profiler_start_region, profiler_end_region
  use file, only : file_t
  use logger, only : neko_log, LOG_SIZE
  use advection, only : advection_t
  use ale_rigid_kinematics, only : ale_config_t, pivot_state_t, &
       point_tracker_t, body_kinematics_t, &
       init_pivot_state, update_pivot_location, &
       compute_body_kinematics_built_in, ab_integrate_point_pos
  use ale_routines_cpu, only : add_kinematics_to_mesh_velocity_cpu, &
       update_ale_mesh_cpu, compute_cheap_dist_v2_cpu
  use ale_routines_device, only : add_kinematics_to_mesh_velocity_device, &
       update_ale_mesh_device, compute_cheap_dist_device
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_HIP, NEKO_BCKND_CUDA
  use mpi_f08, only : MPI_WTIME, MPI_Barrier
  use comm, only : NEKO_COMM
  use registry, only : neko_registry
  use field_series, only : field_series_t
  use import_field_utils, only : import_fields
  use time_state, only : time_state_t
  use fld_file, only : fld_file_t
  use user_intf, only : user_t, user_ale_mesh_velocity_intf, &
       user_ale_base_shapes_intf, user_ale_rigid_kinematics_intf, &
       dummy_user_ale_mesh_velocity, dummy_user_ale_rigid_kinematics, &
       dummy_user_ale_base_shapes
  use math, only : glmin, pi, copy
  use device_math, only : device_glmin, device_copy
  use field_math, only : field_rzero, field_add2, &
       field_cmult
  use scalar_bc_projector, only : scalar_bc_projector_t

  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST, device_sync
  use operators, only : rotate_cyc
  use fld_file_output, only : fld_file_output_t
  use, intrinsic :: iso_c_binding, only : c_associated
  implicit none
  private

  public :: compute_stiffness_ale
  public :: add_kinematics_to_mesh_velocity
  public :: update_ale_mesh
  public :: log_rot_angles
  public :: log_pivot

  type, public :: ale_manager_t
     ! Default
     logical :: active = .false.
     logical :: has_moving_boundary = .false.

     type(bc_list_t) :: bc_list
     type(zero_dirichlet_t) :: bc_moving
     type(zero_dirichlet_t) :: bc_fixed

     type(ale_config_t) :: config

     !> Mesh velocity fields (Registered in neko_registry)
     type(field_t), pointer :: wm_x => null()
     type(field_t), pointer :: wm_y => null()
     type(field_t), pointer :: wm_z => null()

     !> History buffers for mesh velocity components
     type(field_series_t) :: wm_x_lag
     type(field_series_t) :: wm_y_lag
     type(field_series_t) :: wm_z_lag

     !> Reference initial grid to calculate the rotation accurately
     type(field_t) :: x_ref, y_ref, z_ref

     !> Pivot states for each body
     type(pivot_state_t), allocatable :: ale_pivot(:)
     type(body_kinematics_t), allocatable :: body_kin(:)

     !> Base shape fields for mesh movement (Laplace solution for each body)
     !> base_shapes(i) is 1.0 on body i and 0.0 on all others.
     type(field_t), allocatable :: base_shapes(:)

     !> Sum of all base shapes. Should be \f$ \in [0, 1] \f$.
     type(field_t) :: phi_total

     real(kind=rp), pointer :: global_pivot_pos(:) => null()
     real(kind=rp), pointer :: global_pivot_vel_lag(:, :) => null()

     ! Basis Vectors for orientation
     real(kind=rp), pointer :: global_basis_pos(:) => null()
     ! Store history for the ghost trackers
     real(kind=rp), pointer :: global_basis_vel_lag(:, :) => null()
     ! Private handles to the ghost trackers (2 per body)
     integer, allocatable :: ghost_handles(:,:)
     ! Rotation matrices
     real(kind=rp), allocatable :: body_rot_matrices(:,:,:)

     type(point_tracker_t), allocatable :: trackers(:)
     integer :: n_trackers = 0

     procedure(user_ale_mesh_velocity_intf), nopass, pointer :: &
          user_ale_mesh_vel => null()
     procedure(user_ale_base_shapes_intf), nopass, pointer :: &
          user_ale_base_shapes => null()
     procedure(user_ale_rigid_kinematics_intf), nopass, pointer :: &
          user_ale_rigid_kinematics => null()

   contains
     procedure, pass(this) :: init => ale_manager_init
     procedure, pass(this) :: free => ale_manager_free
     procedure, pass(this) :: mesh_preview
     procedure, pass(this) :: solve_base_mesh_displacement
     procedure, pass(this) :: advance_mesh
     procedure, pass(this) :: update_mesh_velocity
     procedure, pass(this) :: set_pivot_restart
     procedure, pass(this) :: sync_chkp
     procedure, pass(this) :: request_tracker
     procedure, pass(this) :: get_tracker_pos
     procedure, pass(this) :: compute_rotation_matrix
     procedure, pass(this) :: prep_checkpoint => set_pivot_basis_for_checkpoint
     procedure, pass(this) :: ghost_tracker_coord_step
     procedure, pass(this) :: log_rot_angles
     procedure, pass(this) :: log_pivot
     procedure, pass(this) :: register_checkpoint_fields
  end type ale_manager_t

  type(ale_manager_t), public, pointer :: neko_ale => null()

contains

  !> Initialize ALE Manager
  !> Sets up solver, registers fields, solves for base shape, etc.
  subroutine ale_manager_init(this, coef, json, user, chkp)
    class(ale_manager_t), intent(inout), target :: this
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: json
    type(user_t), intent(in) :: user
    type(chkp_t), intent(inout) :: chkp
    type(json_file) :: body_sub, bc_subdict
    type(json_file) :: precon_params
    type(time_state_t) :: t_init
    integer, allocatable :: zone_indices(:)
    integer :: time_order
    integer :: n_moving_zones
    integer :: z, tmp_int, ksp_max_iter
    integer, allocatable :: moving_zone_ids(:)
    integer :: i, j, k, n_bcs, n, n_bodies
    real(kind=rp), allocatable :: tmp_vec(:)
    real(kind=rp) :: tmp_val, abstol
    character(len=128) :: log_buf
    character(len=256) :: log_buf_l
    character(len=:), allocatable :: bc_type
    character(len=:), allocatable :: tmp_str
    character(len=:), allocatable :: ksp_solver
    character(len=:), allocatable :: precon_type
    logical :: tmp_logical, oifs
    logical :: moving_
    logical :: found_zone
    logical :: has_user_rigid_kin, has_user_mesh_vel
    logical :: has_builtin_osc, has_builtin_rot, is_rot_active
    logical :: res_monitor, import_base_shapes

    if (json%valid_path('case.fluid.ale')) then
       call json_get(json, 'case.fluid.ale.enabled', this%active)
    end if
    call json_get_or_default(json, 'case.numerics.oifs', oifs, .false.)

    if (.not. this%active) then
       neko_ale => null()
       return
    else if (this%active) then
       ! force all elements as deformed when mesh changes.
       call coef%msh%all_deformed()

       if (NEKO_BCKND_DEVICE .eq. 1) then
          if ((.not. (NEKO_BCKND_HIP .eq. 1)) .and. &
               (.not. (NEKO_BCKND_CUDA .eq. 1))) then
             call neko_error("ALE currently " // &
                  "supported only with HIP or CUDA backend.")
          end if
       end if
       if (oifs) then
          call neko_error("ALE not currently supported with OIFS.")
       end if
       if (json%valid_path('case.checkpoint_format')) then
          call json_get(json, 'case.checkpoint_format', tmp_str)
          if (trim(tmp_str) /= 'chkp') then
             call neko_error("ALE is not supported with the '" // &
                  trim(tmp_str) // &
                  "' checkpoint format. Please use 'chkp'.")
          end if
       end if
       neko_ale => this
    end if

    call neko_log%section("ALE Initialization")
    call neko_log%message(" ")

    if (NEKO_BCKND_HIP .eq. 1) then
       call neko_log%message("Initializing ALE " // &
            "with device backend (HIP).")
    else if (NEKO_BCKND_CUDA .eq. 1) then
       call neko_log%message("Initializing ALE " // &
            "with device backend (CUDA).")
    else
       call neko_log%message("Initializing ALE " // &
            "with CPU backend.")
    end if

    tmp_logical = .false.
    n = coef%dof%size()

    call this%x_ref%init(coef%dof, "x_ref")
    call this%y_ref%init(coef%dof, "y_ref")
    call this%z_ref%init(coef%dof, "z_ref")

    call copy(this%x_ref%x, coef%dof%x, n)
    call copy(this%y_ref%x, coef%dof%y, n)
    call copy(this%z_ref%x, coef%dof%z, n)

    ! Sync to device
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%x_ref%copy_from(HOST_TO_DEVICE, .false.)
       call this%y_ref%copy_from(HOST_TO_DEVICE, .false.)
       call this%z_ref%copy_from(HOST_TO_DEVICE, .true.)
    end if

    ! Set user function pointers.
    this%user_ale_mesh_vel => user%ale_mesh_velocity
    this%user_ale_base_shapes => user%ale_base_shapes
    this%user_ale_rigid_kinematics => user%ale_rigid_kinematics

    ! Check user association states
    has_user_rigid_kin = .not. associated(this%user_ale_rigid_kinematics, &
         dummy_user_ale_rigid_kinematics)
    has_user_mesh_vel = .not. associated(this%user_ale_mesh_vel, &
         dummy_user_ale_mesh_velocity)

    ! Enable B history (Blag, Blaglag)
    call coef%enable_B_history()
    call json_get(json, 'case.numerics.time_order', time_order)

    ! Stuff for zone_id checks
    n_moving_zones = 0
    if (allocated(moving_zone_ids)) deallocate(moving_zone_ids)
    allocate(moving_zone_ids(0))

    ! Register mesh velocity fields
    call neko_registry%add_field(coef%dof, 'wm_x')
    call neko_registry%add_field(coef%dof, 'wm_y')
    call neko_registry%add_field(coef%dof, 'wm_z')
    this%wm_x => neko_registry%get_field('wm_x')
    this%wm_y => neko_registry%get_field('wm_y')
    this%wm_z => neko_registry%get_field('wm_z')

    call get_ale_solver_params_json(this, json, ksp_solver, precon_type, &
         precon_params, abstol, ksp_max_iter, res_monitor, import_base_shapes)

    ! Mark BCs
    call this%bc_moving%init_from_components(coef)
    call this%bc_fixed%init_from_components(coef)

    if (json%valid_path('case.fluid.boundary_conditions')) then
       call json%info('case.fluid.boundary_conditions', n_children = n_bcs)

       do i = 1, n_bcs
          call json_extract_item(json, 'case.fluid.boundary_conditions', &
               i, bc_subdict)

          if (allocated(bc_type)) deallocate(bc_type)
          call json_get(bc_subdict, 'type', bc_type)

          if (allocated(zone_indices)) deallocate(zone_indices)
          call json_get(bc_subdict, 'zone_indices', zone_indices)

          moving_ = .false.
          if (trim(bc_type) .eq. 'no_slip') then
             call json_get_or_default(bc_subdict, 'moving', moving_, .false.)
          end if

          if (moving_) then
             do j = 1, size(zone_indices)
                ! we append unique moving zone ids for future checks
                call append_unique_int(moving_zone_ids, n_moving_zones, &
                     zone_indices(j))
                call this%bc_moving%mark_zone(coef%msh%labeled_zones(&
                     zone_indices(j)))
             end do
             this%has_moving_boundary = .true.
          else
             do j = 1, size(zone_indices)
                call this%bc_fixed%mark_zone(coef%msh%labeled_zones(&
                     zone_indices(j)))
             end do
          end if
       end do
    end if

    call this%bc_moving%finalize()
    call this%bc_fixed%finalize()
    call this%bc_list%init()
    call this%bc_list%append(this%bc_moving)
    call this%bc_list%append(this%bc_fixed)

    ! Mesh Stiffness
    if (json%valid_path('case.fluid.ale.solver.mesh_stiffness.type')) then
       call json%get('case.fluid.ale.solver.mesh_stiffness.type', tmp_str)
       this%config%stiffness_type = tmp_str
       if (.not. (trim(tmp_str) .eq. 'built-in')) then
          call neko_error("ALE: stiffness_type must be 'built-in'")
       end if
    end if

    if ( associated(this%user_ale_base_shapes, &
         dummy_user_ale_base_shapes) .and. (.not. import_base_shapes)) then
       call neko_log%message('Solver Type       : (' // &
            trim(ksp_solver) // ', ' // trim(precon_type) // ')')
       write(log_buf, '(A,ES13.6)') 'Abs tol           :', abstol
       call neko_log%message(log_buf)
       call neko_log%message('Mesh Stiffness    : ' // &
            trim(this%config%stiffness_type))
    end if
    call neko_log%message(' ')

    ! Bodies
    if (json%valid_path('case.fluid.ale.bodies')) then
       call json%info('case.fluid.ale.bodies', n_children = n_bodies)
       this%config%nbodies = n_bodies
       allocate(this%config%bodies(n_bodies))
       allocate(this%ale_pivot(n_bodies))
       allocate(this%body_kin(n_bodies))
       allocate(this%base_shapes(n_bodies))
       allocate(this%global_pivot_pos(3 * this%config%nbodies))
       allocate(this%global_pivot_vel_lag(3 * this%config%nbodies, 3))
       allocate(this%global_basis_pos(6 * this%config%nbodies))
       allocate(this%ghost_handles(2, this%config%nbodies))
       allocate(this%global_basis_vel_lag(6 * this%config%nbodies, 3))
       allocate(this%body_rot_matrices(3, 3, this%config%nbodies))

       this%global_pivot_pos = 0.0_rp
       this%global_pivot_vel_lag = 0.0_rp
       this%global_basis_pos = 0.0_rp
       this%global_basis_vel_lag = 0.0_rp
       this%body_rot_matrices = 0.0_rp

       do i = 1, n_bodies
          this%body_rot_matrices(1, 1, i) = 1.0_rp
          this%body_rot_matrices(2, 2, i) = 1.0_rp
          this%body_rot_matrices(3, 3, i) = 1.0_rp
       end do

       do i = 1, n_bodies
          call json_extract_item(json, 'case.fluid.ale.bodies', i, body_sub)
          this%config%bodies(i)%id = i

          if (body_sub%valid_path('name')) then
             call json_get(body_sub, 'name', tmp_str)
             this%config%bodies(i)%name = tmp_str
          else
             write(this%config%bodies(i)%name, '(A,I0)') 'body_', i
          end if

          if (body_sub%valid_path('zone_indices')) then
             call json_get(body_sub, 'zone_indices', zone_indices)
             this%config%bodies(i)%zone_indices = zone_indices
          else
             call neko_error("ALE: body " // &
                  trim(this%config%bodies(i)%name) // &
                  " must have 'zone_indices'")
          end if

          ! Oscillation
          this%config%bodies(i)%osc_amp = 0.0_rp
          this%config%bodies(i)%osc_freq = 0.0_rp
          if (body_sub%valid_path('oscillation')) then
             call json_get(body_sub, 'oscillation.amplitude', tmp_vec, &
                  expected_size = 3)
             this%config%bodies(i)%osc_amp = tmp_vec
             call json_get(body_sub, 'oscillation.frequency', tmp_vec, &
                  expected_size = 3)
             this%config%bodies(i)%osc_freq = tmp_vec
          end if

          ! Rotation
          if (body_sub%valid_path('rotation')) then
             ! Check if pivot exists.
             if (.not. body_sub%valid_path('pivot')) then
                call neko_error("ale.bodies.pivot is missing " // &
                     "from the case file.")
             end if

             call json_get(body_sub, 'rotation.type', tmp_str)
             this%config%bodies(i)%rotation_type = tmp_str

             select case (trim(tmp_str))
             case ('harmonic')
                call json_get(body_sub, 'rotation.amplitude_deg', tmp_vec, &
                     expected_size = 3)
                this%config%bodies(i)%rot_amp_degree = tmp_vec

                call json_get(body_sub, 'rotation.frequency', tmp_vec, &
                     expected_size = 3)
                this%config%bodies(i)%rot_freq = tmp_vec


             case ('ramp')
                call json_get(body_sub, 'rotation.ramp_t0', tmp_vec, &
                     expected_size = 3)
                this%config%bodies(i)%ramp_t0 = tmp_vec

                call json_get(body_sub, 'rotation.ramp_omega0', tmp_vec, &
                     expected_size = 3)
                this%config%bodies(i)%ramp_omega0 = tmp_vec


             case ('smooth_step')
                call json_get_or_default(body_sub, 'rotation.axis', &
                     tmp_int, 3)
                if (tmp_int .ge. 1 .and. tmp_int .le. 3) then
                   this%config%bodies(i)%rotation_axis = tmp_int
                else
                   call neko_error("ALE: rotation.axis must be (integer) " // &
                        "1 -> x, 2 -> y, or 3 -> z")
                end if
                call json_get(body_sub, 'rotation.step_control_times', &
                     tmp_vec, expected_size = 4)
                this%config%bodies(i)%step_control_times = tmp_vec

                call json_get(body_sub, 'rotation.target_angle_deg', tmp_val)
                this%config%bodies(i)%target_rot_angle_deg = tmp_val

             case default
                call neko_error("ALE: rotation.type must be 'harmonic', " // &
                     "'ramp', or 'smooth_step'")
             end select
          end if

          ! Rotation Center
          if (body_sub%valid_path('pivot')) then
             call json_get_or_default(body_sub, 'pivot.type', tmp_str, &
                  'relative')
             this%config%bodies(i)%rotation_center_type = tmp_str
             call json_get(body_sub, 'pivot.value', tmp_vec, expected_size = 3)
             this%config%bodies(i)%rot_center = tmp_vec


             tmp_str = this%config%bodies(i)%rotation_center_type
             if (trim(tmp_str) /= 'relative' .and. &
                  trim(tmp_str) /= 'relative_sin') then
                call neko_error("ALE: pivot.type must be " // &
                     "'relative', or 'relative_sin'.")
             end if
          end if

          ! Stiff Geom
          if (body_sub%valid_path('stiff_geom')) then
             call json_get(body_sub, 'stiff_geom.type', tmp_str)
             this%config%bodies(i)%stiff_geom%type = tmp_str
             call json_get(body_sub, 'stiff_geom.gain', &
                  this%config%bodies(i)%stiff_geom%gain)
             call json_get(body_sub, 'stiff_geom.decay_profile', tmp_str)
             this%config%bodies(i)%stiff_geom%decay_profile = tmp_str

             select case (trim(this%config%bodies(i)%stiff_geom%decay_profile))
             case ('gaussian')
                call json_get_or_default(body_sub, &
                     'stiff_geom.cutoff_coef', &
                     this%config%bodies(i)%stiff_geom%cutoff_coef, 9.0_rp)
             case ('tanh')
                call json_get_or_default(body_sub, &
                     'stiff_geom.cutoff_coef', &
                     this%config%bodies(i)%stiff_geom%cutoff_coef, 3.5_rp)
             case default
                call neko_error("ALE: Invalid stiff_geom.decay_profile: " // &
                     trim(this%config%bodies(i)%stiff_geom%decay_profile))
             end select

             select case (trim(this%config%bodies(i)%stiff_geom%type))
             case ('cylinder', 'sphere')
                call json_get(body_sub, 'stiff_geom.center', tmp_vec, &
                     expected_size = 3)
                this%config%bodies(i)%stiff_geom%center = tmp_vec

                call json_get(body_sub, 'stiff_geom.radius', &
                     this%config%bodies(i)%stiff_geom%radius)
             case ('cheap_dist')
                call json_get(body_sub, 'stiff_geom.stiff_dist', &
                     this%config%bodies(i)%stiff_geom%stiff_dist)
             case ('box')
                call neko_error("ALE: stiff_geom.type 'box' not yet" // &
                     " implemented.")
             case default
                call neko_error("ALE: Invalid stiff_geom.type: " // &
                     trim(this%config%bodies(i)%stiff_geom%type))
             end select
          elseif (import_base_shapes) then
             ! do nothing.
          else
             call neko_error("ALE: Body '" // &
                  trim(this%config%bodies(i)%name) // &
                  "' must have 'stiff_geom' definition.")
          end if

          ! Initialize the pivots.
          call init_pivot_state(this%ale_pivot(i), this%config%bodies(i))

          call this%base_shapes(i)%init(coef%dof, &
               "phi_" // trim(this%config%bodies(i)%name))
          call field_rzero(this%base_shapes(i))

          ! Create Ghost Trackers for numerically forming the rotation matrix
          ! of each body.
          ! Basis X (Pivot + 1.0 in X)
          this%ghost_handles(1, i) = this%request_tracker( &
               this%config%bodies(i)%rot_center + [1.0_rp, 0.0_rp, 0.0_rp], &
               this%config%bodies(i)%id)
          ! Basis Y (Pivot + 1.0 in Y)
          this%ghost_handles(2, i) = this%request_tracker( &
               this%config%bodies(i)%rot_center + [0.0_rp, 1.0_rp, 0.0_rp], &
               this%config%bodies(i)%id)

          call neko_log%message('Registered Body : ' // &
               trim(this%config%bodies(i)%name))

          ! Logging Stiff Body
          call neko_log%message(' ')
          if (associated(this%user_ale_base_shapes, &
               dummy_user_ale_base_shapes) .and. &
               (.not. import_base_shapes)) then
             write(log_buf, '(A,A)') '   Stiff Type    : ', &
                  trim(this%config%bodies(i)%stiff_geom%type)
             call neko_log%message(log_buf)
             write(log_buf, '(A,ES18.11,A,A,A,ES10.4)') '    Gain      : ', &
                  this%config%bodies(i)%stiff_geom%gain, ' | Profile: ', &
                  trim(this%config%bodies(i)%stiff_geom%decay_profile), &
                  ' | Cutoff: ', this%config%bodies(i)%stiff_geom%cutoff_coef
             call neko_log%message(log_buf)
             select case (trim(this%config%bodies(i)%stiff_geom%type))
             case ('cylinder', 'sphere')
                write(log_buf, '(A,3(ES23.15,1X))') '    Center       :', &
                     this%config%bodies(i)%stiff_geom%center
                call neko_log%message(log_buf)
                write(log_buf, '(A,ES23.15)') '    Radius       :', &
                     this%config%bodies(i)%stiff_geom%radius
                call neko_log%message(log_buf)
             case ('cheap_dist')
                write(log_buf, '(A,ES23.15)') '    Stiff Dist:', &
                     this%config%bodies(i)%stiff_geom%stiff_dist
                call neko_log%message(log_buf)
             end select
          end if
          call neko_log%message(' ')

          ! Logging Oscillation
          has_builtin_osc = any(abs(this%config%bodies(i)%osc_amp) .gt. 0.0_rp)

          if (has_builtin_osc) then
             if (has_user_rigid_kin .or. has_user_mesh_vel) then
                call neko_log%message('   Oscillation    : ' // &
                     'X(t) = Amp*sin(2*pi*Freq*t) + User')
                write(log_buf, '(A,3(ES18.11,1X))') '    Amp       :', &
                     this%config%bodies(i)%osc_amp
                call neko_log%message(log_buf)
                write(log_buf, '(A,3(ES18.11,1X))') '    Freq      :', &
                     this%config%bodies(i)%osc_freq
                call neko_log%message(log_buf)
             else
                call neko_log%message('   Oscillation    : ' // &
                     'X(t) = Amp*sin(2*pi*Freq*t)')
                write(log_buf, '(A,3(ES18.11,1X))') '    Amp       :', &
                     this%config%bodies(i)%osc_amp
                call neko_log%message(log_buf)
                write(log_buf, '(A,3(ES18.11,1X))') '    Freq      :', &
                     this%config%bodies(i)%osc_freq
                call neko_log%message(log_buf)
             end if
          else
             if (has_user_rigid_kin .or. has_user_mesh_vel) then
                call neko_log%message('   Oscillation    : User-defined')
             else
                call neko_log%message('   Oscillation    : None')
             end if
          end if
          call neko_log%message(' ')

          ! Logging Rotation
          has_builtin_rot = (trim(this%config%bodies(i)%rotation_type) &
               /= 'user')

          if (trim(this%config%bodies(i)%rotation_type) .eq. 'user') then

             call neko_log%message('   Rotation Type: User-defined')

          elseif (has_builtin_rot) then

             ! Check parameters active
             is_rot_active = .false.
             select case (trim(this%config%bodies(i)%rotation_type))
             case ('harmonic')
                is_rot_active = any(abs(this%config%bodies(i)%rot_amp_degree) &
                     .gt. 0.0_rp)
             case ('ramp')
                is_rot_active = any(abs(this%config%bodies(i)%ramp_omega0) &
                     .gt. 0.0_rp)
             case ('smooth_step')
                is_rot_active = &
                     (abs(this%config%bodies(i)%target_rot_angle_deg) &
                     .gt. 0.0_rp)
             end select

             if (is_rot_active) then
                ! Harmonic
                if (trim(this%config%bodies(i)%rotation_type) &
                     .eq. 'harmonic') then
                   if (has_user_rigid_kin .or. has_user_mesh_vel) then
                      call neko_log%message('   Rotation     : ' // &
                           'Theta(t) = Amp*sin(2*pi*Freq*t) + User')
                   else
                      call neko_log%message('   Rotation     : ' // &
                           'Theta(t) = Amp*sin(2*pi*Freq*t)')
                   end if
                   write(log_buf, '(A,3(ES18.11,1X))') '    Amp (deg) :', &
                        this%config%bodies(i)%rot_amp_degree
                   call neko_log%message(log_buf)
                   write(log_buf, '(A,3(ES18.11,1X))') '    Freq      :', &
                        this%config%bodies(i)%rot_freq
                   call neko_log%message(log_buf)

                   ! Ramp
                elseif (trim(this%config%bodies(i)%rotation_type) &
                     .eq. 'ramp') then
                   if (has_user_rigid_kin .or. has_user_mesh_vel) then
                      call neko_log%message('   Rotation     : ' // &
                           'Omega(t) = Omega0*(1 - exp(-4.6*t/t0)) + User')
                   else
                      call neko_log%message('   Rotation     : ' // &
                           'Omega(t) = Omega0*(1 - exp(-4.6*t/t0))')
                   end if
                   write(log_buf, '(A,3(ES18.11,1X))') '    Omega0    :', &
                        this%config%bodies(i)%ramp_omega0
                   call neko_log%message(log_buf)
                   write(log_buf, '(A,3(ES18.11,1X))') '    t0        :', &
                        this%config%bodies(i)%ramp_t0
                   call neko_log%message(log_buf)

                   ! Smooth Step
                elseif (trim(this%config%bodies(i)%rotation_type) &
                     .eq. 'smooth_step') then
                   if (has_user_rigid_kin .or. has_user_mesh_vel) then
                      call neko_log%message('   Rotation     : ' // &
                           'Smooth Step Control + User')
                   else
                      call neko_log%message('   Rotation     : ' // &
                           'Smooth Step Control')
                   end if
                   write(log_buf, '(A,I10)') '    Rotation Axis    :', &
                        this%config%bodies(i)%rotation_axis
                   call neko_log%message(log_buf)
                   write(log_buf, '(A,ES18.11)') '    Target Rot ' // &
                        'Angle (deg)  :', &
                        this%config%bodies(i)%target_rot_angle_deg
                   call neko_log%message(log_buf)
                   write(log_buf, '(A,4(ES18.11,1X))') &
                        '        Control Times [t0, t1, t2, t3]    :', &
                        this%config%bodies(i)%step_control_times
                   call neko_log%message(log_buf)
                end if
             else
                if (has_user_rigid_kin .or. has_user_mesh_vel) then
                   call neko_log%message('   Rotation Type: User-defined')
                else
                   call neko_log%message('   Rotation Type: None')
                end if
             end if

          end if

          ! Logging Pivot
          call neko_log%message(' ')
          call neko_log%message('   Pivot Type    : ' // &
               trim(this%config%bodies(i)%rotation_center_type))

          write(log_buf, '(A,3(ES18.11,1X))') '    Init Pivot:', &
               this%config%bodies(i)%rot_center
          call neko_log%message(log_buf)
          call neko_log%message(' ')

       end do
    else
       call neko_error("ALE: No 'ale bodies' found in case file!")
    end if

    if (this%config%nbodies .gt. 1 .and. (.not. import_base_shapes)) then
       call this%phi_total%init(coef%dof, "phi_total")
       call field_rzero(this%phi_total)
    end if

    ! Check to be sure moving no_slip ids belong to an ALE body
    do i = 1, n_moving_zones
       z = moving_zone_ids(i)
       found_zone = .false.
       j = 1
       do while ((.not. found_zone) .and. (j .le. this%config%nbodies))
          if (any(this%config%bodies(j)%zone_indices .eq. z)) then
             found_zone = .true.
          end if
          j = j + 1
       end do
       if (.not. found_zone) then
          write(log_buf_l, '(A,I0,A)') &
               "ALE: zone index ", z, &
               " has BC no_slip with moving: true, " // &
               "but it is not registered in ALE bodies."
          call neko_error(trim(log_buf_l))
       end if
    end do

    ! Any id registered in ALE bodies must have
    ! no_slip with moving: true in BCs.
    do j = 1, this%config%nbodies
       if (allocated(this%config%bodies(j)%zone_indices)) then
          do i = 1, size(this%config%bodies(j)%zone_indices)
             z = this%config%bodies(j)%zone_indices(i)
             found_zone = .false.
             if (n_moving_zones .gt. 0) then
                if (any(moving_zone_ids(1:n_moving_zones) .eq. z)) then
                   found_zone = .true.
                end if
             end if
             if (.not. found_zone) then
                write(log_buf_l, '(A,I0,A,A)') &
                     "ALE: zone index ", z, &
                     " is registered in ALE bodies, ", &
                     "but the BC is not no_slip with moving: true."
                call neko_error(trim(log_buf_l))
             end if
          end do
       end if
    end do

    ! Check no zone ID is assigned to more than one ALE body.
    do j = 1, this%config%nbodies
       if (allocated(this%config%bodies(j)%zone_indices)) then
          do i = 1, size(this%config%bodies(j)%zone_indices)
             z = this%config%bodies(j)%zone_indices(i)

             do k = j + 1, this%config%nbodies
                if (allocated(this%config%bodies(k)%zone_indices)) then
                   if (any(this%config%bodies(k)%zone_indices .eq. z)) then
                      write(log_buf_l, '(A,I0,A,A,A,A,A)') &
                           "ALE: zone index ", z, &
                           " is assigned to multiple bodies ('", &
                           trim(this%config%bodies(j)%name), "' and '", &
                           trim(this%config%bodies(k)%name), "')."
                      call neko_error(trim(log_buf_l))
                   end if
                end if
             end do

          end do
       end if
    end do

    ! Find the smooth blending function for mesh displacement.
    call this%solve_base_mesh_displacement(coef, json, import_base_shapes, &
         abstol, ksp_solver, ksp_max_iter, &
         precon_type, precon_params, res_monitor)

    ! If we are restarting, we skip this. It will be handled
    ! properly by chkp file.
    if (.not. json%valid_path('case.restart_file')) then
       t_init%t = 0.0_rp
       t_init%tstep = 0
       t_init%dt = 0.0_rp
       call this%update_mesh_velocity(coef, t_init)
    end if

    call this%wm_x_lag%init(this%wm_x, 2)
    call this%wm_y_lag%init(this%wm_y, 2)
    call this%wm_z_lag%init(this%wm_z, 2)

    if (allocated(moving_zone_ids)) deallocate(moving_zone_ids)
    if (allocated(bc_type)) deallocate(bc_type)
    if (allocated(zone_indices)) deallocate(zone_indices)
    if (allocated(ksp_solver)) deallocate(ksp_solver)
    if (allocated(precon_type)) deallocate(precon_type)
    if (allocated(tmp_str)) deallocate(tmp_str)
    if (allocated(tmp_vec)) deallocate(tmp_vec)

    ! Performing mesh_preview.
    call this%mesh_preview(coef, json)

    ! Register checkpoint fields
    call this%register_checkpoint_fields(coef, chkp)

    call neko_log%end_section()
  end subroutine ale_manager_init

  !> Solves the Laplace equation to determine the base shape (phi)
  !> for each body.
  !> It finds a smooth blending function for mesh deformation.
  !> For body i: phi_i = 1 on body i zones, phi_i = 0 on all other boundaries.
  subroutine solve_base_mesh_displacement(this, coef, json, &
       import_base_shapes, abstol, ksp_solver, ksp_max_iter, precon_type, &
       precon_params, res_monitor)
    class(ale_manager_t), intent(inout), target :: this
    class(ax_t), allocatable :: Ax
    class(ksp_t), allocatable :: ksp
    class(pc_t), allocatable :: pc
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: json
    logical, intent(in) :: import_base_shapes
    real(kind=rp), intent(in) :: abstol
    logical, intent(in) :: res_monitor
    character(len=*), intent(in) :: ksp_solver, precon_type
    integer, intent(in) :: ksp_max_iter
    type(json_file), intent(inout) :: precon_params
    type(file_t) :: phi_file
    type(field_t), pointer :: phi_ptr => null()
    type(field_t) :: rhs_field
    type(field_t) :: corr_field
    type(ksp_monitor_t) :: monitor(1)
    real(kind=rp) :: sample_start_time, sample_end_time
    real(kind=rp) :: sample_time
    character(len=LOG_SIZE) :: log_buf
    integer :: n, i, m, k, ierr, body_idx, z_idx
    integer :: j
    real(kind=rp), allocatable :: h1_restore(:, :, :, :)
    real(kind=rp), allocatable :: h2_restore(:, :, :, :)
    type(zero_dirichlet_t) :: bc_active_body
    type(zero_dirichlet_t) :: bc_inactive_body
    type(scalar_bc_projector_t) :: bc_projector
    type(scalar_bc_projector_t) :: bc_projector_zeros_only
    type(json_file) :: body_sub
    character(len=256) :: phi_fname
    character(len=:), allocatable :: tmp_str


    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return
    if (this%config%nbodies .eq. 0) return

    if (import_base_shapes) then
       call neko_log%message(" ")
       call neko_log%message("Importing ALE base shapes" // &
            " (skipping Laplace solve)...")

       do body_idx = 1, this%config%nbodies

          call json_extract_item(json, 'case.fluid.ale.bodies', &
               body_idx, body_sub)

          call json_get(body_sub, 'base_shape_import_file', tmp_str)
          phi_fname = tmp_str

          phi_ptr => this%base_shapes(body_idx)

          ! Load the field
          call import_fields(fname = trim(phi_fname), p = phi_ptr)

          call neko_log%message("   Loaded: " // &
               trim(phi_fname) // &
               " for body: " // &
               trim(this%config%bodies(body_idx)%name))
       end do

       return
    end if

    call neko_log%message(" ")
    call neko_log%message("Starting base mesh motion solve ...")
    n = coef%dof%size()

    call ax_helm_allocator(Ax, type_name = "standard")
    call krylov_solver_factory(ksp, n, ksp_solver, &
         ksp_max_iter, abstol, monitor = res_monitor)
    call ale_precon_factory(pc, ksp, coef, coef%dof, &
         coef%gs_h, this%bc_list, precon_type, precon_params)

    ! Save original h1/h2
    h1_restore = coef%h1
    h2_restore = coef%h2

    call rhs_field%init(coef%dof)
    call corr_field%init(coef%dof)


    ! User Defined Base Shapes (Skip Solver).
    if (.not. associated(this%user_ale_base_shapes, &
         dummy_user_ale_base_shapes)) then
       call neko_log%message("   Using user-defined base shapes " // &
            "(skipping Laplace solve)")

       ! Call User Hook (Populates this%base_shapes)
       call this%user_ale_base_shapes(this%base_shapes)

       ! Compute phi_total (Sum of all user shapes)
       if (this%config%nbodies .gt. 1) then
          call field_rzero(this%phi_total)
          do body_idx = 1, this%config%nbodies
             call field_add2(this%phi_total, this%base_shapes(body_idx), n)
          end do
       end if

       ! Output Shapes
       if (this%config%if_output_phi) then
          ! Individual Bodies
          do body_idx = 1, this%config%nbodies
             call phi_file%init('phi_' // &
                  trim(this%config%bodies(body_idx)%name) // '.fld', &
                  precision = rp)
             select type (ft => phi_file%file_type)
             type is (fld_file_t)
                ft%skip_pressure = .false.
             end select
             call phi_file%write(this%base_shapes(body_idx))
             call phi_file%free()
             call neko_log%message('   phi_' // &
                  trim(this%config%bodies(body_idx)%name) // '.fld saved.')
          end do

          ! Total
          if (this%config%nbodies .gt. 1) then
             call neko_log%message("  phi_total.fld saved.")
             select type (ft => phi_file%file_type)
             type is (fld_file_t)
                ft%skip_pressure = .false.
             end select
             call phi_file%init('phi_total.fld', precision = rp)
             call phi_file%write(this%phi_total)
             call phi_file%free()
          end if
       end if
    else
       ! Standard Laplace Solve (Requires Stiffness)

       ! Compute Stiffness
       call compute_stiffness_ale(coef, this%config)

       ! Output Stiffness if requested (for diagnostic)
       if (this%config%if_output_stiffness) then
          rhs_field%x = coef%h1
          call phi_file%init('stiffness.fld')
          call phi_file%write(rhs_field)
          call phi_file%free()
          call field_rzero(rhs_field)
       end if

       ! Loop over bodies and Solve Laplace
       do body_idx = 1, this%config%nbodies
          call MPI_Barrier(NEKO_COMM, ierr)
          sample_start_time = MPI_WTIME()
          call neko_log%message(" Solving laplace for body: " // &
               trim(this%config%bodies(body_idx)%name))

          call bc_active_body%init_from_components(coef)
          call bc_inactive_body%init_from_components(coef)

          ! Mark zones
          do j = 1, size(this%config%bodies(body_idx)%zone_indices)
             z_idx = this%config%bodies(body_idx)%zone_indices(j)
             call bc_active_body%mark_zone(coef%msh%labeled_zones(z_idx))
          end do

          do i = 1, this%config%nbodies
             if (i /= body_idx) then
                do j = 1, size(this%config%bodies(i)%zone_indices)
                   z_idx = this%config%bodies(i)%zone_indices(j)
                   call bc_inactive_body%mark_zone(&
                        coef%msh%labeled_zones(z_idx))
                end do
             end if
          end do

          call bc_active_body%finalize()
          call bc_inactive_body%finalize()

          ! The Full list for the solver (Freeze everything to 0 correction)
          call bc_projector%mark(this%bc_fixed)
          call bc_projector%mark(bc_active_body)
          call bc_projector%mark(bc_inactive_body)

          ! The "Zeros Only" list for the field (Reset other boundaries)
          call bc_projector_zeros_only%mark(this%bc_fixed)
          call bc_projector_zeros_only%mark(bc_inactive_body)

          call field_rzero(this%base_shapes(body_idx))
          this%base_shapes(body_idx)%x = 0.0_rp
          rhs_field%x = 0.0_rp
          corr_field%x = 0.0_rp
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! phi = phi_corr + phi_lifted!
          ! A*phi_corr = -A*phi_lifted !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Lift BC (Dirichlet = 1.0 on moving body)
          m = bc_active_body%msk(0)
          do i = 1, m
             k = bc_active_body%msk(i)
             this%base_shapes(body_idx)%x(k, 1, 1, 1) = 1.0_rp
          end do

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_memcpy(this%base_shapes(body_idx)%x, &
                  this%base_shapes(body_idx)%x_d, n, HOST_TO_DEVICE, .true.)
          end if

          ! Apply Zeros to others.
          ! This ensures fixed walls and other bodies are 0.0,
          ! even if they share grid with a moving wall.
          call bc_projector_zeros_only%apply(this%base_shapes(body_idx)%x, n)

          ! Compute RHS: RHS = -A * Phi_lifted.
          ! The following is motivated by implementation in Nek5000.
          call Ax%compute(rhs_field%x, this%base_shapes(body_idx)%x, &
               coef, coef%msh, coef%Xh)
          call field_cmult(rhs_field, -1.0_rp)

          ! Here we use the FULL list to apply zero Dirichlet BC
          ! on all boundaries.
          call bc_projector%apply(rhs_field%x, n)
          call coef%gs_h%op(rhs_field, GS_OP_ADD)

          ! Solve
          call field_rzero(corr_field)
          call pc%update()
          monitor(1) = ksp%solve(Ax, corr_field, &
               rhs_field%x, n, coef, bc_projector, coef%gs_h)

          ! phi = phi_lifted + phi_corr
          call field_add2(this%base_shapes(body_idx), corr_field, n)

          ! Update Total Phi
          ! phi_total should be between 0 and 1.
          if (this%config%nbodies .gt. 1) then
             call field_add2(this%phi_total, this%base_shapes(body_idx), n)
          end if

          call MPI_Barrier(NEKO_COMM, ierr)
          sample_end_time = MPI_WTIME()
          sample_time = sample_end_time - sample_start_time
          write(log_buf, '(A, A, A, ES11.4, A)') "   Laplace solve for '", &
               trim(this%config%bodies(body_idx)%name), "' took ", &
               sample_time, " (s)"

          call neko_log%message(log_buf)

          call bc_active_body%free()
          call bc_inactive_body%free()
          call bc_projector%free()
          call bc_projector_zeros_only%free()

          ! We let the host to also have the base_shapes so in
          ! user_ale_mesh_vel
          ! it would be easier in general to use it.
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_memcpy(this%base_shapes(body_idx)%x, &
                  this%base_shapes(body_idx)%x_d, n, DEVICE_TO_HOST, .true.)
          end if

          if (this%config%if_output_phi) then
             call phi_file%init('phi_' // &
                  trim(this%config%bodies(body_idx)%name) // '.fld', &
                  precision = rp)
             select type (ft => phi_file%file_type)
             type is (fld_file_t)
                ft%skip_pressure = .false.
             end select
             call phi_file%write(this%base_shapes(body_idx))
             call phi_file%free()
             call neko_log%message('      phi_' // &
                  trim(this%config%bodies(body_idx)%name) // '.fld saved.')
          end if
       end do

       if (this%config%if_output_phi .and. (this%config%nbodies .gt. 1)) then

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_memcpy(this%phi_total%x, this%phi_total%x_d, n, &
                  DEVICE_TO_HOST, .true.)
          end if

          call neko_log%message("   phi_total.fld saved.")
          call phi_file%init('phi_total.fld', precision = rp)
          select type (ft => phi_file%file_type)
          type is (fld_file_t)
             ft%skip_pressure = .false.
          end select
          call phi_file%write(this%phi_total)
          call phi_file%free()
       end if
    end if

    call rhs_field%free()
    call corr_field%free()
    if (this%config%nbodies > 1) then
       call this%phi_total%free()
    end if

    ! Restore h1/h2 to what they were before
    coef%h1(:,:,:,:) = h1_restore(:,:,:,:)
    coef%h2(:,:,:,:) = h2_restore(:,:,:,:)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(coef%h1, coef%h1_d, n, HOST_TO_DEVICE, .false.)
       call device_memcpy(coef%h2, coef%h2_d, n, HOST_TO_DEVICE, .true.)
    end if

    if (allocated(h1_restore)) deallocate(h1_restore)
    if (allocated(h2_restore)) deallocate(h2_restore)
    if (allocated(Ax)) deallocate(Ax)
    if (allocated(ksp)) then
       call ksp%free()
       deallocate(ksp)
    end if
    if (allocated(pc)) then
       call precon_destroy(pc)
       deallocate(pc)
    end if

  end subroutine solve_base_mesh_displacement

  !> Updates the mesh velocity field based on current time and kinematics
  !> Sums contributions from all bodies: mesh_vel = Sum( V_i * Phi_i )
  subroutine update_mesh_velocity(this, coef, time_s)
    class(ale_manager_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(time_state_t), intent(in) :: time_s
    integer :: i, n
    type(body_kinematics_t) :: current_kin
    real(kind=rp) :: rot_mat(3,3)
    real(kind=rp) :: initial_rot_center(3)

    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return
    call profiler_start_region('ALE add mesh velocity')

    call field_rzero(this%wm_x)
    call field_rzero(this%wm_y)
    call field_rzero(this%wm_z)

    do i = 1, this%config%nbodies
       ! Compute kinematics for built-in motions
       ! "current_kin" will be like solid body kinematics at current time
       call compute_body_kinematics_built_in(current_kin, &
            this%config%bodies(i), time_s)

       ! User modifier (Superposition or Override)
       if (.not. associated(this%user_ale_rigid_kinematics, &
            dummy_user_ale_rigid_kinematics)) then
          call this%user_ale_rigid_kinematics(this%config%bodies(i)%id, &
               time_s, &
               current_kin%vel_trans, &
               current_kin%vel_ang)
       end if

       current_kin%center = this%ale_pivot(i)%pos
       this%ale_pivot(i)%vel = current_kin%vel_trans

       this%body_kin(i)%center = this%ale_pivot(i)%pos
       this%body_kin(i)%vel_trans = current_kin%vel_trans
       this%body_kin(i)%vel_ang = current_kin%vel_ang

       ! Compute rotation matrix at current time
       call this%compute_rotation_matrix(i, time_s)
       rot_mat = this%body_rot_matrices(:,:,i)
       initial_rot_center = this%config%bodies(i)%rot_center

       ! Accumulate contribution from each body and add to mesh velocity
       call add_kinematics_to_mesh_velocity(this%wm_x, this%wm_y, &
            this%wm_z, this%x_ref, this%y_ref, this%z_ref , &
            this%base_shapes(i), coef, current_kin, rot_mat, &
            initial_rot_center)

       ! For checkpointing
       call this%prep_checkpoint(i)
    end do

    ! If user has provided a custom function for mesh velocity.
    ! User mesh velocity will be added to the ale computed mesh velocity.
    ! This routine should not be used for rigid body motions!
    if (.not. associated(this%user_ale_mesh_vel, &
         dummy_user_ale_mesh_velocity)) then
       call this%user_ale_mesh_vel(this%wm_x, this%wm_y, this%wm_z, &
            coef, this%x_ref, this%y_ref, this%z_ref, this%base_shapes, time_s)
    end if

    call profiler_end_region('ALE add mesh velocity')

  end subroutine update_mesh_velocity

  !> Main routine to advance the mesh in time
  subroutine advance_mesh(this, coef, time, nadv)
    class(ale_manager_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    integer :: i

    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return
    call profiler_start_region('ALE update mesh')
    do i = 1, this%config%nbodies
       ! Advance Point Trackers attached to this body.
       ! Can be used for torque calculation (simcomp) at a point distanced
       ! from the body.
       ! or other purposes like tracking movement (user_check).
       call this%ghost_tracker_coord_step(this%body_kin(i), time, nadv, i)
       ! Update Pivot Location if requested
       call update_pivot_location(this%ale_pivot(i), &
            this%ale_pivot(i)%pos, &
            this%ale_pivot(i)%vel, &
            time, &
            nadv, &
            this%config%bodies(i))
    end do

    ! Update lagged B terms (geometry history)
    call coef%update_B_history()

    ! Update mesh coordinates
    call update_ale_mesh(coef, this%wm_x, this%wm_y, this%wm_z, &
         this%wm_x_lag, this%wm_y_lag, this%wm_z_lag, &
         time, nadv, "ab")

    ! Update internal history of mesh velocity.
    call this%wm_x_lag%update()
    call this%wm_y_lag%update()
    call this%wm_z_lag%update()
    call profiler_end_region('ALE update mesh')
  end subroutine advance_mesh

  ! Compute mesh stiffness with per-body gain/decay from stiff_geom.
  subroutine compute_stiffness_ale(coef, params)
    type(coef_t), intent(inout) :: coef
    type(ale_config_t), intent(in) :: params
    integer :: i, n, b, ierr
    integer, allocatable :: cheap_map(:)
    integer :: n_cheap, map_idx
    real(kind=rp) :: x, y, z
    real(kind=rp) :: raw_dist, body_stiff_val, max_added_stiff
    real(kind=rp) :: cx, cy, cz
    real(kind=rp) :: arg, decay, gain, norm_dist
    real(kind=rp) :: sample_start_time, sample_end_time, sample_time
    type(field_t), allocatable :: dist_fields(:)
    character(len=128) :: log_buf

    n = coef%dof%size()

    ! Check how many bodies need cheap_dist and create map
    allocate(cheap_map(params%nbodies))
    cheap_map = 0
    n_cheap = 0

    do b = 1, params%nbodies
       if (trim(params%bodies(b)%stiff_geom%type) .eq. 'cheap_dist') then
          n_cheap = n_cheap + 1
          cheap_map(b) = n_cheap
       end if
    end do

    ! Allocate and Compute cheap_dist only for required bodies
    if (n_cheap > 0) then
       allocate(dist_fields(n_cheap))

       do b = 1, params%nbodies
          map_idx = cheap_map(b)
          if (map_idx .gt. 0) then

             call dist_fields(map_idx)%init(coef%dof, "tmp_cheap_dist")

             call neko_log%message(' ')
             call neko_log%message(" Start: cheap dist calculation " // &
                  "for body '" // trim(params%bodies(b)%name) // "'")

             call MPI_Barrier(NEKO_COMM, ierr)
             sample_start_time = MPI_WTIME()

             if (NEKO_BCKND_DEVICE .eq. 1) then
                call compute_cheap_dist_device(dist_fields(map_idx), coef, &
                     coef%msh, params%bodies(b)%zone_indices, &
                     copy_to_host = .true.)
             else
                call compute_cheap_dist_v2_cpu(dist_fields(map_idx), coef, &
                     coef%msh, params%bodies(b)%zone_indices)
             end if

             call MPI_Barrier(NEKO_COMM, ierr)
             sample_end_time = MPI_WTIME()
             sample_time = sample_end_time - sample_start_time

             write(log_buf, '(A, A, A, ES11.4, A)') "    cheap dist for '", &
                  trim(params%bodies(b)%name), "' took ", sample_time, " (s)"
             call neko_log%message(log_buf)
          end if
       end do
    end if
    call neko_log%message(' ')

    ! Build stiffness field on Host
    select case (trim(params%stiffness_type))
    case ('built-in')

       do concurrent (i = 1:n)
          x = coef%dof%x(i, 1, 1, 1)
          y = coef%dof%y(i, 1, 1, 1)
          z = coef%dof%z(i, 1, 1, 1)

          max_added_stiff = 0.0_rp

          ! Loop over bodies, calculate local contribution
          do b = 1, params%nbodies
             gain = params%bodies(b)%stiff_geom%gain
             if (trim(params%bodies(b)%stiff_geom%type) .eq. 'cheap_dist') then
                decay = params%bodies(b)%stiff_geom%stiff_dist
             else
                decay = params%bodies(b)%stiff_geom%radius
             end if

             ! Geometry Center
             cx = params%bodies(b)%stiff_geom%center(1)
             cy = params%bodies(b)%stiff_geom%center(2)
             cz = params%bodies(b)%stiff_geom%center(3)

             raw_dist = huge(0.0_rp)

             ! Calculate Distance
             select case (trim(params%bodies(b)%stiff_geom%type))
             case ('sphere')
                raw_dist = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)

             case ('cylinder')
                ! Distance to Z-axis centered at (cx, cy)
                raw_dist = sqrt((x - cx)**2 + (y - cy)**2)

             case ('box')
                ! ToDO

             case ('cheap_dist')
                map_idx = cheap_map(b)
                if (map_idx .gt. 0) then
                   raw_dist = dist_fields(map_idx)%x(i, 1, 1, 1)
                end if
             end select

             ! Apply Profile
             body_stiff_val = 0.0_rp
             select case (trim(params%bodies(b)%stiff_geom%decay_profile))
             case ('gaussian')
                ! exp( -(r/decay)^2 )
                arg = -(raw_dist**2) / (decay**2)
                arg = arg * params%bodies(b)%stiff_geom%cutoff_coef
                body_stiff_val = gain * exp(arg)

             case ('tanh')
                ! Tanh profile
                norm_dist = (raw_dist / decay)
                norm_dist = norm_dist * params%bodies(b)%stiff_geom%cutoff_coef
                body_stiff_val = gain * (1.0_rp - tanh(norm_dist))
             end select

             if (body_stiff_val .gt. max_added_stiff) then
                max_added_stiff = body_stiff_val
             end if
          end do

          coef%h1(i, 1, 1, 1) = 1.0_rp + max_added_stiff
          coef%h2(i, 1, 1, 1) = 0.0_rp
       end do

    case default
       call neko_error("ALE Manager: Unknown stiffness type")
    end select

    coef%ifh2 = .false.

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(coef%h1, coef%h1_d, n, HOST_TO_DEVICE, .false.)
       call device_memcpy(coef%h2, coef%h2_d, n, HOST_TO_DEVICE, .true.)
    end if

    if (allocated(dist_fields)) then
       do i = 1, size(dist_fields)
          call dist_fields(i)%free()
       end do
       deallocate(dist_fields)
    end if
    if (allocated(cheap_map)) deallocate(cheap_map)

  end subroutine compute_stiffness_ale

  ! Adds kinematics to mesh velocity.
  subroutine add_kinematics_to_mesh_velocity(wx, wy, wz, &
       x_ref, y_ref, z_ref, phi, coef, kinematics, rot_mat, initial_pivot_loc)
    type(field_t), intent(inout) :: wx, wy, wz
    type(field_t), intent(in) :: x_ref, y_ref, z_ref
    type(field_t), intent(in) :: phi
    type(coef_t), intent(in) :: coef
    type(body_kinematics_t), intent(in) :: kinematics
    real(kind=rp), intent(in) :: initial_pivot_loc(3)
    real(kind=rp), intent(in) :: rot_mat(3,3)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call add_kinematics_to_mesh_velocity_device(wx, wy, wz, &
            x_ref, y_ref, z_ref, &
            phi, coef, kinematics, rot_mat, initial_pivot_loc)
    else
       call add_kinematics_to_mesh_velocity_cpu(wx, wy, wz, &
            x_ref, y_ref, z_ref, &
            phi, coef, kinematics, rot_mat, initial_pivot_loc)
    end if
  end subroutine add_kinematics_to_mesh_velocity

  ! Updates mesh position by integrating mesh velocity in time using AB scheme.
  subroutine update_ale_mesh(c_Xh, wm_x, wm_y, wm_z, wm_x_lag, wm_y_lag, &
       wm_z_lag, time, nadv, scheme_)
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: wm_x, wm_y, wm_z
    type(field_series_t), intent(in) :: wm_x_lag, wm_y_lag, wm_z_lag
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: nadv
    character(len=*), intent(in) :: scheme_
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call update_ale_mesh_device(c_Xh, wm_x, wm_y, wm_z, &
            wm_x_lag, wm_y_lag, wm_z_lag, time, nadv, scheme_)
    else
       call update_ale_mesh_cpu(c_Xh, wm_x, wm_y, wm_z, &
            wm_x_lag, wm_y_lag, wm_z_lag, time, nadv, scheme_)
    end if
  end subroutine update_ale_mesh


  subroutine ale_manager_free(this)
    class(ale_manager_t), intent(inout), target :: this
    integer :: i

    if (.not. this%active) return

    call this%bc_moving%free()
    call this%bc_fixed%free()
    call this%bc_list%free()

    if (allocated(this%base_shapes)) then
       do i = 1, size(this%base_shapes)
          call this%base_shapes(i)%free()
       end do
       deallocate(this%base_shapes)
    end if

    call this%wm_x_lag%free()
    call this%wm_y_lag%free()
    call this%wm_z_lag%free()
    call this%x_ref%free()
    call this%y_ref%free()
    call this%z_ref%free()

    if (allocated(this%ale_pivot)) deallocate(this%ale_pivot)
    if (allocated(this%config%bodies)) deallocate(this%config%bodies)
    if (allocated(this%body_kin)) deallocate(this%body_kin)
    if (associated(this%global_pivot_pos)) deallocate(this%global_pivot_pos)
    if (associated(this%global_pivot_vel_lag)) &
         deallocate(this%global_pivot_vel_lag)
    if (associated(this%global_basis_pos)) deallocate(this%global_basis_pos)
    if (associated(this%global_basis_vel_lag)) &
         deallocate(this%global_basis_vel_lag)
    if (allocated(this%ghost_handles)) deallocate(this%ghost_handles)
    if (allocated(this%body_rot_matrices)) deallocate(this%body_rot_matrices)
    if (allocated(this%trackers)) deallocate(this%trackers)
    if (associated(neko_ale, this)) nullify(neko_ale)

  end subroutine ale_manager_free

  !> Factory for ALE Preconditioner
  subroutine ale_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype, params)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*), intent(in) :: pctype
    type(json_file), intent(inout) :: params
    call precon_allocator(pc, pctype)
    select type (pcp => pc)
    type is (jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (hsmg_t)
       call pcp%init(coef, bclst, params)
    type is (phmg_t)
       call pcp%init(coef, bclst, params)
    end select
    call ksp%set_pc(pc)
  end subroutine ale_precon_factory

  ! Sets the pivot state at restart.
  subroutine set_pivot_restart(this, time_restart)
    class(ale_manager_t), intent(inout) :: this
    real(kind=dp), intent(in) :: time_restart
    type(body_kinematics_t) :: kin_restart
    integer :: i, idx, handle_1, handle_2, offset_base
    type(time_state_t) :: time_state_dummy
    time_state_dummy%t = time_restart

    !if (.not. allocated(this%global_pivot_pos)) return

    do i = 1, this%config%nbodies

       call compute_body_kinematics_built_in(kin_restart, &
            this%config%bodies(i), time_state_dummy)

       ! User Modifier (Superposition or Override)
       if (.not. associated(this%user_ale_rigid_kinematics, &
            dummy_user_ale_rigid_kinematics)) then
          call this%user_ale_rigid_kinematics(this%config%bodies(i)%id, &
               time_state_dummy, &
               kin_restart%vel_trans, &
               kin_restart%vel_ang)
       end if

       this%ale_pivot(i)%vel = kin_restart%vel_trans

       idx = (i - 1) * 3
       ! Restore Position
       this%ale_pivot(i)%pos(1:3) = this%global_pivot_pos(idx + 1:idx + 3)
       this%body_kin(i)%center = this%ale_pivot(i)%pos
       this%body_kin(i)%vel_trans = kin_restart%vel_trans
       this%body_kin(i)%vel_ang = kin_restart%vel_ang

       ! Restore Velocity History
       this%ale_pivot(i)%vel_lag(1:3, 1:3) = &
            this%global_pivot_vel_lag(idx + 1:idx + 3, :)


       offset_base = (i-1)*6
       handle_1 = this%ghost_handles(1, i)
       handle_2 = this%ghost_handles(2, i)

       if ((handle_1 .gt. 0) .and. (handle_1 .le. this%n_trackers)) then
          this%trackers(handle_1)%pos = &
               this%global_basis_pos(offset_base + 1 : offset_base + 3)

          ! Restore velocity history for ghost-x
          this%trackers(handle_1)%vel_lag = &
               this%global_basis_vel_lag(offset_base + 1 : offset_base + 3, :)
       end if

       if ((handle_2 .gt. 0) .and. (handle_2 .le. this%n_trackers)) then
          this%trackers(handle_2)%pos = &
               this%global_basis_pos(offset_base + 4 : offset_base + 6)

          ! Restore velocity history for ghost-y
          this%trackers(handle_2)%vel_lag = &
               this%global_basis_vel_lag(offset_base + 4 : offset_base + 6, :)
       end if
    end do
  end subroutine set_pivot_restart

  ! Restores the current coef and related metrics
  ! and the pivot states at restart.
  subroutine sync_chkp(this, coef, Xh, adv, chkp, gs_Xh)
    class(ale_manager_t), intent(inout) :: this
    class(advection_t), intent(inout) :: adv
    type(coef_t), intent(inout) :: coef
    type(space_t), intent(inout) :: Xh
    type(chkp_t), intent(in) :: chkp
    type(gs_t), intent(inout) :: gs_Xh
    integer :: i, j, n

    ! Return if ALE is not active.
    if (.not. this%active) return

    if (allocated(chkp%previous_mesh%elements)) then
       call neko_error("ALE restart failed: " // &
            "The current mesh has a different number " // &
            "of elements than the checkpoint.")
    end if

    ! Restarting from a different polynomial order
    if (chkp%previous_Xh%lx .ne. Xh%lx) then
       n = coef%dof%size()
       associate(wm_x => this%wm_x, wm_y => this%wm_y, wm_z => this%wm_z)
         do concurrent (j = 1:n)
            ! Mesh Velocity
            wm_x%x(j,1,1,1) = wm_x%x(j,1,1,1) * coef%mult(j,1,1,1)
            wm_y%x(j,1,1,1) = wm_y%x(j,1,1,1) * coef%mult(j,1,1,1)
            wm_z%x(j,1,1,1) = wm_z%x(j,1,1,1) * coef%mult(j,1,1,1)
         end do
       end associate

       do i = 1, this%wm_x_lag%size()
          do concurrent (j = 1:n)
             this%wm_x_lag%lf(i)%x(j,1,1,1) = &
                  this%wm_x_lag%lf(i)%x(j,1,1,1) * coef%mult(j,1,1,1)
             this%wm_y_lag%lf(i)%x(j,1,1,1) = &
                  this%wm_y_lag%lf(i)%x(j,1,1,1) * coef%mult(j,1,1,1)
             this%wm_z_lag%lf(i)%x(j,1,1,1) = &
                  this%wm_z_lag%lf(i)%x(j,1,1,1) * coef%mult(j,1,1,1)
          end do
       end do
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%wm_x%copy_from(HOST_TO_DEVICE, sync = .false.)
       call this%wm_y%copy_from(HOST_TO_DEVICE, sync = .false.)
       call this%wm_z%copy_from(HOST_TO_DEVICE, sync = .false.)

       call this%wm_x_lag%lf(1)%copy_from(HOST_TO_DEVICE, &
            sync = .false.)
       call this%wm_x_lag%lf(2)%copy_from(HOST_TO_DEVICE, &
            sync = .false.)

       call this%wm_y_lag%lf(1)%copy_from(HOST_TO_DEVICE, &
            sync = .false.)
       call this%wm_y_lag%lf(2)%copy_from(HOST_TO_DEVICE, &
            sync = .false.)

       call this%wm_z_lag%lf(1)%copy_from(HOST_TO_DEVICE, &
            sync = .false.)
       call this%wm_z_lag%lf(2)%copy_from(HOST_TO_DEVICE, &
            sync = .false.)

       if (c_associated(coef%dof%x_d)) then
          call device_memcpy(coef%dof%x, coef%dof%x_d, &
               size(coef%dof%x), HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(coef%dof%y, coef%dof%y_d, &
               size(coef%dof%y), HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(coef%dof%z, coef%dof%z_d, &
               size(coef%dof%z), HOST_TO_DEVICE, sync = .false.)
       end if

       if (c_associated(coef%Blag_d)) then
          call device_memcpy(coef%Blag, coef%Blag_d, size(coef%Blag), &
               HOST_TO_DEVICE, sync = .false.)
       end if

       if (c_associated(coef%Blaglag_d)) then
          call device_memcpy(coef%Blaglag, coef%Blaglag_d, &
               size(coef%Blaglag), HOST_TO_DEVICE, sync = .false.)
       end if
       call device_sync()
    end if

    ! Restarting from a different polynomial order
    if (chkp%previous_Xh%lx .ne. Xh%lx) then
       call rotate_cyc(this%wm_x%x, this%wm_y%x, this%wm_z%x, 1, coef)
       call gs_Xh%op(this%wm_x, GS_OP_ADD)
       call gs_Xh%op(this%wm_y, GS_OP_ADD)
       call gs_Xh%op(this%wm_z, GS_OP_ADD)
       call rotate_cyc(this%wm_x%x, this%wm_y%x, this%wm_z%x, 0, coef)

       do i = 1, this%wm_x_lag%size()
          call rotate_cyc(this%wm_x_lag%lf(i)%x, this%wm_y_lag%lf(i)%x, &
               this%wm_z_lag%lf(i)%x, 1, coef)
          call gs_Xh%op(this%wm_x_lag%lf(i), GS_OP_ADD)
          call gs_Xh%op(this%wm_y_lag%lf(i), GS_OP_ADD)
          call gs_Xh%op(this%wm_z_lag%lf(i), GS_OP_ADD)
          call rotate_cyc(this%wm_x_lag%lf(i)%x, this%wm_y_lag%lf(i)%x, &
               this%wm_z_lag%lf(i)%x, 0, coef)
       end do
    end if


    call this%set_pivot_restart(chkp%t)
    call coef%recompute_metrics()

    ! If polynomial order changes during restart, we use current's mesh
    ! mass matrix for Blag and Blaglag. This will introduce some error,
    ! but maybe better than
    ! not restarting at all. Otherwise we need to save lagged mesh coordinates
    ! as well in order to be more accurate.
    if (chkp%previous_Xh%lx .ne. Xh%lx) then
       coef%Blag = coef%B
       coef%Blaglag = coef%B
       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (c_associated(coef%Blag_d)) then
             call device_memcpy(coef%Blag, coef%Blag_d, n, &
                  HOST_TO_DEVICE, sync = .false.)
          end if
          if (c_associated(coef%Blaglag_d)) then
             call device_memcpy(coef%Blaglag, coef%Blaglag_d, n, &
                  HOST_TO_DEVICE, sync = .false.)
          end if
          call device_sync()
       end if
    end if

    call adv%recompute_metrics(coef, .true.)
  end subroutine sync_chkp

  subroutine set_pivot_basis_for_checkpoint(this, body_idx)
    class(ale_manager_t), intent(inout) :: this
    integer, intent(in) :: body_idx
    integer :: idx, offset_base, h1, h2

    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return

    idx = (body_idx - 1) * 3
    this%global_pivot_pos(idx + 1:idx + 3) = this%ale_pivot(body_idx)%pos(1:3)
    this%global_pivot_vel_lag(idx + 1:idx + 3, :) = &
         this%ale_pivot(body_idx)%vel_lag(1:3, 1:3)

    h1 = this%ghost_handles(1, body_idx)
    h2 = this%ghost_handles(2, body_idx)

    offset_base = (body_idx-1)*6

    ! Save Positions
    this%global_basis_pos(offset_base + 1 : offset_base + 3) = &
         this%get_tracker_pos(h1)
    this%global_basis_pos(offset_base + 4 : offset_base + 6) = &
         this%get_tracker_pos(h2)

    ! Ghost-x history
    this%global_basis_vel_lag(offset_base + 1 : offset_base + 3, :) = &
         this%trackers(h1)%vel_lag

    ! Ghost-y history
    this%global_basis_vel_lag(offset_base + 4 : offset_base + 6, :) = &
         this%trackers(h2)%vel_lag
  end subroutine set_pivot_basis_for_checkpoint

  ! Append val to arr if not already present.
  subroutine append_unique_int(arr, n, val)
    integer, allocatable, intent(inout) :: arr(:)
    integer, intent(inout) :: n
    integer, intent(in) :: val
    integer, allocatable :: tmp(:)
    integer :: k

    do k = 1, n
       if (arr(k) .eq. val) return
    end do

    allocate(tmp(n + 1))
    if (n .gt. 0) tmp(1:n) = arr(1:n)
    tmp(n + 1) = val

    if (allocated(arr)) deallocate(arr)
    call move_alloc(tmp, arr)

    n = n + 1
  end subroutine append_unique_int

  !> Performs a preview of the mesh motion to verify quality/topology
  subroutine mesh_preview(this, coef, json)
    class(ale_manager_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: json
    type(fld_file_output_t) :: fout
    type(field_t) :: dummy_field
    type(time_state_t) :: t_state
    type(file_t) :: out_file
    real(kind=rp) :: t_start
    real(kind=rp) :: t_end
    real(kind=rp) :: dt
    real(kind=rp) :: min_jac
    integer :: output_freq
    integer :: step, n_steps
    integer :: nadv, nadv_sim
    integer :: n
    logical :: mesh_preview_active
    character(len=128) :: log_buf

    mesh_preview_active = .false.

    if (json%valid_path('case.fluid.ale.mesh_preview.enabled')) then
       call json%get('case.fluid.ale.mesh_preview.enabled', &
            mesh_preview_active)
    end if

    if (.not. mesh_preview_active) return

    call json_get_or_default(json, 'case.fluid.ale.mesh_preview.start_time', &
         t_start, 0.0_rp)
    call json_get(json, 'case.fluid.ale.mesh_preview.end_time', &
         t_end)
    call json_get(json, 'case.fluid.ale.mesh_preview.dt', &
         dt)
    call json_get(json, &
         'case.fluid.ale.mesh_preview.output_freq', &
         output_freq)

    call neko_log%section("ALE Mesh Preview")
    call neko_log%message("Executing mesh motion preview...")

    n_steps = int((t_end - t_start) / dt)
    call json_get(json, 'case.numerics.time_order', nadv_sim)

    write(log_buf, '(A, ES23.15)') '  Start Time : ', t_start
    call neko_log%message(log_buf)
    write(log_buf, '(A, ES23.15)') '  End Time   : ', t_end
    call neko_log%message(log_buf)
    write(log_buf, '(A, ES23.15)') '  dt         : ', dt
    call neko_log%message(log_buf)
    write(log_buf, '(A, I0)') '  Num Steps  :   ', n_steps
    call neko_log%message(log_buf)
    write(log_buf, '(A, I0)') '  Output Freq:   ', output_freq
    call neko_log%message(log_buf)
    call neko_log%message('')

    ! Setup dummy field for output
    call dummy_field%init(coef%dof, "mesh_preview")
    call field_rzero(dummy_field)

    call fout%init(rp, "mesh_preview", 1)
    call fout%fields%assign_to_field(1, dummy_field)
    select type (ft => fout%file_%file_type)
    type is (fld_file_t)
       ft%write_mesh = .true.
       ft%skip_pressure = .false.
    end select


    ! Time Loop Setup
    step = 0
    nadv = 1
    t_state%t = t_start
    t_state%dt = dt
    t_state%tstep = 0
    t_state%dtlag = dt
    n = coef%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       min_jac = device_glmin(coef%jac_d, n)
    else
       min_jac = glmin(coef%jac, n)
    end if

    call sync_mesh_preview_step(coef, dummy_field)
    call fout%sample(t_state%t)

    write(log_buf, '(A,I0, A,ES23.15, A,ES18.11)') &
         "Initial Mesh and Mass matrix saved!  Step: ", step, " | Time:", &
         t_state%t, " | Min Jac: ", min_jac

    call neko_log%message(trim(log_buf))
    call this%update_mesh_velocity(coef, t_state)

    do step = 1, n_steps
       t_state%tstep = step
       t_state%t = t_start + (step * dt)
       nadv = min(step, nadv_sim)

       call this%advance_mesh(coef, t_state, nadv)
       call coef%recompute_metrics()


       if (NEKO_BCKND_DEVICE .eq. 1) then
          min_jac = device_glmin(coef%jac_d, n)
       else
          min_jac = glmin(coef%jac, n)
       end if

       if (min_jac .le. 0.0_rp) then
          write(log_buf, '(A, ES18.11, A, ES23.15)') &
               "Negative Jacobian detected (", min_jac, ") at t = ", &
               t_state%t
          call neko_log%message(log_buf)

          call sync_mesh_preview_step(coef, dummy_field)
          call fout%sample(t_state%t)

          write(log_buf, '(A,I0, A,ES23.15, A,ES18.11)') &
               "Mesh and Mass matrix saved!  Step: ", step, " | Time:", &
               t_state%t, " | Min Jac:", min_jac
          call neko_log%message(trim(log_buf))

          call neko_error("ALE Mesh Preview Aborted: Negative Jacobian found.")
       end if

       if (mod(step, output_freq) .eq. 0) then

          call sync_mesh_preview_step(coef, dummy_field)
          call fout%sample(t_state%t)

          write(log_buf, '(A,I0, A,ES23.15, A,ES18.11)') &
               "Mesh and Mass matrix saved!  Step: ", step, " | Time:", &
               t_state%t, " | Min Jac:", min_jac
          call neko_log%message(trim(log_buf))

       end if

       call this%update_mesh_velocity(coef, t_state)

    end do

    call dummy_field%free()
    call fout%free()
    call neko_log%end_section()
    call neko_log%message("Mesh preview complete.")
    call neko_error("ALE Mesh Preview Finished Successfully.")

  end subroutine mesh_preview

  subroutine sync_mesh_preview_step(coef, dummy_field)
    type(fld_file_output_t) :: fout
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: dummy_field
    integer :: n

    n = coef%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(dummy_field%x_d, coef%B_d, n)
    else
       call copy(dummy_field%x, coef%B, n)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(mesh => coef%dof)
         call device_memcpy(mesh%x, mesh%x_d, mesh%size(), &
              DEVICE_TO_HOST, sync = .false.)
         call device_memcpy(mesh%y, mesh%y_d, mesh%size(), &
              DEVICE_TO_HOST, sync = .false.)
         call device_memcpy(mesh%z, mesh%z_d, mesh%size(), &
              DEVICE_TO_HOST, sync = .false.)
       end associate
    end if

  end subroutine sync_mesh_preview_step

  ! Asign a tracker point to a body. The tracker moves with body's
  ! rigid motion.
  function request_tracker(this, initial_pos, body_id) result(handle)
    class(ale_manager_t), intent(inout) :: this
    real(kind=rp), intent(in) :: initial_pos(3)
    integer, intent(in) :: body_id
    integer :: handle
    type(point_tracker_t), allocatable :: tmp(:)

    handle = -100
    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return

    if (.not. allocated(this%trackers)) then
       allocate(this%trackers(30))
       this%n_trackers = 0
    elseif (this%n_trackers .ge. size(this%trackers)) then
       allocate(tmp(size(this%trackers) + 30))
       tmp(1:size(this%trackers)) = this%trackers
       deallocate(this%trackers)
       call move_alloc(tmp, this%trackers)
    end if
    this%n_trackers = this%n_trackers + 1
    handle = this%n_trackers

    this%trackers(handle)%pos = initial_pos
    this%trackers(handle)%body_id = body_id
    this%trackers(handle)%vel_lag = this%ale_pivot(body_id)%vel_lag
  end function request_tracker

  function get_tracker_pos(this, handle) result(pos)
    class(ale_manager_t), intent(in) :: this
    integer, intent(in) :: handle
    real(kind=rp) :: pos(3)

    if (handle .gt. 0 .and. handle .le. this%n_trackers) then
       pos = this%trackers(handle)%pos
    else
       pos = 0.0_rp
    end if
  end function get_tracker_pos


  !> Computes Rotation Matrix
  subroutine compute_rotation_matrix(this, body_idx, time)
    class(ale_manager_t), intent(inout) :: this
    integer, intent(in) :: body_idx
    type(time_state_t), intent(in) :: time
    integer :: h_x, h_y
    real(kind=rp) :: P(3), Gx(3), Gy(3)
    real(kind=rp) :: u(3), v(3), w(3), v_temp(3)

    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return

    ! Get Points
    h_x = this%ghost_handles(1, body_idx)
    h_y = this%ghost_handles(2, body_idx)

    P = this%ale_pivot(body_idx)%pos
    Gx = this%get_tracker_pos(h_x)
    Gy = this%get_tracker_pos(h_y)

    ! Construct u (New X-axis)
    u = Gx - P
    u = u / sqrt(sum(u**2))

    ! Construct w via cross product (Z-axis)
    v_temp = Gy - P
    w(1) = u(2)*v_temp(3) - u(3)*v_temp(2)
    w(2) = u(3)*v_temp(1) - u(1)*v_temp(3)
    w(3) = u(1)*v_temp(2) - u(2)*v_temp(1)
    w = w / sqrt(sum(w**2))

    ! Construct v via orthogonalization (Y-axis)
    v(1) = w(2)*u(3) - w(3)*u(2)
    v(2) = w(3)*u(1) - w(1)*u(3)
    v(3) = w(1)*u(2) - w(2)*u(1)

    this%body_rot_matrices(:, 1, body_idx) = u
    this%body_rot_matrices(:, 2, body_idx) = v
    this%body_rot_matrices(:, 3, body_idx) = w

  end subroutine compute_rotation_matrix


  !> Logs rotation angles for all or selected bodies.
  !> can be called in user%compute.
  !> eg: call neko_ale%log_rot_angles(time, body_idxs)
  subroutine log_rot_angles(this, time, body_idxs)
    class(ale_manager_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    integer, optional, intent(in) :: body_idxs(:)

    integer :: i, idx, n_log
    real(kind=rp) :: roll_deg, pitch_deg, yaw_deg
    real(kind=rp) :: R(3,3)
    character(len=256) :: log_buf
    real(kind=rp), parameter :: rad_to_deg = 180.0_rp / pi

    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return

    if (present(body_idxs)) then
       n_log = size(body_idxs)
    else
       n_log = this%config%nbodies
    end if

    call neko_log%message(" ")
    call neko_log%message("---------Rotation log---------")
    call neko_log%message("variable, time step, time, body, " // &
         "x_val, y_val, z_val")

    ! If body_idxs is provided, only log those. Otherwise, log all.
    do i = 1, n_log

       if (present(body_idxs)) then
          idx = body_idxs(i)
       else
          idx = i
       end if

       R = this%body_rot_matrices(:, :, idx)

       ! Angles
       yaw_deg = atan2(R(2,1), R(1,1)) * rad_to_deg
       pitch_deg = atan2(-R(3,1), sqrt(R(3,2)**2 + R(3,3)**2)) * rad_to_deg
       roll_deg = atan2(R(3,2), R(3,3)) * rad_to_deg

       ! Log Rotation Angles (Roll, Pitch, Yaw) -> (X, Y, Z)
       write(log_buf, '(A, I0, A, ES13.6, A, A, A, 3(ES17.10, :, 2X))') &
            "Total_Rot_deg    ", time%tstep, "  ", time%t, "  ", &
            trim(this%config%bodies(idx)%name), "  ", &
            roll_deg, pitch_deg, yaw_deg
       call neko_log%message(trim(log_buf))

    end do

  end subroutine log_rot_angles

  !> Logs pivot positions for all or selected bodies.
  !> can be called in user%compute.
  !> eg: call neko_ale%log_pivot(time, body_idxs)
  subroutine log_pivot(this, time, body_idxs)
    class(ale_manager_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    integer, optional, intent(in) :: body_idxs(:)
    integer :: i, idx, n_log
    real(kind=rp) :: pivot_pos(3), pivot_vel(3)
    character(len=256) :: log_buf

    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return

    if (present(body_idxs)) then
       n_log = size(body_idxs)
    else
       n_log = this%config%nbodies
    end if

    call neko_log%message(" ")
    call neko_log%message("----------Pivot Log-----------")
    call neko_log%message("variable, time step, time, body, " // &
         "x_val, y_val, z_val")

    ! If body_idxs is provided, only log those. Otherwise, log all.
    do i = 1, n_log

       if (present(body_idxs)) then
          idx = body_idxs(i)
       else
          idx = i
       end if

       pivot_pos = this%ale_pivot(idx)%pos
       pivot_vel = this%ale_pivot(idx)%vel

       ! Pivot Position
       write(log_buf, '(A, I0, A, ES13.6, A, A, A, 3(ES17.10, :, 2X))') &
            "Total_Pivot_pos  ", time%tstep, "  ", time%t, "  ", &
            trim(this%config%bodies(idx)%name), "  ", &
            this%ale_pivot(idx)%pos
       call neko_log%message(trim(log_buf))

       ! Pivot Velocity
       write(log_buf, '(A, I0, A, ES13.6, A, A, A, 3(ES17.10, :, 2X))') &
            "Total_Pivot_vel  ", time%tstep, "  ", time%t, "  ", &
            trim(this%config%bodies(idx)%name), "  ", &
            this%ale_pivot(idx)%vel
       call neko_log%message(trim(log_buf))
    end do

  end subroutine log_pivot

  subroutine ghost_tracker_coord_step(this, kin_object, time_s, nadv, body_idx)
    class(ale_manager_t), intent(inout) :: this
    type(body_kinematics_t), intent(in) :: kin_object
    type(time_state_t), intent(in) :: time_s
    integer, intent(in) :: nadv
    integer, intent(in) :: body_idx
    integer :: t
    real(kind=rp) :: p_vel(3), rel_pos(3), v_tan(3)

    if (.not. this%active) return
    if (.not. this%has_moving_boundary) return

    if (allocated(this%trackers)) then
       do t = 1, this%n_trackers
          if (this%trackers(t)%body_id .eq. &
               this%config%bodies(body_idx)%id) then
             if (t .eq. this%ghost_handles(1, body_idx) .or. &
                  t .eq. this%ghost_handles(2, body_idx)) then

                ! Calculate the Arm vector (r) at current step
                rel_pos = this%trackers(t)%pos - kin_object%center

                ! Calculate tangential velocity (Omega \cross r)
                v_tan(1) = kin_object%vel_ang(2) * rel_pos(3) - &
                     kin_object%vel_ang(3) * rel_pos(2)
                v_tan(2) = kin_object%vel_ang(3) * rel_pos(1) - &
                     kin_object%vel_ang(1) * rel_pos(3)
                v_tan(3) = kin_object%vel_ang(1) * rel_pos(2) - &
                     kin_object%vel_ang(2) * rel_pos(1)

                ! Total velocity
                p_vel = kin_object%vel_trans + v_tan

                if (time_s%tstep .gt. 0) then
                   call ab_integrate_point_pos(this%trackers(t)%pos, &
                        this%trackers(t)%vel_lag, p_vel, time_s, nadv)
                end if

             end if

          end if
       end do

    end if
  end subroutine ghost_tracker_coord_step

  subroutine get_ale_solver_params_json(this, json, ksp_solver, precon_type, &
       precon_params, abstol, ksp_max_iter, res_monitor, import_base_shapes)
    class(ale_manager_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable, intent(inout) :: ksp_solver
    character(len=:), allocatable, intent(inout) :: precon_type
    type(json_file), intent(inout) :: precon_params
    real(kind=rp), intent(out) :: abstol
    integer, intent(out) :: ksp_max_iter
    logical, intent(out) :: res_monitor
    logical, intent(out) :: import_base_shapes
    logical :: tmp_logical
    character(len=:), allocatable :: tmp_str

    if (allocated(ksp_solver)) deallocate(ksp_solver)
    if (allocated(precon_type)) deallocate(precon_type)

    call json_get_or_default(json, &
         'case.fluid.ale.solver.import_base_shape', &
         import_base_shapes, .false.)

    call json_get_or_default(json, 'case.fluid.ale.solver.type', &
         ksp_solver, 'cg')

    call json_get_or_default(json, &
         'case.fluid.ale.solver.preconditioner.type', precon_type, 'jacobi')

    if (json%valid_path('case.fluid.ale.solver.preconditioner')) then
       call json_get(json, 'case.fluid.ale.solver.preconditioner', &
            precon_params)
    end if

    call json_get_or_default(json, &
         'case.fluid.ale.solver.absolute_tolerance', abstol, 1.0e-10_rp)
    call json_get_or_default(json, 'case.fluid.ale.solver.monitor', &
         res_monitor, .false.)
    call json_get_or_default(json, 'case.fluid.ale.solver.max_iterations', &
         ksp_max_iter, 10000)

    if (json%valid_path('case.fluid.ale.solver.output_base_shape')) then
       call json%get('case.fluid.ale.solver.output_base_shape', tmp_logical)
       this%config%if_output_phi = tmp_logical
    end if
    if (json%valid_path('case.fluid.ale.solver.output_stiffness')) then
       call json%get('case.fluid.ale.solver.output_stiffness', tmp_logical)
       this%config%if_output_stiffness = tmp_logical
    end if

    ! Mesh Stiffness
    if (json%valid_path('case.fluid.ale.solver.mesh_stiffness.type')) then
       call json%get('case.fluid.ale.solver.mesh_stiffness.type', tmp_str)
       this%config%stiffness_type = tmp_str
       if (.not. (trim(tmp_str) .eq. 'built-in')) then
          call neko_error("ALE: stiffness_type must be 'built-in'")
       end if
    end if
  end subroutine get_ale_solver_params_json

  ! Register ALE fields for checkpointing.
  subroutine register_checkpoint_fields(this, coef, checkpoint)
    class(ale_manager_t), intent(inout), target :: this
    type(coef_t), intent(inout) :: coef
    type(chkp_t), intent(inout) :: checkpoint
    integer :: i

    if (.not. this%active) return

    ! Add checkpoint data for ALE.
    call checkpoint%add_ale(coef%dof%x, coef%dof%y, &
         coef%dof%z, coef%dof%x_d, coef%dof%y_d, &
         coef%dof%z_d, &
         coef%Blag, coef%Blaglag, coef%Blag_d, coef%Blaglag_d, &
         this%wm_x, this%wm_y, this%wm_z, &
         this%wm_x_lag, this%wm_y_lag, &
         this%wm_z_lag, &
         this%global_pivot_pos, &
         this%global_pivot_vel_lag, &
         this%global_basis_pos, &
         this%global_basis_vel_lag)

  end subroutine register_checkpoint_fields
end module ale_manager
