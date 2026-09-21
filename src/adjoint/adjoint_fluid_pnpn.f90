!> @file adjoint_fluid_pnpn.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Adjoint Pn/Pn formulation.
module adjoint_fluid_pnpn
  use, intrinsic :: iso_fortran_env, only: error_unit
  use coefs, only: coef_t
  use symmetry_aligned, only: symmetry_aligned_t
  use registry, only: neko_registry
  use logger, only: neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use num_types, only: rp
  use krylov, only: ksp_monitor_t
  use adjoint_pnpn_residual, only: adjoint_pnpn_prs_res_t, &
       adjoint_pnpn_vel_res_t, adjoint_pnpn_prs_res_factory, &
       adjoint_pnpn_vel_res_factory
  use rhs_maker, only: rhs_maker_sumab_t, rhs_maker_bdf_t, rhs_maker_ext_t, &
       rhs_maker_oifs_t, rhs_maker_sumab_fctry, rhs_maker_bdf_fctry, &
       rhs_maker_ext_fctry, rhs_maker_oifs_fctry
  use adjoint_fluid_scheme_incompressible, only: &
       adjoint_fluid_scheme_incompressible_t
  use device_mathops, only: device_opadd2cm
  use fluid_aux, only: fluid_step_info
  use scratch_registry, only: neko_scratch_registry
  use projection, only: projection_t
  use projection_vel, only: projection_vel_t
  use device, only: device_memcpy, HOST_TO_DEVICE, device_event_sync, &
       glb_cmd_event
  use advection_adjoint, only: advection_adjoint_t, advection_adjoint_factory
  use profiler, only: profiler_start_region, profiler_end_region
  use json_module, only: json_file, json_core, json_value
  use json_utils, only: json_get, json_get_or_default, json_extract_item
  use ax_product, only: ax_t, ax_helm_allocator
  use field, only: field_t
  use dirichlet, only: dirichlet_t
  use shear_stress, only: shear_stress_t
  use wall_model_bc, only: wall_model_bc_t
  use facet_normal, only: facet_normal_t
  use non_normal_aligned, only: non_normal_aligned_t
  use checkpoint, only: chkp_t
  use mesh, only: mesh_t
  use user_intf, only: user_t
  use time_step_controller, only: time_step_controller_t
  use gs_ops, only: GS_OP_ADD
  use neko_config, only: NEKO_BCKND_DEVICE
  use mathops, only: opadd2cm
  use scalar_bc_projector, only: scalar_bc_projector_t
  use vector_bc_projector, only: vector_bc_projector_t, &
       segregated_vector_bc_projector_t
  use zero_dirichlet, only: zero_dirichlet_t
  use utils, only: neko_error
  use field_math, only: field_add2, field_copy, &
       field_add2s2
  use bc, only: bc_t, BC_DIRICHLET
  use file, only: file_t
  use operators, only: ortho
  use opr_device, only: device_ortho
  use inflow, only: inflow_t
  use field_dirichlet, only: field_dirichlet_t
  use blasius, only: blasius_t
  use field_dirichlet_vector, only: field_dirichlet_vector_t
  use dong_outflow, only: dong_outflow_t
  use time_state, only: time_state_t
  use vector, only: vector_t
  use device_math, only: device_vlsc3, device_cmult, device_col2
  use math, only: vlsc3, cmult, col2
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only: mpi_sum, mpi_max, mpi_allreduce, MPI_INTEGER, &
       MPI_LOGICAL, MPI_LOR
  use operators, only : opgrad, curl, grad
  use normal_vec_bcs, only: normal_vec_bcs_t

  implicit none
  private

  type, public, extends(adjoint_fluid_scheme_incompressible_t) :: &
       adjoint_fluid_pnpn_t

     !> The right-hand sides in the linear solves.
     type(field_t) :: p_res, u_res, v_res, w_res

     !> The unknowns in the linear solves, i.e. the solution increments with
     !! respect to the previous time-step.
     type(field_t) :: dp, du, dv, dw

     !
     ! Implicit operators, i.e. the left-hand-side of the Helmholz problem.
     !

     ! Coupled Helmholz operator for velocity
     class(ax_t), allocatable :: Ax_vel
     ! Helmholz operator for pressure
     class(ax_t), allocatable :: Ax_prs

     !
     ! Projections for solver speed-up
     !

     !> Pressure projection
     type(projection_t) :: proj_prs
     type(projection_vel_t) :: proj_vel

     !
     ! Special Karniadakis scheme boundary conditions in the pressure equation
     !

     !> Surface term in pressure rhs. Masks all strong velocity bcs.
     type(facet_normal_t) :: bc_prs_surface

     !> Surface term in pressure rhs. Masks symmetry bcs.
     type(facet_normal_t) :: bc_sym_surface

     !> Surface term in pressure rhs. Masks symmetry bcs.
     type(normal_vec_bcs_t) :: bc_curl_curl

     !> Projector for velocity increment constraints.
     class(vector_bc_projector_t), allocatable :: bcs_vel_projector
     !> Projector for pressure increment constraints.
     type(scalar_bc_projector_t) :: bcs_prs_projector


     ! Checker for wether we have a strong pressure bc. If not, the pressure
     ! is demeaned at every time step.
     logical :: prs_dirichlet = .false.


     ! The advection operator.
     class(advection_adjoint_t), allocatable :: adv

     ! Time OIFS interpolation scheme for advection.
     logical :: oifs

     ! Time variables
     type(field_t) :: abx1, aby1, abz1
     type(field_t) :: abx2, aby2, abz2

     ! Advection terms for the oifs method
     type(field_t) :: advx, advy, advz

     !> Pressure residual equation for computing `p_res`.
     class(adjoint_pnpn_prs_res_t), allocatable :: prs_res

     !> Velocity residual equation for computing `u_res`, `v_res`, `w_res`.
     class(adjoint_pnpn_vel_res_t), allocatable :: vel_res

     !> Summation of AB/BDF contributions
     class(rhs_maker_sumab_t), allocatable :: sumab

     !> Contributions to kth order extrapolation scheme
     class(rhs_maker_ext_t), allocatable :: makeabf

     !> Contributions to F from lagged BD terms
     class(rhs_maker_bdf_t), allocatable :: makebdf

     !> Contributions to the RHS from the OIFS method
     class(rhs_maker_oifs_t), allocatable :: makeoifs

     !  !> Adjust flow volume
     !  type(fluid_volflow_t) :: vol_flow

     !> Whether to use the full formulation of the viscous stress term
     logical :: full_stress_formulation = .false.

     ! ======================================================================= !
     ! Addressable attributes

     real(kind=rp) :: norm_scaling !< Constant for the norm of the velocity.
     real(kind=rp) :: norm_target !< Target norm for the velocity field.
     real(kind=rp) :: norm_tolerance !< Tolerance for when to rescale the flow.

     ! ======================================================================= !
     ! Definition of shorthands and local variables

     !> Norm of the base field.
     real(kind=rp) :: norm_l2_base

     !> Upper limit for the norm
     real(kind=rp) :: norm_l2_upper
     !> Lower limit for the norm
     real(kind=rp) :: norm_l2_lower

     !> Output file
     type(file_t) :: file_output

   contains
     !> Constructor.
     procedure, pass(this) :: init => adjoint_fluid_pnpn_init
     !> Destructor.
     procedure, pass(this) :: free => adjoint_fluid_pnpn_free
     !> Perform a single time-step of the scheme.
     procedure, pass(this) :: step => adjoint_fluid_pnpn_step
     !> Restart from a previous solution.
     procedure, pass(this) :: restart => adjoint_fluid_pnpn_restart
     !> Set up boundary conditions.
     procedure, pass(this) :: setup_bcs => adjoint_fluid_pnpn_setup_bcs
     !> Write a field with boundary condition specifications.
     procedure, pass(this) :: write_boundary_conditions => &
          adjoint_fluid_pnpn_write_boundary_conditions

     !> Compute the power_iterations field.
     procedure, public, pass(this) :: PW_compute_ => power_iterations_compute

  end type adjoint_fluid_pnpn_t

  !> Boundary condition factory for pressure.
  !! @details Will mark a mesh zone for the bc and finalize.
  !! @param[inout] object The object to be allocated.
  !! @param[in] scheme The `scalar_pnpn` scheme.
  !! @param[inout] json JSON object for initializing the bc.
  !! @param[in] coef SEM coefficients.
  !! @param[in] user The user interface.
  interface pressure_bc_factory
     module subroutine pressure_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(adjoint_fluid_pnpn_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), target, intent(in) :: coef
       type(user_t), target, intent(in) :: user
     end subroutine pressure_bc_factory
  end interface pressure_bc_factory

  !> Boundary condition factory for velocity
  !! @details Will mark a mesh zone for the bc and finalize.
  !! @param[inout] object The object to be allocated.
  !! @param[in] scheme The `scalar_pnpn` scheme.
  !! @param[inout] json JSON object for initializing the bc.
  !! @param[in] coef SEM coefficients.
  !! @param[in] user The user interface.
  interface velocity_bc_factory
     module subroutine velocity_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(adjoint_fluid_pnpn_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), target, intent(in) :: coef
       type(user_t), target, intent(in) :: user
     end subroutine velocity_bc_factory
  end interface velocity_bc_factory

contains

  subroutine adjoint_fluid_pnpn_init(this, msh, lx, params, user, chkp)
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    character(len=15), parameter :: scheme = 'Adjoint (Pn/Pn)'
    real(kind=rp) :: abs_tol
    character(len=LOG_SIZE) :: log_buf
    integer :: integer_val, solver_maxiter
    character(len=:), allocatable :: solver_type, precon_type
    logical :: monitor, found
    logical :: advection
    type(json_file) :: numerics_params, precon_params

    ! Temporary field pointers
    character(len=:), allocatable :: file_name
    character(len=256) :: header_line

    call this%free()

    ! Initialize base class
    call this%init_base(msh, lx, params, scheme, user, .true.)

    ! Add pressure field to the registry. For this scheme it is in the same
    ! Xh as the velocity
    call neko_registry%add_field(this%dm_Xh, 'p_adj')
    this%p_adj => neko_registry%get_field('p_adj')

    !
    ! Select governing equations via associated residual and Ax types
    !

    call json_get(params, 'case.numerics.time_order', integer_val)
    allocate(this%ext_bdf)
    call this%ext_bdf%init(integer_val)

    call json_get_or_default(params, "case.fluid.full_stress_formulation", &
         this%full_stress_formulation, .false.)

    if (this%full_stress_formulation .eqv. .true.) then
       call neko_error( &
            "Full stress formulation is not supported in the adjoint module.")
       !  ! Setup backend dependent Ax routines
       !  call ax_helm_allocator(this%Ax_vel, type_name = "full")

       !  ! Setup backend dependent prs residual routines
       !  call pnpn_prs_res_stress_factory(this%prs_res)

       !  ! Setup backend dependent vel residual routines
       !  call pnpn_vel_res_stress_factory(this%vel_res)
    else
       ! Setup backend dependent Ax routines
       call ax_helm_allocator(this%Ax_vel, type_name = "standard")

       ! Setup backend dependent prs residual routines
       call adjoint_pnpn_prs_res_factory(this%prs_res)

       ! Setup backend dependent vel residual routines
       call adjoint_pnpn_vel_res_factory(this%vel_res)

       allocate(segregated_vector_bc_projector_t :: this%bcs_vel_projector)
       call this%bcs_vel_projector%init(this%c_Xh)
    end if

    if (params%valid_path('case.fluid.nut_field')) then
       if (this%full_stress_formulation .eqv. .false.) then
          call neko_error("You need to set full_stress_formulation to " // &
               "true for the fluid to have a spatially varying " // &
               "viscocity field.")
       end if
       call json_get(params, 'case.fluid.nut_field', this%nut_field_name)
    else
       this%nut_field_name = ""
    end if

    ! Setup Ax for the pressure
    call ax_helm_allocator(this%Ax_prs, type_name = "standard")


    ! Setup backend dependent summation of AB/BDF
    call rhs_maker_sumab_fctry(this%sumab)

    ! Setup backend dependent summation of extrapolation scheme
    call rhs_maker_ext_fctry(this%makeabf)

    ! Setup backend depenent contributions to F from lagged BD terms
    call rhs_maker_bdf_fctry(this%makebdf)

    ! Setup backend dependent summations of the OIFS method
    call rhs_maker_oifs_fctry(this%makeoifs)

    ! Initialize variables specific to this plan
    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call this%p_res%init(dm_Xh, "p_res")
      call this%u_res%init(dm_Xh, "u_res")
      call this%v_res%init(dm_Xh, "v_res")
      call this%w_res%init(dm_Xh, "w_res")
      call this%abx1%init(dm_Xh, "abx1")
      call this%aby1%init(dm_Xh, "aby1")
      call this%abz1%init(dm_Xh, "abz1")
      call this%abx2%init(dm_Xh, "abx2")
      call this%aby2%init(dm_Xh, "aby2")
      call this%abz2%init(dm_Xh, "abz2")
      call this%advx%init(dm_Xh, "advx")
      call this%advy%init(dm_Xh, "advy")
      call this%advz%init(dm_Xh, "advz")
    end associate

    call this%du%init(this%dm_Xh, 'du')
    call this%dv%init(this%dm_Xh, 'dv')
    call this%dw%init(this%dm_Xh, 'dw')
    call this%dp%init(this%dm_Xh, 'dp')

    ! Set up boundary conditions
    call this%setup_bcs(user, params)

    ! Check if we need to output boundaries
    call json_get_or_default(params, 'case.output_boundary', found, .false.)
    if (found) call this%write_boundary_conditions()

    call this%proj_prs%init(this%dm_Xh%size(), this%pr_projection_dim, &
         this%pr_projection_activ_step)

    call this%proj_vel%init(this%dm_Xh%size(), this%vel_projection_dim, &
         this%vel_projection_activ_step)

    ! Determine the time-interpolation scheme
    call json_get_or_default(params, 'case.numerics.oifs', this%oifs, .false.)
    if (params%valid_path('case.fluid.flow_rate_force')) then
       call neko_error("Flow rate forcing not available for adjoint_fluid_pnpn")
    end if

    ! Setup pressure solver
    call neko_log%section("Pressure solver")

    call json_get_or_default(params, &
         'case.fluid.pressure_solver.max_iterations', &
         solver_maxiter, 800)
    call json_get(params, 'case.fluid.pressure_solver.type', solver_type)
    call json_get(params, 'case.fluid.pressure_solver.preconditioner.type', &
         precon_type)
    call json_get(params, &
         'case.fluid.pressure_solver.preconditioner', precon_params)
    call json_get(params, 'case.fluid.pressure_solver.absolute_tolerance', &
         abs_tol)
    call json_get_or_default(params, 'case.fluid.pressure_solver.monitor', &
         monitor, .false.)
    call neko_log%message('Type       : ('// trim(solver_type) // &
         ', ' // trim(precon_type) // ')')
    write(log_buf, '(A,ES13.6)') 'Abs tol    :', abs_tol
    call neko_log%message(log_buf)

    call this%solver_factory(this%ksp_prs, this%dm_Xh%size(), &
         solver_type, solver_maxiter, abs_tol, monitor)
    call this%precon_factory_(this%pc_prs, this%ksp_prs, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs_prs, &
         precon_type, precon_params)
    call neko_log%end_section()

    ! Initialize the advection factory
    call neko_log%section("Advection factory")
    call json_get_or_default(params, 'case.fluid.advection', advection, .true.)
    call json_get(params, 'case.numerics', numerics_params)
    call advection_adjoint_factory(this%adv, numerics_params, this%c_Xh, &
         this%ulag, this%vlag, this%wlag, &
         chkp%dtlag, chkp%tlag, this%ext_bdf, &
         .not. advection)
    ! Should be in init_base maybe?
    this%chkp => chkp
    ! This is probably scheme specific
    ! Should not be init really, but more like, add fluid or something...
    call this%chkp%add_fluid(this%u_adj, this%v_adj, this%w_adj, this%p_adj)

    this%chkp%abx1 => this%abx1
    this%chkp%abx2 => this%abx2
    this%chkp%aby1 => this%aby1
    this%chkp%aby2 => this%aby2
    this%chkp%abz1 => this%abz1
    this%chkp%abz2 => this%abz2
    call this%chkp%add_lag(this%ulag, this%vlag, this%wlag)

    call neko_log%end_section()

    ! ------------------------------------------------------------------------ !
    ! Handling the rescaling and baseflow

    ! Read the norm scaling from the json file
    call json_get_or_default(params, 'norm_scaling', &
         this%norm_scaling, 0.5_rp)

    ! The baseflow is the solution to the forward.
    ! Userdefined baseflows can be invoked via setting initial conditions
    ! call neko_registry%add_field(this%dm_Xh, 'u')
    ! call neko_registry%add_field(this%dm_Xh, 'v')
    ! call neko_registry%add_field(this%dm_Xh, 'w')
    ! call neko_registry%add_field(this%dm_Xh, 'p')
    this%u_b => neko_registry%get_field('u')
    this%v_b => neko_registry%get_field('v')
    this%w_b => neko_registry%get_field('w')
    this%p_b => neko_registry%get_field('p')

    ! Read the json file
    call json_get_or_default(params, 'norm_target', &
         this%norm_target, -1.0_rp)
    call json_get_or_default(params, 'norm_tolerance', &
         this%norm_tolerance, 10.0_rp)

    ! Build the header
    call json_get_or_default(params, 'output_file', &
         file_name, 'power_iterations.csv')
    call this%file_output%init(trim(file_name))
    write(header_line, '(A)') 'Time, Norm, Scaling'
    call this%file_output%set_header(header_line)

  end subroutine adjoint_fluid_pnpn_init

  subroutine adjoint_fluid_pnpn_restart(this, chkp)
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    real(kind=rp) :: dtlag(10), tlag(10)
    integer :: i, n

    dtlag = chkp%dtlag
    tlag = chkp%tlag

    n = this%u_adj%dof%size()
    if (allocated(this%chkp%previous_mesh%elements) .or. &
         chkp%previous_Xh%lx .ne. this%Xh%lx) then
       associate(u => this%u_adj, v => this%v_adj, w => this%w_adj, &
            p => this%p_adj, c_Xh => this%c_Xh, ulag => this%ulag, &
            vlag => this%vlag, wlag => this%wlag)
         call col2(u%x, c_Xh%mult, n)
         call col2(v%x, c_Xh%mult, n)
         call col2(w%x, c_Xh%mult, n)
         call col2(p%x, c_Xh%mult, n)
         do i = 1, this%ulag%size()
            call col2(ulag%lf(i)%x, c_Xh%mult, n)
            call col2(vlag%lf(i)%x, c_Xh%mult, n)
            call col2(wlag%lf(i)%x, c_Xh%mult, n)
         end do
       end associate
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(u => this%u_adj, v => this%v_adj, w => this%w_adj, &
            ulag => this%ulag, vlag => this%vlag, wlag => this%wlag,&
            p => this%p_adj)
         call device_memcpy(u%x, u%x_d, u%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(v%x, v%x_d, v%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(w%x, w%x_d, w%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(p%x, p%x_d, p%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(ulag%lf(1)%x, ulag%lf(1)%x_d, &
              u%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(ulag%lf(2)%x, ulag%lf(2)%x_d, &
              u%dof%size(), HOST_TO_DEVICE, sync = .false.)

         call device_memcpy(vlag%lf(1)%x, vlag%lf(1)%x_d, &
              v%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(vlag%lf(2)%x, vlag%lf(2)%x_d, &
              v%dof%size(), HOST_TO_DEVICE, sync = .false.)

         call device_memcpy(wlag%lf(1)%x, wlag%lf(1)%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(wlag%lf(2)%x, wlag%lf(2)%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abx1%x, this%abx1%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abx2%x, this%abx2%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%aby1%x, this%aby1%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%aby2%x, this%aby2%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abz1%x, this%abz1%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%abz2%x, this%abz2%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%advx%x, this%advx%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%advy%x, this%advy%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(this%advz%x, this%advz%x_d, &
              w%dof%size(), HOST_TO_DEVICE, sync = .false.)
       end associate
    end if
    ! Make sure that continuity is maintained (important for interpolation)
    ! Do not do this for lagged rhs
    ! (derivatives are not necessairly coninous across elements)

    if (allocated(this%chkp%previous_mesh%elements) &
         .or. this%chkp%previous_Xh%lx .ne. this%Xh%lx) then
       call this%gs_Xh%op(this%u_adj, GS_OP_ADD)
       call this%gs_Xh%op(this%v_adj, GS_OP_ADD)
       call this%gs_Xh%op(this%w_adj, GS_OP_ADD)
       call this%gs_Xh%op(this%p_adj, GS_OP_ADD)

       do i = 1, this%ulag%size()
          call this%gs_Xh%op(this%ulag%lf(i), GS_OP_ADD)
          call this%gs_Xh%op(this%vlag%lf(i), GS_OP_ADD)
          call this%gs_Xh%op(this%wlag%lf(i), GS_OP_ADD)
       end do
    end if

  end subroutine adjoint_fluid_pnpn_restart

  subroutine adjoint_fluid_pnpn_free(this)
    class(adjoint_fluid_pnpn_t), intent(inout) :: this

    !Deallocate velocity and pressure fields
    call this%scheme_free()

    call this%bc_prs_surface%free()
    call this%bc_sym_surface%free()
    call this%bc_curl_curl%free()
    if (allocated(this%bcs_vel_projector)) then
       call this%bcs_vel_projector%free()
       deallocate(this%bcs_vel_projector)
    end if
    call this%bcs_prs_projector%free()
    call this%proj_prs%free()
    call this%proj_vel%free()

    call this%p_res%free()
    call this%u_res%free()
    call this%v_res%free()
    call this%w_res%free()

    call this%du%free()
    call this%dv%free()
    call this%dw%free()
    call this%dp%free()

    call this%abx1%free()
    call this%aby1%free()
    call this%abz1%free()

    call this%abx2%free()
    call this%aby2%free()
    call this%abz2%free()

    call this%advx%free()
    call this%advy%free()
    call this%advz%free()

    if (allocated(this%Ax_vel)) then
       deallocate(this%Ax_vel)
    end if

    if (allocated(this%Ax_prs)) then
       deallocate(this%Ax_prs)
    end if

    if (allocated(this%prs_res)) then
       deallocate(this%prs_res)
    end if

    if (allocated(this%vel_res)) then
       deallocate(this%vel_res)
    end if

    if (allocated(this%sumab)) then
       deallocate(this%sumab)
    end if

    if (allocated(this%makeabf)) then
       deallocate(this%makeabf)
    end if

    if (allocated(this%makebdf)) then
       deallocate(this%makebdf)
    end if

    if (allocated(this%makeoifs)) then
       deallocate(this%makeoifs)
    end if

    if (allocated(this%ext_bdf)) then
       deallocate(this%ext_bdf)
    end if

    ! call this%vol_flow%free()

  end subroutine adjoint_fluid_pnpn_free

  !> Advance fluid simulation in time.
  !! @param this The fluid simulation object.
  !! @param time The time state object.
  !! @param dt_controller timestep controller
  subroutine adjoint_fluid_pnpn_step(this, time, dt_controller)
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    ! number of degrees of freedom
    integer :: n
    ! Solver results monitors (pressure + 3 velocity)
    type(ksp_monitor_t) :: ksp_results(4)
    type(field_t), pointer :: dx_p_adj, dy_p_adj, dz_p_adj, nx1, nx2, nx3, &
         work1, work2
    integer :: temp_indices(3)
    integer :: cc_indices(8)
    real(kind=rp) :: rho_val, mu_val

    if (this%freeze) return

    n = this%dm_Xh%size()

    call profiler_start_region('Adjoint')
    associate(u => this%u_adj, v => this%v_adj, w => this%w_adj, &
         p => this%p_adj, &
         u_e => this%u_adj_e, v_e => this%v_adj_e, w_e => this%w_adj_e, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_b => this%u_b, v_b => this%v_b, w_b => this%w_b, &
         u_res => this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax_vel => this%Ax_vel, Ax_prs => this%Ax_prs, &
         Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         msh => this%msh, prs_res => this%prs_res, &
         source_term => this%source_term, vel_res => this%vel_res, &
         sumab => this%sumab, &
         makeabf => this%makeabf, makebdf => this%makebdf, &
         vel_projection_dim => this%vel_projection_dim, &
         pr_projection_dim => this%pr_projection_dim, &
         oifs => this%oifs, &
         rho => this%rho, mu => this%mu, &
         f_x => this%f_adj_x, f_y => this%f_adj_y, f_z => this%f_adj_z, &
         t => time%t, tstep => time%tstep, dt => time%dt, &
         ext_bdf => this%ext_bdf, event => glb_cmd_event)

      ! Extrapolate the velocity if it's not done in nut_field estimation
      call sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
           ulag, vlag, wlag, ext_bdf%advection_coeffs%x, ext_bdf%nadv)

      ! Compute the source terms
      call this%source_term%compute(time)

      ! Add Neumann bc contributions to the RHS
      call this%bcs_vel%apply_vector(f_x%x, f_y%x, f_z%x, &
           this%dm_Xh%size(), time, strong = .false.)

      if (oifs) then
         call neko_error("OIFS not implemented for adjoint")

      else
         ! Add the advection operators to the right-hand-side.
         call this%adv%compute_adjoint(u, v, w, u_b, v_b, w_b, &
              f_x, f_y, f_z, &
              Xh, this%c_Xh, dm_Xh%size())

         ! At this point the RHS contains the sum of the advection operator and
         ! additional source terms, evaluated using the velocity field from the
         ! previous time-step. Now, this value is used in the explicit time
         ! scheme to advance both terms in time.
         call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
              this%abx2, this%aby2, this%abz2, &
              f_x%x, f_y%x, f_z%x, &
              rho%x(1,1,1,1), ext_bdf%advection_coeffs%x, n)

         ! Add the RHS contributions coming from the BDF scheme.
         call makebdf%compute_fluid(ulag, vlag, wlag, f_x%x, f_y%x, f_z%x, &
              u, v, w, c_Xh%B, c_Xh%Blag, c_Xh%Blaglag, rho%x(1,1,1,1), dt, &
              ext_bdf%diffusion_coeffs%x, ext_bdf%ndiff, n)
      end if

      call ulag%update()
      call vlag%update()
      call wlag%update()

      call this%bc_apply_vel(time, strong = .true.)
      call this%bc_apply_prs(time)

      ! Now we need the surface contribution of the curl curl BC.(explicit in p)
      call neko_scratch_registry%request_field(dx_p_adj, cc_indices(1), .false.)
      call neko_scratch_registry%request_field(dy_p_adj, cc_indices(2), .false.)
      call neko_scratch_registry%request_field(dz_p_adj, cc_indices(3), .false.)

      ! Note: zero interior
      call neko_scratch_registry%request_field(nx1, cc_indices(4), .true.)
      call neko_scratch_registry%request_field(nx2, cc_indices(5), .true.)
      call neko_scratch_registry%request_field(nx3, cc_indices(6), .true.)

      call neko_scratch_registry%request_field(work1, cc_indices(7), .false.)
      call neko_scratch_registry%request_field(work2, cc_indices(8), .false.)

      ! gradient of adjoint pressure (explicit)
      call grad(dx_p_adj%x, dy_p_adj%x, dz_p_adj%x, this%p_adj%x, c_Xh)

      ! Now we compute the n x grad(p) (this include 2D weights)
      call this%bc_curl_curl%apply_n_cross(nx1%x, nx2%x, nx3%x, dx_p_adj%x, &
           dy_p_adj%x, dz_p_adj%x, dx_p_adj%size())

      ! Now we need curl on the test function, note that transpose of curl is
      ! negative curl
      ! reuse dx_p_adj etc as fx, fy, fz etc
      call curl(dx_p_adj, dy_p_adj, dz_p_adj, nx1, nx2, nx3, work1, work2, c_Xh)

      ! Forward does gsop on the residual (which has the pressure gradient)
      call gs_Xh%op(dx_p_adj, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(dy_p_adj, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(dz_p_adj, GS_OP_ADD, event)
      call device_event_sync(event)

      ! multiplcity
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(dx_p_adj%x_d, c_Xh%mult_d, dx_p_adj%size())
         call device_col2(dy_p_adj%x_d, c_Xh%mult_d, dx_p_adj%size())
         call device_col2(dz_p_adj%x_d, c_Xh%mult_d, dx_p_adj%size())
      else
         call col2(dx_p_adj%x, c_Xh%mult, dx_p_adj%size())
         call col2(dy_p_adj%x, c_Xh%mult, dx_p_adj%size())
         call col2(dz_p_adj%x, c_Xh%mult, dx_p_adj%size())
      end if

      rho_val = rho%x(1,1,1,1)
      mu_val = mu%x(1,1,1,1)
      call field_add2s2(f_x, dx_p_adj, -mu_val / rho_val)
      call field_add2s2(f_y, dy_p_adj, -mu_val / rho_val)
      call field_add2s2(f_z, dz_p_adj, -mu_val / rho_val)

      call neko_scratch_registry%relinquish_field(cc_indices)

      ! Update material properties if necessary
      call this%update_material_properties(time)

      ! Compute intermediate velocity residual.

      call profiler_start_region('Adjoint_velocity_residual')

      call vel_res%compute(Ax_vel, u, v, w, &
           u_res, v_res, w_res, &
           p, &
           f_x, f_y, f_z, &
           c_Xh, msh, Xh, &
           mu, rho, ext_bdf%diffusion_coeffs%x(1), &
           dt, dm_Xh%size())

      call gs_Xh%op(u_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(v_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(w_res, GS_OP_ADD, event)
      call device_event_sync(event)

      ! Set residual to zero at strong velocity boundaries.
      call this%bcs_vel_projector%apply(u_res%x, v_res%x, w_res%x, n)

      call profiler_end_region('Adjoint_velocity_residual')

      call this%proj_vel%pre_solving(u_res%x, v_res%x, w_res%x, &
           tstep, c_Xh, n, dt_controller, 'Velocity')

      call this%pc_vel%update()

      call profiler_start_region("Adjoint_velocity_solve")
      ksp_results(1:3) = this%ksp_vel%solve_coupled(Ax_vel, du, dv, dw, &
           u_res%x, v_res%x, w_res%x, n, c_Xh, &
           this%bcs_vel_projector, gs_Xh, &
           this%ksp_vel%max_iter)

      call this%proj_vel%post_solving(du%x, dv%x, dw%x, Ax_vel, c_Xh, &
           this%bcs_vel_projector, gs_Xh, n, tstep, &
           dt_controller)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opadd2cm(u%x_d, v%x_d, w%x_d, &
              du%x_d, dv%x_d, dw%x_d, 1.0_rp, n, msh%gdim)
      else
         call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, 1.0_rp, n, msh%gdim)
      end if
      call profiler_end_region("Adjoint_velocity_solve")

      !------------------------------------------------------------------------!
      ! now the RHS of our pressure eqn is the new adjoint velocity
      ! be careful with the order of the gsops here, we will handle this in the
      ! residual calculation. So we enter WITHOUT a mass matrix.
      call field_copy(f_x, u)
      call field_copy(f_y, v)
      call field_copy(f_z, w)
      !------------------------------------------------------------------------!
      call profiler_start_region('Adjoint_pressure_residual')

      call prs_res%compute(p, p_res, &
           u, v, w, &
           f_x, f_y, f_z, &
           c_Xh, gs_Xh, &
           this%bc_prs_surface, this%bc_sym_surface, &
           Ax_prs, ext_bdf%diffusion_coeffs%x(1), dt, &
           mu, rho, event)

      ! De-mean the pressure residual when no strong pressure boundaries present
      if (.not. this%prs_dirichlet .and. NEKO_BCKND_DEVICE .eq. 1) then
         call device_ortho(p_res%x_d, this%glb_n_points, n)
      else if (.not. this%prs_dirichlet) then
         call ortho(p_res%x, this%glb_n_points, n)
      end if

      call gs_Xh%op(p_res, GS_OP_ADD, event)
      call device_event_sync(event)

      ! Set the residual to zero at strong pressure boundaries.
      call this%bcs_prs_projector%apply(p_res%x, p%dof%size())


      call profiler_end_region('Adjoint_pressure_residual')


      call this%proj_prs%pre_solving(p_res%x, tstep, c_Xh, n, dt_controller, &
           'Pressure')

      call this%pc_prs%update()

      call profiler_start_region('Adjoint_pressure_solve')

      ! Solve for the pressure increment.
      ksp_results(4) = &
           this%ksp_prs%solve(Ax_prs, dp, p_res%x, n, c_Xh, &
           this%bcs_prs_projector, gs_Xh)


      call profiler_end_region('Adjoint_pressure_solve')

      call this%proj_prs%post_solving(dp%x, Ax_prs, c_Xh, &
           this%bcs_prs_projector, gs_Xh, n, tstep, dt_controller)

      ! Update the pressure with the increment. Demean if necessary.
      call field_add2(p, dp, n)
      if (.not. this%prs_dirichlet .and. NEKO_BCKND_DEVICE .eq. 1) then
         call device_ortho(p%x_d, this%glb_n_points, n)
      else if (.not. this%prs_dirichlet) then
         call ortho(p%x, this%glb_n_points, n)
      end if

      ksp_results(4)%name = 'Adjoint Pressure'
      ksp_results(1)%name = 'Adjoint Velocity U'
      ksp_results(2)%name = 'Adjoint Velocity V'
      ksp_results(3)%name = 'Adjoint Velocity W'

      if (this%forced_flow_rate) then
         call neko_error('Forced flow rate is not implemented for the adjoint')

      end if

      !------------------------------------------------------------------------!
      ! correct the velocity with the pressure
      call neko_scratch_registry%request_field(dx_p_adj, temp_indices(1), &
           .false.)
      call neko_scratch_registry%request_field(dy_p_adj, temp_indices(2), &
           .false.)
      call neko_scratch_registry%request_field(dz_p_adj, temp_indices(3), &
           .false.)

      ! gradient of adjoint pressure (explicit)
      call opgrad(dx_p_adj%x, dy_p_adj%x, dz_p_adj%x, this%p_adj%x, c_Xh)

      ! they gsop the residual (which has the pressure gradient)
      call gs_Xh%op(dx_p_adj, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(dy_p_adj, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(dz_p_adj, GS_OP_ADD, event)
      call device_event_sync(event)

      ! divide by mass matrix
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(dx_p_adj%x_d, c_Xh%Binv_d, dx_p_adj%size())
         call device_col2(dy_p_adj%x_d, c_Xh%Binv_d, dx_p_adj%size())
         call device_col2(dz_p_adj%x_d, c_Xh%Binv_d, dx_p_adj%size())
      else
         ! NOTE. This term comes from the handling of the pressure RHS, which
         ! DOES include the multiplicity in the op.
         call col2(dx_p_adj%x, c_Xh%Binv, dx_p_adj%size())
         call col2(dy_p_adj%x, c_Xh%Binv, dx_p_adj%size())
         call col2(dz_p_adj%x, c_Xh%Binv, dx_p_adj%size())
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opadd2cm(u%x_d, v%x_d, w%x_d, dx_p_adj%x_d, &
              dy_p_adj%x_d, dz_p_adj%x_d, -1.0_rp, n, msh%gdim)
      else
         call opadd2cm(u%x, v%x, w%x, dx_p_adj%x, dy_p_adj%x, dz_p_adj%x, &
              -1.0_rp, n, msh%gdim)
      end if

      call neko_scratch_registry%relinquish_field(temp_indices)
      !------------------------------------------------------------------------!

      call fluid_step_info(time, ksp_results, &
           this%full_stress_formulation, this%strict_convergence)

    end associate
    call profiler_end_region('Adjoint')

  end subroutine adjoint_fluid_pnpn_step

  !> Sets up the boundary condition for the scheme.
  !! @param this The fluid simulation object.
  !! @param user The user interface.
  !! @param params The json file containing the parameters.
  subroutine adjoint_fluid_pnpn_setup_bcs(this, user, params)
    use mpi_f08, only: MPI_IN_PLACE
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    type(user_t), target, intent(in) :: user
    type(json_file), intent(inout) :: params
    integer :: i, n_bcs, j, zone_size, global_zone_size, ierr
    class(bc_t), pointer :: bc_i
    type(json_core) :: core
    type(json_value), pointer :: bc_object
    type(json_file) :: bc_subdict
    logical :: found
    ! Monitor which boundary zones have been marked
    logical, allocatable :: marked_zones(:)
    integer, allocatable :: zone_indices(:)
    character(len=:), allocatable :: json_key

    ! Special PnPn boundary conditions for pressure
    call this%bc_prs_surface%init_from_components(this%c_Xh)
    call this%bc_sym_surface%init_from_components(this%c_Xh)
    call this%bc_curl_curl%init_from_components(this%c_Xh)

    json_key = 'case.adjoint_fluid.boundary_conditions'

    ! Populate bcs_vel and bcs_prs based on the case file
    if (params%valid_path(json_key)) then
       call params%info(json_key, n_children = n_bcs)
       call params%get_core(core)
       call params%get(json_key, bc_object, found)

       !
       ! Velocity bcs
       !
       call this%bcs_vel%init(n_bcs)

       allocate(marked_zones(size(this%msh%labeled_zones)))
       marked_zones = .false.

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)

          call json_get(bc_subdict, "zone_indices", zone_indices)

          ! Check that we are not trying to assing a bc to zone, for which one
          ! has already been assigned and that the zone has more than 0 size
          ! in the mesh.
          do j = 1, size(zone_indices)
             zone_size = this%msh%labeled_zones(zone_indices(j))%size
             call MPI_Allreduce(zone_size, global_zone_size, 1, &
                  MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

             if (global_zone_size .eq. 0) then
                write(error_unit, '(A, A, I0, A, A, I0, A)') "*** ERROR ***: ",&
                     "Zone index ", zone_indices(j), &
                     " is invalid as this zone has 0 size, meaning it does ", &
                     "not in the mesh. Check adjoint boundary condition ", &
                     i, "."
                error stop
             end if

             if (marked_zones(zone_indices(j)) .eqv. .true.) then
                write(error_unit, '(A, A, I0, A, A, A, A)') "*** ERROR ***: ", &
                     "Zone with index ", zone_indices(j), &
                     " has already been assigned a boundary condition. ", &
                     "Please check your boundary_conditions entry for the ", &
                     "adjoint and make sure that each zone index appears ", &
                     "only in a single boundary condition."
                error stop
             else
                marked_zones(zone_indices(j)) = .true.
             end if
          end do

          bc_i => null()
          call velocity_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for velocity in particular,
          ! so we check.
          if (associated(bc_i)) then

             select type (bc_i)
             type is (symmetry_aligned_t)
                call this%bcs_vel_projector%mark(bc_i%bc_x, component = 'x')
                call this%bcs_vel_projector%mark(bc_i%bc_y, component = 'y')
                call this%bcs_vel_projector%mark(bc_i%bc_z, component = 'z')
                call this%bcs_vel%append(bc_i)
                call this%bc_sym_surface%mark_facets(bc_i%marked_facet)
             type is (non_normal_aligned_t)
                call this%bcs_vel_projector%mark(bc_i%bc_x, component = 'x')
                call this%bcs_vel_projector%mark(bc_i%bc_y, component = 'y')
                call this%bcs_vel_projector%mark(bc_i%bc_z, component = 'z')
                call this%bcs_vel%append(bc_i)
             type is (shear_stress_t)
                call neko_error("The shear_stress boundary condition " // &
                     "requires the full stress formulation to be enabled.")
             type is (wall_model_bc_t)
                call neko_error("The wall_model boundary condition " // &
                     "requires the full stress formulation to be enabled.")
             class default

                ! Additionally we mark the special PnPn pressure  bc.
                if (bc_i%bc_type .eq. BC_DIRICHLET) then
                   call this%bc_prs_surface%mark_labeled_zones( &
                        bc_i%zone_indices)
                   call this%bcs_vel_projector%mark(bc_i, component = 'x')
                   call this%bcs_vel_projector%mark(bc_i, component = 'y')
                   call this%bcs_vel_projector%mark(bc_i, component = 'z')
                end if

                ! add all BCs to curl curl
                call this%bc_curl_curl%mark_facets(bc_i%marked_facet)

                call this%bcs_vel%append(bc_i)
             end select
          end if
       end do

       ! Make sure all labeled zones with non-zero size have been marked
       do i = 1, size(this%msh%labeled_zones)
          if ((this%msh%labeled_zones(i)%size .gt. 0) .and. &
               (marked_zones(i) .eqv. .false.)) then
             write(error_unit, '(A, A, I0)') "*** ERROR ***: ", &
                  "No adjoint boundary condition assigned to zone ", i
             error stop
          end if
       end do

       !
       ! Pressure bcs
       !
       call this%bcs_prs%init(n_bcs)

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)
          bc_i => null()
          call pressure_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for pressure in particular,
          ! so we check.
          if (associated(bc_i)) then
             call this%bcs_prs%append(bc_i)

             ! Mark strong pressure bcs in the projector to force zero change.
             if (bc_i%bc_type .eq. BC_DIRICHLET) then
                call this%bcs_prs_projector%mark(bc_i)
             end if

          end if

       end do
    else
       ! Check that there are no labeled zones, i.e. all are periodic.
       do i = 1, size(this%msh%labeled_zones)
          if (this%msh%labeled_zones(i)%size .gt. 0) then
             call neko_error("No boundary_conditions entry in the case file!")
          end if
       end do

       call this%bcs_vel%init()
       call this%bcs_prs%init()

    end if

    call this%bc_prs_surface%finalize()
    call this%bc_sym_surface%finalize()
    call this%bc_curl_curl%finalize()
    call this%bcs_vel_projector%finalize(rebuild_mask = .true.)

    ! If we have no strong pressure bcs, we will demean the pressure
    this%prs_dirichlet = this%bcs_prs_projector%dof_mask%is_set()
    call MPI_Allreduce(MPI_IN_PLACE, this%prs_dirichlet, 1, &
         MPI_LOGICAL, MPI_LOR, NEKO_COMM)

  end subroutine adjoint_fluid_pnpn_setup_bcs

  !> Write a field with boundary condition specifications
  subroutine adjoint_fluid_pnpn_write_boundary_conditions(this)
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    type(dirichlet_t) :: bdry_mask
    type(field_t), pointer :: bdry_field
    type(file_t) :: bdry_file
    integer :: temp_index, i
    class(bc_t), pointer :: bci
    character(len=LOG_SIZE) :: log_buf

    call neko_log%section("Adjoint boundary conditions")
    write(log_buf, '(A)') &
         'Marking using integer keys in boundary_adjoint0.f00000'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') 'Condition-value pairs: '
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  no_slip                         = 1'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  velocity_value                  = 2'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  outflow, normal_outflow (+dong) = 3'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  symmetry                        = 4'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  user_velocity_pointwise         = 5'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  periodic                        = 6'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  user_velocity                   = 7'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  user_pressure                   = 8'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  shear_stress                    = 9'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  wall_modelling                  = 10'
    call neko_log%message(log_buf)
    write(log_buf, '(A)') '  blasius_profile                 = 11'
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call neko_scratch_registry%request_field(bdry_field, temp_index, .true.)



    call bdry_mask%init_from_components(this%c_Xh, 6.0_rp)
    call bdry_mask%mark_zone(this%msh%periodic)
    call bdry_mask%finalize()
    call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
    call bdry_mask%free()

    do i = 1, this%bcs_prs%size()
       bci => this%bcs_prs%get(i)
       select type (bc => bci)
       type is (zero_dirichlet_t)
          call bdry_mask%init_from_components(this%c_Xh, 3.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (dong_outflow_t)
          call bdry_mask%init_from_components(this%c_Xh, 3.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (field_dirichlet_t)
          call bdry_mask%init_from_components(this%c_Xh, 8.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       end select
    end do

    do i = 1, this%bcs_vel%size()
       bci => this%bcs_vel%get(i)
       select type (bc => bci)
       type is (zero_dirichlet_t)
          call bdry_mask%init_from_components(this%c_Xh, 1.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (inflow_t)
          call bdry_mask%init_from_components(this%c_Xh, 2.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (symmetry_aligned_t)
          call bdry_mask%init_from_components(this%c_Xh, 4.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (field_dirichlet_vector_t)
          call bdry_mask%init_from_components(this%c_Xh, 7.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (shear_stress_t)
          call bdry_mask%init_from_components(this%c_Xh, 9.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (wall_model_bc_t)
          call bdry_mask%init_from_components(this%c_Xh, 10.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       type is (blasius_t)
          call bdry_mask%init_from_components(this%c_Xh, 11.0_rp)
          call bdry_mask%mark_facets(bci%marked_facet)
          call bdry_mask%finalize()
          call bdry_mask%apply_scalar(bdry_field%x, this%dm_Xh%size())
          call bdry_mask%free()
       end select
    end do


    call bdry_file%init('boundary_adjoint.fld')
    call bdry_file%write(bdry_field)

    call neko_scratch_registry%relinquish_field(temp_index)
  end subroutine adjoint_fluid_pnpn_write_boundary_conditions

  ! End of section to verify
  ! ========================================================================== !

  subroutine rescale_fluid(fluid_data, scale)
    !> Fluid data
    class(adjoint_fluid_pnpn_t), intent(inout) :: fluid_data
    !> Scaling factor
    real(kind=rp), intent(in) :: scale

    ! Local variables
    integer :: i

    ! Scale the velocity fields
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(fluid_data%u_adj%x_d, scale, fluid_data%u_adj%size())
       call device_cmult(fluid_data%v_adj%x_d, scale, fluid_data%v_adj%size())
       call device_cmult(fluid_data%w_adj%x_d, scale, fluid_data%w_adj%size())
    else
       call cmult(fluid_data%u_adj%x, scale, fluid_data%u_adj%size())
       call cmult(fluid_data%v_adj%x, scale, fluid_data%v_adj%size())
       call cmult(fluid_data%w_adj%x, scale, fluid_data%w_adj%size())
    end if

    ! Scale the right hand sides
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(fluid_data%f_adj_x%x_d, scale, &
            fluid_data%f_adj_x%size())
       call device_cmult(fluid_data%f_adj_y%x_d, scale, &
            fluid_data%f_adj_y%size())
       call device_cmult(fluid_data%f_adj_z%x_d, scale, &
            fluid_data%f_adj_z%size())
       ! HARRY
       ! maybe the abx's too
       call device_cmult(fluid_data%abx1%x_d, scale, fluid_data%abx1%size())
       call device_cmult(fluid_data%aby1%x_d, scale, fluid_data%aby1%size())
       call device_cmult(fluid_data%abz1%x_d, scale, fluid_data%abz1%size())
       call device_cmult(fluid_data%abx2%x_d, scale, fluid_data%abx2%size())
       call device_cmult(fluid_data%aby2%x_d, scale, fluid_data%aby2%size())
       call device_cmult(fluid_data%abz2%x_d, scale, fluid_data%abz2%size())

    else
       call cmult(fluid_data%f_adj_x%x, scale, fluid_data%f_adj_x%size())
       call cmult(fluid_data%f_adj_y%x, scale, fluid_data%f_adj_y%size())
       call cmult(fluid_data%f_adj_z%x, scale, fluid_data%f_adj_z%size())

       call cmult(fluid_data%abx1%x, scale, fluid_data%abx1%size())
       call cmult(fluid_data%aby1%x, scale, fluid_data%aby1%size())
       call cmult(fluid_data%abz1%x, scale, fluid_data%abz1%size())

       call cmult(fluid_data%abx2%x, scale, fluid_data%abx2%size())
       call cmult(fluid_data%aby2%x, scale, fluid_data%aby2%size())
       call cmult(fluid_data%abz2%x, scale, fluid_data%abz2%size())
    end if

    ! Scale the lag terms
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, fluid_data%ulag%size()
          call device_cmult(fluid_data%ulag%lf(i)%x_d, &
               scale, fluid_data%ulag%lf(i)%size())
       end do

       do i = 1, fluid_data%vlag%size()
          call device_cmult(fluid_data%vlag%lf(i)%x_d, &
               scale, fluid_data%vlag%lf(i)%size())
       end do

       do i = 1, fluid_data%wlag%size()
          call device_cmult(fluid_data%wlag%lf(i)%x_d, &
               scale, fluid_data%wlag%lf(i)%size())
       end do
    else
       do i = 1, fluid_data%ulag%size()
          call cmult(fluid_data%ulag%lf(i)%x, &
               scale, fluid_data%ulag%lf(i)%size())
       end do

       do i = 1, fluid_data%vlag%size()
          call cmult(fluid_data%vlag%lf(i)%x, &
               scale, fluid_data%vlag%lf(i)%size())
       end do

       do i = 1, fluid_data%wlag%size()
          call cmult(fluid_data%wlag%lf(i)%x, &
               scale, fluid_data%wlag%lf(i)%size())
       end do
    end if

  end subroutine rescale_fluid

  function norm(x, y, z, B, volume, n)
    use mpi_f08, only: MPI_IN_PLACE

    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: x, y, z
    real(kind=rp), dimension(n), intent(in) :: B
    real(kind=rp), intent(in) :: volume

    real(kind=rp) :: norm

    norm = vlsc3(x, x, B, n) + vlsc3(y, y, B, n) + vlsc3(z, z, B, n)

    call mpi_allreduce(MPI_IN_PLACE, norm, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)

    norm = sqrt(norm / volume)
  end function norm

  function device_norm(x_d, y_d, z_d, B_d, volume, n)
    use mpi_f08, only: MPI_IN_PLACE

    type(c_ptr), intent(in) :: x_d, y_d, z_d
    type(c_ptr), intent(in) :: B_d
    real(kind=rp), intent(in) :: volume
    integer, intent(in) :: n

    real(kind=rp) :: device_norm

    device_norm = device_vlsc3(x_d, x_d, B_d, n) + &
         device_vlsc3(y_d, y_d, B_d, n) + &
         device_vlsc3(z_d, z_d, B_d, n)

    call mpi_allreduce(MPI_IN_PLACE, device_norm, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)

    device_norm = sqrt(device_norm / volume)

  end function device_norm

  !> Compute the power_iterations field.
  !! @param this The fluid simulation object.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine power_iterations_compute(this, t, tstep)
    class(adjoint_fluid_pnpn_t), target, intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Local variables
    real(kind=rp) :: scaling_factor
    real(kind=rp) :: norm_l2, norm_l2_base
    character(len=256) :: log_message
    type(vector_t) :: data_line
    integer :: n

    n = this%c_Xh%dof%size()
    if (tstep .eq. 1) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          norm_l2_base = device_norm(this%u_adj%x_d, this%v_adj%x_d, &
               this%w_adj%x_d, &
               this%c_Xh%B_d, this%c_Xh%volume, n)
       else
          norm_l2_base = this%norm_scaling * norm(this%u_adj%x, this%v_adj%x, &
               this%w_adj%x, &
               this%c_Xh%B, this%c_Xh%volume, n)
       end if
       if (this%norm_target .lt. 0.0_rp) then
          this%norm_target = norm_l2_base
       end if

       this%norm_l2_upper = this%norm_tolerance * this%norm_target
       this%norm_l2_lower = this%norm_target / this%norm_tolerance

    end if

    ! Compute the norm of the velocity field and eigenvalue estimate
    if (NEKO_BCKND_DEVICE .eq. 1) then
       norm_l2 = device_norm(this%u_adj%x_d, this%v_adj%x_d, this%w_adj%x_d, &
            this%c_Xh%B_d, this%c_Xh%volume, n)
    else
       norm_l2 = norm(this%u_adj%x, this%v_adj%x, this%w_adj%x, &
            this%c_Xh%B, this%c_Xh%volume, n)
    end if
    norm_l2 = sqrt(this%norm_scaling) * norm_l2
    scaling_factor = 1.0_rp

    ! Rescale the flow if necessary
    if (norm_l2 .gt. this%norm_l2_upper &
         .or. norm_l2 .lt. this%norm_l2_lower) then
       scaling_factor = this%norm_target / norm_l2
       call rescale_fluid(this, scaling_factor)
       norm_l2 = this%norm_target

       if (tstep .eq. 1) then
          scaling_factor = 1.0_rp
       end if
    end if

    ! Log the results
    !call neko_log%section('Power Iterations', lvl = NEKO_LOG_DEBUG)
    call neko_log%section('Power Iterations')

    write (log_message, '(A7,E20.14)') 'Norm: ', norm_l2
    call neko_log%message(log_message, lvl = NEKO_LOG_DEBUG)
    write (log_message, '(A7,E20.14)') 'Scaling: ', scaling_factor
    call neko_log%message(log_message, lvl = NEKO_LOG_DEBUG)

    ! Save to file
    call data_line%init(2)
    data_line%x = [norm_l2, scaling_factor]
    call this%file_output%write(data_line, t)

    !call neko_log%end_section('Power Iterations', lvl = NEKO_LOG_DEBUG)
    call neko_log%end_section('Power Iterations')
  end subroutine power_iterations_compute

end module adjoint_fluid_pnpn
