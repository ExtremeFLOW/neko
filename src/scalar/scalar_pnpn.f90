! Copyright (c) 2022-2024, The Neko Authors
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
!> Contains the `scalar_pnpn_t` type.

module scalar_pnpn
  use num_types, only: rp
  use, intrinsic :: iso_fortran_env, only: error_unit
  use rhs_maker, only : rhs_maker_bdf_t, rhs_maker_ext_t, rhs_maker_oifs_t, &
       rhs_maker_ext_fctry, rhs_maker_bdf_fctry, rhs_maker_oifs_fctry
  use scalar_scheme, only : scalar_scheme_t
  use checkpoint, only : chkp_t
  use field, only : field_t
  use bc_list, only : bc_list_t
  use mesh, only : mesh_t
  use coefs, only : coef_t
  use device, only : HOST_TO_DEVICE, device_memcpy
  use gather_scatter, only : gs_t, GS_OP_ADD
  use scalar_residual, only : scalar_residual_t, scalar_residual_factory
  use ax_product, only : ax_t, ax_helm_factory
  use field_series, only: field_series_t
  use field_registry, only: neko_field_registry
  use facet_normal, only : facet_normal_t
  use krylov, only : ksp_monitor_t
  use device_math, only : device_add2s2, device_col2
  use scalar_aux, only : scalar_step_info
  use time_scheme_controller, only : time_scheme_controller_t
  use projection, only : projection_t
  use math, only : glsc2, col2, add2s2
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use advection, only : advection_t, advection_factory
  use profiler, only : profiler_start_region, profiler_end_region
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use json_module, only : json_file, json_core, json_value
  use user_intf, only : user_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use zero_dirichlet, only : zero_dirichlet_t
  use time_step_controller, only : time_step_controller_t
  use scratch_registry, only : neko_scratch_registry
  use time_state, only : time_state_t
  use bc, only : bc_t
  use comm, only : NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_MAX
  implicit none
  private


  type, public, extends(scalar_scheme_t) :: scalar_pnpn_t

     !> The residual of the transport equation.
     type(field_t) :: s_res

     !> Solution increment.
     type(field_t) :: ds

     !> Helmholz operator.
     class(ax_t), allocatable :: Ax

     !> Solution projection.
     type(projection_t) :: proj_s

     !> Dirichlet conditions for the residual
     !! Collects all the Dirichlet condition facets into one bc and applies 0,
     !! Since the values never change there during the solve.
     type(zero_dirichlet_t) :: bc_res

     !> A bc list for the bc_res. Contains only that, essentially just to wrap
     !! the if statement determining whether to apply on the device or CPU.
     !! Also needed since a bc_list is the type that is sent to, e.g. solvers,
     !! cannot just send `bc_res` on its own.
     type(bc_list_t) :: bclst_ds

     !> Advection operator.
     class(advection_t), allocatable :: adv

     ! Time interpolation scheme
     logical :: oifs

     ! Advection terms for the oifs method
     type(field_t) :: advs

     !> Computes the residual of the equation, i.e. `s_res`.
     class(scalar_residual_t), allocatable :: res

     !> Contributions to kth order extrapolation scheme.
     class(rhs_maker_ext_t), allocatable :: makeext

     !> Contributions to the RHS from lagged BDF terms.
     class(rhs_maker_bdf_t), allocatable :: makebdf

     !> Contributions to the RHS from the OIFS method.
     class(rhs_maker_oifs_t), allocatable :: makeoifs

     !> Lag arrays
     type(field_t) :: abx1, abx2

   contains
     !> Constructor.
     procedure, pass(this) :: init => scalar_pnpn_init
     !> To restart
     procedure, pass(this) :: restart => scalar_pnpn_restart
     !> Destructor.
     procedure, pass(this) :: free => scalar_pnpn_free
     !> Solve for the current timestep.
     procedure, pass(this) :: step => scalar_pnpn_step
     !> Setup the boundary conditions
     procedure, pass(this) :: setup_bcs_ => scalar_pnpn_setup_bcs_
     !> Sync lag field data to registry for checkpointing
  end type scalar_pnpn_t

  interface
     !> Boundary condition factory. Both constructs and initializes the object.
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `scalar_pnpn` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     module subroutine bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(scalar_pnpn_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine bc_factory
  end interface

contains

  !> Constructor.
  !! @param msh The mesh.
  !! @param coef The coefficients.
  !! @param gs The gather-scatter.
  !! @param params The case parameter file in json.
  !! @param user Type with user-defined procedures.
  !! @param chkp Set up checkpoint for restarts.
  !! @param ulag Lag arrays for the x velocity component.
  !! @param vlag Lag arrays for the y velocity component.
  !! @param wlag Lag arrays for the z velocity component.
  !! @param time_scheme The time-integration controller.
  !! @param rho The fluid density.
  subroutine scalar_pnpn_init(this, msh, coef, gs, params, numerics_params, &
       user, chkp, ulag, vlag, wlag, time_scheme, rho)
    class(scalar_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    type(field_t), target, intent(in) :: rho
    integer :: i
    class(bc_t), pointer :: bc_i
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'
    logical :: advection

    call this%free()

    ! Initiliaze base type.
    call this%scheme_init(msh, coef, gs, params, scheme, user, rho)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%ax, full_formulation = .false.)

    ! Setup backend dependent scalar residual routines
    call scalar_residual_factory(this%res)

    ! Setup backend dependent summation of extrapolation scheme
    call rhs_maker_ext_fctry(this%makeext)

    ! Setup backend dependent contributions to F from lagged BD terms
    call rhs_maker_bdf_fctry(this%makebdf)

    ! Setup backend dependent contributions of the OIFS scheme
    call rhs_maker_oifs_fctry(this%makeoifs)

    ! Initialize variables specific to this plan
    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call this%s_res%init(dm_Xh, "s_res")

      call this%abx1%init(dm_Xh, trim(this%name)//"_abx1")
      call neko_field_registry%add_field(dm_Xh, trim(this%name)//"_abx1", ignore_existing = .true.)

      call this%abx2%init(dm_Xh, trim(this%name)//"_abx2")
      call neko_field_registry%add_field(dm_Xh, trim(this%name)//"_abx2", ignore_existing = .true.)

      call this%advs%init(dm_Xh, "advs")

      call this%ds%init(dm_Xh, 'ds')

    end associate

    ! Set up boundary conditions
    call this%setup_bcs_(user)

    ! Initialize dirichlet bcs for scalar residual
    call this%bc_res%init(this%c_Xh, params)
    do i = 1, this%bcs%size()
       if (this%bcs%strong(i)) then
          bc_i => this%bcs%get(i)
          call this%bc_res%mark_facets(bc_i%marked_facet)
       end if
    end do

!    call this%bc_res%mark_zones_from_list('d_s', this%bc_labels)
    call this%bc_res%finalize()

    call this%bclst_ds%init()
    call this%bclst_ds%append(this%bc_res)


    ! Initialize projection space
    call this%proj_s%init(this%dm_Xh%size(), this%projection_dim, &
         this%projection_activ_step)

    ! Determine the time-interpolation scheme
    call json_get_or_default(params, 'case.numerics.oifs', this%oifs, .false.)
    ! Point to case checkpoint
    this%chkp => chkp
    ! Initialize advection factory
    call json_get_or_default(params, 'advection', advection, .true.)

    call advection_factory(this%adv, numerics_params, this%c_Xh, &
         ulag, vlag, wlag, this%chkp%dtlag, &
         this%chkp%tlag, time_scheme, .not. advection, &
         this%slag)
  end subroutine scalar_pnpn_init

  ! Restarts the scalar from a checkpoint
  subroutine scalar_pnpn_restart(this, chkp)
    class(scalar_pnpn_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    real(kind=rp) :: dtlag(10), tlag(10)
    integer :: n
    type(field_t), pointer :: temp_field
    dtlag = chkp%dtlag
    tlag = chkp%tlag

    n = this%s%dof%size()

    ! Lag fields are restored through the checkpoint's fsp mechanism

    call col2(this%s%x, this%c_Xh%mult, n)
    call col2(this%slag%lf(1)%x, this%c_Xh%mult, n)
    call col2(this%slag%lf(2)%x, this%c_Xh%mult, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%s%x, this%s%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%slag%lf(1)%x, this%slag%lf(1)%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%slag%lf(2)%x, this%slag%lf(2)%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%abx1%x, this%abx1%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%abx2%x, this%abx2%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%advs%x, this%advs%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
    end if

    call this%gs_Xh%op(this%s, GS_OP_ADD)
    call this%gs_Xh%op(this%slag%lf(1), GS_OP_ADD)
    call this%gs_Xh%op(this%slag%lf(2), GS_OP_ADD)

  end subroutine scalar_pnpn_restart

  subroutine scalar_pnpn_free(this)
    class(scalar_pnpn_t), intent(inout) :: this

    !Deallocate scalar field
    call this%scheme_free()

    call this%bc_res%free()
    call this%bclst_ds%free()
    call this%proj_s%free()

    call this%s_res%free()

    call this%ds%free()

    call this%abx1%free()
    call this%abx2%free()

    call this%advs%free()

    if (allocated(this%adv)) then
       call this%adv%free()
       deallocate(this%adv)
    end if

    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

    if (allocated(this%res)) then
       deallocate(this%res)
    end if

    if (allocated(this%makeext)) then
       deallocate(this%makeext)
    end if

    if (allocated(this%makebdf)) then
       deallocate(this%makebdf)
    end if

    if (allocated(this%makeoifs)) then
       deallocate(this%makeoifs)
    end if

  end subroutine scalar_pnpn_free

  subroutine scalar_pnpn_step(this, time, ext_bdf, dt_controller, &
       ksp_results)
    class(scalar_pnpn_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(time_step_controller_t), intent(in) :: dt_controller
    type(ksp_monitor_t), intent(inout) :: ksp_results
    ! Number of degrees of freedom
    integer :: n

    n = this%dm_Xh%size()

    call profiler_start_region(trim(this%name), 2)
    associate(u => this%u, v => this%v, w => this%w, s => this%s, &
         cp => this%cp, rho => this%rho, lambda_tot => this%lambda_tot, &
         ds => this%ds, &
         s_res => this%s_res, &
         Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         slag => this%slag, oifs => this%oifs, &
         projection_dim => this%projection_dim, &
         msh => this%msh, res => this%res, makeoifs => this%makeoifs, &
         makeext => this%makeext, makebdf => this%makebdf, &
         t => time%t, tstep => time%tstep, dt => time%dt)

      ! Logs extra information the log level is NEKO_LOG_DEBUG or above.
      call print_debug(this)
      ! Compute the source terms
      call this%source_term%compute(time)

      ! Apply weak boundary conditions, that contribute to the source terms.
      call this%bcs%apply_scalar(this%f_Xh%x, dm_Xh%size(), time, .false.)

      if (oifs) then
         ! Add the advection operators to the right-hans-side.
         call this%adv%compute_scalar(u, v, w, s, this%advs, &
              Xh, this%c_Xh, dm_Xh%size())

         call makeext%compute_scalar(this%abx1, this%abx2, f_Xh%x, &
              rho%x(1,1,1,1), ext_bdf%advection_coeffs, n)

         call makeoifs%compute_scalar(this%advs%x, f_Xh%x, rho%x(1,1,1,1), dt,&
              n)
      else
         ! Add the advection operators to the right-hans-side.
         call this%adv%compute_scalar(u, v, w, s, f_Xh, &
              Xh, this%c_Xh, dm_Xh%size())

         ! At this point the RHS contains the sum of the advection operator,
         ! Neumann boundary sources and additional source terms, evaluated using
         ! the scalar field from the previous time-step. Now, this value is used in
         ! the explicit time scheme to advance these terms in time.
         call makeext%compute_scalar(this%abx1, this%abx2, f_Xh%x, &
              rho%x(1,1,1,1), ext_bdf%advection_coeffs, n)

         ! Add the RHS contributions coming from the BDF scheme.
         call makebdf%compute_scalar(slag, f_Xh%x, s, c_Xh%B, rho%x(1,1,1,1), &
              dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
      end if

      call slag%update()

      !> Apply strong boundary conditions.
      call this%bcs%apply_scalar(this%s%x, this%dm_Xh%size(), time, .true.)

      ! Update material properties if necessary
      call this%update_material_properties(time)

      ! Compute scalar residual.
      call profiler_start_region(trim(this%name) // '_residual', 20)
      call res%compute(Ax, s, s_res, f_Xh, c_Xh, msh, Xh, lambda_tot, &
           rho%x(1,1,1,1)*cp%x(1,1,1,1), ext_bdf%diffusion_coeffs(1), dt, &
           dm_Xh%size())

      call gs_Xh%op(s_res, GS_OP_ADD)

      ! Apply a 0-valued Dirichlet boundary conditions on the ds.
      call this%bclst_ds%apply_scalar(s_res%x, dm_Xh%size())

      call profiler_end_region(trim(this%name) // '_residual', 20)

      call this%proj_s%pre_solving(s_res%x, tstep, c_Xh, n, dt_controller)

      call this%pc%update()
      call profiler_start_region(trim(this%name) // '_solve', 21)
      ksp_results = this%ksp%solve(Ax, ds, s_res%x, n, &
           c_Xh, this%bclst_ds, gs_Xh)
      ksp_results%name = trim(this%name)
      call profiler_end_region(trim(this%name) // '_solve', 21)

      call this%proj_s%post_solving(ds%x, Ax, c_Xh, this%bclst_ds, gs_Xh, &
           n, tstep, dt_controller)

      ! Update the solution
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(s%x_d, ds%x_d, 1.0_rp, n)
      else
         call add2s2(s%x, ds%x, 1.0_rp, n)
      end if

    end associate
    call profiler_end_region(trim(this%name), 2)
  end subroutine scalar_pnpn_step

  subroutine print_debug(this)
    class(scalar_pnpn_t), intent(inout) :: this
    character(len=LOG_SIZE) :: log_buf
    integer :: n

    n = this%dm_Xh%size()

    write(log_buf, '(A, A, E15.7, A, E15.7, A, E15.7)') 'Scalar debug', &
         ' l2norm s', glsc2(this%s%x, this%s%x, n), &
         ' slag1', glsc2(this%slag%lf(1)%x, this%slag%lf(1)%x, n), &
         ' slag2', glsc2(this%slag%lf(2)%x, this%slag%lf(2)%x, n)
    call neko_log%message(log_buf, lvl = NEKO_LOG_DEBUG)
    write(log_buf, '(A, A, E15.7, A, E15.7)') 'Scalar debug2', &
         ' l2norm abx1', glsc2(this%abx1%x, this%abx1%x, n), &
         ' abx2', glsc2(this%abx2%x, this%abx2%x, n)
    call neko_log%message(log_buf, lvl = NEKO_LOG_DEBUG)
  end subroutine print_debug

  !> Initialize boundary conditions
  !! @param user The user object binding the user-defined routines.
  subroutine scalar_pnpn_setup_bcs_(this, user)
    class(scalar_pnpn_t), intent(inout) :: this
    type(user_t), target, intent(in) :: user
    integer :: i, j, n_bcs, zone_size, global_zone_size, ierr
    type(json_core) :: core
    type(json_value), pointer :: bc_object
    type(json_file) :: bc_subdict
    class(bc_t), pointer :: bc_i
    logical :: found
    ! Monitor which boundary zones have been marked
    logical, allocatable :: marked_zones(:)
    integer, allocatable :: zone_indices(:)

    if (this%params%valid_path('boundary_conditions')) then
       call this%params%info('boundary_conditions', &
            n_children = n_bcs)
       call this%params%get_core(core)
       call this%params%get('boundary_conditions', bc_object, found)

       call this%bcs%init(n_bcs)

       allocate(marked_zones(size(this%msh%labeled_zones)))
       marked_zones = .false.

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)

          ! Check that we are not trying to assing a bc to zone, for which one
          ! has already been assigned and that the zone has more than 0 size
          ! in the mesh.
          call json_get(bc_subdict, "zone_indices", zone_indices)

          do j = 1, size(zone_indices)
             zone_size = this%msh%labeled_zones(zone_indices(j))%size
             call MPI_Allreduce(zone_size, global_zone_size, 1, &
                  MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

             if (global_zone_size .eq. 0) then
                write(error_unit, '(A, A, I0, A, A, I0, A)') "*** ERROR ***: ",&
                     "Zone index ", zone_indices(j), &
                     " is invalid as this zone has 0 size, meaning it ", &
                     "does not exist in the mesh. Check scalar boundary condition ", &
                     i, "."
                error stop
             end if

             if (marked_zones(zone_indices(j)) .eqv. .true.) then
                write(error_unit, '(A, A, I0, A, A, A, A)') "*** ERROR ***: ", &
                     "Zone with index ", zone_indices(j), &
                     " has already been assigned a boundary condition. ", &
                     "Please check your boundary_conditions entry for the ", &
                     "scalar and make sure that each zone index appears only ",&
                     "in a single boundary condition."
                error stop
             else
                marked_zones(zone_indices(j)) = .true.
             end if
          end do

          bc_i => null()

          call bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)
          call this%bcs%append(bc_i)
       end do

       ! Make sure all labeled zones with non-zero size have been marked
       do i = 1, size(this%msh%labeled_zones)
          if ((this%msh%labeled_zones(i)%size .gt. 0) .and. &
               (marked_zones(i) .eqv. .false.)) then
             write(error_unit, '(A, A, I0)') "*** ERROR ***: ", &
                  "No scalar boundary condition assigned to zone ", i
             error stop
          end if
       end do
    else
       ! Check that there are no labeled zones, i.e. all are periodic.
       do i = 1, size(this%msh%labeled_zones)
          if (this%msh%labeled_zones(i)%size .gt. 0) then
             write(error_unit, '(A, A, A)') "*** ERROR ***: ", &
                  "No boundary_conditions entry in the case file for scalar ", &
                  this%s%name
             error stop
          end if
       end do

       ! For a pure periodic case, we still need to initilise the bc lists
       ! to a zero size to avoid issues with apply() in step()
       call this%bcs%init()

    end if
  end subroutine scalar_pnpn_setup_bcs_


end module scalar_pnpn
