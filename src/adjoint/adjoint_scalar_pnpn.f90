!> @file adjoint_scalar_pnpn.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
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
!> Contains the `adjoint_scalar_pnpn_t` type.

module adjoint_scalar_pnpn
  use comm, only: NEKO_COMM
  use utils, only: neko_error
  use num_types, only: rp
  use, intrinsic :: iso_fortran_env, only: error_unit
  use rhs_maker, only : rhs_maker_bdf_t, rhs_maker_ext_t, rhs_maker_oifs_t, &
       rhs_maker_ext_fctry, rhs_maker_bdf_fctry, rhs_maker_oifs_fctry
  use adjoint_scalar_scheme, only : adjoint_scalar_scheme_t
  use checkpoint, only : chkp_t
  use field, only : field_t
  use scalar_bc_projector, only : scalar_bc_projector_t
  use mesh, only : mesh_t
  use coefs, only : coef_t
  use device, only : HOST_TO_DEVICE, device_memcpy
  use gather_scatter, only : gs_t, GS_OP_ADD
  use scalar_residual, only : scalar_residual_t, scalar_residual_factory
  use ax_product, only : ax_t, ax_helm_allocator
  use field_series, only: field_series_t
  use facet_normal, only : facet_normal_t
  use krylov, only : ksp_monitor_t
  use device_math, only : device_add2s2, device_col2
  use time_scheme_controller, only : time_scheme_controller_t
  use projection, only : projection_t
  use math, only : glsc2, col2, add2s2
  use field_math, only : field_col3, field_col2
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use advection_adjoint, only : advection_adjoint_t, advection_adjoint_factory
  use profiler, only : profiler_start_region, profiler_end_region
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use json_module, only : json_file, json_core, json_value
  use user_intf, only : user_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use zero_dirichlet, only : zero_dirichlet_t
  use time_step_controller, only : time_step_controller_t
  use scratch_registry, only : neko_scratch_registry
  use time_state, only : time_state_t
  use bc, only : bc_t, BC_DIRICHLET
  use mpi_f08, only: MPI_INTEGER, MPI_SUM, MPI_MAX
  implicit none
  private


  type, public, extends(adjoint_scalar_scheme_t) :: adjoint_scalar_pnpn_t

     !> The residual of the transport equation.
     type(field_t) :: s_adj_res

     !> Solution increment.
     type(field_t) :: ds_adj

     !> Helmholz operator.
     class(ax_t), allocatable :: Ax

     !> Solution projection.
     type(projection_t) :: proj_s

     !> Dirichlet conditions for the residual
     !! Collects all the Dirichlet condition facets into one bc and applies 0,
     !! Since the values never change there during the solve.
     type(zero_dirichlet_t) :: bc_res

     !> Projector for the adjoint scalar increment constraints.
     type(scalar_bc_projector_t) :: bc_projector

     !> Advection operator.
     class(advection_adjoint_t), allocatable :: adv

     ! Time interpolation scheme
     logical :: oifs

     ! Advection terms for the oifs method
     type(field_t) :: advs

     !> Computes the residual of the equation, i.e. `s_adj_res`.
     class(scalar_residual_t), allocatable :: res

     !> Contributions to kth order extrapolation scheme.
     class(rhs_maker_ext_t), allocatable :: makeext

     !> Contributions to the RHS from lagged BDF terms.
     class(rhs_maker_bdf_t), allocatable :: makebdf

     !> Contributions to the RHS from the OIFS method.
     class(rhs_maker_oifs_t), allocatable :: makeoifs

   contains
     !> Constructor.
     procedure, pass(this) :: init => adjoint_scalar_pnpn_init
     !> To restart
     procedure, pass(this) :: restart => adjoint_scalar_pnpn_restart
     !> Destructor.
     procedure, pass(this) :: free => adjoint_scalar_pnpn_free
     !> Solve for the current timestep.
     procedure, pass(this) :: step => adjoint_scalar_pnpn_step
     !> Setup the boundary conditions
     procedure, pass(this) :: setup_bcs_ => adjoint_scalar_pnpn_setup_bcs_
  end type adjoint_scalar_pnpn_t

  !> Boundary condition factory. Both constructs and initializes the object.
  !! @details Will mark a mesh zone for the bc and finalize.
  !! @param[inout] object The object to be allocated.
  !! @param[in] scheme The `adjoint_scalar_pnpn` scheme.
  !! @param[inout] json JSON object for initializing the bc.
  !! @param[in] coef SEM coefficients.
  interface adjoint_bc_factory
     module subroutine adjoint_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(adjoint_scalar_pnpn_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), target, intent(in) :: coef
       type(user_t), target, intent(in) :: user
     end subroutine adjoint_bc_factory
  end interface adjoint_bc_factory

contains

  !> Constructor.
  !! @details initialize the scheme.
  !! @param[inout] this The object.
  !! @param[in] msh The mesh.
  !! @param[in] coef The coefficients of the mesh.
  !! @param[in] gs The gather-scatter.
  !! @param[inout] params_adjoint The snippet of the json for the adjoint scalar.
  !! @param[inout] params_primal The snippet of the json for the primal scalar.
  !! @param[inout] numerics_params The numeric parameters json.
  !! @param[in] user Type with user-defined procedures.
  !! @param[inout] chkp The checkpoint object.
  !! @param[in] ulag Lag arrays for the x velocity component.
  !! @param[in] vlag Lag arrays for the y velocity component.
  !! @param[in] wlag Lag arrays for the z velocity component.
  !! @param[in] time_scheme The time-integration controller.
  !! @param[in] rho The fluid density.
  subroutine adjoint_scalar_pnpn_init(this, msh, coef, gs, params_adjoint, &
       params_primal, numerics_params, user, chkp, ulag, vlag, wlag, &
       time_scheme, rho)
    class(adjoint_scalar_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params_adjoint
    type(json_file), target, intent(inout) :: params_primal
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
    call this%scheme_init(msh, coef, gs, params_adjoint, params_primal, &
         scheme, user, rho)

    ! Setup backend dependent Ax routines
    call ax_helm_allocator(this%ax, type_name = "standard")

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

      call this%s_adj_res%init(dm_Xh, "s_adj_res")

      call this%abx1%init(dm_Xh, "abx1")

      call this%abx2%init(dm_Xh, "abx2")

      call this%advs%init(dm_Xh, "advs")

      call this%ds_adj%init(dm_Xh, 'ds_adj')

    end associate

    ! Set up boundary conditions
    call this%setup_bcs_(user)

    ! Initialize dirichlet bcs for scalar residual
    call this%bc_res%init(this%c_Xh, params_adjoint)
    do i = 1, this%bcs%size()
       if (this%bcs%bc_type(i) .eq. BC_DIRICHLET) then
          bc_i => this%bcs%get(i)
          call this%bc_res%mark_facets(bc_i%marked_facet)
       end if
    end do

!    call this%bc_res%mark_zones_from_list('d_s', this%bc_labels)
    call this%bc_res%finalize()
    call this%bc_projector%mark(this%bc_res)


    ! Initialize projection space
    call this%proj_s%init(this%dm_Xh%size(), this%projection_dim, &
         this%projection_activ_step)

    ! Determine the time-interpolation scheme
    call json_get_or_default(numerics_params, 'oifs', this%oifs, .false.)

    ! Point to case checkpoint
    this%chkp => chkp

    ! Initialize advection factory
    call json_get_or_default(params_adjoint, 'advection', advection, .true.)

    call advection_adjoint_factory(this%adv, numerics_params, this%c_Xh, &
         ulag, vlag, wlag, this%chkp%dtlag, &
         this%chkp%tlag, time_scheme, .not. advection, &
         this%s_adj_lag)
    ! Add lagged term to checkpoint
    ! @todo Init chkp object, note, adding 3 slags

    ! Add scalar info to checkpoint
    ! call this%chkp%add_scalar(this%s)
    ! this%chkp%abs1 => this%abx1
    ! this%chkp%abs2 => this%abx2
    ! this%chkp%slag => this%slag

  end subroutine adjoint_scalar_pnpn_init

  !> I envision the arguments to this func might need to be expanded
  subroutine adjoint_scalar_pnpn_restart(this, chkp)
    class(adjoint_scalar_pnpn_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    real(kind=rp) :: dtlag(10), tlag(10)
    integer :: n
    dtlag = chkp%dtlag
    tlag = chkp%tlag

    n = this%s_adj%dof%size()

    call col2(this%s_adj%x, this%c_Xh%mult, n)
    call col2(this%s_adj_lag%lf(1)%x, this%c_Xh%mult, n)
    call col2(this%s_adj_lag%lf(2)%x, this%c_Xh%mult, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%s_adj%x, this%s_adj%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%s_adj_lag%lf(1)%x, this%s_adj_lag%lf(1)%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%s_adj_lag%lf(2)%x, this%s_adj_lag%lf(2)%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%abx1%x, this%abx1%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%abx2%x, this%abx2%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%advs%x, this%advs%x_d, &
            n, HOST_TO_DEVICE, sync = .false.)
    end if

    call this%gs_Xh%op(this%s_adj, GS_OP_ADD)
    call this%gs_Xh%op(this%s_adj_lag%lf(1), GS_OP_ADD)
    call this%gs_Xh%op(this%s_adj_lag%lf(2), GS_OP_ADD)

  end subroutine adjoint_scalar_pnpn_restart

  subroutine adjoint_scalar_pnpn_free(this)
    class(adjoint_scalar_pnpn_t), intent(inout) :: this

    !Deallocate scalar field
    call this%scheme_free()

    call this%bc_projector%free()
    call this%bc_res%free()
    call this%proj_s%free()

    call this%s_adj_res%free()

    call this%ds_adj%free()

    call this%abx1%free()
    call this%abx2%free()

    call this%advs%free()

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

  end subroutine adjoint_scalar_pnpn_free

  subroutine adjoint_scalar_pnpn_step(this, time, ext_bdf, dt_controller, &
       ksp_results)
    class(adjoint_scalar_pnpn_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(time_step_controller_t), intent(in) :: dt_controller
    type(ksp_monitor_t), intent(inout) :: ksp_results
    type(field_t), pointer :: rho_cp
    integer :: rho_cp_index
    ! Number of degrees of freedom
    integer :: n

    if (this%freeze) return

    n = this%dm_Xh%size()
    call neko_scratch_registry%request_field(rho_cp, rho_cp_index, .false.)

    call profiler_start_region('Adjoint Scalar')
    associate(u => this%u, v => this%v, w => this%w, s_adj => this%s_adj, &
         cp => this%cp, rho => this%rho, lambda => this%lambda, &
         ds_adj => this%ds_adj, &
         s_adj_res => this%s_adj_res, &
         Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         s_adj_lag => this%s_adj_lag, oifs => this%oifs, &
         projection_dim => this%projection_dim, &
         msh => this%msh, res => this%res, makeoifs => this%makeoifs, &
         makeext => this%makeext, makebdf => this%makebdf, &
         t => time%t, tstep => time%tstep, dt => time%dt)

      ! Logs extra information the log level is NEKO_LOG_DEBUG or above.
      call print_debug(this)

      ! Update material properties and their pointwise product.
      ! This MUST happen before rho_cp is used below (by makebdf and by the
      ! source/advection scaling). It used to sit after those uses, which left
      ! makebdf reading an uninitialised scratch field while res%compute used
      ! the correct value -- so the BDF mass terms could not cancel at steady
      ! state and the converged adjoint scalar came out proportional to dt.
      ! Mirrors the ordering in Neko's forward scalar_pnpn.
      call this%update_material_properties(time)
      call field_col3(rho_cp, rho, cp)

      ! Compute the source terms
      call this%source_term%compute(time)

      ! if (oifs) then
      !    call neko_error("oifs not implemented for adjoint scalar")
      ! ! Add the advection operators to the right-hans-side.
      ! call this%adv%compute_scalar(u, v, w, s_adj, this%advs, &
      !                           Xh, this%c_Xh, dm_Xh%size())

      ! call makeext%compute_scalar(this%abx1, this%abx2, f_Xh%x, rho, &
      !                             ext_bdf%advection_coeffs, n)

      ! call makeoifs%compute_scalar(this%advs%x, f_Xh%x, rho, dt, n)
      ! else
      ! Add the advection operators to the right-hans-side.
      call this%adv%compute_adjoint_scalar(u, v, w, s_adj, f_Xh, &
           Xh, this%c_Xh, dm_Xh%size())

      ! Scale the volumetric source and advection terms by rho * cp, then add
      ! the weak boundary fluxes without that scaling -- same split as the
      ! forward scalar. (The weak BC application was previously done before
      ! the advection term and never scaled, so it could not follow this
      ! convention.)
      call field_col2(f_Xh, rho_cp)

      ! Apply weak boundary conditions, that contribute to the source terms.
      call this%bcs%apply_scalar(this%f_Xh%x, dm_Xh%size(), time, .false.)

      ! At this point the RHS contains the sum of the advection operator,
      ! Neumann boundary sources and additional source terms, evaluated using
      ! the scalar field from the previous time-step.
      ! Now, this value is used in the explicit time scheme to advance these
      ! terms in time.
      call makeext%compute_scalar(this%abx1, this%abx2, f_Xh%x, &
           ext_bdf%advection_coeffs%x, n)

      ! Add the RHS contributions coming from the BDF scheme.
      call makebdf%compute_scalar(s_adj_lag, f_Xh%x, s_adj, c_Xh%B, &
           rho_cp, dt, ext_bdf%diffusion_coeffs%x, ext_bdf%ndiff, n)
      ! end if

      call s_adj_lag%update()

      !> Apply strong boundary conditions.
      call this%bcs%apply_scalar(this%s_adj%x, this%dm_Xh%size(), time, &
           .true.)

      ! Compute scalar residual.
      call profiler_start_region('Adjoint_scalar_residual')
      call res%compute(Ax, s_adj, s_adj_res, f_Xh, c_Xh, msh, Xh, &
           lambda, rho_cp, ext_bdf%diffusion_coeffs%x(1), &
           dt, dm_Xh%size())

      call gs_Xh%op(s_adj_res, GS_OP_ADD)


      ! Apply a 0-valued Dirichlet boundary conditions on the ds_adj.
      call this%bc_projector%apply(s_adj_res%x, dm_Xh%size())

      call profiler_end_region('Adjoint_scalar_residual')

      call this%proj_s%pre_solving(s_adj_res%x, tstep, c_Xh, n, dt_controller)

      call this%pc%update()
      call profiler_start_region('Adjoint_scalar_solve')
      ksp_results = this%ksp%solve(Ax, ds_adj, s_adj_res%x, n, &
           c_Xh, this%bc_projector, gs_Xh)
      call profiler_end_region('Adjoint_scalar_solve')

      ksp_results%name = 'Adjoint Scalar'

      call this%proj_s%post_solving(ds_adj%x, Ax, c_Xh, this%bc_projector, &
           gs_Xh, n, tstep, dt_controller)

      ! Update the solution
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(s_adj%x_d, ds_adj%x_d, 1.0_rp, n)
      else
         call add2s2(s_adj%x, ds_adj%x, 1.0_rp, n)
      end if

    end associate
    call neko_scratch_registry%relinquish_field(rho_cp_index)
    call profiler_end_region('Adjoint Scalar')
  end subroutine adjoint_scalar_pnpn_step

  subroutine print_debug(this)
    class(adjoint_scalar_pnpn_t), intent(inout) :: this
    ! character(len=LOG_SIZE) :: log_buf
    integer :: n

    n = this%dm_Xh%size()
    ! TODO come back to this
    !write(log_buf,'(A,A,E15.7,A,E15.7,A,E15.7)') 'Adjoint scalar debug', &
    !   ' l2norm s_adj', glsc2(this%s_adj%x, this%s_adj%x, n), &
    !   ' slag1', glsc2(this%s_adj_lag%lf(1)%x, this%s_adj_lag%lf(1)%x, n), &
    !   ' slag2', glsc2(this%s_adj_lag%lf(2)%x, this%s_adj_lag%lf(2)%x, n)
    !call neko_log%message(log_buf, lvl=NEKO_LOG_DEBUG)
    !write(log_buf,'(A,A,E15.7,A,E15.7)') 'Adjoint scalar debug2', &
    !   ' l2norm abx1', glsc2(this%abx1%x, this%abx1%x, n), &
    !   ' abx2', glsc2(this%abx2%x, this%abx2%x, n)
    !call neko_log%message(log_buf, lvl=NEKO_LOG_DEBUG)
  end subroutine print_debug

  !> Initialize boundary conditions
  !! @param[inout] this The this.
  !! @param user The user object binding the user-defined routines.
  subroutine adjoint_scalar_pnpn_setup_bcs_(this, user)
    class(adjoint_scalar_pnpn_t), target, intent(inout) :: this
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
                     "does not exist in the mesh. Check adjoint scalar BC ", &
                     i, "."
                error stop
             end if

             if (marked_zones(zone_indices(j)) .eqv. .true.) then
                write(error_unit, '(A, A, I0, A, A, A, A)') "*** ERROR ***: ", &
                     "Zone with index ", zone_indices(j), &
                     " has already been assigned a boundary condition. ", &
                     "Please check your boundary_conditions entry for the ", &
                     "adjoint scalar and make sure that each zone index ", &
                     "appears only in a single boundary condition."
                error stop
             else
                marked_zones(zone_indices(j)) = .true.
             end if
          end do

          bc_i => null()

          call adjoint_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)
          call this%bcs%append(bc_i)
       end do

       ! Make sure all labeled zones with non-zero size have been marked
       do i = 1, size(this%msh%labeled_zones)
          if ((this%msh%labeled_zones(i)%size .gt. 0) .and. &
               (marked_zones(i) .eqv. .false.)) then
             write(error_unit, '(A, A, I0)') "*** ERROR ***: ", &
                  "No adjoint scalar boundary condition assigned to zone ", i
             error stop
          end if
       end do
    else
       ! Check that there are no labeled zones, i.e. all are periodic.
       do i = 1, size(this%msh%labeled_zones)
          if (this%msh%labeled_zones(i)%size .gt. 0) then
             write(error_unit, '(A, A, A)') "*** ERROR ***: ", &
                  "No boundary_conditions entry in the case file for " // &
                  " adjoint scalar ", &
                  this%s%name
             error stop
          end if
       end do
    end if
  end subroutine adjoint_scalar_pnpn_setup_bcs_

end module adjoint_scalar_pnpn
