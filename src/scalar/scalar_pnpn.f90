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
!> Containts the scalar_pnpn_t type.

module scalar_pnpn
  use num_types, only: rp
  use rhs_maker, only : rhs_maker_bdf_t, rhs_maker_ext_t, rhs_maker_oifs_t, &
       rhs_maker_ext_fctry, rhs_maker_bdf_fctry, rhs_maker_oifs_fctry
  use scalar_scheme, only : scalar_scheme_t
  use dirichlet, only : dirichlet_t
  use neumann, only : neumann_t
  use field, only : field_t
  use bc, only : bc_list_t, bc_list_init, bc_list_free, bc_list_apply_scalar, &
                 bc_list_add
  use mesh, only : mesh_t
  use checkpoint, only : chkp_t
  use coefs, only : coef_t
  use device, only : HOST_TO_DEVICE, device_memcpy
  use gather_scatter, only : gs_t, GS_OP_ADD
  use scalar_residual, only : scalar_residual_t, scalar_residual_factory
  use ax_product, only : ax_t, ax_helm_factory
  use field_series, only: field_series_t
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
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use user_intf, only : user_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use time_step_controller, only : time_step_controller_t
  use scratch_registry, only: neko_scratch_registry
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
     type(dirichlet_t) :: bc_res

     !> A bc list for the bc_res. Contains only that, essentially just to wrap
     !! the if statement determining whether to apply on the device or CPU.
     type(bc_list_t) :: bclst_ds

     !> Advection operator.
     class(advection_t), allocatable :: adv

     ! Time interpolation scheme
     logical :: oifs

     ! Lag arrays for the RHS.
     type(field_t) :: abx1, abx2

     ! Advection terms for the oifs method
     type(field_t) :: advs

     !> Computes the residual.
     class(scalar_residual_t), allocatable :: res

     !> Contributions to kth order extrapolation scheme.
     class(rhs_maker_ext_t), allocatable :: makeext

     !> Contributions to the RHS from lagged BDF terms.
     class(rhs_maker_bdf_t), allocatable :: makebdf

     !> Contributions to the RHS from the OIFS method
     class(rhs_maker_oifs_t), allocatable :: makeoifs

   contains
     !> Constructor.
     procedure, pass(this) :: init => scalar_pnpn_init
     !> To restart
     procedure, pass(this) :: restart => scalar_pnpn_restart
     !> Destructor.
     procedure, pass(this) :: free => scalar_pnpn_free
     !> Solve for the current timestep.
     procedure, pass(this) :: step => scalar_pnpn_step
  end type scalar_pnpn_t

contains

  !> Constructor.
  !! @param msh The mesh.
  !! @param coef The coefficients.
  !! @param gs The gather-scatter.
  !! @param params The case parameter file in json.
  !! @param user Type with user-defined procedures.
  subroutine scalar_pnpn_init(this, msh, coef, gs, params, user, &
       ulag, vlag, wlag, time_scheme, rho)
    class(scalar_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    type(coef_t), target, intent(inout) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    real(kind=rp), intent(in) :: rho
    integer :: i
    character(len=15), parameter :: scheme = 'Modular (Pn/Pn)'

    call this%free()

    ! Initiliaze base type.
    call this%scheme_init(msh, coef, gs, params, scheme, user, rho)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%ax, full_formulation = .false.)

    ! Setup backend dependent scalar residual routines
    call scalar_residual_factory(this%res)

    ! Setup backend dependent summation of extrapolation scheme
    call rhs_maker_ext_fctry(this%makeext)

    ! Setup backend depenent contributions to F from lagged BD terms
    call rhs_maker_bdf_fctry(this%makebdf)

    ! Setup backend dependent contributions of the OIFS scheme
    call rhs_maker_oifs_fctry(this%makeoifs)

    ! Initialize variables specific to this plan
    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call this%s_res%init(dm_Xh, "s_res")

      call this%abx1%init(dm_Xh, "abx1")

      call this%abx2%init(dm_Xh, "abx2")

      call this%advs%init(dm_Xh, "advs")

      call this%ds%init(dm_Xh, 'ds')

    end associate

    ! Initialize dirichlet bcs for scalar residual
    ! todo: look that this works
    call this%bc_res%init_base(this%c_Xh)
    do i = 1, this%n_dir_bcs
       call this%bc_res%mark_facets(this%dir_bcs(i)%marked_facet)
    end do

    ! Check for user bcs
    if (this%user_bc%msk(0) .gt. 0) then
       call this%bc_res%mark_facets(this%user_bc%marked_facet)
    end if

    call this%bc_res%mark_zones_from_list(msh%labeled_zones, 'd_s', &
                                         this%bc_labels)
    call this%bc_res%finalize()
    call this%bc_res%set_g(0.0_rp)

    call bc_list_init(this%bclst_ds)
    call bc_list_add(this%bclst_ds, this%bc_res)


    ! Intialize projection space
    call this%proj_s%init(this%dm_Xh%size(), this%projection_dim,  &
                            this%projection_activ_step)

    ! Add lagged term to checkpoint
    ! @todo Init chkp object, note, adding 3 slags
    ! call this%chkp%add_lag(this%slag, this%slag, this%slag)

    ! Determine the time-interpolation scheme
    call json_get_or_default(params, 'case.numerics.oifs', this%oifs, .false.)

    ! Initialize advection factory
    call advection_factory(this%adv, params, this%c_Xh, &
                           ulag, vlag, wlag, this%chkp%dtlag, &
                           this%chkp%tlag, time_scheme, this%slag)
  end subroutine scalar_pnpn_init

  !> I envision the arguments to this func might need to be expanded
  subroutine scalar_pnpn_restart(this, dtlag, tlag)
    class(scalar_pnpn_t), target, intent(inout) :: this
    real(kind=rp) :: dtlag(10), tlag(10)
    integer :: n


    n = this%s%dof%size()

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

    call bc_list_free(this%bclst_ds)
    call this%proj_s%free()

    call this%s_res%free()

    call this%ds%free()

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

  end subroutine scalar_pnpn_free

  subroutine scalar_pnpn_step(this, t, tstep, dt, ext_bdf, dt_controller)
    class(scalar_pnpn_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: t
    integer, intent(inout) :: tstep
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    type(time_step_controller_t), intent(in) :: dt_controller
    ! Number of degrees of freedom
    integer :: n
    ! Linear solver results monitor
    type(ksp_monitor_t) :: ksp_results(1)
    character(len=LOG_SIZE) :: log_buf

    n = this%dm_Xh%size()

    call profiler_start_region('Scalar', 2)
    associate(u => this%u, v => this%v, w => this%w, s => this%s, &
         cp => this%cp, rho => this%rho, &
         ds => this%ds, &
         s_res => this%s_res, &
         Ax => this%Ax, f_Xh => this%f_Xh, Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         slag => this%slag, oifs => this%oifs, &
         lambda_field => this%lambda_field, &
         projection_dim => this%projection_dim, &
         msh => this%msh, res => this%res, makeoifs => this%makeoifs, &
         makeext => this%makeext, makebdf => this%makebdf, &
         if_variable_dt => dt_controller%if_variable_dt, &
         dt_last_change => dt_controller%dt_last_change)

      if (neko_log%level_ .ge. NEKO_LOG_DEBUG) then
         write(log_buf, '(A,A,E15.7,A,E15.7,A,E15.7)') 'Scalar debug', &
              ' l2norm s', glsc2(this%s%x, this%s%x, n), &
              ' slag1', glsc2(this%slag%lf(1)%x, this%slag%lf(1)%x, n), &
              ' slag2', glsc2(this%slag%lf(2)%x, this%slag%lf(2)%x, n)
         call neko_log%message(log_buf)
         write(log_buf, '(A,A,E15.7,A,E15.7)') 'Scalar debug2', &
              ' l2norm abx1', glsc2(this%abx1%x, this%abx1%x, n), &
              ' abx2', glsc2(this%abx2%x, this%abx2%x, n)
         call neko_log%message(log_buf)
      end if

      ! Compute the source terms
      call this%source_term%compute(t, tstep)

      ! Compute the grandient jump penalty term
      if (this%if_gradient_jump_penalty .eqv. .true.) then
         call this%gradient_jump_penalty%compute(u, v, w, s)
         call this%gradient_jump_penalty%perform(f_Xh)
      end if

      ! Apply Neumann boundary conditions
      call bc_list_apply_scalar(this%bclst_neumann, this%f_Xh%x, dm_Xh%size())

      if (oifs) then
         ! Add the advection operators to the right-hans-side.
         call this%adv%compute_scalar(u, v, w, s, this%advs, &
                                   Xh, this%c_Xh, dm_Xh%size())

         call makeext%compute_scalar(this%abx1, this%abx2, f_Xh%x, rho, &
                                     ext_bdf%advection_coeffs, n)

         call makeoifs%compute_scalar(this%advs%x, f_Xh%x, rho, dt, n)
      else
         ! Add the advection operators to the right-hans-side.
         call this%adv%compute_scalar(u, v, w, s, f_Xh, &
                                      Xh, this%c_Xh, dm_Xh%size())

         ! At this point the RHS contains the sum of the advection operator,
         ! Neumann boundary sources and additional source terms, evaluated using
         ! the scalar field from the previous time-step. Now, this value is used in
         ! the explicit time scheme to advance these terms in time.
         call makeext%compute_scalar(this%abx1, this%abx2, f_Xh%x, rho, &
                                     ext_bdf%advection_coeffs, n)

         ! Add the RHS contributions coming from the BDF scheme.
         call makebdf%compute_scalar(slag, f_Xh%x, s, c_Xh%B, rho, dt, &
              ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
      end if

      call slag%update()

      !> Apply Dirichlet boundary conditions
      !! We assume that no change of boundary conditions
      !! occurs between elements. i.e. we do not apply gsop here like in Nek5000
      call this%field_dir_bc%update(this%field_dir_bc%field_list, &
           this%field_dirichlet_bcs, this%c_Xh, t, tstep, "scalar")
      call bc_list_apply_scalar(this%bclst_dirichlet, &
           this%s%x, this%dm_Xh%size())


      ! Update material properties if necessary
      call this%update_material_properties()

      ! Compute scalar residual.
      call profiler_start_region('Scalar_residual', 20)
      call res%compute(Ax, s,  s_res, f_Xh, c_Xh, msh, Xh, lambda_field, &
           rho*cp, ext_bdf%diffusion_coeffs(1), dt, dm_Xh%size())

      call gs_Xh%op(s_res, GS_OP_ADD)

      ! Apply a 0-valued Dirichlet boundary conditions on the ds.
      call bc_list_apply_scalar(this%bclst_ds, s_res%x, dm_Xh%size())

      call profiler_end_region('Scalar_residual', 20)

      call this%proj_s%pre_solving(s_res%x, tstep, c_Xh, n, dt_controller)

      call this%pc%update()
      call profiler_start_region('Scalar_solve', 21)
      ksp_results(1) = this%ksp%solve(Ax, ds, s_res%x, n, &
           c_Xh, this%bclst_ds, gs_Xh)
      call profiler_end_region('Scalar_solve', 21)

     call this%proj_s%post_solving(ds%x, Ax, c_Xh, &
                                 this%bclst_ds, gs_Xh, n, tstep, dt_controller)

      ! Update the solution
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(s%x_d, ds%x_d, 1.0_rp, n)
      else
         call add2s2(s%x, ds%x, 1.0_rp, n)
      end if

      call scalar_step_info(tstep, t, dt, ksp_results)

    end associate
    call profiler_end_region('Scalar', 2)
  end subroutine scalar_pnpn_step


end module scalar_pnpn
