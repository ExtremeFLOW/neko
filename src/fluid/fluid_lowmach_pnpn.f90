! Copyright (c) 2026, The Neko Authors
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
!> Low-Mach Pn-Pn fluid scheme.
!!
!! Extends fluid_pnpn_t. Contributes:
!!   * EOS plumbing: reads P0, R_specific, temperature-field name from JSON
!!     and updates the density field rho = P0 / (R_specific * T) each step.
!!   * Owns lowmach_prs_res_t and lowmach_vel_res_t residual objects. When
!!     low-Mach mode is enabled these replace the inherited pnpn residuals
!!     inside step(). A Q_T (thermal divergence) field is allocated and
!!     passed to the residuals.
!!
!! Physics content as of this commit: Q_T is held at zero and the residual
!! bodies reproduce pnpn_res_cpu exactly, so enabling low-Mach still gives
!! bit-identical output while the infrastructure is in place to add the Q_T
!! source and variable-property stress terms in subsequent commits.
module fluid_lowmach_pnpn
  use num_types, only : rp
  use fluid_pnpn, only : fluid_pnpn_t
  use field, only : field_t
  use registry, only : neko_registry
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use logger, only : neko_log, LOG_SIZE
  use mesh, only : mesh_t
  use user_intf, only : user_t
  use checkpoint, only : chkp_t
  use utils, only : neko_error
  use krylov, only : ksp_monitor_t
  use bc, only : bc_t
  use non_normal, only : non_normal_t
  use file, only : file_t
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use profiler, only : profiler_start_region, profiler_end_region
  use field_math, only : field_add2
  use field_series, only : field_series_t
  use operators, only : ortho, rotate_cyc, opgrad, cdtp
  use opr_device, only : device_ortho
  use scratch_registry, only : neko_scratch_registry
  use device, only : device_event_sync, glb_cmd_event, device_memcpy, &
       HOST_TO_DEVICE
  use device_mathops, only : device_opadd2cm
  use mathops, only : opadd2cm
  use fluid_aux, only : fluid_step_info
  use gs_ops, only : GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE
  use lowmach_residual, only : lowmach_prs_res_t, lowmach_vel_res_t, &
       lowmach_prs_res_factory, lowmach_vel_res_factory
  implicit none
  private

  type, public, extends(fluid_pnpn_t) :: fluid_lowmach_pnpn_t
     !> Master switch. When false, step() runs as vanilla pnpn.
     logical :: low_mach_enabled = .false.
     !> Thermodynamic (background) pressure — constant for open domains.
     real(kind=rp) :: P0 = 1.0_rp
     !> Specific gas constant R/M.
     real(kind=rp) :: R_specific = 1.0_rp
     !> Thermal conductivity used in Q_T = div(k grad T) / (rho cp T).
     real(kind=rp) :: k_cond = 1.0_rp
     !> Specific heat at constant pressure used in the same expression.
     real(kind=rp) :: cp = 1.0_rp
     !> Name of the temperature field in the global registry (default "s").
     character(len=:), allocatable :: T_field_name
     !> Pointer to the temperature field, resolved lazily on first step.
     type(field_t), pointer :: T_ptr => null()
     !> Thermal divergence source field. Currently held at zero; subsequent
     !! commits will populate it from the energy equation.
     type(field_t) :: Q_T_field
     !> Low-Mach pressure residual operator.
     class(lowmach_prs_res_t), allocatable :: lm_prs_res
     !> Low-Mach velocity residual operator.
     class(lowmach_vel_res_t), allocatable :: lm_vel_res
   contains
     procedure, pass(this) :: init => fluid_lowmach_pnpn_init
     procedure, pass(this) :: step => fluid_lowmach_pnpn_step
  end type fluid_lowmach_pnpn_t

contains

  !> Initialise the scheme. Delegates to the parent Pn-Pn init for all the
  !! heavy lifting, then reads low-Mach parameters, allocates residuals and
  !! the Q_T field.
  subroutine fluid_lowmach_pnpn_init(this, msh, lx, params, user, chkp)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    character(len=LOG_SIZE) :: log_buf
    character(len=:), allocatable :: tname
    integer :: i, n

    ! Run the standard Pn-Pn init to set up mesh, dofmap, fields, BCs, solvers.
    call this%fluid_pnpn_t%init(msh, lx, params, user, chkp)

    ! Low-Mach parameters — all under case.fluid.low_mach, all optional.
    call json_get_or_default(params, 'case.fluid.low_mach.enabled', &
         this%low_mach_enabled, .false.)
    call json_get_or_default(params, 'case.fluid.low_mach.P0', &
         this%P0, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.R_specific', &
         this%R_specific, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.k_conductivity', &
         this%k_cond, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.cp', &
         this%cp, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.temperature_field', &
         tname, 's')
    this%T_field_name = tname

    if (this%low_mach_enabled .and. NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("fluid_lowmach_pnpn: device backend not yet implemented")
    end if

    ! Allocate Q_T field (zero-initialised) and the residual operators.
    call this%Q_T_field%init(this%dm_Xh, 'Q_T')
    n = this%dm_Xh%size()
    do i = 1, n
       this%Q_T_field%x(i,1,1,1) = 0.0_rp
    end do

    call lowmach_prs_res_factory(this%lm_prs_res)
    call lowmach_vel_res_factory(this%lm_vel_res)

    if (this%low_mach_enabled) then
       call neko_log%section('Low-Mach')
       write(log_buf, '(A,E15.7)') 'P0         : ', this%P0
       call neko_log%message(log_buf)
       write(log_buf, '(A,E15.7)') 'R_specific : ', this%R_specific
       call neko_log%message(log_buf)
       write(log_buf, '(A,E15.7)') 'k_cond     : ', this%k_cond
       call neko_log%message(log_buf)
       write(log_buf, '(A,E15.7)') 'cp         : ', this%cp
       call neko_log%message(log_buf)
       write(log_buf, '(A,A)')     'T field    : ', trim(this%T_field_name)
       call neko_log%message(log_buf)
       call neko_log%end_section()
    end if

  end subroutine fluid_lowmach_pnpn_init

  !> Update the density field from the temperature field using the ideal-gas
  !! EOS rho = P0 / (R_specific * T).
  subroutine lowmach_update_density(this)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    integer :: i, n

    if (.not. associated(this%T_ptr)) then
       if (neko_registry%field_exists(this%T_field_name)) then
          this%T_ptr => neko_registry%get_field(this%T_field_name)
       else
          return
       end if
    end if

    n = this%dm_Xh%size()
    do concurrent (i = 1:n)
       this%rho%x(i,1,1,1) = this%P0 / (this%R_specific &
            * this%T_ptr%x(i,1,1,1))
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%rho%x, this%rho%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine lowmach_update_density

  !> Populate Q_T from the current temperature field using
  !!   Q_T = div(k grad T) / (rho * cp * T).
  !! The algebra relies on the identity D_t T = div(k grad T) / (rho cp)
  !! from the energy equation (no heat release), so the advection term in
  !! the material derivative cancels and Q_T becomes a pure spatial
  !! expression.
  subroutine lowmach_update_Q_T(this)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(field_t), pointer :: gx, gy, gz, w1
    integer :: temp_indices(4)
    integer :: i, n

    if (.not. associated(this%T_ptr)) return

    n = this%dm_Xh%size()

    call neko_scratch_registry%request_field(gx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(gy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(gz, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(w1, temp_indices(4), .false.)

    ! Weak-form gradient of T.
    call opgrad(gx%x, gy%x, gz%x, this%T_ptr%x, this%c_Xh)

    ! Assemble into the strong-form gradient and multiply by k.
    call this%gs_Xh%op(gx, GS_OP_ADD)
    call this%gs_Xh%op(gy, GS_OP_ADD)
    call this%gs_Xh%op(gz, GS_OP_ADD)
    do concurrent (i = 1:n)
       gx%x(i,1,1,1) = gx%x(i,1,1,1) * this%c_Xh%Binv(i,1,1,1) * this%k_cond
       gy%x(i,1,1,1) = gy%x(i,1,1,1) * this%c_Xh%Binv(i,1,1,1) * this%k_cond
       gz%x(i,1,1,1) = gz%x(i,1,1,1) * this%c_Xh%Binv(i,1,1,1) * this%k_cond
    end do

    ! Weak divergence via cdtp, accumulating into Q_T.
    call cdtp(this%Q_T_field%x, gx%x, this%c_Xh%drdx, this%c_Xh%dsdx, &
         this%c_Xh%dtdx, this%c_Xh)
    call cdtp(w1%x, gy%x, this%c_Xh%drdy, this%c_Xh%dsdy, &
         this%c_Xh%dtdy, this%c_Xh)
    do concurrent (i = 1:n)
       this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) &
            + w1%x(i,1,1,1)
    end do
    call cdtp(w1%x, gz%x, this%c_Xh%drdz, this%c_Xh%dsdz, &
         this%c_Xh%dtdz, this%c_Xh)
    do concurrent (i = 1:n)
       this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) &
            + w1%x(i,1,1,1)
    end do

    ! Assemble, convert to strong form, and divide by (rho * cp * T).
    call this%gs_Xh%op(this%Q_T_field, GS_OP_ADD)
    do concurrent (i = 1:n)
       this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) &
            * this%c_Xh%Binv(i,1,1,1) &
            / (this%rho%x(i,1,1,1) * this%cp * this%T_ptr%x(i,1,1,1))
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lowmach_update_Q_T

  !> Variable-density counterpart of rhs_maker_ext_cpu. Applies the AB
  !! extrapolation of the momentum forcing history and multiplies the
  !! result by the spatially-varying rho field instead of a scalar.
  subroutine lm_rhs_ext_field_rho(fx_lag, fy_lag, fz_lag, &
       fx_laglag, fy_laglag, fz_laglag, fx, fy, fz, rho, ext_coeffs, n)
    type(field_t), intent(inout) :: fx_lag, fy_lag, fz_lag
    type(field_t), intent(inout) :: fx_laglag, fy_laglag, fz_laglag
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: ext_coeffs(4)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: fx(n), fy(n), fz(n)
    type(field_t), pointer :: temp1, temp2, temp3
    integer :: temp_indices(3)
    integer :: i

    call neko_scratch_registry%request_field(temp1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(temp2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(temp3, temp_indices(3), .false.)

    do concurrent (i = 1:n)
       temp1%x(i,1,1,1) = ext_coeffs(2) * fx_lag%x(i,1,1,1) + &
            ext_coeffs(3) * fx_laglag%x(i,1,1,1)
       temp2%x(i,1,1,1) = ext_coeffs(2) * fy_lag%x(i,1,1,1) + &
            ext_coeffs(3) * fy_laglag%x(i,1,1,1)
       temp3%x(i,1,1,1) = ext_coeffs(2) * fz_lag%x(i,1,1,1) + &
            ext_coeffs(3) * fz_laglag%x(i,1,1,1)
    end do

    do concurrent (i = 1:n)
       fx_laglag%x(i,1,1,1) = fx_lag%x(i,1,1,1)
       fy_laglag%x(i,1,1,1) = fy_lag%x(i,1,1,1)
       fz_laglag%x(i,1,1,1) = fz_lag%x(i,1,1,1)
       fx_lag%x(i,1,1,1) = fx(i)
       fy_lag%x(i,1,1,1) = fy(i)
       fz_lag%x(i,1,1,1) = fz(i)
    end do

    do concurrent (i = 1:n)
       fx(i) = (ext_coeffs(1) * fx(i) + temp1%x(i,1,1,1)) * rho%x(i,1,1,1)
       fy(i) = (ext_coeffs(1) * fy(i) + temp2%x(i,1,1,1)) * rho%x(i,1,1,1)
       fz(i) = (ext_coeffs(1) * fz(i) + temp3%x(i,1,1,1)) * rho%x(i,1,1,1)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lm_rhs_ext_field_rho

  !> Variable-density counterpart of rhs_maker_bdf_cpu. Assembles the BDF
  !! lagged-velocity contribution to the momentum RHS weighting each
  !! pointwise term by rho(x)/dt.
  subroutine lm_rhs_bdf_field_rho(ulag, vlag, wlag, bfx, bfy, bfz, &
       u, v, w, B, rho, dt, bd, nbd, n)
    integer, intent(in) :: n, nbd
    type(field_t), intent(in) :: u, v, w
    type(field_series_t), intent(in) :: ulag, vlag, wlag
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(inout) :: bfx(n), bfy(n), bfz(n)
    real(kind=rp), intent(in) :: B(n)
    real(kind=rp), intent(in) :: dt, bd(4)
    type(field_t), pointer :: tb1, tb2, tb3
    integer :: temp_indices(3)
    integer :: i, ilag

    call neko_scratch_registry%request_field(tb1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(tb2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(tb3, temp_indices(3), .false.)

    do concurrent (i = 1:n)
       tb1%x(i,1,1,1) = u%x(i,1,1,1) * B(i) * bd(2)
       tb2%x(i,1,1,1) = v%x(i,1,1,1) * B(i) * bd(2)
       tb3%x(i,1,1,1) = w%x(i,1,1,1) * B(i) * bd(2)
    end do

    do ilag = 2, nbd
       do concurrent (i = 1:n)
          tb1%x(i,1,1,1) = tb1%x(i,1,1,1) + &
               (ulag%lf(ilag-1)%x(i,1,1,1) * B(i) * bd(ilag+1))
          tb2%x(i,1,1,1) = tb2%x(i,1,1,1) + &
               (vlag%lf(ilag-1)%x(i,1,1,1) * B(i) * bd(ilag+1))
          tb3%x(i,1,1,1) = tb3%x(i,1,1,1) + &
               (wlag%lf(ilag-1)%x(i,1,1,1) * B(i) * bd(ilag+1))
       end do
    end do

    do concurrent (i = 1:n)
       bfx(i) = bfx(i) + tb1%x(i,1,1,1) * (rho%x(i,1,1,1) / dt)
       bfy(i) = bfy(i) + tb2%x(i,1,1,1) * (rho%x(i,1,1,1) / dt)
       bfz(i) = bfz(i) + tb3%x(i,1,1,1) * (rho%x(i,1,1,1) / dt)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lm_rhs_bdf_field_rho

  !> Variable-density counterpart of rhs_maker_oifs_cpu. OIFS branch of
  !! the momentum RHS assembly with rho as a spatial field.
  subroutine lm_rhs_oifs_field_rho(phi_x, phi_y, phi_z, &
       bf_x, bf_y, bf_z, rho, dt, n)
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: bf_x(n), bf_y(n), bf_z(n)
    real(kind=rp), intent(inout) :: phi_x(n), phi_y(n), phi_z(n)
    integer :: i

    do concurrent (i = 1:n)
       bf_x(i) = bf_x(i) + phi_x(i) * (rho%x(i,1,1,1) / dt)
       bf_y(i) = bf_y(i) + phi_y(i) * (rho%x(i,1,1,1) / dt)
       bf_z(i) = bf_z(i) + phi_z(i) * (rho%x(i,1,1,1) / dt)
    end do

  end subroutine lm_rhs_oifs_field_rho

  !> Advance the solution by one time step. When low_mach_enabled is false,
  !! runs the inherited Pn-Pn path bit-for-bit. When true, updates rho from
  !! the EOS and routes the pressure/velocity residuals through the new
  !! lowmach_residual operators (which currently still ignore Q_T and match
  !! pnpn_res_cpu exactly).
  subroutine fluid_lowmach_pnpn_step(this, time, dt_controller)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    integer :: n
    type(ksp_monitor_t) :: ksp_results(4)
    type(file_t) :: dump_file
    class(bc_t), pointer :: bc_i
    type(non_normal_t), pointer :: bc_j

    if (this%freeze) return

    n = this%dm_Xh%size()

    if (this%low_mach_enabled) then
       call lowmach_update_density(this)
       call lowmach_update_Q_T(this)
    end if

    call profiler_start_region('Fluid', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_res => this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax_vel => this%Ax_vel, Ax_prs => this%Ax_prs, &
         Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         msh => this%msh, prs_res => this%prs_res, &
         source_term => this%source_term, vel_res => this%vel_res, &
         sumab => this%sumab, makeoifs => this%makeoifs, &
         makeabf => this%makeabf, makebdf => this%makebdf, &
         vel_projection_dim => this%vel_projection_dim, &
         pr_projection_dim => this%pr_projection_dim, &
         oifs => this%oifs, &
         rho => this%rho, mu_tot => this%mu_tot, &
         f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
         t => time%t, tstep => time%tstep, dt => time%dt, &
         ext_bdf => this%ext_bdf, event => glb_cmd_event)

      call sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
           ulag, vlag, wlag, ext_bdf%advection_coeffs%x, ext_bdf%nadv)

      call this%source_term%compute(time)

      call this%bcs_vel%apply_vector(f_x%x, f_y%x, f_z%x, &
           this%dm_Xh%size(), time, strong = .false.)

      if (oifs) then
         call this%adv%compute(u, v, w, &
              this%advx, this%advy, this%advz, &
              Xh, this%c_Xh, dm_Xh%size(), dt)

         if (this%low_mach_enabled) then
            call lm_rhs_ext_field_rho(this%abx1, this%aby1, this%abz1, &
                 this%abx2, this%aby2, this%abz2, &
                 f_x%x, f_y%x, f_z%x, &
                 rho, ext_bdf%advection_coeffs%x, n)
            call lm_rhs_oifs_field_rho(this%advx%x, this%advy%x, this%advz%x, &
                 f_x%x, f_y%x, f_z%x, rho, dt, n)
         else
            call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
                 this%abx2, this%aby2, this%abz2, &
                 f_x%x, f_y%x, f_z%x, &
                 rho%x(1,1,1,1), ext_bdf%advection_coeffs%x, n)

            call makeoifs%compute_fluid(this%advx%x, this%advy%x, this%advz%x, &
                 f_x%x, f_y%x, f_z%x, &
                 rho%x(1,1,1,1), dt, n)
         end if
      else
         call this%adv%compute(u, v, w, &
              f_x, f_y, f_z, &
              Xh, this%c_Xh, dm_Xh%size())

         if (this%low_mach_enabled) then
            call lm_rhs_ext_field_rho(this%abx1, this%aby1, this%abz1, &
                 this%abx2, this%aby2, this%abz2, &
                 f_x%x, f_y%x, f_z%x, &
                 rho, ext_bdf%advection_coeffs%x, n)
            call lm_rhs_bdf_field_rho(ulag, vlag, wlag, &
                 f_x%x, f_y%x, f_z%x, u, v, w, c_Xh%B, rho, dt, &
                 ext_bdf%diffusion_coeffs%x, ext_bdf%ndiff, n)
         else
            call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
                 this%abx2, this%aby2, this%abz2, &
                 f_x%x, f_y%x, f_z%x, &
                 rho%x(1,1,1,1), ext_bdf%advection_coeffs%x, n)

            call makebdf%compute_fluid(ulag, vlag, wlag, f_x%x, f_y%x, f_z%x, &
                 u, v, w, c_Xh%B, rho%x(1,1,1,1), dt, &
                 ext_bdf%diffusion_coeffs%x, ext_bdf%ndiff, n)
         end if
      end if

      call ulag%update()
      call vlag%update()
      call wlag%update()

      call this%bc_apply_vel(time, strong = .true.)
      call this%bc_apply_prs(time)

      call this%update_material_properties(time)

      call profiler_start_region('Pressure_residual', 18)

      if (this%low_mach_enabled) then
         call this%lm_prs_res%compute(p, p_res, &
              u, v, w, &
              u_e, v_e, w_e, &
              f_x, f_y, f_z, &
              c_Xh, gs_Xh, &
              this%bc_prs_surface, this%bc_sym_surface, &
              Ax_prs, ext_bdf%diffusion_coeffs%x(1), dt, &
              mu_tot, rho, this%Q_T_field, event)
      else
         call prs_res%compute(p, p_res, &
              u, v, w, &
              u_e, v_e, w_e, &
              f_x, f_y, f_z, &
              c_Xh, gs_Xh, &
              this%bc_prs_surface, this%bc_sym_surface, &
              Ax_prs, ext_bdf%diffusion_coeffs%x(1), dt, &
              mu_tot, rho, event)
      end if

      if (.not. this%prs_dirichlet .and. NEKO_BCKND_DEVICE .eq. 1) then
         call device_ortho(p_res%x_d, this%glb_n_points, n)
      else if (.not. this%prs_dirichlet) then
         call ortho(p_res%x, this%glb_n_points, n)
      end if

      call gs_Xh%op(p_res, GS_OP_ADD, event)
      call device_event_sync(event)

      call this%bclst_dp%apply_scalar(p_res%x, p%dof%size(), time)

      call profiler_end_region('Pressure_residual', 18)

      call this%proj_prs%pre_solving(p_res%x, tstep, c_Xh, n, dt_controller, &
           Ax=Ax_prs, gs_h=gs_Xh, bclst=this%bclst_dp, string='Pressure')

      call this%pc_prs%update()

      call profiler_start_region('Pressure_solve', 3)

      ksp_results(1) = &
           this%ksp_prs%solve(Ax_prs, dp, p_res%x, n, c_Xh, &
           this%bclst_dp, gs_Xh)
      ksp_results(1)%name = 'Pressure'

      call profiler_end_region('Pressure_solve', 3)

      call this%proj_prs%post_solving(dp%x, Ax_prs, c_Xh, &
           this%bclst_dp, gs_Xh, n, tstep, dt_controller)

      call field_add2(p, dp, n)
      if (.not. this%prs_dirichlet .and. NEKO_BCKND_DEVICE .eq. 1) then
         call device_ortho(p%x_d, this%glb_n_points, n)
      else if (.not. this%prs_dirichlet) then
         call ortho(p%x, this%glb_n_points, n)
      end if

      call profiler_start_region('Velocity_residual', 19)
      if (this%low_mach_enabled) then
         call this%lm_vel_res%compute(Ax_vel, u, v, w, &
              u_res, v_res, w_res, &
              p, &
              f_x, f_y, f_z, &
              c_Xh, msh, Xh, &
              mu_tot, rho, this%Q_T_field, &
              ext_bdf%diffusion_coeffs%x(1), dt, dm_Xh%size())
      else
         call vel_res%compute(Ax_vel, u, v, w, &
              u_res, v_res, w_res, &
              p, &
              f_x, f_y, f_z, &
              c_Xh, msh, Xh, &
              mu_tot, rho, ext_bdf%diffusion_coeffs%x(1), &
              dt, dm_Xh%size())
      end if

      call rotate_cyc(u_res%x, v_res%x, w_res%x, 1, c_Xh)
      call gs_Xh%op(u_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(v_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call gs_Xh%op(w_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call rotate_cyc(u_res%x, v_res%x, w_res%x, 0, c_Xh)

      call this%bclst_vel_res%apply(u_res, v_res, w_res, time)

      call profiler_end_region('Velocity_residual', 19)

      call this%proj_vel%pre_solving(u_res%x, v_res%x, w_res%x, &
           tstep, c_Xh, n, dt_controller, 'Velocity')

      call this%pc_vel%update()

      call profiler_start_region("Velocity_solve", 4)
      ksp_results(2:4) = this%ksp_vel%solve_coupled(Ax_vel, du, dv, dw, &
           u_res%x, v_res%x, w_res%x, n, c_Xh, &
           this%bclst_du, this%bclst_dv, this%bclst_dw, gs_Xh, &
           this%ksp_vel%max_iter)
      call profiler_end_region("Velocity_solve", 4)
      if (this%full_stress_formulation) then
         ksp_results(2)%name = 'Momentum'
      else
         ksp_results(2)%name = 'X-Velocity'
         ksp_results(3)%name = 'Y-Velocity'
         ksp_results(4)%name = 'Z-Velocity'
      end if

      call this%proj_vel%post_solving(du%x, dv%x, dw%x, Ax_vel, c_Xh, &
           this%bclst_du, this%bclst_dv, this%bclst_dw, gs_Xh, n, tstep, &
           dt_controller)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_opadd2cm(u%x_d, v%x_d, w%x_d, &
              du%x_d, dv%x_d, dw%x_d, 1.0_rp, n, msh%gdim)
      else
         call opadd2cm(u%x, v%x, w%x, du%x, dv%x, dw%x, 1.0_rp, n, msh%gdim)
      end if

      if (this%forced_flow_rate) then
         call this%vol_flow%adjust( u, v, w, p, u_res, v_res, w_res, p_res, &
              c_Xh, gs_Xh, ext_bdf, rho%x(1,1,1,1), mu_tot, &
              dt, time, this%bclst_dp, this%bclst_du, this%bclst_dv, &
              this%bclst_dw, this%bclst_vel_res, Ax_vel, Ax_prs, this%ksp_prs, &
              this%ksp_vel, this%pc_prs, this%pc_vel, this%ksp_prs%max_iter, &
              this%ksp_vel%max_iter)
      end if

      call fluid_step_info(time, ksp_results, &
           this%full_stress_formulation, this%strict_convergence, &
           this%allow_stabilization)

    end associate
    call profiler_end_region('Fluid', 1)
  end subroutine fluid_lowmach_pnpn_step

end module fluid_lowmach_pnpn
