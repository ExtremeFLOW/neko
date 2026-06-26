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
!! Physics content: Q_T = div(k grad T)/(rho cp T) is computed from the
!! temperature field and (i) drives the continuity constraint div u = Q_T via
!! a source in the pressure residual, and (ii) appears in the momentum
!! residual through the dilatation stress -(2/3) grad(mu Q_T). The deviatoric
!! stress div[mu(grad u + grad u^T)] is handled by the coupled full-stress
!! operator, so low_mach requires full_stress_formulation = true.
!!
!! Both the pressure and velocity residuals carry the variable-property viscous
!! stress and the Q_T terms. Known follow-up: there is no device backend yet.
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
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use profiler, only : profiler_start_region, profiler_end_region
  use field_math, only : field_add2, field_copy
  use math, only : glsc2, glmin, glmax, rzero
  use, intrinsic :: iso_c_binding, only : c_ptr
  use operators, only : ortho, rotate_cyc, opgrad, cdtp, dudxyz
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
     !> Thermodynamic (background) pressure — constant for open domains.
     real(kind=rp) :: P0 = 1.0_rp
     !> Specific gas constant R/M.
     real(kind=rp) :: R_specific = 1.0_rp
     !> Temperature clamp for the EOS (guards rho against unphysical T).
     real(kind=rp) :: T_eos_min = tiny(1.0_rp)
     real(kind=rp) :: T_eos_max = huge(1.0_rp)
     !> Thermal conductivity used in Q_T = div(k grad T) / (rho cp T).
     real(kind=rp) :: k_cond = 1.0_rp
     !> Specific heat at constant pressure used in the same expression.
     real(kind=rp) :: cp = 1.0_rp
     !> Optional built-in temperature-dependent property model. One of
     !! 'none' (default; properties come from the user material_properties
     !! hook), 'sutherland', or 'power_law'. When set, mu(T) and lambda(T) are
     !! filled from the law each step, so a user hook is not required.
     !!   power_law:   q(T) = q_ref (T/T_ref)**exponent
     !!   sutherland:  q(T) = q_ref (T/T_ref)**1.5 (T_ref + S_q)/(T + S_q)
     character(len=:), allocatable :: property_model
     real(kind=rp) :: prop_mu_ref = 1.0_rp
     real(kind=rp) :: prop_lambda_ref = 1.0_rp
     real(kind=rp) :: prop_T_ref = 1.0_rp
     real(kind=rp) :: prop_S_mu = 0.0_rp
     real(kind=rp) :: prop_S_lambda = 0.0_rp
     real(kind=rp) :: prop_exponent = 0.7_rp
     !> Name of the temperature field in the global registry (default "s").
     character(len=:), allocatable :: T_field_name
     !> Pointer to the temperature field, resolved lazily on first step.
     type(field_t), pointer :: T_ptr => null()
     !> Pointers to the scalar's conductivity (lambda) and cp fields, resolved
     !! lazily. Using these (instead of the constant k_cond/cp) makes Q_T
     !! consistent with the energy equation the scalar actually solves
     !! (matches Nek5000 userqtl_scig, which uses vdiff/vtrans of the T field).
     type(field_t), pointer :: lambda_ptr => null()
     type(field_t), pointer :: cp_ptr => null()
     !> Thermal divergence source field, Q_T = div(k grad T)/(rho cp T),
     !! populated each step from the temperature field by lowmach_update_Q_T.
     type(field_t), pointer :: Q_T_field => null()
     !> Low-Mach pressure residual operator.
     class(lowmach_prs_res_t), allocatable :: lm_prs_res
     !> Low-Mach velocity residual operator.
     class(lowmach_vel_res_t), allocatable :: lm_vel_res
   contains
     procedure, pass(this) :: init => fluid_lowmach_pnpn_init
     procedure, pass(this) :: free => fluid_lowmach_pnpn_free
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

    ! Low-Mach is a variable-density scheme whose momentum stress is assembled
    ! with the coupled full-stress operator, which in turn requires a coupled
    ! velocity solver. Both are mandatory, so force them here rather than make
    ! the user set them (and risk a silently inconsistent configuration). They
    ! must be forced BEFORE the parent init, which builds Ax_vel from these keys.
    call json_force_logical(params, 'case.fluid.full_stress_formulation', .true.)
    call json_force_string(params, 'case.fluid.velocity_solver.type', 'coupled_cg')

    ! Run the standard Pn-Pn init to set up mesh, dofmap, fields, BCs, solvers.
    call this%fluid_pnpn_t%init(msh, lx, params, user, chkp)

    ! Low-Mach parameters — all under case.fluid.low_mach, all optional.
    call json_get_or_default(params, 'case.fluid.low_mach.P0', &
         this%P0, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.R_specific', &
         this%R_specific, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.T_eos_min', &
         this%T_eos_min, tiny(1.0_rp))
    call json_get_or_default(params, 'case.fluid.low_mach.T_eos_max', &
         this%T_eos_max, huge(1.0_rp))
    call json_get_or_default(params, 'case.fluid.low_mach.k_conductivity', &
         this%k_cond, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.cp', &
         this%cp, 1.0_rp)
    call json_get_or_default(params, 'case.fluid.low_mach.temperature_field', &
         tname, 's')
    this%T_field_name = tname

    ! Optional built-in property model (mu(T), lambda(T)). Default 'none' keeps
    ! the user material_properties hook as the source of properties.
    call json_get_or_default(params, &
         'case.fluid.low_mach.property_model.type', this%property_model, 'none')
    call json_get_or_default(params, &
         'case.fluid.low_mach.property_model.mu_ref', this%prop_mu_ref, 1.0_rp)
    call json_get_or_default(params, &
         'case.fluid.low_mach.property_model.lambda_ref', &
         this%prop_lambda_ref, 1.0_rp)
    call json_get_or_default(params, &
         'case.fluid.low_mach.property_model.T_ref', this%prop_T_ref, 1.0_rp)
    call json_get_or_default(params, &
         'case.fluid.low_mach.property_model.S_mu', this%prop_S_mu, 0.0_rp)
    call json_get_or_default(params, &
         'case.fluid.low_mach.property_model.S_lambda', &
         this%prop_S_lambda, 0.0_rp)
    call json_get_or_default(params, &
         'case.fluid.low_mach.property_model.exponent', &
         this%prop_exponent, 0.7_rp)

    if (this%property_model .ne. 'none' .and. &
         this%property_model .ne. 'sutherland' .and. &
         this%property_model .ne. 'power_law') then
       call neko_error("fluid_lowmach_pnpn: unknown property_model.type '" // &
            this%property_model // "' (use none, sutherland or power_law)")
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("fluid_lowmach_pnpn: device backend not yet implemented")
    end if

    ! Allocate Q_T field (zero-initialised) and the residual operators.
    ! Register it so it is readable by user hooks / the field writer (e.g. for
    ! the lowMach_test QTL validation).
    call neko_registry%add_field(this%dm_Xh, 'Q_T', ignore_existing = .true.)
    this%Q_T_field => neko_registry%get_field('Q_T')
    n = this%dm_Xh%size()
    do i = 1, n
       this%Q_T_field%x(i,1,1,1) = 0.0_rp
    end do

    call lowmach_prs_res_factory(this%lm_prs_res)
    call lowmach_vel_res_factory(this%lm_vel_res)

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
    write(log_buf, '(A,A)')     'Property   : ', trim(this%property_model)
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine fluid_lowmach_pnpn_init

  !> Force a logical JSON key to a value. json_file%update replaces the value
  !! if the key exists and creates it (with any missing parent objects) via
  !! add_by_path if it does not, so this overrides whatever the user set.
  subroutine json_force_logical(params, path, val)
    type(json_file), intent(inout) :: params
    character(len=*), intent(in) :: path
    logical, intent(in) :: val
    logical :: found
    call params%update(path, val, found)
  end subroutine json_force_logical

  !> Force a string JSON key to a value (see json_force_logical).
  subroutine json_force_string(params, path, val)
    type(json_file), intent(inout) :: params
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: val
    logical :: found
    call params%update(path, val, found)
  end subroutine json_force_string

  !> Free the low-Mach members (the Q_T field and the residual operators),
  !! then delegate to the parent for everything inherited from fluid_pnpn.
  subroutine fluid_lowmach_pnpn_free(this)
    class(fluid_lowmach_pnpn_t), intent(inout) :: this

    nullify(this%Q_T_field)   ! owned by neko_registry

    if (allocated(this%lm_prs_res)) deallocate(this%lm_prs_res)
    if (allocated(this%lm_vel_res)) deallocate(this%lm_vel_res)

    nullify(this%T_ptr)
    nullify(this%lambda_ptr)
    nullify(this%cp_ptr)

    call this%fluid_pnpn_t%free()

  end subroutine fluid_lowmach_pnpn_free

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

    ! Clamp the temperature used in the EOS to a physical range. Near-wall /
    ! inlet temperature discontinuities can produce Gibbs over/undershoots that
    ! would otherwise drive rho = P0/(R T) to infinity or negative values,
    ! destroying both the physics and the pressure-solve conditioning.
    n = this%dm_Xh%size()
    do concurrent (i = 1:n)
       this%rho%x(i,1,1,1) = this%P0 / (this%R_specific &
            * min(max(this%T_ptr%x(i,1,1,1), this%T_eos_min), this%T_eos_max))
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%rho%x, this%rho%x_d, n, &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine lowmach_update_density

  !> Fill mu(T) and lambda(T) from the built-in property model (Sutherland or
  !! power law) when one is configured, so the user does not have to supply a
  !! material_properties hook. Writes the molecular fields (this%mu and the
  !! scalar's "<T>_lambda") and their "_tot" counterparts so the values are
  !! immediately consistent for the Q_T assembly and the Helmholtz solves this
  !! step. A no-op (returns early) when property_model == 'none', preserving the
  !! user-hook path bit-for-bit.
  subroutine lowmach_update_properties(this)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(field_t), pointer :: lam, lam_tot
    integer :: i, n
    real(kind=rp) :: Tval, fac_mu, fac_lam

    if (this%property_model .eq. 'none') return

    if (.not. associated(this%T_ptr)) then
       if (neko_registry%field_exists(this%T_field_name)) then
          this%T_ptr => neko_registry%get_field(this%T_field_name)
       else
          return
       end if
    end if

    ! Resolve the scalar's base/total conductivity fields. They are absent if
    ! no scalar is configured, in which case only mu is updated.
    nullify(lam)
    nullify(lam_tot)
    if (neko_registry%field_exists(this%T_field_name // "_lambda")) &
         lam => neko_registry%get_field(this%T_field_name // "_lambda")
    if (neko_registry%field_exists(this%T_field_name // "_lambda_tot")) &
         lam_tot => neko_registry%get_field(this%T_field_name // "_lambda_tot")

    n = this%dm_Xh%size()
    do i = 1, n
       Tval = min(max(this%T_ptr%x(i,1,1,1), this%T_eos_min), this%T_eos_max)
       if (this%property_model .eq. 'sutherland') then
          fac_mu = (Tval / this%prop_T_ref)**1.5_rp &
               * (this%prop_T_ref + this%prop_S_mu) / (Tval + this%prop_S_mu)
          fac_lam = (Tval / this%prop_T_ref)**1.5_rp &
               * (this%prop_T_ref + this%prop_S_lambda) &
               / (Tval + this%prop_S_lambda)
       else  ! power_law
          fac_mu = (Tval / this%prop_T_ref)**this%prop_exponent
          fac_lam = fac_mu
       end if
       this%mu%x(i,1,1,1) = this%prop_mu_ref * fac_mu
       if (associated(lam)) lam%x(i,1,1,1) = this%prop_lambda_ref * fac_lam
    end do

    ! Sync the *_tot fields so Q_T (which reads lambda_tot) and the velocity /
    ! energy Helmholtz operators see the fresh values within this step.
    call field_copy(this%mu_tot, this%mu)
    if (associated(lam) .and. associated(lam_tot)) call field_copy(lam_tot, lam)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%mu%x, this%mu%x_d, n, HOST_TO_DEVICE, &
            sync = .false.)
    end if

  end subroutine lowmach_update_properties

  !> Populate Q_T = div(lambda grad T) / (rho cp T), the low-Mach thermal
  !! divergence (div u = Q_T) for an ideal gas with constant background
  !! pressure (no volumetric heat release, open domain so dP0/dt = 0).
  !!
  !! This mirrors Nek5000 userqtl_scig (plan4.f): grad T (weak) -> assemble ->
  !! 1/B (strong grad) -> multiply by the conductivity field lambda -> weak
  !! divergence -> assemble -> 1/B (strong) -> divide by (rho cp T). The
  !! crucial point is that lambda and cp here are the SAME spatially varying
  !! fields the scalar (energy) equation uses (Nek5000's vdiff/vtrans of the
  !! temperature field), so the imposed velocity divergence is consistent with
  !! the temperature the scalar actually transports. Using a constant k/cp
  !! instead makes div u = Q_T inconsistent with the energy balance and was the
  !! cause of the heated-channel blow-up.
  subroutine lowmach_update_Q_T(this)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(field_t), pointer :: gx, gy, gz, w1, qdot
    integer :: temp_indices(4)
    integer :: i, n
    logical :: have_lambda, have_cp

    if (.not. associated(this%T_ptr)) return

    ! Resolve the scalar's conductivity / cp fields lazily. These are created
    ! by the scalar scheme as "<temperature>_lambda(_tot)" and "<temperature>_cp"
    ! and filled (e.g. by a Sutherland material_properties hook) with lambda(T).
    if (.not. associated(this%lambda_ptr)) then
       if (neko_registry%field_exists(this%T_field_name // "_lambda_tot")) then
          this%lambda_ptr => &
               neko_registry%get_field(this%T_field_name // "_lambda_tot")
       else if (neko_registry%field_exists(this%T_field_name // "_lambda")) then
          this%lambda_ptr => &
               neko_registry%get_field(this%T_field_name // "_lambda")
       end if
    end if
    if (.not. associated(this%cp_ptr)) then
       if (neko_registry%field_exists(this%T_field_name // "_cp")) then
          this%cp_ptr => neko_registry%get_field(this%T_field_name // "_cp")
       end if
    end if
    have_lambda = associated(this%lambda_ptr)
    have_cp = associated(this%cp_ptr)

    n = this%dm_Xh%size()

    call neko_scratch_registry%request_field(gx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(gy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(gz, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(w1, temp_indices(4), .false.)

    ! Conduction div(lambda grad T) via dudxyz (strong collocation derivative),
    ! the same primitive Neko's div() operator uses. cdtp (D^T, weak/transpose)
    ! does NOT compose with Binv into div(lambda grad T) (it gave O(1e2) garbage
    ! even for T = 2 + x^2 whose Laplacian is the constant 2); and the weak
    ! ax_helm route, while consistent with the energy operator, drops the
    ! boundary-flux term and grows Q_T near the all-Dirichlet boundary. dudxyz is
    ! the accurate, bounded choice (Q_T err ~6e-3 on the manufactured tanh).
    call dudxyz(gx%x, this%T_ptr%x, this%c_Xh%drdx, this%c_Xh%dsdx, &
         this%c_Xh%dtdx, this%c_Xh)
    call dudxyz(gy%x, this%T_ptr%x, this%c_Xh%drdy, this%c_Xh%dsdy, &
         this%c_Xh%dtdy, this%c_Xh)
    call dudxyz(gz%x, this%T_ptr%x, this%c_Xh%drdz, this%c_Xh%dsdz, &
         this%c_Xh%dtdz, this%c_Xh)

    ! Make the gradient continuous (gs_op ADD then average via mult = 1/count)
    ! and multiply by the conductivity lambda(T) (else constant k_cond).
    call this%gs_Xh%op(gx, GS_OP_ADD)
    call this%gs_Xh%op(gy, GS_OP_ADD)
    call this%gs_Xh%op(gz, GS_OP_ADD)
    if (have_lambda) then
       do concurrent (i = 1:n)
          gx%x(i,1,1,1) = gx%x(i,1,1,1) * this%c_Xh%mult(i,1,1,1) &
               * this%lambda_ptr%x(i,1,1,1)
          gy%x(i,1,1,1) = gy%x(i,1,1,1) * this%c_Xh%mult(i,1,1,1) &
               * this%lambda_ptr%x(i,1,1,1)
          gz%x(i,1,1,1) = gz%x(i,1,1,1) * this%c_Xh%mult(i,1,1,1) &
               * this%lambda_ptr%x(i,1,1,1)
       end do
    else
       do concurrent (i = 1:n)
          gx%x(i,1,1,1) = gx%x(i,1,1,1) * this%c_Xh%mult(i,1,1,1) * this%k_cond
          gy%x(i,1,1,1) = gy%x(i,1,1,1) * this%c_Xh%mult(i,1,1,1) * this%k_cond
          gz%x(i,1,1,1) = gz%x(i,1,1,1) * this%c_Xh%mult(i,1,1,1) * this%k_cond
       end do
    end if

    ! Strong-form divergence div(lambda grad T) via dudxyz (matches Neko's div()).
    call dudxyz(this%Q_T_field%x, gx%x, this%c_Xh%drdx, this%c_Xh%dsdx, &
         this%c_Xh%dtdx, this%c_Xh)
    call dudxyz(w1%x, gy%x, this%c_Xh%drdy, this%c_Xh%dsdy, &
         this%c_Xh%dtdy, this%c_Xh)
    do concurrent (i = 1:n)
       this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) + w1%x(i,1,1,1)
    end do
    call dudxyz(w1%x, gz%x, this%c_Xh%drdz, this%c_Xh%dsdz, &
         this%c_Xh%dtdz, this%c_Xh)
    do concurrent (i = 1:n)
       this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) + w1%x(i,1,1,1)
    end do

    call this%gs_Xh%op(this%Q_T_field, GS_OP_ADD)
    do concurrent (i = 1:n)
       this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) &
            * this%c_Xh%mult(i,1,1,1)
    end do

    ! Add the volumetric heat source q (STRONG form, point values) so that
    !   Q_T_field := div(lambda grad T) + q,
    ! matching Nek5000 userqtl_scig, whose QTL uses the full energy-equation RHS
    ! (including the source). q is read from the OPTIONAL registry field
    ! "<temperature>_qdot" (created/filled by the user, e.g. a manufactured
    ! source or chemical heat release); when absent, behaviour is the previous
    ! conduction-only Q_T, so existing cases are bit-for-bit unchanged.
    if (neko_registry%field_exists(this%T_field_name // "_qdot")) then
       qdot => neko_registry%get_field(this%T_field_name // "_qdot")
       do concurrent (i = 1:n)
          this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) &
               + qdot%x(i,1,1,1)
       end do
    end if

    ! Divide by (rho cp T) to obtain Q_T. The temperature used here MUST be the
    ! SAME clamped value used to form rho = P0/(R*T_clamp) in
    ! lowmach_update_density, otherwise rho*T no longer equals the EOS constant
    ! P0/R. With the same clamp, rho*T_clamp = P0/R holds identically.
    if (have_cp) then
       do concurrent (i = 1:n)
          this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) &
               / (this%rho%x(i,1,1,1) * this%cp_ptr%x(i,1,1,1) &
                  * min(max(this%T_ptr%x(i,1,1,1), this%T_eos_min), &
                        this%T_eos_max))
       end do
    else
       do concurrent (i = 1:n)
          this%Q_T_field%x(i,1,1,1) = this%Q_T_field%x(i,1,1,1) &
               / (this%rho%x(i,1,1,1) * this%cp &
                  * min(max(this%T_ptr%x(i,1,1,1), this%T_eos_min), &
                        this%T_eos_max))
       end do
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine lowmach_update_Q_T

  !> Advance the solution by one time step. Updates rho from the EOS and Q_T
  !! from the temperature field, then routes the pressure/velocity residuals
  !! through the lowmach_residual operators: Q_T enters the pressure residual
  !! as the div u = Q_T source and the velocity residual as the
  !! -(2/3) grad(mu Q_T) dilatation stress. The variable-density momentum RHS
  !! is assembled by the standard rhs makers, which take the rho field.
  subroutine fluid_lowmach_pnpn_step(this, time, dt_controller)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    integer :: n
    type(ksp_monitor_t) :: ksp_results(4)

    if (this%freeze) return

    n = this%dm_Xh%size()

    ! Update the density field from the EOS, fill mu(T)/lambda(T) from the
    ! built-in property model (if configured), then the thermal divergence Q_T
    ! from the temperature field; all feed the variable-density residuals below.
    call lowmach_update_density(this)
    call lowmach_update_properties(this)
    call lowmach_update_Q_T(this)

    call profiler_start_region('Fluid', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         u_e => this%u_e, v_e => this%v_e, w_e => this%w_e, &
         du => this%du, dv => this%dv, dw => this%dw, dp => this%dp, &
         u_res => this%u_res, v_res => this%v_res, w_res => this%w_res, &
         p_res => this%p_res, Ax_vel => this%Ax_vel, Ax_prs => this%Ax_prs, &
         Xh => this%Xh, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         ulag => this%ulag, vlag => this%vlag, wlag => this%wlag, &
         msh => this%msh, source_term => this%source_term, &
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

      ! The momentum RHS is weighted by the spatially varying density field
      ! rho%x; the standard rhs makers consume the full density array.
      if (oifs) then
         call this%adv%compute(u, v, w, &
              this%advx, this%advy, this%advz, &
              Xh, this%c_Xh, dm_Xh%size(), dt)

         call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
              this%abx2, this%aby2, this%abz2, &
              f_x%x, f_y%x, f_z%x, &
              rho%x, ext_bdf%advection_coeffs%x, n)

         call makeoifs%compute_fluid(this%advx%x, this%advy%x, this%advz%x, &
              f_x%x, f_y%x, f_z%x, &
              rho%x, dt, n)
      else
         call this%adv%compute(u, v, w, &
              f_x, f_y, f_z, &
              Xh, this%c_Xh, dm_Xh%size())

         call makeabf%compute_fluid(this%abx1, this%aby1, this%abz1,&
              this%abx2, this%aby2, this%abz2, &
              f_x%x, f_y%x, f_z%x, &
              rho%x, ext_bdf%advection_coeffs%x, n)

         call makebdf%compute_fluid(ulag, vlag, wlag, f_x%x, f_y%x, f_z%x, &
              u, v, w, c_Xh%B, rho%x, dt, &
              ext_bdf%diffusion_coeffs%x, ext_bdf%ndiff, n)
      end if

      call ulag%update()
      call vlag%update()
      call wlag%update()

      call this%bc_apply_vel(time, strong = .true.)
      call this%bc_apply_prs(time)

      call this%update_material_properties(time)

      call profiler_start_region('Pressure_solve', 3)

      ! Variable-density pressure solve: flexible GMRES on the true 1/rho(x)
      ! operator preconditioned by a constant-coefficient hsmg, exactly like
      ! Nek5000's low-Mach pressure step (see lowmach_pressure_solve).
      call lowmach_pressure_solve(this, time, dt_controller, &
           ext_bdf%diffusion_coeffs%x(1), dt, n, ksp_results(1))

      call profiler_end_region('Pressure_solve', 3)

      call profiler_start_region('Velocity_residual', 19)
      call this%lm_vel_res%compute(Ax_vel, u, v, w, &
           u_res, v_res, w_res, &
           p, &
           f_x, f_y, f_z, &
           c_Xh, msh, Xh, &
           mu_tot, rho, this%Q_T_field, &
           ext_bdf%diffusion_coeffs%x(1), dt, dm_Xh%size())

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

  !> Pressure solve for the variable-density (low-Mach) case -- a direct
  !! replica of Nek5000's low-Mach pressure step (plan4.f / hmh_gmres).
  !!
  !! Nek5000 solves the TRUE variable-density operator A = div((1/rho) grad p)
  !! with a flexible, residual-projected GMRES, preconditioned by a
  !! CONSTANT-coefficient geometric multigrid (its h1mg sets mg_h1 = const and
  !! the Schwarz/FDM smoother is geometry-only). The Krylov method absorbs the
  !! spatially varying 1/rho; the preconditioner only accelerates.
  !!
  !! We reproduce that here WITHOUT touching the shared hsmg or the
  !! incompressible solver: pc_prs%update is called while c_Xh%h1 holds a single
  !! constant cref, so hsmg only ever sees a uniform coefficient (exactly as in
  !! the constant-density case); c_Xh%h1 is then restored to 1/rho(x) for the
  !! GMRES operator. One GMRES pass per step -- no defect-correction loop.
  subroutine lowmach_pressure_solve(this, time, dt_controller, bd, dt, n, &
       ksp_result)
    class(fluid_lowmach_pnpn_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    real(kind=rp), intent(in) :: bd, dt
    integer, intent(in) :: n
    type(ksp_monitor_t), intent(out) :: ksp_result
    type(c_ptr) :: event
    real(kind=rp) :: cref, hmin, hmax
    integer :: i, nl

    associate (c_Xh => this%c_Xh, gs_Xh => this%gs_Xh, Ax_prs => this%Ax_prs, &
         p => this%p, p_res => this%p_res, dp => this%dp, rho => this%rho)

      event = glb_cmd_event
      ! Mutable copy: projection pre/post_solving declare their size arg inout.
      nl = n

      ! 1. True variable-density pressure residual. lm_prs_res sets
      !    c_Xh%h1 = 1/rho(x), i.e. the operator A = div((1/rho) grad p),
      !    matching Nek5000 plan4 (h1 = 1/vtrans, h2 = 0).
      call this%lm_prs_res%compute(p, p_res, this%u, this%v, this%w, &
           this%u_e, this%v_e, this%w_e, this%f_x, this%f_y, this%f_z, &
           c_Xh, gs_Xh, this%bc_prs_surface, this%bc_sym_surface, &
           Ax_prs, bd, dt, this%mu_tot, rho, this%Q_T_field, event)

      if (.not. this%prs_dirichlet) call ortho(p_res%x, this%glb_n_points, n)
      call gs_Xh%op(p_res, GS_OP_ADD, event)
      call device_event_sync(event)
      call this%bclst_dp%apply_scalar(p_res%x, n, time)

      ! 2. A-conjugate residual projection against the TRUE variable operator
      !    (Nek5000 project1); c_Xh%h1 is still 1/rho here.
      call this%proj_prs%pre_solving(p_res%x, time%tstep, c_Xh, nl, &
           dt_controller, Ax = Ax_prs, gs_h = gs_Xh, bclst = this%bclst_dp, &
           string = 'Pressure')

      ! 3. Refresh the preconditioner on a CONSTANT-coefficient operator
      !    (Nek5000 h1mg uses mg_h1 = const). cref = midpoint of the 1/rho
      !    range. hsmg only ever sees a uniform h1 here, so neither pc_hsmg
      !    nor the incompressible solver is affected.
      hmin = glmin(c_Xh%h1, n)
      hmax = glmax(c_Xh%h1, n)
      cref = 0.5_rp * (hmin + hmax)
      do concurrent (i = 1:n)
         c_Xh%h1(i,1,1,1) = cref
      end do
      c_Xh%ifh2 = .false.
      call this%pc_prs%update()

      ! 4. Restore the TRUE variable operator and solve it with ONE flexible
      !    GMRES pass (Nek5000 hmh_gmres). The constant-coefficient hsmg only
      !    preconditions; GMRES handles 1/rho(x).
      do concurrent (i = 1:n)
         c_Xh%h1(i,1,1,1) = 1.0_rp / rho%x(i,1,1,1)
      end do
      call rzero(dp%x, n)
      ksp_result = this%ksp_prs%solve(Ax_prs, dp, p_res%x, n, c_Xh, &
           this%bclst_dp, gs_Xh)
      ksp_result%name = 'Pressure'

      ! 5. Reconstruct/store the projection basis (Nek5000 project2), update p.
      call this%proj_prs%post_solving(dp%x, Ax_prs, c_Xh, this%bclst_dp, &
           gs_Xh, nl, time%tstep, dt_controller)
      call field_add2(p, dp, n)
      if (.not. this%prs_dirichlet) call ortho(p%x, this%glb_n_points, n)

    end associate
  end subroutine lowmach_pressure_solve

end module fluid_lowmach_pnpn
