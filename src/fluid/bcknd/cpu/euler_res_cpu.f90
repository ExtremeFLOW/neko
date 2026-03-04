! Copyright (c) 2025, The Neko Authors
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
!> This module implements CPU-based residual calculations for the Euler equations.
!? It handles the time advancement of primitive variables using Runge-Kutta methods
!? and evaluates the right-hand side terms of the Euler equations including artificial viscosity.
module euler_res_cpu
  use euler_residual, only : euler_rhs_t
  use viscous_flux, only : VISCOUS_FLUX_MONOLITHIC, VISCOUS_FLUX_NAVIER_STOKES
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use operators, only : div, rotate_cyc, opgrad, cdtp
  use math, only : vdot3, cmult, copy, col2, add2
  use gs_ops, only : GS_OP_ADD
  use scratch_registry, only : neko_scratch_registry
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  implicit none
  private

  type, public, extends(euler_rhs_t) :: euler_res_cpu_t
   contains
     procedure, nopass :: step => advance_primitive_variables_cpu
     procedure, nopass :: evaluate_rhs => evaluate_rhs_cpu
  end type euler_res_cpu_t

  !> Module variables to store viscous flux parameters (set by factory)
  integer, public :: euler_res_cpu_viscous_flux_type = VISCOUS_FLUX_MONOLITHIC
  real(kind=rp), public :: euler_res_cpu_gamma = 1.4_rp

contains
  !> Advances the primitive variables (density, momentum, energy)
  !> in time using a Runge-Kutta scheme
  !> @param rho_field Density field
  !> @param m_x,m_y,m_z Momentum components
  !> @param E Total energy
  !> @param p Pressure field
  !> @param u,v,w Velocity components
  !> @param Ax Matrix-vector product operator
  !> @param coef Coefficients for spatial discretization
  !> @param gs Gather-scatter operator for parallel communication
  !> @param h Mesh size field
  !> @param artificial_visc Artificial viscosity field (entropy viscosity, min with low-order)
  !> @param mu Dynamic viscosity field (physical viscosity for momentum)
  !> @param kappa Thermal conductivity field (physical viscosity for energy)
  !> @param rk_scheme Runge-Kutta time integration scheme
  !> @param dt Time step size
  subroutine advance_primitive_variables_cpu(rho_field, m_x, m_y, m_z, &
       E, p, u, v, w, Ax, Ax_navier_stokes, &
       coef, gs, h, artificial_visc, mu, kappa, rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax
    class(Ax_t), intent(inout) :: Ax_navier_stokes
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j, k
    type(field_t), pointer :: k_rho_1, k_rho_2, k_rho_3, k_rho_4, &
         k_m_x_1, k_m_x_2, k_m_x_3, k_m_x_4, &
         k_m_y_1, k_m_y_2, k_m_y_3, k_m_y_4, &
         k_m_z_1, k_m_z_2, k_m_z_3, k_m_z_4, &
         k_E_1, k_E_2, k_E_3, k_E_4, &
         temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E
    integer :: tmp_indices(25)
    type(field_list_t) :: k_rho, k_m_x, k_m_y, k_m_z, k_E

    n = p%dof%size()
    s = rk_scheme%order
    call neko_scratch_registry%request_field(k_rho_1, tmp_indices(1), .true.)
    call neko_scratch_registry%request_field(k_rho_2, tmp_indices(2), .true.)
    call neko_scratch_registry%request_field(k_rho_3, tmp_indices(3), .true.)
    call neko_scratch_registry%request_field(k_rho_4, tmp_indices(4), .true.)
    call neko_scratch_registry%request_field(k_m_x_1, tmp_indices(5), .true.)
    call neko_scratch_registry%request_field(k_m_x_2, tmp_indices(6), .true.)
    call neko_scratch_registry%request_field(k_m_x_3, tmp_indices(7), .true.)
    call neko_scratch_registry%request_field(k_m_x_4, tmp_indices(8), .true.)
    call neko_scratch_registry%request_field(k_m_y_1, tmp_indices(9), .true.)
    call neko_scratch_registry%request_field(k_m_y_2, tmp_indices(10), .true.)
    call neko_scratch_registry%request_field(k_m_y_3, tmp_indices(11), .true.)
    call neko_scratch_registry%request_field(k_m_y_4, tmp_indices(12), .true.)
    call neko_scratch_registry%request_field(k_m_z_1, tmp_indices(13), .true.)
    call neko_scratch_registry%request_field(k_m_z_2, tmp_indices(14), .true.)
    call neko_scratch_registry%request_field(k_m_z_3, tmp_indices(15), .true.)
    call neko_scratch_registry%request_field(k_m_z_4, tmp_indices(16), .true.)
    call neko_scratch_registry%request_field(k_E_1, tmp_indices(17), .true.)
    call neko_scratch_registry%request_field(k_E_2, tmp_indices(18), .true.)
    call neko_scratch_registry%request_field(k_E_3, tmp_indices(19), .true.)
    call neko_scratch_registry%request_field(k_E_4, tmp_indices(20), .true.)
    call neko_scratch_registry%request_field(temp_rho, tmp_indices(21), .false.)
    call neko_scratch_registry%request_field(temp_m_x, tmp_indices(22), .false.)
    call neko_scratch_registry%request_field(temp_m_y, tmp_indices(23), .false.)
    call neko_scratch_registry%request_field(temp_m_z, tmp_indices(24), .false.)
    call neko_scratch_registry%request_field(temp_E, tmp_indices(25), .false.)

    ! Initialize Runge-Kutta stage variables for each conserved quantity
    call k_rho%init(4)
    call k_rho%assign(1, k_rho_1)
    call k_rho%assign(2, k_rho_2)
    call k_rho%assign(3, k_rho_3)
    call k_rho%assign(4, k_rho_4)
    call k_m_x%init(4)
    call k_m_x%assign(1, k_m_x_1)
    call k_m_x%assign(2, k_m_x_2)
    call k_m_x%assign(3, k_m_x_3)
    call k_m_x%assign(4, k_m_x_4)
    call k_m_y%init(4)
    call k_m_y%assign(1, k_m_y_1)
    call k_m_y%assign(2, k_m_y_2)
    call k_m_y%assign(3, k_m_y_3)
    call k_m_y%assign(4, k_m_y_4)
    call k_m_z%init(4)
    call k_m_z%assign(1, k_m_z_1)
    call k_m_z%assign(2, k_m_z_2)
    call k_m_z%assign(3, k_m_z_3)
    call k_m_z%assign(4, k_m_z_4)
    call k_E%init(4)
    call k_E%assign(1, k_E_1)
    call k_E%assign(2, k_E_2)
    call k_E%assign(3, k_E_3)
    call k_E%assign(4, k_E_4)

    ! Loop over Runge-Kutta stages
    do i = 1, s
       ! Copy current solution state to temporary arrays for this RK stage
       do concurrent (k = 1:n)
          temp_rho%x(k,1,1,1) = rho_field%x(k,1,1,1)
          temp_m_x%x(k,1,1,1) = m_x%x(k,1,1,1)
          temp_m_y%x(k,1,1,1) = m_y%x(k,1,1,1)
          temp_m_z%x(k,1,1,1) = m_z%x(k,1,1,1)
          temp_E%x(k,1,1,1) = E%x(k,1,1,1)
       end do

       ! Accumulate previous stage contributions using RK coefficients
       do j = 1, i-1
          do concurrent (k = 1:n)
             temp_rho%x(k,1,1,1) = temp_rho%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * &
                  k_rho%items(j)%ptr%x(k,1,1,1)
             temp_m_x%x(k,1,1,1) = temp_m_x%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * &
                  k_m_x%items(j)%ptr%x(k,1,1,1)
             temp_m_y%x(k,1,1,1) = temp_m_y%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * &
                  k_m_y%items(j)%ptr%x(k,1,1,1)
             temp_m_z%x(k,1,1,1) = temp_m_z%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * &
                  k_m_z%items(j)%ptr%x(k,1,1,1)
             temp_E%x(k,1,1,1) = temp_E%x(k,1,1,1) &
                  + dt * rk_scheme%coeffs_A(i, j) * &
                  k_E%items(j)%ptr%x(k,1,1,1)
          end do
       end do

       ! Evaluate RHS terms for current stage using intermediate solution values
       call evaluate_rhs_cpu(k_rho%items(i)%ptr, k_m_x%items(i)%ptr, &
            k_m_y%items(i)%ptr, k_m_z%items(i)%ptr, &
            k_E%items(i)%ptr, &
            temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
            p, u, v, w, Ax, Ax_navier_stokes, &
            coef, gs, h, artificial_visc, mu, kappa)
    end do

    ! Update the solution
    do i = 1, s
       do concurrent (k = 1:n)
          rho_field%x(k,1,1,1) = rho_field%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_rho%items(i)%ptr%x(k,1,1,1)
          m_x%x(k,1,1,1) = m_x%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_m_x%items(i)%ptr%x(k,1,1,1)
          m_y%x(k,1,1,1) = m_y%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_m_y%items(i)%ptr%x(k,1,1,1)
          m_z%x(k,1,1,1) = m_z%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_m_z%items(i)%ptr%x(k,1,1,1)
          E%x(k,1,1,1) = E%x(k,1,1,1) &
               + dt * rk_scheme%coeffs_b(i) * k_E%items(i)%ptr%x(k,1,1,1)
       end do
    end do

    call neko_scratch_registry%relinquish_field(tmp_indices)

  end subroutine advance_primitive_variables_cpu

  !> Evaluates the right-hand side of the Euler equations.
  !> For Navier-Stokes viscous flux: all viscous terms (stress tensor,
  !> thermal diffusion, density diffusion) are included directly in the
  !> fluxes and go through div(). No Ax operator is used.
  !> For monolithic viscous flux: inviscid fluxes go through div(), and
  !> diffusion is applied separately via the Ax (Laplacian) operator.
  !> @param rhs_rho_field Output array for density RHS terms
  !> @param rhs_m_x Output array for x-momentum RHS terms
  !> @param rhs_m_y Output array for y-momentum RHS terms
  !> @param rhs_m_z Output array for z-momentum RHS terms
  !> @param rhs_E Output array for energy RHS terms
  !> @param rho_field Input density field
  !> @param m_x Input x-momentum field
  !> @param m_y Input y-momentum field
  !> @param m_z Input z-momentum field
  !> @param E Input total energy field
  !> @param p Input pressure field
  !> @param u Input x-velocity field
  !> @param v Input y-velocity field
  !> @param w Input z-velocity field
  !> @param Ax Matrix-vector product operator (used for monolithic path)
  !> @param Ax_navier_stokes Ax operator for NS (unused in flux-only path)
  !> @param coef Spatial discretization coefficients
  !> @param gs Gather-scatter operator for parallel communication
  !> @param h Mesh size field
  !> @param artificial_visc Artificial viscosity field (entropy viscosity)
  !> @param mu Dynamic viscosity field (physical viscosity for momentum)
  !> @param kappa Thermal conductivity field (for energy diffusion)
  subroutine evaluate_rhs_cpu(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
       rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, Ax_navier_stokes, &
       coef, gs, h, artificial_visc, mu, kappa)
    type(field_t), intent(inout) :: rhs_rho_field, &
         rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax
    class(Ax_t), intent(inout) :: Ax_navier_stokes
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: i, n, n_tmp
    ! Scratch fields for flux components (shared by both paths)
    type(field_t), pointer :: f_x, f_y, f_z
    ! Scratch fields for velocity gradients (NS path only)
    type(field_t), pointer :: s_dudx, s_dudy, s_dudz, &
         s_dvdx, s_dvdy, s_dvdz, s_dwdx, s_dwdy, s_dwdz
    ! Scratch fields for viscous diffusion (both paths)
    type(field_t), pointer :: visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    ! Scratch field for cdtp accumulation (NS path only)
    type(field_t), pointer :: cdtp_tmp
    integer :: tmp_indices(18)
    real(kind=rp) :: mu_eff, kappa_eff
    real(kind=rp) :: txx, txy, txz, tyy, tyz, tzz
    real(kind=rp) :: u_local, v_local, w_local

    n = coef%dof%size()

    ! Common scratch fields for flux vector components
    call neko_scratch_registry%request_field(f_x, tmp_indices(1), .false.)
    call neko_scratch_registry%request_field(f_y, tmp_indices(2), .false.)
    call neko_scratch_registry%request_field(f_z, tmp_indices(3), .false.)
    n_tmp = 3

    if (euler_res_cpu_viscous_flux_type == VISCOUS_FLUX_NAVIER_STOKES) then
       ! ================================================================
       ! NAVIER-STOKES PATH (mixed formulation)
       ! Phase 1: Inviscid fluxes via div (strong form) + gs + mult
       !          (same formulation as monolithic, known stable).
       ! Phase 2: Viscous fluxes (stress tensor, thermal diffusion,
       !          density diffusion) via cdtp (weak form) + gs + Binv.
       !          No Ax operator is used for the viscous part.
       ! ================================================================

       ! --- Phase 1: Inviscid fluxes (strong form, same as monolithic) ---

       ! Density: div(m)
       call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, coef)

       ! x-momentum: div([m_x*m_x/rho + p, m_x*m_y/rho, m_x*m_z/rho])
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = m_x%x(i,1,1,1) * m_x%x(i,1,1,1) / &
               rho_field%x(i,1,1,1) + p%x(i,1,1,1)
          f_y%x(i,1,1,1) = m_x%x(i,1,1,1) * m_y%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_z%x(i,1,1,1) = m_x%x(i,1,1,1) * m_z%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
       end do
       call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, coef)

       ! y-momentum: div([m_y*m_x/rho, m_y*m_y/rho + p, m_y*m_z/rho])
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = m_y%x(i,1,1,1) * m_x%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_y%x(i,1,1,1) = m_y%x(i,1,1,1) * m_y%x(i,1,1,1) / &
               rho_field%x(i,1,1,1) + p%x(i,1,1,1)
          f_z%x(i,1,1,1) = m_y%x(i,1,1,1) * m_z%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
       end do
       call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, coef)

       ! z-momentum: div([m_z*m_x/rho, m_z*m_y/rho, m_z*m_z/rho + p])
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = m_z%x(i,1,1,1) * m_x%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_y%x(i,1,1,1) = m_z%x(i,1,1,1) * m_y%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_z%x(i,1,1,1) = m_z%x(i,1,1,1) * m_z%x(i,1,1,1) / &
               rho_field%x(i,1,1,1) + p%x(i,1,1,1)
       end do
       call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, coef)

       ! Energy: div([(E+p)*m/rho])
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
               * m_x%x(i,1,1,1) / rho_field%x(i,1,1,1)
          f_y%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
               * m_y%x(i,1,1,1) / rho_field%x(i,1,1,1)
          f_z%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
               * m_z%x(i,1,1,1) / rho_field%x(i,1,1,1)
       end do
       call div(rhs_E%x, f_x%x, f_y%x, f_z%x, coef)

       ! Gather-scatter + mult + negate (inviscid part)
       call gs%op(rhs_rho_field, GS_OP_ADD)
       call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 1, coef)
       call gs%op(rhs_m_x, GS_OP_ADD)
       call gs%op(rhs_m_y, GS_OP_ADD)
       call gs%op(rhs_m_z, GS_OP_ADD)
       call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 0, coef)
       call gs%op(rhs_E, GS_OP_ADD)
       do concurrent (i = 1:n)
          rhs_rho_field%x(i,1,1,1) = -rhs_rho_field%x(i,1,1,1) &
               * coef%mult(i,1,1,1)
          rhs_m_x%x(i,1,1,1) = -rhs_m_x%x(i,1,1,1) * coef%mult(i,1,1,1)
          rhs_m_y%x(i,1,1,1) = -rhs_m_y%x(i,1,1,1) * coef%mult(i,1,1,1)
          rhs_m_z%x(i,1,1,1) = -rhs_m_z%x(i,1,1,1) * coef%mult(i,1,1,1)
          rhs_E%x(i,1,1,1) = -rhs_E%x(i,1,1,1) * coef%mult(i,1,1,1)
       end do

       ! ================================================================
       ! Phase 2+3: Artificial viscosity (Ax) + Physical NS flux (cdtp)
       !
       ! Strategy:
       !   - Ax on ALL conserved variables with h1 = artificial_visc
       !     (shock-capturing, operates on conserved variables)
       !   - cdtp of stress tensor (physical NS, operates on velocity)
       !     added to momentum visc accumulators
       !   - cdtp of (tau.u + kappa*grad_E) added to energy visc
       !   - Single gs + Binv pass on combined result
       ! ================================================================

       ! Scratch: vel gradients (9) + visc accumulators (5) + cdtp temp (1)
       call neko_scratch_registry%request_field(s_dudx, tmp_indices(4), .false.)
       call neko_scratch_registry%request_field(s_dudy, tmp_indices(5), .false.)
       call neko_scratch_registry%request_field(s_dudz, tmp_indices(6), .false.)
       call neko_scratch_registry%request_field(s_dvdx, tmp_indices(7), .false.)
       call neko_scratch_registry%request_field(s_dvdy, tmp_indices(8), .false.)
       call neko_scratch_registry%request_field(s_dvdz, tmp_indices(9), .false.)
       call neko_scratch_registry%request_field(s_dwdx, tmp_indices(10), .false.)
       call neko_scratch_registry%request_field(s_dwdy, tmp_indices(11), .false.)
       call neko_scratch_registry%request_field(s_dwdz, tmp_indices(12), .false.)
       call neko_scratch_registry%request_field(visc_rho, tmp_indices(13), &
            .false.)
       call neko_scratch_registry%request_field(visc_m_x, tmp_indices(14), &
            .false.)
       call neko_scratch_registry%request_field(visc_m_y, tmp_indices(15), &
            .false.)
       call neko_scratch_registry%request_field(visc_m_z, tmp_indices(16), &
            .false.)
       call neko_scratch_registry%request_field(visc_E, tmp_indices(17), &
            .false.)
       call neko_scratch_registry%request_field(cdtp_tmp, tmp_indices(18), &
            .false.)
       n_tmp = 18

       ! Compute substage velocity u_sub = m/rho
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = m_x%x(i,1,1,1) / rho_field%x(i,1,1,1)
          f_y%x(i,1,1,1) = m_y%x(i,1,1,1) / rho_field%x(i,1,1,1)
          f_z%x(i,1,1,1) = m_z%x(i,1,1,1) / rho_field%x(i,1,1,1)
       end do
       call opgrad(s_dudx%x, s_dudy%x, s_dudz%x, f_x%x, coef)
       call opgrad(s_dvdx%x, s_dvdy%x, s_dvdz%x, f_y%x, coef)
       call opgrad(s_dwdx%x, s_dwdy%x, s_dwdz%x, f_z%x, coef)

       ! --- Ax on all conserved variables (artificial viscosity) ---
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = artificial_visc%x(i,1,1,1)
       end do
       call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_x%x, m_x%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_y%x, m_y%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_z%x, m_z%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = 1.0_rp
       end do

       ! --- Physical NS: x-momentum, tau_x = mu*[2*du/dx, du/dy+dv/dx, du/dz+dw/dx] ---
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = mu%x(i,1,1,1) * 2.0_rp * s_dudx%x(i,1,1,1)
          f_y%x(i,1,1,1) = mu%x(i,1,1,1) * (s_dudy%x(i,1,1,1) &
               + s_dvdx%x(i,1,1,1))
          f_z%x(i,1,1,1) = mu%x(i,1,1,1) * (s_dudz%x(i,1,1,1) &
               + s_dwdx%x(i,1,1,1))
       end do
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_m_x%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_m_x%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_m_x%x, cdtp_tmp%x, n)

       ! --- Physical NS: y-momentum, tau_y = mu*[du/dy+dv/dx, 2*dv/dy, dv/dz+dw/dy] ---
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = mu%x(i,1,1,1) * (s_dudy%x(i,1,1,1) &
               + s_dvdx%x(i,1,1,1))
          f_y%x(i,1,1,1) = mu%x(i,1,1,1) * 2.0_rp * s_dvdy%x(i,1,1,1)
          f_z%x(i,1,1,1) = mu%x(i,1,1,1) * (s_dvdz%x(i,1,1,1) &
               + s_dwdy%x(i,1,1,1))
       end do
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_m_y%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_m_y%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_m_y%x, cdtp_tmp%x, n)

       ! --- Physical NS: z-momentum, tau_z = mu*[du/dz+dw/dx, dv/dz+dw/dy, 2*dw/dz] ---
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = mu%x(i,1,1,1) * (s_dudz%x(i,1,1,1) &
               + s_dwdx%x(i,1,1,1))
          f_y%x(i,1,1,1) = mu%x(i,1,1,1) * (s_dvdz%x(i,1,1,1) &
               + s_dwdy%x(i,1,1,1))
          f_z%x(i,1,1,1) = mu%x(i,1,1,1) * 2.0_rp * s_dwdz%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_m_z%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_m_z%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_m_z%x, cdtp_tmp%x, n)

       ! --- Physical NS: energy viscous flux ---
       ! F_E_visc = tau.u + kappa*grad(T) (POSITIVE, for rhs -= Binv*cdtp)
       ! where T = p / (rho * (gamma - 1))
       ! Step 1: Compute tau.u (needs all velocity gradients)
       do concurrent (i = 1:n)
          txx = mu%x(i,1,1,1) * 2.0_rp * s_dudx%x(i,1,1,1)
          txy = mu%x(i,1,1,1) * (s_dudy%x(i,1,1,1) + s_dvdx%x(i,1,1,1))
          txz = mu%x(i,1,1,1) * (s_dudz%x(i,1,1,1) + s_dwdx%x(i,1,1,1))
          tyy = mu%x(i,1,1,1) * 2.0_rp * s_dvdy%x(i,1,1,1)
          tyz = mu%x(i,1,1,1) * (s_dvdz%x(i,1,1,1) + s_dwdy%x(i,1,1,1))
          tzz = mu%x(i,1,1,1) * 2.0_rp * s_dwdz%x(i,1,1,1)
          u_local = m_x%x(i,1,1,1) / rho_field%x(i,1,1,1)
          v_local = m_y%x(i,1,1,1) / rho_field%x(i,1,1,1)
          w_local = m_z%x(i,1,1,1) / rho_field%x(i,1,1,1)
          f_x%x(i,1,1,1) = txx * u_local + txy * v_local &
               + txz * w_local
          f_y%x(i,1,1,1) = txy * u_local + tyy * v_local &
               + tyz * w_local
          f_z%x(i,1,1,1) = txz * u_local + tyz * v_local &
               + tzz * w_local
       end do
       ! Step 2: Add kappa*grad(T) where T = p/(rho*(gamma-1))
       ! Compute T locally into s_dvdx (free at this point)
       do concurrent (i = 1:n)
          s_dvdx%x(i,1,1,1) = p%x(i,1,1,1) / &
               (rho_field%x(i,1,1,1) * (euler_res_cpu_gamma - 1.0_rp))
       end do
       call opgrad(s_dudx%x, s_dudy%x, s_dudz%x, s_dvdx%x, coef)
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = f_x%x(i,1,1,1) &
               + kappa%x(i,1,1,1) * s_dudx%x(i,1,1,1)
          f_y%x(i,1,1,1) = f_y%x(i,1,1,1) &
               + kappa%x(i,1,1,1) * s_dudy%x(i,1,1,1)
          f_z%x(i,1,1,1) = f_z%x(i,1,1,1) &
               + kappa%x(i,1,1,1) * s_dudz%x(i,1,1,1)
       end do
       ! cdtp of energy viscous flux, add to visc_E
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_E%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_E%x, cdtp_tmp%x, n)
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_E%x, cdtp_tmp%x, n)

       ! --- Gather-scatter and apply Binv to combined viscous terms ---
       call gs%op(visc_rho, GS_OP_ADD)
       call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 1, coef)
       call gs%op(visc_m_x, GS_OP_ADD)
       call gs%op(visc_m_y, GS_OP_ADD)
       call gs%op(visc_m_z, GS_OP_ADD)
       call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 0, coef)
       call gs%op(visc_E, GS_OP_ADD)
       do concurrent (i = 1:n)
          rhs_rho_field%x(i,1,1,1) = rhs_rho_field%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_rho%x(i,1,1,1)
          rhs_m_x%x(i,1,1,1) = rhs_m_x%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_m_x%x(i,1,1,1)
          rhs_m_y%x(i,1,1,1) = rhs_m_y%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_m_y%x(i,1,1,1)
          rhs_m_z%x(i,1,1,1) = rhs_m_z%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_m_z%x(i,1,1,1)
          rhs_E%x(i,1,1,1) = rhs_E%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_E%x(i,1,1,1)
       end do

    else
       ! ================================================================
       ! MONOLITHIC PATH
       ! Pure inviscid fluxes through div(), then Ax for diffusion.
       ! ================================================================

       ! --- Density flux: F_rho = [m_x, m_y, m_z] ---
       call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, coef)

       ! --- x-momentum flux: F = [m_x*m_x/rho + p, m_x*m_y/rho, m_x*m_z/rho]
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = m_x%x(i,1,1,1) * m_x%x(i,1,1,1) / &
               rho_field%x(i,1,1,1) + p%x(i,1,1,1)
          f_y%x(i,1,1,1) = m_x%x(i,1,1,1) * m_y%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_z%x(i,1,1,1) = m_x%x(i,1,1,1) * m_z%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
       end do
       call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, coef)

       ! --- y-momentum flux: F = [m_y*m_x/rho, m_y*m_y/rho + p, m_y*m_z/rho]
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = m_y%x(i,1,1,1) * m_x%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_y%x(i,1,1,1) = m_y%x(i,1,1,1) * m_y%x(i,1,1,1) / &
               rho_field%x(i,1,1,1) + p%x(i,1,1,1)
          f_z%x(i,1,1,1) = m_y%x(i,1,1,1) * m_z%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
       end do
       call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, coef)

       ! --- z-momentum flux: F = [m_z*m_x/rho, m_z*m_y/rho, m_z*m_z/rho + p]
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = m_z%x(i,1,1,1) * m_x%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_y%x(i,1,1,1) = m_z%x(i,1,1,1) * m_y%x(i,1,1,1) / &
               rho_field%x(i,1,1,1)
          f_z%x(i,1,1,1) = m_z%x(i,1,1,1) * m_z%x(i,1,1,1) / &
               rho_field%x(i,1,1,1) + p%x(i,1,1,1)
       end do
       call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, coef)

       ! --- Energy flux: F = [(E+p)*u, (E+p)*v, (E+p)*w] ---
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) * u%x(i,1,1,1)
          f_y%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) * v%x(i,1,1,1)
          f_z%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) * w%x(i,1,1,1)
       end do
       call div(rhs_E%x, f_x%x, f_y%x, f_z%x, coef)

       ! --- Gather-scatter and multiply ---
       call gs%op(rhs_rho_field, GS_OP_ADD)
       call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 1, coef)
       call gs%op(rhs_m_x, GS_OP_ADD)
       call gs%op(rhs_m_y, GS_OP_ADD)
       call gs%op(rhs_m_z, GS_OP_ADD)
       call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 0, coef)
       call gs%op(rhs_E, GS_OP_ADD)
       do concurrent (i = 1:n)
          rhs_rho_field%x(i,1,1,1) = rhs_rho_field%x(i,1,1,1) &
               * coef%mult(i,1,1,1)
          rhs_m_x%x(i,1,1,1) = rhs_m_x%x(i,1,1,1) * coef%mult(i,1,1,1)
          rhs_m_y%x(i,1,1,1) = rhs_m_y%x(i,1,1,1) * coef%mult(i,1,1,1)
          rhs_m_z%x(i,1,1,1) = rhs_m_z%x(i,1,1,1) * coef%mult(i,1,1,1)
          rhs_E%x(i,1,1,1) = rhs_E%x(i,1,1,1) * coef%mult(i,1,1,1)
       end do

       ! --- Artificial viscosity via Ax (Laplacian) ---
       call neko_scratch_registry%request_field(visc_rho, tmp_indices(4), &
            .false.)
       call neko_scratch_registry%request_field(visc_m_x, tmp_indices(5), &
            .false.)
       call neko_scratch_registry%request_field(visc_m_y, tmp_indices(6), &
            .false.)
       call neko_scratch_registry%request_field(visc_m_z, tmp_indices(7), &
            .false.)
       call neko_scratch_registry%request_field(visc_E, tmp_indices(8), &
            .false.)
       n_tmp = 8

       ! h1 = max(artificial, mu) for density and momentum
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
       end do
       call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_x%x, m_x%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_y%x, m_y%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_z%x, m_z%x, coef, p%msh, p%Xh)

       ! h1 = max(artificial, kappa) for energy
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
       end do
       call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

       ! Reset h1 coefficient back to 1.0
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = 1.0_rp
       end do

       ! Gather-scatter on viscous terms
       call gs%op(visc_rho, GS_OP_ADD)
       call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 1, coef)
       call gs%op(visc_m_x, GS_OP_ADD)
       call gs%op(visc_m_y, GS_OP_ADD)
       call gs%op(visc_m_z, GS_OP_ADD)
       call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 0, coef)
       call gs%op(visc_E, GS_OP_ADD)

       ! rhs = -rhs - Binv * visc_*
       do concurrent (i = 1:n)
          rhs_rho_field%x(i,1,1,1) = -rhs_rho_field%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_rho%x(i,1,1,1)
          rhs_m_x%x(i,1,1,1) = -rhs_m_x%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_m_x%x(i,1,1,1)
          rhs_m_y%x(i,1,1,1) = -rhs_m_y%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_m_y%x(i,1,1,1)
          rhs_m_z%x(i,1,1,1) = -rhs_m_z%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_m_z%x(i,1,1,1)
          rhs_E%x(i,1,1,1) = -rhs_E%x(i,1,1,1) &
               - coef%Binv(i,1,1,1) * visc_E%x(i,1,1,1)
       end do
    end if

    call neko_scratch_registry%relinquish_field(tmp_indices(1:n_tmp))
  end subroutine evaluate_rhs_cpu

end module euler_res_cpu
