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
  use operators, only : div, rotate_cyc, cdtp, dudxyz
  use math, only : vdot3, cmult, copy, col2, add2
  use gs_ops, only : GS_OP_ADD
  use scratch_registry, only : neko_scratch_registry
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  use facet_normal, only : facet_normal_t
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
       E, p, u, v, w, Ax, &
       coef, gs, h, artificial_visc, mu, kappa, bc_visc_surface, &
       rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(facet_normal_t), intent(in) :: bc_visc_surface
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
            p, u, v, w, Ax, &
            coef, gs, h, artificial_visc, mu, kappa, bc_visc_surface)
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
       rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
       coef, gs, h, artificial_visc, mu, kappa, bc_visc_surface)
    type(field_t), intent(inout) :: rhs_rho_field, &
         rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(facet_normal_t), intent(in) :: bc_visc_surface
    integer :: i, n, n_tmp
    ! Scratch fields for flux components (shared by both paths)
    type(field_t), pointer :: f_x, f_y, f_z
    ! Scratch fields for viscous diffusion (both paths)
    type(field_t), pointer :: visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    ! Scratch fields for GP path: density gradients, cdtp workspace
    type(field_t), pointer :: drhodx, drhody, drhodz, cdtp_tmp
    integer :: tmp_indices(12)
    real(kind=rp) :: mu_eff, kappa_eff, half_mu

    n = coef%dof%size()

    ! Common scratch fields for flux vector components
    call neko_scratch_registry%request_field(f_x, tmp_indices(1), .false.)
    call neko_scratch_registry%request_field(f_y, tmp_indices(2), .false.)
    call neko_scratch_registry%request_field(f_z, tmp_indices(3), .false.)
    n_tmp = 3

    if (euler_res_cpu_viscous_flux_type == VISCOUS_FLUX_NAVIER_STOKES) then
       ! ================================================================
       ! GUERMOND-POPOV VISCOUS FLUX PATH
       ! Phase 1: Inviscid fluxes via div (strong form) + gs + mult
       !          (same formulation as monolithic, known stable).
       ! Phase 2: GP viscous regularization with two coefficients:
       !   μ_eff = max(art_visc, μ_physical)
       !   κ_eff = max(art_visc, κ_physical)
       !
       ! GP viscous flux (Guermond & Popov, CMAME 2014):
       !   Density:  κ_eff ∇²ρ
       !   Momentum: μ_eff ∇²m_i + (κ_eff - μ_eff) Σ_j ∂(u_j ∂ρ/∂x_i)/∂x_j
       !   Energy:   κ_eff ∇²E + (μ_eff - κ_eff)(ρ/2) Σ_j ∂(∂|u|²/∂x_j)/∂x_j
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

       ! Energy: div([(E+p)*u, (E+p)*v, (E+p)*w])
       ! Use BC-enforced velocity (u,v,w) to ensure zero energy flux
       ! at no-slip walls, matching the monolithic path.
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
               * u%x(i,1,1,1)
          f_y%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
               * v%x(i,1,1,1)
          f_z%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) &
               * w%x(i,1,1,1)
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
       ! Phase 2: GP viscous terms
       ! ================================================================

       ! --- Allocate scratch fields ---
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
       call neko_scratch_registry%request_field(cdtp_tmp, tmp_indices(9), &
            .false.)
       call neko_scratch_registry%request_field(drhodx, tmp_indices(10), &
            .false.)
       call neko_scratch_registry%request_field(drhody, tmp_indices(11), &
            .false.)
       call neko_scratch_registry%request_field(drhodz, tmp_indices(12), &
            .false.)
       n_tmp = 12

       ! --- Compute density gradients ---
       call dudxyz(drhodx%x, rho_field%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call dudxyz(drhody%x, rho_field%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call dudxyz(drhodz%x, rho_field%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)

       coef%ifh2 = .false.

       ! ---------------------------------------------------------------
       ! GP density: Ax(κ_eff, ρ)
       ! ---------------------------------------------------------------
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = max(artificial_visc%x(i,1,1,1), &
               kappa%x(i,1,1,1))
       end do
       call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)

       ! ---------------------------------------------------------------
       ! GP momentum: Ax(μ_eff, m_i)
       ! ---------------------------------------------------------------
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = max(artificial_visc%x(i,1,1,1), &
               mu%x(i,1,1,1))
       end do
       call Ax%compute(visc_m_x%x, m_x%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_y%x, m_y%x, coef, p%msh, p%Xh)
       call Ax%compute(visc_m_z%x, m_z%x, coef, p%msh, p%Xh)

       ! ---------------------------------------------------------------
       ! GP momentum corrections (symmetric gradient formulation).
       !
       ! The full GP momentum flux is:
       !   F_{ij} = μ ρ S_{ij} + κ (∂_i ρ) u_j
       ! where S_{ij} = (∂_j u_i + ∂_i u_j) / 2.
       !
       ! After subtracting the Ax part (μ ∂_j m_i), the cdtp
       ! correction for component i, direction j is:
       !   corr_{ij} = (μ/2)(∂_i m_j - ∂_j m_i)
       !            - (μ/2) u_i ∂_j ρ
       !            + (κ - μ/2) u_j ∂_i ρ
       !
       ! For j = i (diagonal), this simplifies to (κ - μ) u_i ∂_i ρ.
       ! For j ≠ i (off-diagonal), momentum gradients are needed.
       ! When κ_eff == μ_eff, all corrections vanish.
       ! ---------------------------------------------------------------

       ! ============== x-momentum corrections ==============

       ! j=x (diagonal): (κ - μ) u ∂ρ/∂x
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_x%x(i,1,1,1) = (kappa_eff - mu_eff) &
               * u%x(i,1,1,1) * drhodx%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_m_x%x, cdtp_tmp%x, n)

       ! j=y: (μ/2)(∂m_y/∂x - ∂m_x/∂y) - (μ/2) u ∂ρ/∂y
       !     + (κ - μ/2) v ∂ρ/∂x
       call dudxyz(f_y%x, m_y%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
       call dudxyz(cdtp_tmp%x, m_x%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          half_mu = 0.5_rp * mu_eff
          f_y%x(i,1,1,1) = half_mu &
               * (f_y%x(i,1,1,1) - cdtp_tmp%x(i,1,1,1)) &
               - half_mu * u%x(i,1,1,1) * drhody%x(i,1,1,1) &
               + (kappa_eff - half_mu) &
               * v%x(i,1,1,1) * drhodx%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_m_x%x, cdtp_tmp%x, n)

       ! j=z: (μ/2)(∂m_z/∂x - ∂m_x/∂z) - (μ/2) u ∂ρ/∂z
       !     + (κ - μ/2) w ∂ρ/∂x
       call dudxyz(f_z%x, m_z%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
       call dudxyz(cdtp_tmp%x, m_x%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          half_mu = 0.5_rp * mu_eff
          f_z%x(i,1,1,1) = half_mu &
               * (f_z%x(i,1,1,1) - cdtp_tmp%x(i,1,1,1)) &
               - half_mu * u%x(i,1,1,1) * drhodz%x(i,1,1,1) &
               + (kappa_eff - half_mu) &
               * w%x(i,1,1,1) * drhodx%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_m_x%x, cdtp_tmp%x, n)

       ! ============== y-momentum corrections ==============

       ! j=x: (μ/2)(∂m_x/∂y - ∂m_y/∂x) - (μ/2) v ∂ρ/∂x
       !     + (κ - μ/2) u ∂ρ/∂y
       call dudxyz(f_x%x, m_x%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
       call dudxyz(cdtp_tmp%x, m_y%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          half_mu = 0.5_rp * mu_eff
          f_x%x(i,1,1,1) = half_mu &
               * (f_x%x(i,1,1,1) - cdtp_tmp%x(i,1,1,1)) &
               - half_mu * v%x(i,1,1,1) * drhodx%x(i,1,1,1) &
               + (kappa_eff - half_mu) &
               * u%x(i,1,1,1) * drhody%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_m_y%x, cdtp_tmp%x, n)

       ! j=y (diagonal): (κ - μ) v ∂ρ/∂y
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_y%x(i,1,1,1) = (kappa_eff - mu_eff) &
               * v%x(i,1,1,1) * drhody%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_m_y%x, cdtp_tmp%x, n)

       ! j=z: (μ/2)(∂m_z/∂y - ∂m_y/∂z) - (μ/2) v ∂ρ/∂z
       !     + (κ - μ/2) w ∂ρ/∂y
       call dudxyz(f_z%x, m_z%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
       call dudxyz(cdtp_tmp%x, m_y%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          half_mu = 0.5_rp * mu_eff
          f_z%x(i,1,1,1) = half_mu &
               * (f_z%x(i,1,1,1) - cdtp_tmp%x(i,1,1,1)) &
               - half_mu * v%x(i,1,1,1) * drhodz%x(i,1,1,1) &
               + (kappa_eff - half_mu) &
               * w%x(i,1,1,1) * drhody%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_m_y%x, cdtp_tmp%x, n)

       ! ============== z-momentum corrections ==============

       ! j=x: (μ/2)(∂m_x/∂z - ∂m_z/∂x) - (μ/2) w ∂ρ/∂x
       !     + (κ - μ/2) u ∂ρ/∂z
       call dudxyz(f_x%x, m_x%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
       call dudxyz(cdtp_tmp%x, m_z%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          half_mu = 0.5_rp * mu_eff
          f_x%x(i,1,1,1) = half_mu &
               * (f_x%x(i,1,1,1) - cdtp_tmp%x(i,1,1,1)) &
               - half_mu * w%x(i,1,1,1) * drhodx%x(i,1,1,1) &
               + (kappa_eff - half_mu) &
               * u%x(i,1,1,1) * drhodz%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_m_z%x, cdtp_tmp%x, n)

       ! j=y: (μ/2)(∂m_y/∂z - ∂m_z/∂y) - (μ/2) w ∂ρ/∂y
       !     + (κ - μ/2) v ∂ρ/∂z
       call dudxyz(f_y%x, m_y%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
       call dudxyz(cdtp_tmp%x, m_z%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          half_mu = 0.5_rp * mu_eff
          f_y%x(i,1,1,1) = half_mu &
               * (f_y%x(i,1,1,1) - cdtp_tmp%x(i,1,1,1)) &
               - half_mu * w%x(i,1,1,1) * drhody%x(i,1,1,1) &
               + (kappa_eff - half_mu) &
               * v%x(i,1,1,1) * drhodz%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_m_z%x, cdtp_tmp%x, n)

       ! j=z (diagonal): (κ - μ) w ∂ρ/∂z
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_z%x(i,1,1,1) = (kappa_eff - mu_eff) &
               * w%x(i,1,1,1) * drhodz%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_m_z%x, cdtp_tmp%x, n)

       ! ---------------------------------------------------------------
       ! GP energy: Ax(κ_eff, E)
       ! ---------------------------------------------------------------
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = max(artificial_visc%x(i,1,1,1), &
               kappa%x(i,1,1,1))
       end do
       call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

       ! ---------------------------------------------------------------
       ! GP energy correction (symmetric gradient formulation).
       !
       ! The GP energy flux is:
       !   F_j^E = κ ∂_j(ρe) + (|u|²/2) κ ∂_j ρ + μ ρ Σ_k S_{jk} u_k
       !         = κ ∂_j E + (μ/2 - κ) ρ Σ_k u_k ∂_j u_k
       !                   + (μ/2) ρ Σ_k u_k ∂_k u_j
       !
       ! After subtracting Ax(κ, E), the cdtp correction is:
       !   f_j = (μ/2 - κ) ρ Σ_k u_k ∂_j u_k + (μ/2) ρ Σ_k u_k ∂_k u_j
       !
       ! Combined coefficient per gradient for each j:
       !   ∂_j u_j × u_j: (μ - κ)     [both terms contribute]
       !   ∂_j u_k × u_k: (μ/2 - κ)   [Term 1 only, k≠j]
       !   ∂_k u_j × u_k: μ/2         [Term 2 only, k≠j]
       !
       ! When μ_eff == κ_eff, all corrections vanish.
       ! ---------------------------------------------------------------

       ! --- j = x direction ---
       ! ∂u/∂x: coefficient (μ - κ) × u [shared by both terms]
       call dudxyz(cdtp_tmp%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_x%x(i,1,1,1) = (mu_eff - kappa_eff) &
               * u%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂v/∂x: coefficient (μ/2 - κ) × v [Term 1 only]
       call dudxyz(cdtp_tmp%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_x%x(i,1,1,1) = f_x%x(i,1,1,1) &
               + (0.5_rp * mu_eff - kappa_eff) &
               * v%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂w/∂x: coefficient (μ/2 - κ) × w [Term 1 only]
       call dudxyz(cdtp_tmp%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_x%x(i,1,1,1) = f_x%x(i,1,1,1) &
               + (0.5_rp * mu_eff - kappa_eff) &
               * w%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂u/∂y: coefficient (μ/2) × v [Term 2 only]
       call dudxyz(cdtp_tmp%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          f_x%x(i,1,1,1) = f_x%x(i,1,1,1) &
               + 0.5_rp * mu_eff &
               * v%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂u/∂z: coefficient (μ/2) × w [Term 2 only]
       call dudxyz(cdtp_tmp%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          f_x%x(i,1,1,1) = f_x%x(i,1,1,1) &
               + 0.5_rp * mu_eff &
               * w%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! Multiply accumulated sum by ρ
       do concurrent (i = 1:n)
          f_x%x(i,1,1,1) = rho_field%x(i,1,1,1) * f_x%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_x%x, coef%drdx, coef%dsdx, &
            coef%dtdx, coef)
       call add2(visc_E%x, cdtp_tmp%x, n)

       ! --- j = y direction ---
       ! ∂u/∂y: coefficient (μ/2 - κ) × u [Term 1 only]
       call dudxyz(cdtp_tmp%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_y%x(i,1,1,1) = (0.5_rp * mu_eff - kappa_eff) &
               * u%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂v/∂y: coefficient (μ - κ) × v [both terms]
       call dudxyz(cdtp_tmp%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_y%x(i,1,1,1) = f_y%x(i,1,1,1) &
               + (mu_eff - kappa_eff) &
               * v%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂w/∂y: coefficient (μ/2 - κ) × w [Term 1 only]
       call dudxyz(cdtp_tmp%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_y%x(i,1,1,1) = f_y%x(i,1,1,1) &
               + (0.5_rp * mu_eff - kappa_eff) &
               * w%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂v/∂x: coefficient (μ/2) × u [Term 2 only]
       call dudxyz(cdtp_tmp%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          f_y%x(i,1,1,1) = f_y%x(i,1,1,1) &
               + 0.5_rp * mu_eff &
               * u%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂v/∂z: coefficient (μ/2) × w [Term 2 only]
       call dudxyz(cdtp_tmp%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          f_y%x(i,1,1,1) = f_y%x(i,1,1,1) &
               + 0.5_rp * mu_eff &
               * w%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! Multiply accumulated sum by ρ
       do concurrent (i = 1:n)
          f_y%x(i,1,1,1) = rho_field%x(i,1,1,1) * f_y%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_y%x, coef%drdy, coef%dsdy, &
            coef%dtdy, coef)
       call add2(visc_E%x, cdtp_tmp%x, n)

       ! --- j = z direction ---
       ! ∂u/∂z: coefficient (μ/2 - κ) × u [Term 1 only]
       call dudxyz(cdtp_tmp%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_z%x(i,1,1,1) = (0.5_rp * mu_eff - kappa_eff) &
               * u%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂v/∂z: coefficient (μ/2 - κ) × v [Term 1 only]
       call dudxyz(cdtp_tmp%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_z%x(i,1,1,1) = f_z%x(i,1,1,1) &
               + (0.5_rp * mu_eff - kappa_eff) &
               * v%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂w/∂z: coefficient (μ - κ) × w [both terms]
       call dudxyz(cdtp_tmp%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          kappa_eff = max(artificial_visc%x(i,1,1,1), kappa%x(i,1,1,1))
          f_z%x(i,1,1,1) = f_z%x(i,1,1,1) &
               + (mu_eff - kappa_eff) &
               * w%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂w/∂x: coefficient (μ/2) × u [Term 2 only]
       call dudxyz(cdtp_tmp%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          f_z%x(i,1,1,1) = f_z%x(i,1,1,1) &
               + 0.5_rp * mu_eff &
               * u%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! ∂w/∂y: coefficient (μ/2) × v [Term 2 only]
       call dudxyz(cdtp_tmp%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
       do concurrent (i = 1:n)
          mu_eff = max(artificial_visc%x(i,1,1,1), mu%x(i,1,1,1))
          f_z%x(i,1,1,1) = f_z%x(i,1,1,1) &
               + 0.5_rp * mu_eff &
               * v%x(i,1,1,1) * cdtp_tmp%x(i,1,1,1)
       end do
       ! Multiply accumulated sum by ρ
       do concurrent (i = 1:n)
          f_z%x(i,1,1,1) = rho_field%x(i,1,1,1) * f_z%x(i,1,1,1)
       end do
       call cdtp(cdtp_tmp%x, f_z%x, coef%drdz, coef%dsdz, &
            coef%dtdz, coef)
       call add2(visc_E%x, cdtp_tmp%x, n)

       ! Reset h1 coefficient back to 1.0
       do concurrent (i = 1:n)
          coef%h1(i,1,1,1) = 1.0_rp
       end do

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
