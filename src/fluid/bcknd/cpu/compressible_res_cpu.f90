! Copyright (c) 2025-2026, The Neko Authors
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
!> This implements CPU-based residual calculations for the compressible equations.
!! It handles the time advancement of primitive variables using
!! Runge-Kutta methods and evaluates the right-hand side terms of the
!! compressible equations including artificial viscosity.
module compressible_res_cpu
  use compressible_residual, only : compressible_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use operators, only : grad, opgrad, rotate_cyc
  use math, only : copy
  use gs_ops, only : GS_OP_ADD
  use scratch_registry, only : neko_scratch_registry
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  use compressible_ops_cpu, only : compressible_ops_cpu_update_uvw, &
       compressible_ops_cpu_update_mxyz_p_ruvw
  use bc_list, only : bc_list_t
  use time_state, only : time_state_t
  implicit none
  private

  type, public, extends(compressible_rhs_t) :: compressible_res_cpu_t
   contains
     procedure, nopass :: step => advance_primitive_variables_cpu
     procedure, nopass :: evaluate_rhs => evaluate_rhs_cpu
  end type compressible_res_cpu_t

  !> Whether physical Navier-Stokes fluxes are active for the current step.
  !! Set by `compressible_rhs_set_physical_flux`, see `compressible_residual`.
  logical, public :: compressible_res_cpu_add_physical_flux = .false.
  !> Whether physical viscous stress is active for the current step.
  !! Set by `compressible_rhs_set_physical_flux`, see `compressible_residual`.
  logical, public :: compressible_res_cpu_add_physical_stress = .false.
  !> Module variable to store thermodynamic parameter set by factory.
  real(kind=rp), public :: compressible_res_cpu_gamma = 1.4_rp


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
  !> @param artificial_visc Artificial viscosity field (entropy viscosity, min with low-order)
  !> @param mu Dynamic viscosity field (physical viscosity for momentum)
  !> @param kappa Thermal conductivity field (physical viscosity for energy)
  !> @param bcs_vel Velocity boundary conditions
  !> @param time Current time state
  !> @param rk_scheme Runge-Kutta time integration scheme
  !> @param dt Time step size
  subroutine advance_primitive_variables_cpu(rho_field, m_x, m_y, m_z, &
       E, p, u, v, w, Ax, &
       Ax_stress, coef, gs, artificial_visc, mu, kappa, bcs_vel, time, &
       rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: bcs_vel
    type(time_state_t), intent(in) :: time
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j, l, ntmp
    logical :: first
    type(field_t), pointer :: k_rho_i, k_m_x_i, k_m_y_i, k_m_z_i, k_E_i, &
         temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
         temp_p, temp_u, temp_v, temp_w, temp_ruvw
    integer :: tmp_indices(30)
    type(field_list_t) :: k_rho, k_m_x, k_m_y, k_m_z, k_E
    ! These contiguous pointers are necessary to ensure that
    ! the Fujitsu compiler properly vectorizes the loop nests
    real(kind=rp), contiguous, pointer :: k_rho_ptr(:,:,:,:)
    real(kind=rp), contiguous, pointer :: k_m_x_ptr(:,:,:,:)
    real(kind=rp), contiguous, pointer :: k_m_y_ptr(:,:,:,:)
    real(kind=rp), contiguous, pointer :: k_m_z_ptr(:,:,:,:)
    real(kind=rp), contiguous, pointer :: k_E_ptr(:,:,:,:)

    n = p%dof%size()
    s = rk_scheme%order

    ! One stage derivative per conserved variable and stage. They are
    ! fully written by evaluate_rhs before they are read, so the scratch
    ! fields are requested without clearing them.
    call k_rho%init(s)
    call k_m_x%init(s)
    call k_m_y%init(s)
    call k_m_z%init(s)
    call k_E%init(s)
    do i = 1, s
       call neko_scratch_registry%request_field(k_rho_i, tmp_indices(5*i-4), &
            .false.)
       call neko_scratch_registry%request_field(k_m_x_i, tmp_indices(5*i-3), &
            .false.)
       call neko_scratch_registry%request_field(k_m_y_i, tmp_indices(5*i-2), &
            .false.)
       call neko_scratch_registry%request_field(k_m_z_i, tmp_indices(5*i-1), &
            .false.)
       call neko_scratch_registry%request_field(k_E_i, tmp_indices(5*i), &
            .false.)
       call k_rho%assign(i, k_rho_i)
       call k_m_x%assign(i, k_m_x_i)
       call k_m_y%assign(i, k_m_y_i)
       call k_m_z%assign(i, k_m_z_i)
       call k_E%assign(i, k_E_i)
    end do
    ntmp = 5 * s
    call neko_scratch_registry%request_field(temp_rho, tmp_indices(ntmp+1), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_x, tmp_indices(ntmp+2), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_y, tmp_indices(ntmp+3), &
         .false.)
    call neko_scratch_registry%request_field(temp_m_z, tmp_indices(ntmp+4), &
         .false.)
    call neko_scratch_registry%request_field(temp_E, tmp_indices(ntmp+5), &
         .false.)
    call neko_scratch_registry%request_field(temp_p, tmp_indices(ntmp+6), &
         .false.)
    call neko_scratch_registry%request_field(temp_u, tmp_indices(ntmp+7), &
         .false.)
    call neko_scratch_registry%request_field(temp_v, tmp_indices(ntmp+8), &
         .false.)
    call neko_scratch_registry%request_field(temp_w, tmp_indices(ntmp+9), &
         .false.)
    call neko_scratch_registry%request_field(temp_ruvw, tmp_indices(ntmp+10), &
         .false.)
    ntmp = ntmp + 10

    ! Loop over Runge-Kutta stages. One parallel region per stage forms the
    ! stage state temp = U + dt * sum_j A(i,j) k_j. Zero entries of the
    ! Butcher tableau are skipped (for the classical RK4 only the
    ! sub-diagonal is non-zero, which halves the sweeps), and the first
    ! non-zero term is fused with the copy of U, saving one pass over all
    ! ten arrays. The zero test is exact on purpose: the tableau is filled
    ! with literal zeros. Each thread sets its private k_*_ptr aliases and
    ! its private copy of first; the implicit barrier at the end of each
    ! !$omp do simd separates the sweeps.
    do i = 1, s
       !$omp parallel private(j, l, first, k_rho_ptr) &
       !$omp& private(k_m_x_ptr, k_m_y_ptr, k_m_z_ptr, k_E_ptr)
       first = .true.
       do j = 1, i-1
          if (rk_scheme%coeffs_A(i, j) .eq. 0.0_rp) cycle
          k_rho_ptr => k_rho%items(j)%ptr%x
          k_m_x_ptr => k_m_x%items(j)%ptr%x
          k_m_y_ptr => k_m_y%items(j)%ptr%x
          k_m_z_ptr => k_m_z%items(j)%ptr%x
          k_E_ptr => k_E%items(j)%ptr%x

          if (first) then
             ! Copy of the solution state fused with the first contribution
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ CONCURRENT
             !DIR$ IVDEP
             !GCC$ ivdep
             !$omp do simd
             do l = 1, n
                temp_rho%x(l,1,1,1) = rho_field%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_rho_ptr(l,1,1,1)
                temp_m_x%x(l,1,1,1) = m_x%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_m_x_ptr(l,1,1,1)
                temp_m_y%x(l,1,1,1) = m_y%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_m_y_ptr(l,1,1,1)
                temp_m_z%x(l,1,1,1) = m_z%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_m_z_ptr(l,1,1,1)
                temp_E%x(l,1,1,1) = E%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_E_ptr(l,1,1,1)
             end do
             !$omp end do simd
             first = .false.
          else
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ CONCURRENT
             !DIR$ IVDEP
             !GCC$ ivdep
             !$omp do simd
             do l = 1, n
                temp_rho%x(l,1,1,1) = temp_rho%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_rho_ptr(l,1,1,1)
                temp_m_x%x(l,1,1,1) = temp_m_x%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_m_x_ptr(l,1,1,1)
                temp_m_y%x(l,1,1,1) = temp_m_y%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_m_y_ptr(l,1,1,1)
                temp_m_z%x(l,1,1,1) = temp_m_z%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_m_z_ptr(l,1,1,1)
                temp_E%x(l,1,1,1) = temp_E%x(l,1,1,1) + &
                     dt * rk_scheme%coeffs_A(i, j) * k_E_ptr(l,1,1,1)
             end do
             !$omp end do simd
          end if
       end do

       ! No stage contribution at all (always the case for the first stage):
       ! plain copy of the solution state
       if (first) then
          !OCL NORECURRENCE, NOVREC, NOALIAS
          !DIR$ CONCURRENT
          !DIR$ IVDEP
          !GCC$ ivdep
          !$omp do simd
          do l = 1, n
             temp_rho%x(l,1,1,1) = rho_field%x(l,1,1,1)
             temp_m_x%x(l,1,1,1) = m_x%x(l,1,1,1)
             temp_m_y%x(l,1,1,1) = m_y%x(l,1,1,1)
             temp_m_z%x(l,1,1,1) = m_z%x(l,1,1,1)
             temp_E%x(l,1,1,1) = E%x(l,1,1,1)
          end do
          !$omp end do simd
       end if
       !$omp end parallel

       ! Evaluate RHS terms using primitive variables from the RK stage state.
       call compressible_ops_cpu_update_uvw(temp_u%x, temp_v%x, temp_w%x, &
            temp_m_x%x, temp_m_y%x, temp_m_z%x, temp_rho%x, n)
       if (compressible_res_cpu_add_physical_stress) then
          call bcs_vel%apply_vector(temp_u%x, temp_v%x, temp_w%x, n, time, &
               strong = .true.)
       end if
       call compressible_ops_cpu_update_mxyz_p_ruvw(temp_m_x%x, temp_m_y%x, &
            temp_m_z%x, temp_p%x, temp_ruvw%x, temp_u%x, temp_v%x, temp_w%x, &
            temp_E%x, temp_rho%x, compressible_res_cpu_gamma, n)

       call evaluate_rhs_cpu(k_rho%items(i)%ptr, k_m_x%items(i)%ptr, &
            k_m_y%items(i)%ptr, k_m_z%items(i)%ptr, &
            k_E%items(i)%ptr, &
            temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
            temp_p, temp_u, temp_v, temp_w, Ax, &
            Ax_stress, coef, gs, artificial_visc, mu, kappa)
    end do

    ! Update the solution. Single parallel region covers all s stages.
    !$omp parallel default(shared) private(i, l, k_rho_ptr) &
    !$omp& private(k_m_x_ptr, k_m_y_ptr, k_m_z_ptr, k_E_ptr)
    do i = 1, s
       k_rho_ptr => k_rho%items(i)%ptr%x
       k_m_x_ptr => k_m_x%items(i)%ptr%x
       k_m_y_ptr => k_m_y%items(i)%ptr%x
       k_m_z_ptr => k_m_z%items(i)%ptr%x
       k_E_ptr => k_E%items(i)%ptr%x

       !OCL NORECURRENCE, NOVREC, NOALIAS
       !DIR$ CONCURRENT
       !DIR$ IVDEP
       !GCC$ ivdep
       !$omp do simd
       do l = 1, n
          rho_field%x(l,1,1,1) = rho_field%x(l,1,1,1) + &
               dt * rk_scheme%coeffs_b(i) * k_rho_ptr(l,1,1,1)
          m_x%x(l,1,1,1) = m_x%x(l,1,1,1) + &
               dt * rk_scheme%coeffs_b(i) * k_m_x_ptr(l,1,1,1)
          m_y%x(l,1,1,1) = m_y%x(l,1,1,1) + &
               dt * rk_scheme%coeffs_b(i) * k_m_y_ptr(l,1,1,1)
          m_z%x(l,1,1,1) = m_z%x(l,1,1,1) + &
               dt * rk_scheme%coeffs_b(i) * k_m_z_ptr(l,1,1,1)
          E%x(l,1,1,1) = E%x(l,1,1,1) + &
               dt * rk_scheme%coeffs_b(i) * k_E_ptr(l,1,1,1)
       end do
       !$omp end do simd
    end do
    !$omp end parallel

    call neko_scratch_registry%relinquish_field(tmp_indices(1:ntmp))

  end subroutine advance_primitive_variables_cpu

  !> Evaluates the right-hand side of the compressible equations.
  !> Inviscid terms are evaluated through div(). Stabilization is always the
  !> existing scalar Laplacian. In navier-stokes mode, physical mu and kappa
  !> add a separate compressible Navier-Stokes viscous flux.
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
  !> @param Ax Matrix-vector product operator for Laplacian terms
  !> @param coef Spatial discretization coefficients
  !> @param gs Gather-scatter operator for parallel communication
  !> @param artificial_visc Artificial viscosity field
  !> @param mu Dynamic viscosity field for physical momentum diffusion
  !> @param kappa Thermal conductivity field for physical energy diffusion
  subroutine evaluate_rhs_cpu(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
       rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
       Ax_stress, coef, gs, artificial_visc, mu, kappa)
    type(field_t), intent(inout) :: rhs_rho_field, &
         rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: i, n
    type(field_t), pointer :: visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    integer :: tmp_indices(5)

    n = coef%dof%size()

    call neko_scratch_registry%request_field(visc_rho, tmp_indices(1), .false.)
    call neko_scratch_registry%request_field(visc_m_x, tmp_indices(2), .false.)
    call neko_scratch_registry%request_field(visc_m_y, tmp_indices(3), .false.)
    call neko_scratch_registry%request_field(visc_m_z, tmp_indices(4), .false.)
    call neko_scratch_registry%request_field(visc_E, tmp_indices(5), .false.)

    !> Inviscid part, rhs = div(F) for all five conserved variables in one
    !! element-local sweep: the mass flux m, the momentum flux
    !! rho u u^T + p I from m and rho, and the energy flux (E + p) u from
    !! the stage velocity.
    call compressible_res_inviscid_div_cpu(rhs_rho_field%x, rhs_m_x%x, &
         rhs_m_y%x, rhs_m_z%x, rhs_E%x, rho_field%x, m_x%x, m_y%x, m_z%x, &
         E%x, p%x, u%x, v%x, w%x, coef)

    ! Density, momentum and energy are all stabilized by the same
    ! artificial viscosity, which is h1 for the Laplacians below.
    call copy(coef%h1, artificial_visc%x, n)

    coef%ifh2 = .false.

    ! Calculate artificial diffusion with variable viscosity. The momentum
    ! components share the operator, so they use the fused vector apply
    ! (geometric factors and h1 are loaded once for all three components).
    call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
    call Ax%compute_vector(visc_m_x%x, visc_m_y%x, visc_m_z%x, &
         m_x%x, m_y%x, m_z%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

    if (compressible_res_cpu_add_physical_flux) then
       call add_navier_stokes_flux_cpu(visc_m_x, visc_m_y, visc_m_z, visc_E, &
            rho_field, p, u, v, w, mu, kappa, Ax, Ax_stress, coef)
    end if

    ! Combine the inviscid divergence and the viscous Laplacian before a
    ! single gather-scatter per field. The divergence is pointwise and wants
    ! the multiplicity average, the Laplacian is a weak form and wants the
    ! inverse of the assembled mass matrix. Both weights are continuous,
    ! i.e. identical on every node of a gather-scatter class, so
    !   gs(a) * mult + Binv * gs(b) = gs(a * mult + Binv * b)
    ! and one exchange replaces two. The sweep also moves div to the rhs
    ! and resets h1 to 1; nothing between here and the exchange reads h1.
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       rhs_rho_field%x(i,1,1,1) = -(rhs_rho_field%x(i,1,1,1) * &
            coef%mult(i,1,1,1) + coef%Binv(i,1,1,1) * visc_rho%x(i,1,1,1))
       rhs_m_x%x(i,1,1,1) = -(rhs_m_x%x(i,1,1,1) * coef%mult(i,1,1,1) + &
            coef%Binv(i,1,1,1) * visc_m_x%x(i,1,1,1))
       rhs_m_y%x(i,1,1,1) = -(rhs_m_y%x(i,1,1,1) * coef%mult(i,1,1,1) + &
            coef%Binv(i,1,1,1) * visc_m_y%x(i,1,1,1))
       rhs_m_z%x(i,1,1,1) = -(rhs_m_z%x(i,1,1,1) * coef%mult(i,1,1,1) + &
            coef%Binv(i,1,1,1) * visc_m_z%x(i,1,1,1))
       rhs_E%x(i,1,1,1) = -(rhs_E%x(i,1,1,1) * coef%mult(i,1,1,1) + &
            coef%Binv(i,1,1,1) * visc_E%x(i,1,1,1))
       coef%h1(i,1,1,1) = 1.0_rp
    end do
    !$omp end parallel do simd

    call gs%op(rhs_rho_field, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 1, coef)
    call gs%op(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, n, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 0, coef)
    call gs%op(rhs_E, GS_OP_ADD)

    call neko_scratch_registry%relinquish_field(tmp_indices)
  end subroutine evaluate_rhs_cpu

  !> Inviscid part of the compressible residual: the divergence of the
  !! density, momentum and energy fluxes, rhs = div(F), in one element-local
  !! sweep.
  !!
  !! The fluxes are formed element by element and differentiated in place,
  !! so no flux component round-trips through a global field and the
  !! geometric factors are read once per element instead of once per
  !! divergence. The momentum flux tensor is symmetric: its three
  !! off-diagonal components are differentiated once and feed two
  !! divergences each, 12 fluxes and 36 local derivatives per point instead
  !! of 15 and 45. The discrete operator is the one div() computes,
  !! jacinv * sum_i (D_r F_i drdx_i + D_s F_i dsdx_i + D_t F_i dtdx_i), with
  !! the x, y and z terms accumulated in that order.
  subroutine compressible_res_inviscid_div_cpu(rhs_rho, rhs_m_x, rhs_m_y, &
       rhs_m_z, rhs_E, rho, m_x, m_y, m_z, E, p, u, v, w, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(inout) :: &
         rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    real(kind=rp), dimension(coef%Xh%lxyz, coef%msh%nelv), intent(in) :: &
         rho, m_x, m_y, m_z, E, p, u, v, w
    real(kind=rp), allocatable :: f(:), fr(:), fs(:), ft(:)
    integer :: el

    associate(Xh => coef%Xh)
      ! One set of element work arrays per thread for the whole sweep
      !$omp parallel private(f, fr, fs, ft)
      allocate(f(Xh%lxyz), fr(Xh%lxyz), fs(Xh%lxyz), ft(Xh%lxyz))
      !$omp do
      do el = 1, coef%msh%nelv
         call inviscid_div_element(Xh%lx, Xh%lxyz, rhs_rho(1,el), &
              rhs_m_x(1,el), rhs_m_y(1,el), rhs_m_z(1,el), rhs_E(1,el), &
              rho(1,el), m_x(1,el), m_y(1,el), m_z(1,el), E(1,el), p(1,el), &
              u(1,el), v(1,el), w(1,el), &
              coef%drdx(1,1,1,el), coef%dsdx(1,1,1,el), coef%dtdx(1,1,1,el), &
              coef%drdy(1,1,1,el), coef%dsdy(1,1,1,el), coef%dtdy(1,1,1,el), &
              coef%drdz(1,1,1,el), coef%dsdz(1,1,1,el), coef%dtdz(1,1,1,el), &
              coef%jacinv(1,1,1,el), Xh%dx, Xh%dy, Xh%dz, f, fr, fs, ft)
      end do
      !$omp end do
      deallocate(f, fr, fs, ft)
      !$omp end parallel
    end associate

  end subroutine compressible_res_inviscid_div_cpu

  !> Fused inviscid divergence of one element, see
  !! compressible_res_inviscid_div_cpu. The 12 fluxes are taken in the
  !! order that sums every divergence x, then y, then z: the mass flux m,
  !! the six distinct components of the symmetric momentum flux
  !! rho u u^T + p I, and the energy flux (E + p) u from the stage velocity.
  !! An off-diagonal momentum flux F_ab is the b component of the flux of
  !! m_a and the a component of the flux of m_b.
  subroutine inviscid_div_element(lx, n, rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, &
       rhs_E, rho, m_x, m_y, m_z, E, p, u, v, w, drdx, dsdx, dtdx, &
       drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv, dx, dy, dz, f, fr, fs, ft)
    integer, intent(in) :: lx, n
    real(kind=rp), dimension(n), intent(inout) :: rhs_rho, rhs_m_x, &
         rhs_m_y, rhs_m_z, rhs_E
    real(kind=rp), dimension(n), intent(in) :: rho, m_x, m_y, m_z, E, p, &
         u, v, w
    real(kind=rp), dimension(n), intent(in) :: drdx, dsdx, dtdx, &
         drdy, dsdy, dtdy, drdz, dsdz, dtdz, jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(n), intent(inout) :: f, fr, fs, ft

    ! rho: F = m
    call local_derivatives(lx, m_x, dx, dy, dz, fr, fs, ft)
    call div_term(n, rhs_rho, fr, fs, ft, drdx, dsdx, dtdx, jacinv)
    call local_derivatives(lx, m_y, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_rho, fr, fs, ft, drdy, dsdy, dtdy, jacinv)
    call local_derivatives(lx, m_z, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_rho, fr, fs, ft, drdz, dsdz, dtdz, jacinv)

    ! m: F_xx, F_xy, F_xz
    call flux_quad_p(n, f, m_x, m_x, rho, p)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term(n, rhs_m_x, fr, fs, ft, drdx, dsdx, dtdx, jacinv)
    call flux_quad(n, f, m_x, m_y, rho)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_m_x, fr, fs, ft, drdy, dsdy, dtdy, jacinv)
    call div_term(n, rhs_m_y, fr, fs, ft, drdx, dsdx, dtdx, jacinv)
    call flux_quad(n, f, m_x, m_z, rho)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_m_x, fr, fs, ft, drdz, dsdz, dtdz, jacinv)
    call div_term(n, rhs_m_z, fr, fs, ft, drdx, dsdx, dtdx, jacinv)

    ! m: F_yy, F_yz
    call flux_quad_p(n, f, m_y, m_y, rho, p)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_m_y, fr, fs, ft, drdy, dsdy, dtdy, jacinv)
    call flux_quad(n, f, m_y, m_z, rho)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_m_y, fr, fs, ft, drdz, dsdz, dtdz, jacinv)
    call div_term_add(n, rhs_m_z, fr, fs, ft, drdy, dsdy, dtdy, jacinv)

    ! m: F_zz
    call flux_quad_p(n, f, m_z, m_z, rho, p)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_m_z, fr, fs, ft, drdz, dsdz, dtdz, jacinv)

    ! E: F = (E + p) u
    call flux_energy(n, f, E, p, u)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term(n, rhs_E, fr, fs, ft, drdx, dsdx, dtdx, jacinv)
    call flux_energy(n, f, E, p, v)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_E, fr, fs, ft, drdy, dsdy, dtdy, jacinv)
    call flux_energy(n, f, E, p, w)
    call local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    call div_term_add(n, rhs_E, fr, fs, ft, drdz, dsdz, dtdz, jacinv)

  end subroutine inviscid_div_element

  !> Momentum flux component f = a b / rho
  subroutine flux_quad(n, f, a, b, rho)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    real(kind=rp), dimension(n), intent(in) :: a, b, rho
    integer :: i

    do i = 1, n
       f(i) = a(i) * b(i) / rho(i)
    end do

  end subroutine flux_quad

  !> Diagonal momentum flux component f = a b / rho + p
  subroutine flux_quad_p(n, f, a, b, rho, p)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    real(kind=rp), dimension(n), intent(in) :: a, b, rho, p
    integer :: i

    do i = 1, n
       f(i) = a(i) * b(i) / rho(i) + p(i)
    end do

  end subroutine flux_quad_p

  !> Energy flux component f = (E + p) u
  subroutine flux_energy(n, f, E, p, u)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    real(kind=rp), dimension(n), intent(in) :: E, p, u
    integer :: i

    do i = 1, n
       f(i) = (E(i) + p(i)) * u(i)
    end do

  end subroutine flux_energy

  !> One direction of a divergence from the local derivatives of a flux,
  !! out = ((fr dr + fs ds) + ft dt) jacinv, as dudxyz forms it.
  subroutine div_term(n, out, fr, fs, ft, dr, ds, dt, jacinv)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: out
    real(kind=rp), dimension(n), intent(in) :: fr, fs, ft, dr, ds, dt, jacinv
    integer :: i

    do i = 1, n
       out(i) = ((fr(i) * dr(i) + fs(i) * ds(i)) + ft(i) * dt(i)) * jacinv(i)
    end do

  end subroutine div_term

  !> As div_term, added to the divergence accumulated so far.
  subroutine div_term_add(n, out, fr, fs, ft, dr, ds, dt, jacinv)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: out
    real(kind=rp), dimension(n), intent(in) :: fr, fs, ft, dr, ds, dt, jacinv
    integer :: i

    do i = 1, n
       out(i) = out(i) + &
            ((fr(i) * dr(i) + fs(i) * ds(i)) + ft(i) * dt(i)) * jacinv(i)
    end do

  end subroutine div_term_add

  !> Local derivatives of an element field, fr = D_r f, fs = D_s f and
  !! ft = D_t f, in one sweep. The first term starts each sum and the
  !! others accumulate in order, so the sums match those of the unrolled
  !! kernels of cpu_dudxyz. Specialised on lx like the other CPU kernels,
  !! so that the short inner loops unroll; the generic version serves the
  !! remaining orders.
  subroutine local_derivatives(lx, f, dx, dy, dz, fr, fs, ft)
    integer, intent(in) :: lx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft

    select case (lx)
    case (14)
       call local_derivatives_lx14(f, dx, dy, dz, fr, fs, ft)
    case (13)
       call local_derivatives_lx13(f, dx, dy, dz, fr, fs, ft)
    case (12)
       call local_derivatives_lx12(f, dx, dy, dz, fr, fs, ft)
    case (11)
       call local_derivatives_lx11(f, dx, dy, dz, fr, fs, ft)
    case (10)
       call local_derivatives_lx10(f, dx, dy, dz, fr, fs, ft)
    case (9)
       call local_derivatives_lx9(f, dx, dy, dz, fr, fs, ft)
    case (8)
       call local_derivatives_lx8(f, dx, dy, dz, fr, fs, ft)
    case (7)
       call local_derivatives_lx7(f, dx, dy, dz, fr, fs, ft)
    case (6)
       call local_derivatives_lx6(f, dx, dy, dz, fr, fs, ft)
    case (5)
       call local_derivatives_lx5(f, dx, dy, dz, fr, fs, ft)
    case (4)
       call local_derivatives_lx4(f, dx, dy, dz, fr, fs, ft)
    case (3)
       call local_derivatives_lx3(f, dx, dy, dz, fr, fs, ft)
    case (2)
       call local_derivatives_lx2(f, dx, dy, dz, fr, fs, ft)
    case default
       call local_derivatives_lx(lx, f, dx, dy, dz, fr, fs, ft)
    end select

  end subroutine local_derivatives

  subroutine local_derivatives_lx(lx, f, dx, dy, dz, fr, fs, ft)
    integer, intent(in) :: lx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx

  subroutine local_derivatives_lx14(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 14
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx14

  subroutine local_derivatives_lx13(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 13
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx13

  subroutine local_derivatives_lx12(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 12
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx12

  subroutine local_derivatives_lx11(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 11
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx11

  subroutine local_derivatives_lx10(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 10
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx10

  subroutine local_derivatives_lx9(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 9
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx9

  subroutine local_derivatives_lx8(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 8
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx8

  subroutine local_derivatives_lx7(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 7
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx7

  subroutine local_derivatives_lx6(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 6
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx6

  subroutine local_derivatives_lx5(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 5
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx5

  subroutine local_derivatives_lx4(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 4
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx4

  subroutine local_derivatives_lx3(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 3
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx3

  subroutine local_derivatives_lx2(f, dx, dy, dz, fr, fs, ft)
    integer, parameter :: lx = 2
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: f
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: fr, fs, ft
    integer :: i, j, k, l

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             fr(i,j,k) = dx(i,1) * f(1,j,k)
             fs(i,j,k) = dy(j,1) * f(i,1,k)
             ft(i,j,k) = dz(k,1) * f(i,j,1)
          end do
          do l = 2, lx
             do i = 1, lx
                fr(i,j,k) = fr(i,j,k) + dx(i,l) * f(l,j,k)
                fs(i,j,k) = fs(i,j,k) + dy(j,l) * f(i,l,k)
                ft(i,j,k) = ft(i,j,k) + dz(k,l) * f(i,j,l)
             end do
          end do
       end do
    end do

  end subroutine local_derivatives_lx2

  !> Add the physical Navier-Stokes flux contribution to the viscous residual.
  !! @param visc_m_x Viscous residual for x-momentum.
  !! @param visc_m_y Viscous residual for y-momentum.
  !! @param visc_m_z Viscous residual for z-momentum.
  !! @param visc_E Viscous residual for total energy.
  !! @param rho_field Density field.
  !! @param p Pressure field.
  !! @param u X-velocity field.
  !! @param v Y-velocity field.
  !! @param w Z-velocity field.
  !! @param mu Dynamic viscosity field.
  !! @param kappa Thermal conductivity field.
  !! @param Ax Matrix-vector product operator for scalar diffusion.
  !! @param Ax_stress Matrix-vector product operator for full stress.
  !! @param coef Spatial discretization coefficients.
  subroutine add_navier_stokes_flux_cpu(visc_m_x, visc_m_y, visc_m_z, visc_E, &
       rho_field, p, u, v, w, mu, kappa, Ax, Ax_stress, coef)
    type(field_t), intent(inout) :: visc_m_x, visc_m_y, visc_m_z, visc_E
    type(field_t), intent(in) :: rho_field
    type(field_t), intent(in) :: p, u, v, w, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(field_t), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, &
         dwdx, dwdy, dwdz, tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, &
         tau_zz, f_x, f_y, f_z, div_flux, dissipation
    integer :: tmp_indices(20)
    integer :: i, n
    real(kind=rp) :: div_u, two_thirds

    n = coef%dof%size()
    two_thirds = 2.0_rp / 3.0_rp

    call neko_scratch_registry%request_field(dudx, tmp_indices(1), .false.)
    call neko_scratch_registry%request_field(dudy, tmp_indices(2), .false.)
    call neko_scratch_registry%request_field(dudz, tmp_indices(3), .false.)
    call neko_scratch_registry%request_field(dvdx, tmp_indices(4), .false.)
    call neko_scratch_registry%request_field(dvdy, tmp_indices(5), .false.)
    call neko_scratch_registry%request_field(dvdz, tmp_indices(6), .false.)
    call neko_scratch_registry%request_field(dwdx, tmp_indices(7), .false.)
    call neko_scratch_registry%request_field(dwdy, tmp_indices(8), .false.)
    call neko_scratch_registry%request_field(dwdz, tmp_indices(9), .false.)
    call neko_scratch_registry%request_field(tau_xx, tmp_indices(10), .false.)
    call neko_scratch_registry%request_field(tau_xy, tmp_indices(11), .false.)
    call neko_scratch_registry%request_field(tau_xz, tmp_indices(12), .false.)
    call neko_scratch_registry%request_field(tau_yy, tmp_indices(13), .false.)
    call neko_scratch_registry%request_field(tau_yz, tmp_indices(14), .false.)
    call neko_scratch_registry%request_field(tau_zz, tmp_indices(15), .false.)
    call neko_scratch_registry%request_field(f_x, tmp_indices(16), .false.)
    call neko_scratch_registry%request_field(f_y, tmp_indices(17), .false.)
    call neko_scratch_registry%request_field(f_z, tmp_indices(18), .false.)
    call neko_scratch_registry%request_field(div_flux, tmp_indices(19), .false.)
    call neko_scratch_registry%request_field(dissipation, tmp_indices(20), &
         .false.)

    call grad(dudx%x, dudy%x, dudz%x, u%x, coef)
    call grad(dvdx%x, dvdy%x, dvdz%x, v%x, coef)
    call grad(dwdx%x, dwdy%x, dwdz%x, w%x, coef)

    ! Viscous stress tensor, dilatational flux, viscous dissipation and
    ! h1 = mu for the stress Laplacian, fused into a single sweep.
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd private(div_u)
    do i = 1, n
       div_u = dudx%x(i,1,1,1) + dvdy%x(i,1,1,1) + dwdz%x(i,1,1,1)
       div_flux%x(i,1,1,1) = mu%x(i,1,1,1) * div_u
       coef%h1(i,1,1,1) = mu%x(i,1,1,1)
       tau_xx%x(i,1,1,1) = mu%x(i,1,1,1) * &
            (2.0_rp * dudx%x(i,1,1,1) - two_thirds * div_u)
       tau_yy%x(i,1,1,1) = mu%x(i,1,1,1) * &
            (2.0_rp * dvdy%x(i,1,1,1) - two_thirds * div_u)
       tau_zz%x(i,1,1,1) = mu%x(i,1,1,1) * &
            (2.0_rp * dwdz%x(i,1,1,1) - two_thirds * div_u)
       tau_xy%x(i,1,1,1) = mu%x(i,1,1,1) * &
            (dudy%x(i,1,1,1) + dvdx%x(i,1,1,1))
       tau_xz%x(i,1,1,1) = mu%x(i,1,1,1) * &
            (dudz%x(i,1,1,1) + dwdx%x(i,1,1,1))
       tau_yz%x(i,1,1,1) = mu%x(i,1,1,1) * &
            (dvdz%x(i,1,1,1) + dwdy%x(i,1,1,1))
       dissipation%x(i,1,1,1) = &
            tau_xx%x(i,1,1,1) * dudx%x(i,1,1,1) &
            + tau_xy%x(i,1,1,1) * (dudy%x(i,1,1,1) + dvdx%x(i,1,1,1)) &
            + tau_xz%x(i,1,1,1) * (dudz%x(i,1,1,1) + dwdx%x(i,1,1,1)) &
            + tau_yy%x(i,1,1,1) * dvdy%x(i,1,1,1) &
            + tau_yz%x(i,1,1,1) * (dvdz%x(i,1,1,1) + dwdy%x(i,1,1,1)) &
            + tau_zz%x(i,1,1,1) * dwdz%x(i,1,1,1)
    end do
    !$omp end parallel do simd

    call Ax_stress%compute_vector(f_x%x, f_y%x, f_z%x, u%x, v%x, w%x, coef, &
         p%msh, p%Xh)
    call opgrad(dudx%x, dudy%x, dudz%x, div_flux%x, coef)

    ! Accumulate the stress contribution into the viscous residual, add
    ! the dissipation to the energy residual, compute the internal energy
    ! (temperature) for the heat flux and set h1 = kappa for its Laplacian.
    ! Fused: one fork-join instead of four.
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x%x(i,1,1,1) = f_x%x(i,1,1,1) &
            - two_thirds * dudx%x(i,1,1,1)
       f_y%x(i,1,1,1) = f_y%x(i,1,1,1) &
            - two_thirds * dudy%x(i,1,1,1)
       f_z%x(i,1,1,1) = f_z%x(i,1,1,1) &
            - two_thirds * dudz%x(i,1,1,1)
       visc_m_x%x(i,1,1,1) = visc_m_x%x(i,1,1,1) + f_x%x(i,1,1,1)
       visc_m_y%x(i,1,1,1) = visc_m_y%x(i,1,1,1) + f_y%x(i,1,1,1)
       visc_m_z%x(i,1,1,1) = visc_m_z%x(i,1,1,1) + f_z%x(i,1,1,1)
       visc_E%x(i,1,1,1) = visc_E%x(i,1,1,1) &
            + u%x(i,1,1,1) * f_x%x(i,1,1,1) &
            + v%x(i,1,1,1) * f_y%x(i,1,1,1) &
            + w%x(i,1,1,1) * f_z%x(i,1,1,1) &
            - coef%B(i,1,1,1) * dissipation%x(i,1,1,1)
       div_flux%x(i,1,1,1) = p%x(i,1,1,1) / &
            (rho_field%x(i,1,1,1) * (compressible_res_cpu_gamma - 1.0_rp))
       coef%h1(i,1,1,1) = kappa%x(i,1,1,1)
    end do
    !$omp end parallel do simd

    call Ax%compute(dudx%x, div_flux%x, coef, p%msh, p%Xh)

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       visc_E%x(i,1,1,1) = visc_E%x(i,1,1,1) + dudx%x(i,1,1,1)
    end do
    !$omp end parallel do simd

    call neko_scratch_registry%relinquish_field(tmp_indices)

  end subroutine add_navier_stokes_flux_cpu

end module compressible_res_cpu
