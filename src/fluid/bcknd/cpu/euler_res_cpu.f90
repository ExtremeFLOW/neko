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
!> This module implements CPU-based residual calculations for the Euler equations.
!? It handles the time advancement of primitive variables using Runge-Kutta methods
!? and evaluates the right-hand side terms of the Euler equations including artificial viscosity.
module euler_res_cpu
  use euler_residual, only : euler_rhs_t
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use num_types, only : rp
  use operators, only : div, rotate_cyc
  use gs_ops, only : GS_OP_ADD
  use scratch_registry, only : neko_scratch_registry
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use field_list, only : field_list_t
  implicit none
  private

  type, public, extends(euler_rhs_t) :: euler_res_cpu_t
   contains
     procedure, nopass :: step => advance_primitive_variables_cpu
     procedure, nopass :: evaluate_rhs_cpu
  end type euler_res_cpu_t

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
  !> @param effective_visc Effective artificial viscosity field
  !> @param rk_scheme Runge-Kutta time integration scheme
  !> @param dt Time step size
  subroutine advance_primitive_variables_cpu(rho_field, m_x, m_y, m_z, &
       E, p, u, v, w, Ax, &
       coef, gs, h, effective_visc, rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, effective_visc
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j
    type(field_t), pointer :: k_rho_1, k_rho_2, k_rho_3, k_rho_4, &
         k_m_x_1, k_m_x_2, k_m_x_3, k_m_x_4, &
         k_m_y_1, k_m_y_2, k_m_y_3, k_m_y_4, &
         k_m_z_1, k_m_z_2, k_m_z_3, k_m_z_4, &
         k_E_1, k_E_2, k_E_3, k_E_4, &
         temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E
    integer :: tmp_indices(25)
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
       call euler_res_cpu_copy5(temp_rho%x, temp_m_x%x, temp_m_y%x, &
            temp_m_z%x, temp_E%x, rho_field%x, m_x%x, m_y%x, m_z%x, &
            E%x, n)

       ! Accumulate previous stage contributions using RK coefficients
       do j = 1, i-1
          k_rho_ptr => k_rho%items(j)%ptr%x
          k_m_x_ptr => k_m_x%items(j)%ptr%x
          k_m_y_ptr => k_m_y%items(j)%ptr%x
          k_m_z_ptr => k_m_z%items(j)%ptr%x
          k_E_ptr => k_E%items(j)%ptr%x
          call euler_res_cpu_axpy5(temp_rho%x, temp_m_x%x, temp_m_y%x, &
               temp_m_z%x, temp_E%x, k_rho_ptr, k_m_x_ptr, k_m_y_ptr, &
               k_m_z_ptr, k_E_ptr, dt * rk_scheme%coeffs_A(i, j), n)
       end do

       ! Evaluate RHS terms for current stage using intermediate solution values
       call evaluate_rhs_cpu(k_rho%items(i)%ptr, k_m_x%items(i)%ptr, &
            k_m_y%items(i)%ptr, k_m_z%items(i)%ptr, &
            k_E%items(i)%ptr, &
            temp_rho, temp_m_x, temp_m_y, temp_m_z, temp_E, &
            p, u, v, w, Ax, &
            coef, gs, h, effective_visc)
    end do

    ! Update the solution
    do i = 1, s
       k_rho_ptr => k_rho%items(i)%ptr%x
       k_m_x_ptr => k_m_x%items(i)%ptr%x
       k_m_y_ptr => k_m_y%items(i)%ptr%x
       k_m_z_ptr => k_m_z%items(i)%ptr%x
       k_E_ptr => k_E%items(i)%ptr%x
       call euler_res_cpu_axpy5(rho_field%x, m_x%x, m_y%x, m_z%x, E%x, &
            k_rho_ptr, k_m_x_ptr, k_m_y_ptr, k_m_z_ptr, k_E_ptr, &
            dt * rk_scheme%coeffs_b(i), n)
    end do

    call neko_scratch_registry%relinquish_field(tmp_indices)

  end subroutine advance_primitive_variables_cpu

  !> Evaluates the right-hand side of the Euler equations
  !> Computes fluxes, divergence terms, and artificial viscosity
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
  !> @param Ax Matrix-vector product operator
  !> @param coef Spatial discretization coefficients
  !> @param gs Gather-scatter operator for parallel communication
  !> @param h Mesh size field
  !> @param effective_visc Effective artificial viscosity field
  subroutine evaluate_rhs_cpu(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
       rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
       coef, gs, h, effective_visc)
    type(field_t), intent(inout) :: rhs_rho_field, &
         rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, effective_visc
    class(Ax_t), intent(inout) :: Ax
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n
    type(field_t), pointer :: f_x, f_y, f_z, &
         visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    integer :: tmp_indices(8)

    n = coef%dof%size()
    call neko_scratch_registry%request_field(f_x, tmp_indices(1), .false.)
    call neko_scratch_registry%request_field(f_y, tmp_indices(2), .false.)
    call neko_scratch_registry%request_field(f_z, tmp_indices(3), .false.)

    !> rho = rho - dt * div(m)
    ! Compute density flux divergence
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, coef)

    !> m = m - dt * div(rho * u * u^T + p*I)
    ! Compute momentum flux divergences
    ! m_x
    call euler_res_cpu_flux_mx(f_x%x, f_y%x, f_z%x, m_x%x, m_y%x, &
         m_z%x, rho_field%x, p%x, n)
    call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_y
    call euler_res_cpu_flux_my(f_x%x, f_y%x, f_z%x, m_x%x, m_y%x, &
         m_z%x, rho_field%x, p%x, n)
    call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_z
    call euler_res_cpu_flux_mz(f_x%x, f_y%x, f_z%x, m_x%x, m_y%x, &
         m_z%x, rho_field%x, p%x, n)
    call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, coef)

    !> E = E - dt * div(u * (E + p))
    ! Compute energy flux divergence
    call euler_res_cpu_flux_E(f_x%x, f_y%x, f_z%x, E%x, p%x, u%x, v%x, &
         w%x, n)
    call div(rhs_E%x, f_x%x, f_y%x, f_z%x, coef)

    ! gs
    call gs%op(rhs_rho_field, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 1, coef)
    call gs%op(rhs_m_x, GS_OP_ADD)
    call gs%op(rhs_m_y, GS_OP_ADD)
    call gs%op(rhs_m_z, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 0, coef)
    call gs%op(rhs_E, GS_OP_ADD)
    call euler_res_cpu_scale5(rhs_rho_field%x, rhs_m_x%x, rhs_m_y%x, &
         rhs_m_z%x, rhs_E%x, coef%mult, rhs_E%dof%size())

    call neko_scratch_registry%request_field(visc_rho, tmp_indices(4), .false.)
    call neko_scratch_registry%request_field(visc_m_x, tmp_indices(5), .false.)
    call neko_scratch_registry%request_field(visc_m_y, tmp_indices(6), .false.)
    call neko_scratch_registry%request_field(visc_m_z, tmp_indices(7), .false.)
    call neko_scratch_registry%request_field(visc_E, tmp_indices(8), .false.)

    ! Set h1 coefficient to the effective viscosity for the Laplacian operator
    call euler_res_cpu_set_h1(coef%h1, effective_visc%x, n)

    ! Calculate artificial diffusion with variable viscosity
    call Ax%compute(visc_rho%x, rho_field%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_x%x, m_x%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_y%x, m_y%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_m_z%x, m_z%x, coef, p%msh, p%Xh)
    call Ax%compute(visc_E%x, E%x, coef, p%msh, p%Xh)

    ! Reset h1 coefficient back to 1.0 for other operations
    call euler_res_cpu_reset_h1(coef%h1, n)

    ! gs
    call gs%op(visc_rho, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 1, coef)
    call gs%op(visc_m_x, GS_OP_ADD)
    call gs%op(visc_m_y, GS_OP_ADD)
    call gs%op(visc_m_z, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 0, coef)
    call gs%op(visc_E, GS_OP_ADD)

    ! Move div to the rhs and apply artificial viscosity
    ! The viscosity coefficient is already included in the Laplacian operator
    call euler_res_cpu_apply_visc(rhs_rho_field%x, rhs_m_x%x, rhs_m_y%x, &
         rhs_m_z%x, rhs_E%x, visc_rho%x, visc_m_x%x, visc_m_y%x, &
         visc_m_z%x, visc_E%x, coef%Binv, n)

    call neko_scratch_registry%relinquish_field(tmp_indices)
  end subroutine evaluate_rhs_cpu

  subroutine euler_res_cpu_copy5(y1, y2, y3, y4, y5, x1, x2, x3, x4, x5, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: y1(n), y2(n), y3(n), y4(n), y5(n)
    real(kind=rp), intent(in) :: x1(n), x2(n), x3(n), x4(n), x5(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       y1(i) = x1(i)
       y2(i) = x2(i)
       y3(i) = x3(i)
       y4(i) = x4(i)
       y5(i) = x5(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_copy5

  subroutine euler_res_cpu_axpy5(y1, y2, y3, y4, y5, x1, x2, x3, x4, x5, a, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: a
    real(kind=rp), intent(inout) :: y1(n), y2(n), y3(n), y4(n), y5(n)
    real(kind=rp), intent(in) :: x1(n), x2(n), x3(n), x4(n), x5(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       y1(i) = y1(i) + a * x1(i)
       y2(i) = y2(i) + a * x2(i)
       y3(i) = y3(i) + a * x3(i)
       y4(i) = y4(i) + a * x4(i)
       y5(i) = y5(i) + a * x5(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_axpy5

  subroutine euler_res_cpu_flux_mx(f_x, f_y, f_z, m_x, m_y, m_z, rho, p, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: f_x(n), f_y(n), f_z(n)
    real(kind=rp), intent(in) :: m_x(n), m_y(n), m_z(n), rho(n), p(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x(i) = m_x(i) * m_x(i) / rho(i) + p(i)
       f_y(i) = m_x(i) * m_y(i) / rho(i)
       f_z(i) = m_x(i) * m_z(i) / rho(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_flux_mx

  subroutine euler_res_cpu_flux_my(f_x, f_y, f_z, m_x, m_y, m_z, rho, p, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: f_x(n), f_y(n), f_z(n)
    real(kind=rp), intent(in) :: m_x(n), m_y(n), m_z(n), rho(n), p(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x(i) = m_y(i) * m_x(i) / rho(i)
       f_y(i) = m_y(i) * m_y(i) / rho(i) + p(i)
       f_z(i) = m_y(i) * m_z(i) / rho(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_flux_my

  subroutine euler_res_cpu_flux_mz(f_x, f_y, f_z, m_x, m_y, m_z, rho, p, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: f_x(n), f_y(n), f_z(n)
    real(kind=rp), intent(in) :: m_x(n), m_y(n), m_z(n), rho(n), p(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x(i) = m_z(i) * m_x(i) / rho(i)
       f_y(i) = m_z(i) * m_y(i) / rho(i)
       f_z(i) = m_z(i) * m_z(i) / rho(i) + p(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_flux_mz

  subroutine euler_res_cpu_flux_E(f_x, f_y, f_z, E, p, u, v, w, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: f_x(n), f_y(n), f_z(n)
    real(kind=rp), intent(in) :: E(n), p(n), u(n), v(n), w(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x(i) = (E(i) + p(i)) * u(i)
       f_y(i) = (E(i) + p(i)) * v(i)
       f_z(i) = (E(i) + p(i)) * w(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_flux_E

  subroutine euler_res_cpu_scale5(y1, y2, y3, y4, y5, mult, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: y1(n), y2(n), y3(n), y4(n), y5(n)
    real(kind=rp), intent(in) :: mult(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       y1(i) = y1(i) * mult(i)
       y2(i) = y2(i) * mult(i)
       y3(i) = y3(i) * mult(i)
       y4(i) = y4(i) * mult(i)
       y5(i) = y5(i) * mult(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_scale5

  subroutine euler_res_cpu_set_h1(h1, effective_visc, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: h1(n)
    real(kind=rp), intent(in) :: effective_visc(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       h1(i) = effective_visc(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_set_h1

  subroutine euler_res_cpu_reset_h1(h1, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: h1(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       h1(i) = 1.0_rp
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_reset_h1

  subroutine euler_res_cpu_apply_visc(y1, y2, y3, y4, y5, v1, v2, v3, v4, &
       v5, Binv, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: y1(n), y2(n), y3(n), y4(n), y5(n)
    real(kind=rp), intent(in) :: v1(n), v2(n), v3(n), v4(n), v5(n), Binv(n)
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       y1(i) = -y1(i) - Binv(i) * v1(i)
       y2(i) = -y2(i) - Binv(i) * v2(i)
       y3(i) = -y3(i) - Binv(i) * v3(i)
       y4(i) = -y4(i) - Binv(i) * v4(i)
       y5(i) = -y5(i) - Binv(i) * v5(i)
    end do
    !$omp end parallel do simd
  end subroutine euler_res_cpu_apply_visc

end module euler_res_cpu
