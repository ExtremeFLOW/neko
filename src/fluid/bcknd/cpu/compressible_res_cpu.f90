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
  use operators, only : div, grad, opgrad, rotate_cyc
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
  logical :: compressible_res_cpu_add_physical_flux = .false.
  !> Whether physical viscous stress is active for the current step.
  logical :: compressible_res_cpu_add_physical_stress = .false.
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
  !> @param h Mesh size field
  !> @param artificial_visc Artificial viscosity field (entropy viscosity, min with low-order)
  !> @param mu Dynamic viscosity field (physical viscosity for momentum)
  !> @param kappa Thermal conductivity field (physical viscosity for energy)
  !> @param bcs_vel Velocity boundary conditions
  !> @param time Current time state
  !> @param rk_scheme Runge-Kutta time integration scheme
  !> @param dt Time step size
  subroutine advance_primitive_variables_cpu(rho_field, m_x, m_y, m_z, &
       E, p, u, v, w, Ax, &
       Ax_stress, coef, gs, h, artificial_visc, mu, kappa, bcs_vel, time, &
       rk_scheme, dt)
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: bcs_vel
    type(time_state_t), intent(in) :: time
    class(runge_kutta_time_scheme_t), intent(in) :: rk_scheme
    real(kind=rp), intent(in) :: dt
    integer :: n, s, i, j, l
    type(field_t), pointer :: k_rho_1, k_rho_2, k_rho_3, k_rho_4, &
         k_m_x_1, k_m_x_2, k_m_x_3, k_m_x_4, &
         k_m_y_1, k_m_y_2, k_m_y_3, k_m_y_4, &
         k_m_z_1, k_m_z_2, k_m_z_3, k_m_z_4, &
         k_E_1, k_E_2, k_E_3, k_E_4, &
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
    call neko_scratch_registry%request_field(temp_p, tmp_indices(26), .false.)
    call neko_scratch_registry%request_field(temp_u, tmp_indices(27), .false.)
    call neko_scratch_registry%request_field(temp_v, tmp_indices(28), .false.)
    call neko_scratch_registry%request_field(temp_w, tmp_indices(29), .false.)
    call neko_scratch_registry%request_field(temp_ruvw, tmp_indices(30), &
         .false.)

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

    compressible_res_cpu_add_physical_flux = &
         any(mu%x .ne. 0.0_rp) .or. any(kappa%x .ne. 0.0_rp)
    compressible_res_cpu_add_physical_stress = any(mu%x .ne. 0.0_rp)

    ! Loop over Runge-Kutta stages. One parallel region per stage covers
    ! both the initial copy and all (i-1) accumulation sweeps.
    do i = 1, s
       !$omp parallel private(j, l, k_rho_ptr)  &
       !$omp& private(k_m_x_ptr, k_m_y_ptr, k_m_z_ptr, k_E_ptr)

       ! Copy current solution state to temporary arrays for this RK stage
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

       ! Accumulate previous stage contributions using RK coefficients.
       ! Each thread independently sets its private k_*_ptr aliases before
       ! the worksharing loop; the implicit barrier at end of !$omp do simd
       ! synchronises threads between j iterations.
       do j = 1, i-1
          k_rho_ptr => k_rho%items(j)%ptr%x
          k_m_x_ptr => k_m_x%items(j)%ptr%x
          k_m_y_ptr => k_m_y%items(j)%ptr%x
          k_m_z_ptr => k_m_z%items(j)%ptr%x
          k_E_ptr => k_E%items(j)%ptr%x

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
       end do
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
            Ax_stress, coef, gs, h, artificial_visc, mu, kappa)
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

    call neko_scratch_registry%relinquish_field(tmp_indices)

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
  !> @param h Mesh size field
  !> @param artificial_visc Artificial viscosity field
  !> @param mu Dynamic viscosity field for physical momentum diffusion
  !> @param kappa Thermal conductivity field for physical energy diffusion
  subroutine evaluate_rhs_cpu(rhs_rho_field, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, &
       rho_field, m_x, m_y, m_z, E, p, u, v, w, Ax, &
       Ax_stress, coef, gs, h, artificial_visc, mu, kappa)
    type(field_t), intent(inout) :: rhs_rho_field, &
         rhs_m_x, rhs_m_y, rhs_m_z, rhs_E
    type(field_t), intent(inout) :: rho_field, m_x, m_y, m_z, E
    type(field_t), intent(in) :: p, u, v, w, h, artificial_visc, mu, kappa
    class(Ax_t), intent(inout) :: Ax, Ax_stress
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: i, n
    type(field_t), pointer :: f_x, f_y, f_z
    type(field_t), pointer :: visc_rho, visc_m_x, visc_m_y, visc_m_z, visc_E
    integer :: tmp_indices(8)

    n = coef%dof%size()

    call neko_scratch_registry%request_field(f_x, tmp_indices(1), .false.)
    call neko_scratch_registry%request_field(f_y, tmp_indices(2), .false.)
    call neko_scratch_registry%request_field(f_z, tmp_indices(3), .false.)

    ! Hoisted from below so the registry calls are kept outside the
    ! following !$omp parallel do simd regions.
    ! fused into one !$omp parallel do simd (registry calls are not
    ! thread-safe so they cannot sit between two !$omp do constructs).
    call neko_scratch_registry%request_field(visc_rho, tmp_indices(4), .false.)
    call neko_scratch_registry%request_field(visc_m_x, tmp_indices(5), .false.)
    call neko_scratch_registry%request_field(visc_m_y, tmp_indices(6), .false.)
    call neko_scratch_registry%request_field(visc_m_z, tmp_indices(7), .false.)
    call neko_scratch_registry%request_field(visc_E, tmp_indices(8), .false.)

    !> rho = rho - dt * div(m)
    ! Compute density flux divergence
    call div(rhs_rho_field%x, m_x%x, m_y%x, m_z%x, coef)

    !> m = m - dt * div(rho * u * u^T + p*I)
    ! Compute momentum flux divergences
    ! m_x
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x%x(i,1,1,1) = m_x%x(i,1,1,1) * m_x%x(i,1,1,1) / &
            rho_field%x(i,1,1,1) + p%x(i,1,1,1)
       f_y%x(i,1,1,1) = m_x%x(i,1,1,1) * m_y%x(i,1,1,1) / &
            rho_field%x(i,1,1,1)
       f_z%x(i,1,1,1) = m_x%x(i,1,1,1) * m_z%x(i,1,1,1) / &
            rho_field%x(i,1,1,1)
    end do
    !$omp end parallel do simd
    call div(rhs_m_x%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_y
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x%x(i,1,1,1) = m_y%x(i,1,1,1) * m_x%x(i,1,1,1) / &
            rho_field%x(i,1,1,1)
       f_y%x(i,1,1,1) = m_y%x(i,1,1,1) * m_y%x(i,1,1,1) / &
            rho_field%x(i,1,1,1) + p%x(i,1,1,1)
       f_z%x(i,1,1,1) = m_y%x(i,1,1,1) * m_z%x(i,1,1,1) / &
            rho_field%x(i,1,1,1)
    end do
    !$omp end parallel do simd
    call div(rhs_m_y%x, f_x%x, f_y%x, f_z%x, coef)
    ! m_z
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x%x(i,1,1,1) = m_z%x(i,1,1,1) * m_x%x(i,1,1,1) / &
            rho_field%x(i,1,1,1)
       f_y%x(i,1,1,1) = m_z%x(i,1,1,1) * m_y%x(i,1,1,1) / &
            rho_field%x(i,1,1,1)
       f_z%x(i,1,1,1) = m_z%x(i,1,1,1) * m_z%x(i,1,1,1) / &
            rho_field%x(i,1,1,1) + p%x(i,1,1,1)
    end do
    !$omp end parallel do simd
    call div(rhs_m_z%x, f_x%x, f_y%x, f_z%x, coef)

    !> E = E - dt * div(u * (E + p))
    ! Compute energy flux divergence
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       f_x%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) * &
            u%x(i,1,1,1)
       f_y%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) * &
            v%x(i,1,1,1)
       f_z%x(i,1,1,1) = (E%x(i,1,1,1) + p%x(i,1,1,1)) * &
            w%x(i,1,1,1)
    end do
    !$omp end parallel do simd
    call div(rhs_E%x, f_x%x, f_y%x, f_z%x, coef)

    call gs%op(rhs_rho_field, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 1, coef)
    call gs%op(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, n, GS_OP_ADD)
    call rotate_cyc(rhs_m_x%x, rhs_m_y%x, rhs_m_z%x, 0, coef)
    call gs%op(rhs_E, GS_OP_ADD)

    ! Apply multiplicity to the inviscid RHS and set h1 to the artificial
    ! viscosity for the Laplacian (density, momentum and energy are all
    ! stabilized by the same artificial viscosity). Fused: one fork-join
    ! instead of two, and coef%mult / artificial_visc share L1.
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       rhs_rho_field%x(i,1,1,1) = rhs_rho_field%x(i,1,1,1) * &
            coef%mult(i,1,1,1)
       rhs_m_x%x(i,1,1,1) = rhs_m_x%x(i,1,1,1) * coef%mult(i,1,1,1)
       rhs_m_y%x(i,1,1,1) = rhs_m_y%x(i,1,1,1) * coef%mult(i,1,1,1)
       rhs_m_z%x(i,1,1,1) = rhs_m_z%x(i,1,1,1) * coef%mult(i,1,1,1)
       rhs_E%x(i,1,1,1) = rhs_E%x(i,1,1,1) * coef%mult(i,1,1,1)
       coef%h1(i,1,1,1) = artificial_visc%x(i,1,1,1)
    end do
    !$omp end parallel do simd

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

    ! gs. h1=1.0 reset is deferred into the final accumulation below; safe
    ! because nothing here (gs%op, rotate_cyc) reads coef%h1.
    call gs%op(visc_rho, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 1, coef)
    call gs%op(visc_m_x%x, visc_m_y%x, visc_m_z%x, n, GS_OP_ADD)
    call rotate_cyc(visc_m_x%x, visc_m_y%x, visc_m_z%x, 0, coef)
    call gs%op(visc_E, GS_OP_ADD)

    ! Move div to the rhs, apply artificial viscosity, and reset h1 to 1.
    ! The viscosity coefficient is already included in the Laplacian operator.
    ! Fused: one fork-join instead of two, and coef%Binv stays in L1.
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       rhs_rho_field%x(i,1,1,1) = -rhs_rho_field%x(i,1,1,1) - &
            coef%Binv(i,1,1,1) * visc_rho%x(i,1,1,1)
       rhs_m_x%x(i,1,1,1) = -rhs_m_x%x(i,1,1,1) - &
            coef%Binv(i,1,1,1) * visc_m_x%x(i,1,1,1)
       rhs_m_y%x(i,1,1,1) = -rhs_m_y%x(i,1,1,1) - &
            coef%Binv(i,1,1,1) * visc_m_y%x(i,1,1,1)
       rhs_m_z%x(i,1,1,1) = -rhs_m_z%x(i,1,1,1) - &
            coef%Binv(i,1,1,1) * visc_m_z%x(i,1,1,1)
       rhs_E%x(i,1,1,1) = -rhs_E%x(i,1,1,1) - &
            coef%Binv(i,1,1,1) * visc_E%x(i,1,1,1)
       coef%h1(i,1,1,1) = 1.0_rp
    end do
    !$omp end parallel do simd

    call neko_scratch_registry%relinquish_field(tmp_indices)
  end subroutine evaluate_rhs_cpu

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
