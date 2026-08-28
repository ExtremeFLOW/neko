! Copyright (c) 2021-2026, The Neko Authors
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
!> Provides a CPU implementation of the BiCGStab method.
module bicgstab
  use num_types, only: rp, xp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use math, only : glsc3, copy, NEKO_EPS, add2s2, p_update
  use utils, only : neko_error
  use comm, only : NEKO_COMM, MPI_EXTRA_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_SUM
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  !> CPU implementation of the right-preconditioned BiCGStab method.
  !!
  !! The solver starts from a zero initial guess and uses the preconditioner
  !! provided by [ksp_t](#krylov::ksp_t). The coupled interface solves the
  !! three components independently and does not apply a coupled operator.
  type, public, extends(ksp_t) :: bicgstab_t
     !> Search direction \f$p\f$.
     real(kind=rp), allocatable :: p(:)
     !> Preconditioned search direction \f$\hat{p} = M^{-1}p\f$.
     real(kind=rp), allocatable :: p_hat(:)
     !> Residual \f$r\f$.
     real(kind=rp), allocatable :: r(:)
     !> Intermediate residual \f$s = r - \alpha v\f$.
     real(kind=rp), allocatable :: s(:)
     !> Preconditioned intermediate residual \f$\hat{s} = M^{-1}s\f$.
     real(kind=rp), allocatable :: s_hat(:)
     !> Operator action \f$t = A\hat{s}\f$.
     real(kind=rp), allocatable :: t(:)
     !> Operator action \f$v = A\hat{p}\f$.
     real(kind=rp), allocatable :: v(:)
   contains
     !> Initialise a CPU BiCGStab solver.
     procedure, pass(this) :: init => bicgstab_init
     !> Free a CPU BiCGStab solver.
     procedure, pass(this) :: free => bicgstab_free
     !> Solve a linear system with the CPU BiCGStab method.
     procedure, pass(this) :: solve => bicgstab_solve
     !> Solve three independent systems with the CPU BiCGStab method.
     procedure, pass(this) :: solve_coupled => bicgstab_solve_coupled
  end type bicgstab_t

contains

  !> Initialise a CPU BiCGStab solver.
  !! @param n Number of degrees of freedom.
  !! @param max_iter Maximum number of iterations.
  !! @param M Optional preconditioner. An identity preconditioner is used if
  !! absent.
  !! @param rel_tol Optional relative convergence tolerance.
  !! @param abs_tol Optional absolute convergence tolerance.
  !! @param monitor Optional switch for logging the residual at each iteration.
  subroutine bicgstab_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(bicgstab_t), target, intent(inout) :: this
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor


    call this%free()

    allocate(this%p(n))
    allocate(this%p_hat(n))
    allocate(this%r(n))
    allocate(this%s(n))
    allocate(this%s_hat(n))
    allocate(this%t(n))
    allocate(this%v(n))
    if (present(M)) then
       this%M => M
    end if

    if (present(rel_tol) .and. present(abs_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol)
    else if (present(monitor) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, monitor = monitor)
    else if (present(rel_tol)) then
       call this%ksp_init(max_iter, rel_tol = rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol)
    else if (present(monitor)) then
       call this%ksp_init(max_iter, monitor = monitor)
    else
       call this%ksp_init(max_iter)
    end if

  end subroutine bicgstab_init

  !> Free a CPU BiCGStab solver.
  subroutine bicgstab_free(this)
    class(bicgstab_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%v)) then
       deallocate(this%v)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%t)) then
       deallocate(this%t)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if

    if (allocated(this%p_hat)) then
       deallocate(this%p_hat)
    end if

    if (allocated(this%s)) then
       deallocate(this%s)
    end if

    if (allocated(this%s_hat)) then
       deallocate(this%s_hat)
    end if

    nullify(this%M)


  end subroutine bicgstab_free

  !> Solve a linear system with the CPU BiCGStab method.
  !!
  !! The initial guess in `x` is discarded and the iteration starts from zero.
  !! @param Ax Linear operator.
  !! @param x Solution field.
  !! @param f Right-hand side.
  !! @param n Number of degrees of freedom.
  !! @param coef Spectral element coefficients and multiplicity weights.
  !! @param blst Boundary conditions applied to the operator result.
  !! @param gs_h Gather-scatter handle used to assemble the operator result.
  !! @param niter Optional maximum number of iterations, overriding the
  !! configured value.
  !! @return Convergence information for the solve.
  function bicgstab_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) &
       result(ksp_results)
    class(bicgstab_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, i, ierr
    real(kind=rp) :: rnorm, rtr, norm_fac, gamma
    real(kind=rp) :: r_norm, s_norm, shadow_norm, t_norm, v_norm
    ! s^T s, f^T v, v^T v, s^T t, t^T t
    real(kind=rp) :: sts, ftv, vtv, stt, ttt
    real(kind=rp) :: beta, alpha, omega, rho_1, rho_2
    ! Extra-precision accumulator for the fused residual reductions
    real(kind=xp) :: res_sum

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(r => this%r, t => this%t, s => this%s, v => this%v, &
         p => this%p, s_hat => this%s_hat, p_hat => this%p_hat)

      res_sum = 0.0_xp
      !$omp parallel do reduction(+:res_sum)
      do i = 1, n
         x%x(i,1,1,1) = 0.0_rp
         r(i) = f(i)
         res_sum = res_sum + (r(i) * coef%mult(i,1,1,1) * r(i))
      end do
      !$omp end parallel do

      call MPI_Allreduce(MPI_IN_PLACE, res_sum, 1, &
           MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      rtr = res_sum

      ! The implementation deliberately starts from x = 0, so f is both the
      ! initial residual and the fixed shadow residual used by BiCGStab.
      r_norm = bicgstab_sqrt(rtr, 'initial residual')
      shadow_norm = r_norm
      rnorm = r_norm * norm_fac
      gamma = rnorm * this%rel_tol
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0

      ! Avoid entering the recurrence if the zero initial guess already meets
      ! either stopping criterion. Apart from saving work, this prevents a
      ! small right-hand side from being mistaken for a rho breakdown.
      if (r_norm .le. 0.0_rp .or. rnorm .lt. this%abs_tol .or. &
           rnorm .lt. gamma) then
         ksp_results%converged = .true.
         return
      end if

      call this%monitor_start('BiCGStab')
      do iter = 1, max_iter

         rho_1 = glsc3(f, coef%mult, r, n)

         ! BiCGStab breaks down if the residual becomes orthogonal to the
         ! shadow residual. Test the inner product relative to the norms of
         ! both vectors so that uniformly scaling the system does not change
         ! the decision.
         call bicgstab_check_inner_product(rho_1, shadow_norm, r_norm, &
              'rho inner product')

         if (iter .eq. 1) then
            call copy(p, r, n)
         else
            beta = (rho_1 / rho_2) * (alpha / omega)
            if (.not. ieee_is_finite(beta)) then
               call neko_error('BiCGStab failure: non-finite beta')
            end if
            call p_update(p, r, v, beta, omega, n)
         end if

         call this%M%solve(p_hat, p, n)
         call Ax%compute(v, p_hat, coef, x%msh, x%Xh)
         call gs_h%op(v, n, GS_OP_ADD)
         call blst%apply(v, n)

         ! The alpha denominator is another BiCG breakdown point. Computing it
         ! together with ||v|| permits a scale-aware orthogonality check without
         ! an additional global synchronisation.
         call bicgstab_product_and_norm(ftv, vtv, f, v, coef%mult, n)
         v_norm = bicgstab_sqrt(vtv, 'operator result v')
         call bicgstab_check_inner_product(ftv, shadow_norm, &
              v_norm, 'alpha denominator')
         alpha = rho_1 / ftv
         if (.not. ieee_is_finite(alpha)) then
            call neko_error('BiCGStab failure: non-finite alpha')
         end if

         res_sum = 0.0_xp
         !$omp parallel do reduction(+:res_sum)
         do i = 1, n
            s(i) = r(i) - alpha * v(i)
            res_sum = res_sum + s(i) * coef%mult(i,1,1,1) * s(i)
         end do
         !$omp end parallel do

         call MPI_Allreduce(MPI_IN_PLACE, res_sum, 1, &
              MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
         sts = res_sum

         s_norm = bicgstab_sqrt(sts, 'intermediate residual')
         rnorm = s_norm * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            call add2s2(x%x, p_hat, alpha, n)
            call this%monitor_iter(iter, rnorm)
            exit
         end if

         call this%M%solve(s_hat, s, n)
         call Ax%compute(t, s_hat, coef, x%msh, x%Xh)
         call gs_h%op(t, n, GS_OP_ADD)
         call blst%apply(t, n)

         call bicgstab_product_and_norm(stt, ttt, s, t, coef%mult, n)
         t_norm = bicgstab_sqrt(ttt, 'operator result t')
         if (t_norm .le. 0.0_rp) then
            call neko_error('BiCGStab breakdown: zero omega denominator')
         end if
         if (.not. ieee_is_finite(stt)) then
            call neko_error('BiCGStab failure: non-finite omega numerator')
         end if
         omega = stt / ttt
         if (.not. ieee_is_finite(omega)) then
            call neko_error('BiCGStab failure: non-finite omega')
         end if

         res_sum = 0.0_xp
         !$omp parallel do reduction(+:res_sum)
         do i = 1, n
            x%x(i,1,1,1) = x%x(i,1,1,1) + alpha * p_hat(i) + omega * s_hat(i)
            r(i) = s(i) - omega * t(i)
            res_sum = res_sum + r(i) * coef%mult(i,1,1,1) * r(i)
         end do
         !$omp end parallel do

         call MPI_Allreduce(MPI_IN_PLACE, res_sum, 1, &
              MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
         rtr = res_sum

         r_norm = bicgstab_sqrt(rtr, 'recursive residual')
         rnorm = r_norm * norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            exit
         end if

         ! A negative omega is valid. Breakdown occurs only when the local
         ! minimal-residual step stagnates, i.e. when t and s are numerically
         ! orthogonal.
         call bicgstab_check_inner_product(stt, t_norm, s_norm, &
              'omega numerator')
         rho_2 = rho_1

      end do
    end associate
    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)
  end function bicgstab_solve

  !> Check an inner product for a BiCGStab breakdown.
  !!
  !! The comparison is normalised by both vector norms. Dividing by the larger
  !! norm first avoids overflow in their product and preserves scale invariance.
  !! @param inner_product Weighted inner product of the two vectors.
  !! @param norm_a Norm of the first vector.
  !! @param norm_b Norm of the second vector.
  !! @param quantity Name of the algorithmic quantity for error reporting.
  subroutine bicgstab_check_inner_product(inner_product, norm_a, norm_b, &
       quantity)
    real(kind=rp), intent(in) :: inner_product
    real(kind=rp), intent(in) :: norm_a
    real(kind=rp), intent(in) :: norm_b
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: large_norm, small_norm

    if (.not. ieee_is_finite(inner_product)) then
       call neko_error('BiCGStab failure: non-finite ' // trim(quantity))
    end if

    large_norm = max(norm_a, norm_b)
    small_norm = min(norm_a, norm_b)
    if (large_norm .le. 0.0_rp .or. &
         abs(inner_product) / large_norm .le. NEKO_EPS * small_norm) then
       call neko_error('BiCGStab breakdown: near-zero ' // trim(quantity))
    end if

  end subroutine bicgstab_check_inner_product

  !> Compute a weighted inner product and squared norm in one reduction.
  !!
  !! Combining the two values avoids adding a global synchronisation solely
  !! for the scale used by the breakdown checks.
  !! @param product Weighted inner product of `a` and `b`.
  !! @param norm_squared Weighted squared norm of `b`.
  !! @param a First vector in the inner product.
  !! @param b Second vector in the inner product and vector whose norm is taken.
  !! @param mult Spectral element multiplicity weights.
  !! @param n Number of vector entries.
  subroutine bicgstab_product_and_norm(product, norm_squared, a, b, mult, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: product
    real(kind=rp), intent(out) :: norm_squared
    real(kind=rp), dimension(n), intent(in) :: a
    real(kind=rp), dimension(n), intent(in) :: b
    real(kind=rp), dimension(n), intent(in) :: mult
    real(kind=xp) :: product_sum, norm_sum
    real(kind=xp) :: reductions(2)
    integer :: i, ierr

    product_sum = 0.0_xp
    norm_sum = 0.0_xp
    !$omp parallel do reduction(+:product_sum,norm_sum)
    do i = 1, n
       product_sum = product_sum + a(i) * mult(i) * b(i)
       norm_sum = norm_sum + b(i) * mult(i) * b(i)
    end do
    !$omp end parallel do

    reductions(1) = product_sum
    reductions(2) = norm_sum
    call MPI_Allreduce(MPI_IN_PLACE, reductions, 2, MPI_EXTRA_PRECISION, &
         MPI_SUM, NEKO_COMM, ierr)
    product = reductions(1)
    norm_squared = reductions(2)

  end subroutine bicgstab_product_and_norm

  !> Return the square root of a valid squared norm.
  !! @param value Weighted squared norm.
  !! @param quantity Name of the vector for error reporting.
  !! @return The non-negative square root.
  function bicgstab_sqrt(value, quantity) result(root)
    real(kind=rp), intent(in) :: value
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: root

    if (.not. ieee_is_finite(value) .or. value .lt. 0.0_rp) then
       call neko_error('BiCGStab failure: invalid ' // trim(quantity) // &
            ' norm')
    end if
    root = sqrt(value)

  end function bicgstab_sqrt

  !> Solve three independent systems with the CPU BiCGStab method.
  !!
  !! This routine sequentially invokes the scalar solver for each component.
  !! @param Ax Linear operator.
  !! @param x Solution field for the first component.
  !! @param y Solution field for the second component.
  !! @param z Solution field for the third component.
  !! @param fx Right-hand side for the first component.
  !! @param fy Right-hand side for the second component.
  !! @param fz Right-hand side for the third component.
  !! @param n Number of degrees of freedom per component.
  !! @param coef Spectral element coefficients and multiplicity weights.
  !! @param blstx Boundary conditions for the first component.
  !! @param blsty Boundary conditions for the second component.
  !! @param blstz Boundary conditions for the third component.
  !! @param gs_h Gather-scatter handle used to assemble operator results.
  !! @param niter Optional maximum number of iterations, overriding the
  !! configured value.
  !! @return Convergence information for each component.
  function bicgstab_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(bicgstab_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: fx
    real(kind=rp), dimension(n), intent(in) :: fy
    real(kind=rp), dimension(n), intent(in) :: fz
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx
    type(bc_list_t), intent(inout) :: blsty
    type(bc_list_t), intent(inout) :: blstz
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer, optional, intent(in) :: niter

    ksp_results(1) = this%solve(Ax, x, fx, n, coef, blstx, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, blsty, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, blstz, gs_h, niter)

  end function bicgstab_solve_coupled

end module bicgstab
