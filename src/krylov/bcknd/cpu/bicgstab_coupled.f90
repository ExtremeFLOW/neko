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
!> Provides a coupled CPU implementation of the BiCGStab method.
module bicgstab_cpld
  use num_types, only : rp, xp
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use scalar_bc_projector, only : scalar_bc_projector_t
  use vector_bc_projector, only : vector_bc_projector_t
  use host_array, only : host_array_t
  use scratch_registry, only : neko_scratch_registry
  use comm, only : NEKO_COMM, MPI_EXTRA_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_SUM
  use operators, only : rotate_cyc
  use math, only : NEKO_EPS
  use utils, only : neko_error
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  !> Coupled right-preconditioned CPU BiCGStab method.
  !!
  !! The method uses a single Krylov recurrence and combined inner products
  !! over all components. The scalar preconditioner is applied to each
  !! component independently.
  type, public, extends(ksp_t) :: bicgstab_cpld_t
   contains
     !> Initialise a coupled CPU BiCGStab solver.
     procedure, pass(this) :: init => bicgstab_cpld_init
     !> Free a coupled CPU BiCGStab solver.
     procedure, pass(this) :: free => bicgstab_cpld_free
     !> Reject scalar solves, which are not supported by this type.
     procedure, pass(this) :: solve => bicgstab_cpld_solve_scalar
     !> Solve a three-component coupled system with BiCGStab.
     procedure, pass(this) :: solve_coupled => bicgstab_cpld_solve
  end type bicgstab_cpld_t

contains

  !> Initialise a coupled CPU BiCGStab solver.
  !! @param n Number of degrees of freedom per component.
  !! @param max_iter Maximum number of iterations.
  !! @param M Optional preconditioner. An identity preconditioner is used if
  !! absent.
  !! @param rel_tol Optional relative convergence tolerance.
  !! @param abs_tol Optional absolute convergence tolerance.
  !! @param monitor Optional switch for logging the residual at each iteration.
  subroutine bicgstab_cpld_init(this, n, max_iter, M, rel_tol, abs_tol, &
       monitor)
    class(bicgstab_cpld_t), target, intent(inout) :: this
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(in), target :: M
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()

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

  end subroutine bicgstab_cpld_init

  !> Free a coupled CPU BiCGStab solver.
  subroutine bicgstab_cpld_free(this)
    class(bicgstab_cpld_t), intent(inout) :: this

    call this%ksp_free()

    nullify(this%M)

  end subroutine bicgstab_cpld_free

  !> Reject a scalar solve with a coupled BiCGStab solver.
  !! @param Ax Linear operator.
  !! @param x Solution field.
  !! @param f Right-hand side.
  !! @param n Number of degrees of freedom.
  !! @param coef Spectral element coefficients and multiplicity weights.
  !! @param bc_projector Projector for Dirichlet boundary nodes.
  !! @param gs_h Gather-scatter handle used to assemble the operator result.
  !! @param niter Optional maximum number of iterations.
  !! @return Unused convergence information.
  function bicgstab_cpld_solve_scalar(this, Ax, x, f, n, coef, &
       bc_projector, gs_h, niter) result(ksp_results)
    class(bicgstab_cpld_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    class(scalar_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    type(ksp_monitor_t) :: ksp_results

    call neko_error('Coupled BiCGStab is only defined for coupled solves')

    ksp_results%res_start = 0.0_rp
    ksp_results%res_final = 0.0_rp
    ksp_results%iter = 0
    ksp_results%converged = .false.

  end function bicgstab_cpld_solve_scalar

  !> Solve a three-component coupled system with the CPU BiCGStab method.
  !!
  !! The initial guesses are discarded. All inner products combine the three
  !! components, so the method advances one Krylov recurrence for the complete
  !! coupled system.
  !! @param Ax Coupled linear operator.
  !! @param x Solution field for the first component.
  !! @param y Solution field for the second component.
  !! @param z Solution field for the third component.
  !! @param fx Right-hand side for the first component.
  !! @param fy Right-hand side for the second component.
  !! @param fz Right-hand side for the third component.
  !! @param n Number of degrees of freedom per component.
  !! @param coef Spectral element coefficients and multiplicity weights.
  !! @param bc_projector Projector for vector Dirichlet boundary nodes.
  !! @param gs_h Gather-scatter handle used to assemble operator results.
  !! @param niter Optional maximum number of iterations, overriding the
  !! configured value.
  !! @return Identical combined convergence information for all components.
  function bicgstab_cpld_solve(this, Ax, x, y, z, fx, fy, fz, n, coef, &
       bc_projector, gs_h, niter) result(ksp_results)
    class(bicgstab_cpld_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: fx
    real(kind=rp), dimension(n), intent(in) :: fy
    real(kind=rp), dimension(n), intent(in) :: fz
    type(coef_t), intent(inout) :: coef
    class(vector_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer :: i, iter, max_iter, ierr
    real(kind=rp) :: alpha, beta, omega, rho_1, rho_2
    real(kind=rp) :: rnorm, norm_fac, gamma
    real(kind=rp) :: r_norm, s_norm, shadow_norm, t_norm, v_norm
    ! r^T r, s^T s, f^T v, v^T v, s^T t, t^T t
    real(kind=rp) :: rtr, sts, ftv, vtv, stt, ttt
    real(kind=xp) :: norm_sum

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    block
      type(host_array_t), pointer :: p_tmp, p_hat_tmp, r_tmp
      type(host_array_t), pointer :: s_hat_tmp, t_tmp, v_tmp
      real(kind=rp), pointer :: p(:, :), p_hat(:, :), r(:, :)
      real(kind=rp), pointer :: s_hat(:, :), t(:, :), v(:, :)
      integer :: temp_indices(6)

      call neko_scratch_registry%request_host_array(p_tmp, temp_indices(1), &
           3 * n, .false.)
      call neko_scratch_registry%request_host_array(p_hat_tmp, &
           temp_indices(2), 3 * n, .false.)
      call neko_scratch_registry%request_host_array(r_tmp, temp_indices(3), &
           3 * n, .false.)
      call neko_scratch_registry%request_host_array(s_hat_tmp, &
           temp_indices(4), 3 * n, .false.)
      call neko_scratch_registry%request_host_array(t_tmp, temp_indices(5), &
           3 * n, .false.)
      call neko_scratch_registry%request_host_array(v_tmp, temp_indices(6), &
           3 * n, .false.)

      p(1:n, 1:3) => p_tmp%x
      p_hat(1:n, 1:3) => p_hat_tmp%x
      r(1:n, 1:3) => r_tmp%x
      s_hat(1:n, 1:3) => s_hat_tmp%x
      t(1:n, 1:3) => t_tmp%x
      v(1:n, 1:3) => v_tmp%x

      ! BiCGStab starts from zero. The right-hand side is consequently both
      ! the initial residual and the fixed shadow residual.
      norm_sum = 0.0_xp
      !$omp parallel do reduction(+:norm_sum)
      do i = 1, n
         x%x(i, 1, 1, 1) = 0.0_rp
         y%x(i, 1, 1, 1) = 0.0_rp
         z%x(i, 1, 1, 1) = 0.0_rp
         r(i, 1) = fx(i)
         r(i, 2) = fy(i)
         r(i, 3) = fz(i)
         norm_sum = norm_sum + coef%mult(i, 1, 1, 1) * &
              (r(i, 1)**2 + r(i, 2)**2 + r(i, 3)**2)
      end do
      !$omp end parallel do
      call MPI_Allreduce(MPI_IN_PLACE, norm_sum, 1, MPI_EXTRA_PRECISION, &
           MPI_SUM, NEKO_COMM, ierr)

      rtr = norm_sum
      r_norm = bicgstab_cpld_sqrt(rtr, 'initial residual')
      shadow_norm = r_norm
      rnorm = r_norm * norm_fac
      gamma = rnorm * this%rel_tol
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0

      ! Besides saving an iteration, the early check prevents a converged
      ! right-hand side from being reported as a rho breakdown.
      if (r_norm .le. 0.0_rp .or. rnorm .lt. this%abs_tol .or. &
           rnorm .lt. gamma) then
         ksp_results%converged = .true.
         call neko_scratch_registry%relinquish_host_array(temp_indices)
         return
      end if

      call this%monitor_start('Coupled BiCGStab')
      do iter = 1, max_iter
         call bicgstab_cpld_product(rho_1, fx, fy, fz, r(:, 1), r(:, 2), &
              r(:, 3), coef%mult, n)

         ! The combined product-space norms make the breakdown decision
         ! invariant under a uniform scaling of the complete vector system.
         call bicgstab_cpld_check_inner_product(rho_1, shadow_norm, r_norm, &
              'rho inner product')

         if (iter .eq. 1) then
            !$omp parallel do
            do i = 1, n
               p(i, 1) = r(i, 1)
               p(i, 2) = r(i, 2)
               p(i, 3) = r(i, 3)
            end do
            !$omp end parallel do
         else
            beta = (rho_1 / rho_2) * (alpha / omega)
            if (.not. ieee_is_finite(beta)) then
               call neko_error('Coupled BiCGStab failure: non-finite beta')
            end if

            !$omp parallel do
            do i = 1, n
               p(i, 1) = r(i, 1) + &
                    beta * (p(i, 1) - omega * v(i, 1))
               p(i, 2) = r(i, 2) + &
                    beta * (p(i, 2) - omega * v(i, 2))
               p(i, 3) = r(i, 3) + &
                    beta * (p(i, 3) - omega * v(i, 3))
            end do
            !$omp end parallel do
         end if

         ! The preconditioner interface is scalar, so apply the same
         ! preconditioner separately to all three components.
         call this%M%solve(p_hat(:, 1), p(:, 1), n)
         call this%M%solve(p_hat(:, 2), p(:, 2), n)
         call this%M%solve(p_hat(:, 3), p(:, 3), n)

         call Ax%compute_vector(v(:, 1), v(:, 2), v(:, 3), p_hat(:, 1), &
              p_hat(:, 2), p_hat(:, 3), coef, x%msh, x%Xh)
         call bicgstab_cpld_assemble(v, n, coef, bc_projector, gs_h)

         ! Compute f^T v and v^T v in a single global reduction. The latter
         ! supplies the scale for the alpha-denominator breakdown check.
         call bicgstab_cpld_product_and_norm(ftv, vtv, fx, fy, fz, &
              v(:, 1), v(:, 2), v(:, 3), coef%mult, n)
         v_norm = bicgstab_cpld_sqrt(vtv, 'operator result v')
         call bicgstab_cpld_check_inner_product(ftv, shadow_norm, v_norm, &
              'alpha denominator')
         alpha = rho_1 / ftv
         if (.not. ieee_is_finite(alpha)) then
            call neko_error('Coupled BiCGStab failure: non-finite alpha')
         end if

         ! Store the intermediate residual in r and compute its combined norm
         ! in the same pass. The previous residual is no longer needed after
         ! p has been formed.
         norm_sum = 0.0_xp
         !$omp parallel do reduction(+:norm_sum)
         do i = 1, n
            r(i, 1) = r(i, 1) - alpha * v(i, 1)
            r(i, 2) = r(i, 2) - alpha * v(i, 2)
            r(i, 3) = r(i, 3) - alpha * v(i, 3)
            norm_sum = norm_sum + coef%mult(i, 1, 1, 1) * &
                 (r(i, 1)**2 + r(i, 2)**2 + r(i, 3)**2)
         end do
         !$omp end parallel do
         call MPI_Allreduce(MPI_IN_PLACE, norm_sum, 1, MPI_EXTRA_PRECISION, &
              MPI_SUM, NEKO_COMM, ierr)

         sts = norm_sum
         s_norm = bicgstab_cpld_sqrt(sts, 'intermediate residual')
         rnorm = s_norm * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            !$omp parallel do
            do i = 1, n
               x%x(i, 1, 1, 1) = x%x(i, 1, 1, 1) + alpha * p_hat(i, 1)
               y%x(i, 1, 1, 1) = y%x(i, 1, 1, 1) + alpha * p_hat(i, 2)
               z%x(i, 1, 1, 1) = z%x(i, 1, 1, 1) + alpha * p_hat(i, 3)
            end do
            !$omp end parallel do
            call this%monitor_iter(iter, rnorm)
            exit
         end if

         call this%M%solve(s_hat(:, 1), r(:, 1), n)
         call this%M%solve(s_hat(:, 2), r(:, 2), n)
         call this%M%solve(s_hat(:, 3), r(:, 3), n)

         call Ax%compute_vector(t(:, 1), t(:, 2), t(:, 3), s_hat(:, 1), &
              s_hat(:, 2), s_hat(:, 3), coef, x%msh, x%Xh)
         call bicgstab_cpld_assemble(t, n, coef, bc_projector, gs_h)

         ! The numerator and denominator of omega share one reduction.
         call bicgstab_cpld_product_and_norm(stt, ttt, r(:, 1), r(:, 2), &
              r(:, 3), t(:, 1), t(:, 2), t(:, 3), coef%mult, n)
         t_norm = bicgstab_cpld_sqrt(ttt, 'operator result t')
         if (t_norm .le. 0.0_rp) then
            call neko_error(&
                 'Coupled BiCGStab breakdown: zero omega denominator')
         end if
         if (.not. ieee_is_finite(stt)) then
            call neko_error(&
                 'Coupled BiCGStab failure: non-finite omega numerator')
         end if
         omega = stt / ttt
         if (.not. ieee_is_finite(omega)) then
            call neko_error('Coupled BiCGStab failure: non-finite omega')
         end if

         ! Update the solution, recursive residual, and residual norm together
         ! to avoid another traversal of all three component vectors.
         norm_sum = 0.0_xp
         !$omp parallel do reduction(+:norm_sum)
         do i = 1, n
            x%x(i, 1, 1, 1) = x%x(i, 1, 1, 1) + alpha * p_hat(i, 1) + &
                 omega * s_hat(i, 1)
            y%x(i, 1, 1, 1) = y%x(i, 1, 1, 1) + alpha * p_hat(i, 2) + &
                 omega * s_hat(i, 2)
            z%x(i, 1, 1, 1) = z%x(i, 1, 1, 1) + alpha * p_hat(i, 3) + &
                 omega * s_hat(i, 3)
            r(i, 1) = r(i, 1) - omega * t(i, 1)
            r(i, 2) = r(i, 2) - omega * t(i, 2)
            r(i, 3) = r(i, 3) - omega * t(i, 3)
            norm_sum = norm_sum + coef%mult(i, 1, 1, 1) * &
                 (r(i, 1)**2 + r(i, 2)**2 + r(i, 3)**2)
         end do
         !$omp end parallel do
         call MPI_Allreduce(MPI_IN_PLACE, norm_sum, 1, MPI_EXTRA_PRECISION, &
              MPI_SUM, NEKO_COMM, ierr)

         rtr = norm_sum
         r_norm = bicgstab_cpld_sqrt(rtr, 'recursive residual')
         rnorm = r_norm * norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            exit
         end if

         ! A negative omega is valid. Only numerical orthogonality between s
         ! and t constitutes a breakdown of the minimal-residual step.
         call bicgstab_cpld_check_inner_product(stt, s_norm, t_norm, &
              'omega numerator')
         rho_2 = rho_1
      end do

      call this%monitor_stop()
      ksp_results%res_final = rnorm
      ksp_results%iter = iter
      ksp_results%converged = this%is_converged(iter, rnorm)
      call neko_scratch_registry%relinquish_host_array(temp_indices)
    end block

  end function bicgstab_cpld_solve

  !> Check an inner product for a coupled BiCGStab breakdown.
  !! @param inner_product Combined inner product of the two vectors.
  !! @param norm_a Combined norm of the first vector.
  !! @param norm_b Combined norm of the second vector.
  !! @param quantity Name of the algorithmic quantity for error reporting.
  subroutine bicgstab_cpld_check_inner_product(inner_product, norm_a, norm_b, &
       quantity)
    real(kind=rp), intent(in) :: inner_product
    real(kind=rp), intent(in) :: norm_a
    real(kind=rp), intent(in) :: norm_b
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: large_norm, small_norm

    if (.not. ieee_is_finite(inner_product)) then
       call neko_error('Coupled BiCGStab failure: non-finite ' // &
            trim(quantity))
    end if

    ! Divide by the larger norm first to avoid overflow in their product.
    large_norm = max(norm_a, norm_b)
    small_norm = min(norm_a, norm_b)
    if (large_norm .le. 0.0_rp .or. &
         abs(inner_product) / large_norm .le. NEKO_EPS * small_norm) then
       call neko_error('Coupled BiCGStab breakdown: near-zero ' // &
            trim(quantity))
    end if

  end subroutine bicgstab_cpld_check_inner_product

  !> Return the square root of a valid combined squared norm.
  !! @param value Combined weighted squared norm.
  !! @param quantity Name of the vector for error reporting.
  !! @return The non-negative square root.
  function bicgstab_cpld_sqrt(value, quantity) result(root)
    real(kind=rp), intent(in) :: value
    character(len=*), intent(in) :: quantity
    real(kind=rp) :: root

    if (.not. ieee_is_finite(value) .or. value .lt. 0.0_rp) then
       call neko_error('Coupled BiCGStab failure: invalid ' // &
            trim(quantity) // ' norm')
    end if
    root = sqrt(value)

  end function bicgstab_cpld_sqrt

  !> Assemble a coupled operator result and project its boundary data.
  !! @param vector Three-component operator result.
  !! @param n Number of degrees of freedom per component.
  !! @param coef Spectral element coefficients used for cyclic rotation.
  !! @param bc_projector Projector for vector Dirichlet boundary nodes.
  !! @param gs_h Gather-scatter handle used to assemble the vector.
  subroutine bicgstab_cpld_assemble(vector, n, coef, bc_projector, gs_h)
    integer, intent(in) :: n
    real(kind=rp), dimension(n, 3), intent(inout) :: vector
    type(coef_t), intent(inout) :: coef
    class(vector_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h

    ! Cyclic faces must be expressed in their common coordinate system while
    ! the three components are gathered, then rotated back to physical space.
    call rotate_cyc(vector(:, 1), vector(:, 2), vector(:, 3), 1, coef)
    call gs_h%op(vector(:, 1), vector(:, 2), vector(:, 3), n, GS_OP_ADD)
    call rotate_cyc(vector(:, 1), vector(:, 2), vector(:, 3), 0, coef)

    call bc_projector%apply(vector(:, 1), vector(:, 2), vector(:, 3), n)

  end subroutine bicgstab_cpld_assemble

  !> Compute a coupled weighted inner product in one global reduction.
  !! @param product Combined weighted inner product.
  !! @param ax First component of the first vector.
  !! @param ay Second component of the first vector.
  !! @param az Third component of the first vector.
  !! @param bx First component of the second vector.
  !! @param by Second component of the second vector.
  !! @param bz Third component of the second vector.
  !! @param mult Spectral element multiplicity weights.
  !! @param n Number of degrees of freedom per component.
  subroutine bicgstab_cpld_product(product, ax, ay, az, bx, by, bz, mult, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: product
    real(kind=rp), dimension(n), intent(in) :: ax, ay, az
    real(kind=rp), dimension(n), intent(in) :: bx, by, bz
    real(kind=rp), dimension(n), intent(in) :: mult
    real(kind=xp) :: product_sum
    integer :: i, ierr

    product_sum = 0.0_xp
    !$omp parallel do reduction(+:product_sum)
    do i = 1, n
       product_sum = product_sum + mult(i) * &
            (ax(i) * bx(i) + ay(i) * by(i) + az(i) * bz(i))
    end do
    !$omp end parallel do

    call MPI_Allreduce(MPI_IN_PLACE, product_sum, 1, MPI_EXTRA_PRECISION, &
         MPI_SUM, NEKO_COMM, ierr)
    product = product_sum

  end subroutine bicgstab_cpld_product

  !> Compute a coupled inner product and squared norm in one reduction.
  !! @param product Combined weighted inner product of `a` and `b`.
  !! @param norm_squared Combined weighted squared norm of `b`.
  !! @param ax First component of `a`.
  !! @param ay Second component of `a`.
  !! @param az Third component of `a`.
  !! @param bx First component of `b`.
  !! @param by Second component of `b`.
  !! @param bz Third component of `b`.
  !! @param mult Spectral element multiplicity weights.
  !! @param n Number of degrees of freedom per component.
  subroutine bicgstab_cpld_product_and_norm(product, norm_squared, ax, ay, &
       az, bx, by, bz, mult, n)
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: product
    real(kind=rp), intent(out) :: norm_squared
    real(kind=rp), dimension(n), intent(in) :: ax, ay, az
    real(kind=rp), dimension(n), intent(in) :: bx, by, bz
    real(kind=rp), dimension(n), intent(in) :: mult
    real(kind=xp) :: product_sum, norm_sum
    real(kind=xp) :: reductions(2)
    integer :: i, ierr

    product_sum = 0.0_xp
    norm_sum = 0.0_xp
    !$omp parallel do reduction(+:product_sum,norm_sum)
    do i = 1, n
       product_sum = product_sum + mult(i) * &
            (ax(i) * bx(i) + ay(i) * by(i) + az(i) * bz(i))
       norm_sum = norm_sum + mult(i) * &
            (bx(i)**2 + by(i)**2 + bz(i)**2)
    end do
    !$omp end parallel do

    reductions(1) = product_sum
    reductions(2) = norm_sum
    call MPI_Allreduce(MPI_IN_PLACE, reductions, 2, MPI_EXTRA_PRECISION, &
         MPI_SUM, NEKO_COMM, ierr)
    product = reductions(1)
    norm_squared = reductions(2)

  end subroutine bicgstab_cpld_product_and_norm

end module bicgstab_cpld
