! Copyright (c) 2020-2026, The Neko Authors
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
!> Defines various Conjugate Gradient methods
module cg
  use neko_config, only : NEKO_BLK_SIZE
  use num_types, only : rp, xp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use scalar_bc_projector, only : scalar_bc_projector_t
  use vector_bc_projector, only : vector_bc_projector_t, &
       vector_bc_projector_components
  use host_array, only : host_array_t
  use scratch_registry, only : neko_scratch_registry
  use math, only : glsc3, abscmp
  use comm, only : MPI_EXTRA_PRECISION, MPI_REAL_PRECISION, NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_SUM
  implicit none
  private

  integer, parameter :: CG_P_SPACE = 7

  !> CPU implementation of the preconditioned conjugate gradient method.
  !!
  !! Workspace pointers are associated with host arrays from the scratch
  !! registry only for the duration of a solve.
  type, public, extends(ksp_t) :: cg_t
     !> Operator action \f$w = A p\f$.
     real(kind=rp), pointer :: w(:) => null()
     !> Residual \f$r = f - A x\f$.
     real(kind=rp), pointer :: r(:) => null()
     !> Rolling space of search directions \f$p\f$.
     real(kind=rp), pointer :: p(:,:) => null()
     !> Preconditioned residual \f$z = M^{-1} r\f$.
     real(kind=rp), pointer :: z(:) => null()
     !> Step lengths associated with the stored search directions.
     real(kind=rp), pointer :: alpha(:) => null()
   contains
     !> Initialise a CPU PCG solver.
     procedure, pass(this) :: init => cg_init
     !> Free a CPU PCG solver.
     procedure, pass(this) :: free => cg_free
     !> Solve a linear system with the CPU PCG method.
     procedure, pass(this) :: solve => cg_solve
     !> Solve three independent systems with the CPU PCG method.
     procedure, pass(this) :: solve_coupled => cg_solve_coupled
  end type cg_t

contains

  !> Initialise a CPU PCG solver.
  subroutine cg_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(cg_t), intent(inout), target :: this
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
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

  end subroutine cg_init

  !> Free a CPU PCG solver.
  subroutine cg_free(this)
    class(cg_t), intent(inout) :: this

    call this%ksp_free()

    nullify(this%w)
    nullify(this%r)
    nullify(this%p)
    nullify(this%z)
    nullify(this%alpha)

    nullify(this%M)

  end subroutine cg_free

  !> Solve a linear system with the CPU PCG method.
  function cg_solve(this, Ax, x, f, n, coef, bc_projector, gs_h, niter) &
       result(ksp_results)
    class(cg_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    class(scalar_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, i, j, k, p_cur, p_prev, ierr
    real(kind=rp) :: rnorm, rtr, rtz2, rtz1, x_plus(NEKO_BLK_SIZE)
    real(kind=rp) :: beta, pap, norm_fac, tmp
    type(host_array_t), pointer :: w_tmp, r_tmp, p_tmp, z_tmp, alpha_tmp
    integer :: temp_indices(5)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    call neko_scratch_registry%request_host_array(w_tmp, temp_indices(1), &
         n, .false.)
    call neko_scratch_registry%request_host_array(r_tmp, temp_indices(2), &
         n, .false.)
    call neko_scratch_registry%request_host_array(p_tmp, temp_indices(3), &
         n * CG_P_SPACE, .false.)
    call neko_scratch_registry%request_host_array(z_tmp, temp_indices(4), &
         n, .false.)
    call neko_scratch_registry%request_host_array(alpha_tmp, &
         temp_indices(5), CG_P_SPACE, .false.)

    this%w => w_tmp%x
    this%r => r_tmp%x
    this%p(1:n, 1:CG_P_SPACE) => p_tmp%x
    this%z => z_tmp%x
    this%alpha => alpha_tmp%x

    associate(w => this%w, r => this%r, p => this%p, &
         z => this%z, alpha => this%alpha)

      rtz1 = 1.0_rp
      rtr = 0.0_rp
      !$omp parallel do reduction(+:rtr)
      do i = 1, n
         x%x(i,1,1,1) = 0.0_rp
         p(i, CG_P_SPACE) = 0.0_rp
         r(i) = f(i)
         rtr = rtr + (r(i) * coef%mult(i,1,1,1) * r(i))
      end do
      !$omp end parallel do

      call MPI_Allreduce(MPI_IN_PLACE, rtr, 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      rnorm = sqrt(rtr) * norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if (abscmp(rnorm, 0.0_rp)) then
         ksp_results%converged = .true.
         nullify(this%w, this%r, this%p, this%z, this%alpha)
         call neko_scratch_registry%relinquish_host_array(temp_indices)
         return
      end if

      p_prev = CG_P_SPACE
      p_cur = 1
      call this%monitor_start('CG')
      do iter = 1, max_iter
         call this%M%solve(z, r, n)
         rtz2 = rtz1
         rtz1 = glsc3(r, coef%mult, z, n)

         beta = rtz1 / rtz2
         if (iter .eq. 1) beta = 0.0_rp
         !$omp parallel do
         do i = 1, n
            p(i, p_cur) = z(i) + beta * p(i, p_prev)
         end do
         !$omp end parallel do

         call Ax%compute(w, p(1, p_cur), coef, x%msh, x%Xh)
         call gs_h%op(w, n, GS_OP_ADD)
         call bc_projector%apply(w, n)

         pap = glsc3(w, coef%mult, p(1, p_cur), n)

         alpha(p_cur) = rtz1 / pap
         call second_cg_part(rtr, r, coef%mult, w, alpha(p_cur), n)
         rnorm = sqrt(rtr) * norm_fac
         call this%monitor_iter(iter, rnorm)

         if ((p_cur .eq. CG_P_SPACE) .or. &
              (rnorm .lt. this%abs_tol) .or. iter .eq. max_iter) then
            !$omp parallel do private(j, k, x_plus, tmp)
            do i = 0, n-1, NEKO_BLK_SIZE
               if (i + NEKO_BLK_SIZE .le. n) then
                  !$omp simd
                  do k = 1, NEKO_BLK_SIZE
                     x_plus(k) = 0.0_rp
                  end do
                  do j = 1, p_cur
                     !$omp simd
                     do k = 1, NEKO_BLK_SIZE
                        x_plus(k) = x_plus(k) + alpha(j) * p(i+k,j)
                     end do
                  end do
                  !$omp simd
                  do k = 1, NEKO_BLK_SIZE
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                  end do
               else
                  do k = 1, n - i
                     tmp = 0.0_rp
                     do j = 1, p_cur
                        tmp = tmp + alpha(j) * p(i+k,j)
                     end do
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + tmp
                  end do
               end if
            end do
            !$omp end parallel do
            p_prev = p_cur
            p_cur = 1
            if (rnorm .lt. this%abs_tol) exit
         else
            p_prev = p_cur
            p_cur = p_cur + 1
         end if
      end do
    end associate
    nullify(this%w, this%r, this%p, this%z, this%alpha)
    call neko_scratch_registry%relinquish_host_array(temp_indices)
    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)
  end function cg_solve

  subroutine second_cg_part(rtr, r, mult, w, alpha, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: r(n), rtr
    real(kind=xp) :: tmp
    real(kind=rp), intent(in) ::mult(n), w(n), alpha
    integer :: i, ierr

    tmp = 0.0_xp
    !$omp parallel do reduction(+:tmp)
    do i = 1, n
       r(i) = r(i) - alpha*w(i)
       tmp = tmp + r(i) * r(i) * mult(i)
    end do
    !$omp end parallel do
    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    rtr = tmp

  end subroutine second_cg_part

  !> Solve three independent systems with the CPU PCG method.
  function cg_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, bc_projector, gs_h, niter) result(ksp_results)
    class(cg_t), intent(inout) :: this
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
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer, optional, intent(in) :: niter
    type(scalar_bc_projector_t), pointer :: bc_x, bc_y, bc_z

    call vector_bc_projector_components(bc_projector, bc_x, bc_y, bc_z)
    ksp_results(1) = this%solve(Ax, x, fx, n, coef, bc_x, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, bc_y, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, bc_z, gs_h, niter)

  end function cg_solve_coupled

end module cg
