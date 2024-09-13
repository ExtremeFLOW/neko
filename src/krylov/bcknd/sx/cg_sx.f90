! Copyright (c) 2021, The Neko Authors
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
module cg_sx
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, add2s1, abscmp
  implicit none
  private

  !> Standard preconditioned conjugate gradient method (SX version)
  type, public, extends(ksp_t) :: sx_cg_t
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: z(:)
   contains
     procedure, pass(this) :: init => sx_cg_init
     procedure, pass(this) :: free => sx_cg_free
     procedure, pass(this) :: solve => sx_cg_solve
     procedure, pass(this) :: solve_coupled => sx_cg_solve_coupled
  end type sx_cg_t

contains

  !> Initialise a standard PCG solver
  subroutine sx_cg_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(sx_cg_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()

    allocate(this%w(n))
    allocate(this%r(n))
    allocate(this%p(n))
    allocate(this%z(n))

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

  end subroutine sx_cg_init

  !> Deallocate a standard PCG solver
  subroutine sx_cg_free(this)
    class(sx_cg_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if

    if (allocated(this%z)) then
       deallocate(this%z)
    end if

    nullify(this%M)

  end subroutine sx_cg_free

  !> Standard PCG solve
  function sx_cg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(sx_cg_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    real(kind=rp), parameter :: one = 1.0
    real(kind=rp), parameter :: zero = 0.0
    integer :: i, iter, max_iter
    real(kind=rp) :: rnorm, rtr, rtr0, rtz2, rtz1
    real(kind=rp) :: beta, pap, alpha, alphm, norm_fac

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = one / sqrt(coef%volume)

    rtz1 = one
    do i = 1, n
       x%x(i,1,1,1) = 0.0_rp
       this%p(i) = 0.0_rp
       this%r(i) = f(i)
    end do

    rtr = glsc3(this%r, coef%mult, this%r, n)
    rnorm = sqrt(rtr)*norm_fac
    ksp_results%res_start = rnorm
    ksp_results%res_final = rnorm
    ksp_results%iter = 0
    if(abscmp(rnorm, zero)) return
    call this%monitor_start('CG')
    do iter = 1, max_iter
       call this%M%solve(this%z, this%r, n)
       rtz2 = rtz1
       rtz1 = glsc3(this%r, coef%mult, this%z, n)

       beta = rtz1 / rtz2
       if (iter .eq. 1) beta = zero
       call add2s1(this%p, this%z, beta, n)

       call Ax%compute(this%w, this%p, coef, x%msh, x%Xh)
       call gs_h%op(this%w, n, GS_OP_ADD)
       call bc_list_apply(blst, this%w, n)

       pap = glsc3(this%w, coef%mult, this%p, n)

       alpha = rtz1 / pap
       alphm = -alpha
       do i = 1, n
          x%x(i,1,1,1) = x%x(i,1,1,1) + alpha * this%p(i)
          this%r(i) = this%r(i) + alphm * this%w(i)
       end do

       rtr = glsc3(this%r, coef%mult, this%r, n)
       if (iter .eq. 1) rtr0 = rtr
       rnorm = sqrt(rtr) * norm_fac
       call this%monitor_iter(iter, rnorm)
       if (rnorm .lt. this%abs_tol) then
          exit
       end if
    end do
    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function sx_cg_solve

  !> Standard PCG coupled solve
  function sx_cg_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(sx_cg_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: fx
    real(kind=rp), dimension(n), intent(inout) :: fy
    real(kind=rp), dimension(n), intent(inout) :: fz
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx
    type(bc_list_t), intent(inout) :: blsty
    type(bc_list_t), intent(inout) :: blstz
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer, optional, intent(in) :: niter

    ksp_results(1) =  this%solve(Ax, x, fx, n, coef, blstx, gs_h, niter)
    ksp_results(2) =  this%solve(Ax, y, fy, n, coef, blsty, gs_h, niter)
    ksp_results(3) =  this%solve(Ax, z, fz, n, coef, blstz, gs_h, niter)

  end function sx_cg_solve_coupled

end module cg_sx


