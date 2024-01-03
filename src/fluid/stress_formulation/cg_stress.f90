! Copyright (c) 2024, The Neko Authors
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
module cg_stress
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, glsc2, add2s1
  use stress_formulation
  implicit none
  private

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_stress_t
     real(kind=rp), allocatable :: w1(:)
     real(kind=rp), allocatable :: w2(:)
     real(kind=rp), allocatable :: w3(:)
     real(kind=rp), allocatable :: r1(:)
     real(kind=rp), allocatable :: r2(:)
     real(kind=rp), allocatable :: r3(:)
     real(kind=rp), allocatable :: p1(:)
     real(kind=rp), allocatable :: p2(:)
     real(kind=rp), allocatable :: p3(:)
     real(kind=rp), allocatable :: z1(:)
     real(kind=rp), allocatable :: z2(:)
     real(kind=rp), allocatable :: z3(:)
     real(kind=rp), allocatable :: tmp(:)
   contains
     procedure, pass(this) :: init => cg_stress_init
     procedure, pass(this) :: free => cg_stress_free
     procedure, pass(this) :: solve => cg_stress_nop
     procedure, pass(this) :: solve_stress => cg_stress_solve
  end type cg_stress_t

contains

  !> Initialise a standard PCG solver
  subroutine cg_stress_init(this, n, M, rel_tol, abs_tol)
    class(cg_stress_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol

    call this%free()

    allocate(this%w1(n))
    allocate(this%w2(n))
    allocate(this%w3(n))
    allocate(this%r1(n))
    allocate(this%r2(n))
    allocate(this%r3(n))
    allocate(this%p1(n))
    allocate(this%p2(n))
    allocate(this%p3(n))
    allocate(this%z1(n))
    allocate(this%z2(n))
    allocate(this%z3(n))
    allocate(this%tmp(n))

    if (present(M)) then
       this%M => M
    end if

    if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(rel_tol, abs_tol)
    else if (present(rel_tol)) then
       call this%ksp_init(rel_tol=rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(abs_tol=abs_tol)
    else
       call this%ksp_init()
    end if

  end subroutine cg_stress_init

  !> Deallocate a standard PCG solver
  subroutine cg_stress_free(this)
    class(cg_stress_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%w1)) then
       deallocate(this%w1)
    end if

    if (allocated(this%w2)) then
       deallocate(this%w2)
    end if

    if (allocated(this%w3)) then
       deallocate(this%w3)
    end if

    if (allocated(this%r1)) then
       deallocate(this%r1)
    end if

    if (allocated(this%r2)) then
       deallocate(this%r2)
    end if

    if (allocated(this%r3)) then
       deallocate(this%r3)
    end if
    
    if (allocated(this%p1)) then
       deallocate(this%p1)
    end if

    if (allocated(this%p2)) then
       deallocate(this%p2)
    end if

    if (allocated(this%p3)) then
       deallocate(this%p3)
    end if

    if (allocated(this%z1)) then
       deallocate(this%z1)
    end if

    if (allocated(this%z2)) then
       deallocate(this%z2)
    end if

    if (allocated(this%z3)) then
       deallocate(this%z3)
    end if

    if (allocated(this%tmp)) then
       deallocate(this%tmp)
    end if

    nullify(this%M)

  end subroutine cg_stress_free
  
  function cg_stress_nop(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(cg_stress_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
  end function cg_stress_nop
  

  function cg_stress_solve(this, x, y, z, fx, fy, fz, n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(cg_stress_t), intent(inout) :: this
    type(field_t), intent(inout) :: x, y, z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: fx, fy, fz
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx, blsty, blstz
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
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = one / coef%volume

    rtz1 = one
    do i = 1, n
       x%x(i,1,1,1) = 0.0_rp
       y%x(i,1,1,1) = 0.0_rp
       z%x(i,1,1,1) = 0.0_rp
       this%p1(i) = 0.0_rp
       this%p2(i) = 0.0_rp
       this%p3(i) = 0.0_rp
       this%z1(i) = 0.0_rp
       this%z2(i) = 0.0_rp
       this%z3(i) = 0.0_rp
       this%r1(i) = fx(i)
       this%r2(i) = fy(i)
       this%r3(i) = fz(i)
       this%tmp(i) = this%r1(i)**2 + this%r2(i)**2 + this%r3(i)**2
    end do

    rtr = glsc3(this%tmp, coef%mult, coef%binv, n)
    rnorm = sqrt(rtr*norm_fac)
    ksp_results%res_start = rnorm
    ksp_results%res_final = rnorm
    ksp_results%iter = 0
    if(rnorm .eq. zero) return

    do iter = 1, max_iter
       call this%M%solve(this%z1, this%r1, n)
       call this%M%solve(this%z2, this%r2, n)
       call this%M%solve(this%z3, this%r3, n)
       rtz2 = rtz1

       do i = 1, n
          this%tmp(i) = this%z1(i) * this%r1(i) &
                      + this%z2(i) * this%r2(i) &
                      + this%z3(i) * this%r3(i)
       end do
       
       rtz1 = glsc2(this%tmp, coef%mult, n)

       beta = rtz1 / rtz2
       if (iter .eq. 1) beta = zero
       call add2s1(this%p1, this%z1, beta, n)
       call add2s1(this%p2, this%z2, beta, n)
       call add2s1(this%p3, this%z3, beta, n)

       call ax_helm_stress_compute(this%w1, this%w2, this%w3, &
                       this%p1, this%p2, this%p3, coef, x%msh, x%Xh)
       call gs_h%op(this%w1, n, GS_OP_ADD)
       call gs_h%op(this%w2, n, GS_OP_ADD)
       call gs_h%op(this%w3, n, GS_OP_ADD)
       call bc_list_apply(blstx, this%w1, n)
       call bc_list_apply(blsty, this%w2, n)
       call bc_list_apply(blstz, this%w3, n)

       do i = 1, n
          this%tmp(i) = this%w1(i) * this%p1(i) &
                      + this%w2(i) * this%p2(i) &
                      + this%w3(i) * this%p3(i)
       end do

       pap = glsc2(this%tmp, coef%mult, n)

       alpha = rtz1 / pap
       alphm = -alpha
       do i = 1, n
          x%x(i,1,1,1) = x%x(i,1,1,1) + alpha * this%p1(i)
          y%x(i,1,1,1) = y%x(i,1,1,1) + alpha * this%p2(i)
          z%x(i,1,1,1) = z%x(i,1,1,1) + alpha * this%p3(i)
          this%r1(i) = this%r1(i) + alphm * this%w1(i)
          this%r2(i) = this%r2(i) + alphm * this%w2(i)
          this%r3(i) = this%r3(i) + alphm * this%w3(i)
          this%tmp(i) = this%r1(i)**2 + this%r2(i)**2 + this%r3(i)**2
       end do

       rtr = glsc3(this%tmp, coef%mult, coef%binv, n)
       if (iter .eq. 1) rtr0 = rtr
       rnorm = sqrt(rtr * norm_fac)
       if (rnorm .lt. this%abs_tol) then
          exit
       end if
    end do
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function cg_stress_solve

end module cg_stress


