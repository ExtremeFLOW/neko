! Copyright (c) 2020-2021, The Neko Authors
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
!> Chebyshev preconditioner
module cheby
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp
  use field, only : field_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use space, only : space_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, rzero, rone, copy, sub2, cmult2, abscmp, glsc2, add2s1, add2s2
  use comm
  implicit none
  private

  !> Defines a Chebyshev preconditioner
  type, public, extends(ksp_t) :: cheby_t
     class(ax_t), pointer :: ax
     real(kind=rp), allocatable :: d(:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp) :: tha, dlt
     integer :: power_its = 50
     type(gs_t), pointer :: gs_h
     type(coef_t), pointer :: coef
     type(bc_list_t), pointer :: blst
     type(mesh_t), pointer :: msh
     type(space_t), pointer :: Xh
   contains
     procedure, pass(this) :: init => cheby_init
     procedure, pass(this) :: update => cheby_update
     procedure, pass(this) :: solve => cheby_solve
     procedure, pass(this) :: solve_coupled => cheby_solve_coupled
     procedure, pass(this) :: free => cheby_free
  end type cheby_t

contains

  subroutine cheby_init(this, n, Ax, coef, gs_h, blst, msh, Xh)
    class(cheby_t), intent(inout) :: this
    integer, intent(in) :: n
    class(ax_t), intent(inout), target :: Ax
    type(coef_t), intent(inout), target :: coef
    type(gs_t), intent(inout), target :: gs_h
    type(bc_list_t), intent(inout), target :: blst
    type(mesh_t), intent(inout), target :: msh
    type(space_t), intent(inout), target :: Xh

    call this%free()
    this%Ax => Ax
    this%gs_h => gs_h
    this%coef => coef
    this%blst => blst
    this%msh => msh
    this%Xh => Xh
    allocate(this%d(n))
    allocate(this%w(n))
    allocate(this%r(n))
    call cheby_update(this)

  end subroutine cheby_init

  !> Update Chebyshev preconditioner if the geometry G has changed
  subroutine cheby_update(this)
    class(cheby_t), intent(inout) :: this
    real(kind=rp) :: lam, b, a
    real(kind=rp) :: boost = 1.2
    real(kind=rp) :: lam_factor = 30.0
    integer :: i, n
    associate(coef => this%coef, msh=>this%msh, Xh=>this%Xh, gs_h => this%gs_h, &
        blst => this%blst, ax => this%ax, w => this%w, d => this%d)

      n = Xh%lx * Xh%ly * Xh%lz * msh%nelv

      !Power method to get lamba max
      do i = 1, this%power_its
        call ax%compute(w, d, coef, msh, Xh)
        call gs_h%op(w, n, GS_OP_ADD)
        call bc_list_apply(blst, w, n)

        call cmult2(d, w, 1.0/sqrt(glsc2(w, w, n)), n)
      end do

      call ax%compute(w, d, coef, msh, Xh)
      call gs_h%op(w, n, GS_OP_ADD)
      call bc_list_apply(blst, w, n)

      lam = glsc2(d, w, n) / glsc2(d, d, n)
      b = lam * boost
      a = lam / lam_factor
      this%tha = (b+a)/2.0
      this%dlt = (b-a)/2.0
    end associate
  end subroutine cheby_update

  subroutine cheby_free(this)
    class(cheby_t), intent(inout) :: this
    if (allocated(this%d)) then
       deallocate(this%d)
    end if
    nullify(this%gs_h)
    nullify(this%coef)
  end subroutine cheby_free

  !> A chebyshev preconditioner
  function cheby_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(cheby_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=rp) :: a, b

    associate( w => this%w, r => this%r, d => this%d)
      ! calculate residual
      call copy(r, f, n)
      call ax%compute(w, x%x, coef, x%msh, x%Xh)
      call gs_h%op(w, n, GS_OP_ADD)
      call bc_list_apply(blst, w, n)
      call sub2(r, w, n)

      ! First iteration
      call copy(d, r, n)
      a = 2.0 / this%tha
      call add2s2(x%x, d, a, n)! x = x + a*d

      ! Rest of the iterations
      do iter = 2, this%max_iter
        ! calculate residual
        call copy(r, f, n)
        call ax%compute(w, x%x, coef, x%msh, x%Xh)
        call gs_h%op(w, n, GS_OP_ADD)
        call bc_list_apply(blst, w, n)
        call sub2(r, w, n)

        b = (this%dlt * a / 2.0)**2
        a = 1.0 / (this%tha - b)
        call add2s1(d, r, b, n)! d = r + b*d

        call add2s2(x%x, d, a, n)! x = x + a*d
      end do
    end associate
  end function cheby_solve

  !> Standard Chebyshev coupled solve
  function cheby_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(cheby_t), intent(inout) :: this
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

    ksp_results(1) = this%solve(Ax, x, fx, n, coef, blstx, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, blsty, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, blstz, gs_h, niter)

  end function cheby_solve_coupled
end module cheby
