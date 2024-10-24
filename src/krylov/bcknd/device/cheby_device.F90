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
!> Chebyshev preconditioner
module cheby_device
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
  use device_math, only : device_cmult2, device_sub2, &
       device_add2s1, device_add2s2, device_glsc3, device_copy
  use device
  implicit none
  private

  !> Defines a Chebyshev preconditioner
  type, public, extends(ksp_t) :: cheby_device_t
     real(kind=rp), allocatable :: d(:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     type(c_ptr) :: d_d = C_NULL_PTR
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: gs_event = C_NULL_PTR
     real(kind=rp) :: tha, dlt
     integer :: power_its = 150
     logical :: recompute_eigs = .true.
   contains
     procedure, pass(this) :: init => cheby_device_init
     procedure, pass(this) :: free => cheby_device_free
     procedure, pass(this) :: solve => cheby_device_solve
     procedure, pass(this) :: solve_coupled => cheby_device_solve_coupled
  end type cheby_device_t

contains

  !> Initialise a standard solver
  subroutine cheby_device_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(cheby_device_t), intent(inout), target :: this
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()
    allocate(this%d(n))
    allocate(this%w(n))
    allocate(this%r(n))

    call device_map(this%d, this%d_d, n)
    call device_map(this%w, this%w_d, n)
    call device_map(this%r, this%r_d, n)

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

    call device_event_create(this%gs_event, 2)
    
  end subroutine cheby_device_init

  subroutine cheby_device_free(this)
    class(cheby_device_t), intent(inout) :: this

    call this%ksp_free()
    
    if (allocated(this%d)) then
       deallocate(this%d)
    end if

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    nullify(this%M)
    
    if (c_associated(this%d_d)) then
       call device_free(this%d_d)
    end if

    if (c_associated(this%w_d)) then
       call device_free(this%w_d)
    end if

    if (c_associated(this%r_d)) then
       call device_free(this%r_d)
    end if

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if
    
  end subroutine cheby_device_free

  subroutine cheby_device_power(this, Ax, x, n, coef, blst, gs_h)
    class(cheby_device_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp) :: lam, b, a, rn
    real(kind=rp) :: boost = 1.2_rp
    real(kind=rp) :: lam_factor = 30.0_rp
    real(kind=rp) :: wtw, dtw, dtd
    integer :: i

    associate(w => this%w, w_d => this%w_d, d => this%d, d_d => this%d_d)
      
      do i = 1, n
         !TODO: replace with a better way to initialize power method
         call random_number(rn)
         d(i) = rn + 10.0_rp
      end do
      call device_memcpy(d, d_d, n, HOST_TO_DEVICE, sync = .true.)
      
      call gs_h%op(d, n, GS_OP_ADD, this%gs_event)
      call bc_list_apply(blst, d, n)

      !Power method to get lamba max
      do i = 1, this%power_its
        call ax%compute(w, d, coef, x%msh, x%Xh)
        call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
        call bc_list_apply(blst, w, n)

        wtw = device_glsc3(w_d, coef%mult_d, w_d, n)
        call device_cmult2(d_d, w_d, 1.0_rp/sqrt(wtw), n)
        call bc_list_apply(blst, d, n)
      end do

      call ax%compute(w, d, coef, x%msh, x%Xh)
      call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
      call bc_list_apply(blst, w, n)

      dtw = device_glsc3(d_d, coef%mult_d, w_d, n)
      dtd = device_glsc3(d_d, coef%mult_d, d_d, n)
      lam = dtw / dtd
      b = lam * boost
      a = lam / lam_factor
      this%tha = (b+a)/2.0_rp
      this%dlt = (b-a)/2.0_rp

      this%recompute_eigs = .false.
    end associate
  end subroutine cheby_device_power

  !> A chebyshev preconditioner
  function cheby_device_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) &
       result(ksp_results)
    class(cheby_device_t), intent(inout) :: this
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
    real(kind=rp) :: a, b, rtr, rnorm, norm_fac
    type(c_ptr) :: f_d

    f_d = device_get_ptr(f)

    if (this%recompute_eigs) then
       call cheby_device_power(this, Ax, x, n, coef, blst, gs_h)
    end if
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate( w => this%w, r => this%r, d => this%d, &
         w_d => this%w_d, r_d => this%r_d, d_d => this%d_d)
      ! calculate residual
      call device_copy(r_d, f_d, n)
      call ax%compute(w, x%x, coef, x%msh, x%Xh)
      call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
      call bc_list_apply(blst, w, n)
      call device_sub2(r_d, w_d, n)

      rtr = device_glsc3(r_d, coef%mult_d, r_d, n)
      rnorm = sqrt(rtr) * norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0

      ! First iteration
      call this%M%solve(w, r, n)
      call device_copy(d_d, w_d, n)
      a = 2.0_rp / this%tha
      call device_add2s2(x%x_d, d_d, a, n)! x = x + a*d

      ! Rest of the iterations
      do iter = 2, max_iter
        ! calculate residual
        call device_copy(r_d, f_d, n)
        call ax%compute(w, x%x, coef, x%msh, x%Xh)
        call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
        call bc_list_apply(blst, w, n)
        call device_sub2(r_d, w_d, n)

        call this%M%solve(w, r, n)

        if (iter .eq. 2) then
          b = 0.5_rp * (this%dlt * a)**2
        else
          b = (this%dlt * a / 2.0_rp)**2
        end if
        a = 1.0_rp/(this%tha - b/a)
        call device_add2s1(d_d, w_d, b, n)! d = w + b*d

        call device_add2s2(x%x_d, d_d, a, n)! x = x + a*d
      end do

      ! calculate residual
      call device_copy(r_d, f_d, n)
      call ax%compute(w, x%x, coef, x%msh, x%Xh)
      call gs_h%op(w, n, GS_OP_ADD, this%gs_event)
      call bc_list_apply(blst, w, n)
      call device_sub2(r_d, w_d, n)
      rtr = device_glsc3(r_d, coef%mult_d, r_d, n)
      rnorm = sqrt(rtr) * norm_fac
      ksp_results%res_final = rnorm
      ksp_results%iter = iter
    end associate
  end function cheby_device_solve

  !> Standard Cheby_Deviceshev coupled solve
  function cheby_device_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(cheby_device_t), intent(inout) :: this
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

  end function cheby_device_solve_coupled

end module cheby_device


