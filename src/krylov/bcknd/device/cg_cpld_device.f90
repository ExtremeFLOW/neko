! Copyright (c) 2025, The Neko Authors
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
!> Defines a coupled  Conjugate Gradient methods for accelerators
module cg_cpld_device
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use math, only : abscmp
  use device
  use device_math, only : device_rzero, device_copy, device_glsc3, &
       device_add2s1, device_vdot3, device_glsc2
  use device_mathops, only : device_opadd2cm
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  implicit none

  !> Device based coupled preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_cpld_device_t
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


     type(c_ptr) :: w1_d = C_NULL_PTR
     type(c_ptr) :: w2_d = C_NULL_PTR
     type(c_ptr) :: w3_d = C_NULL_PTR

     type(c_ptr) :: r1_d = C_NULL_PTR
     type(c_ptr) :: r2_d = C_NULL_PTR
     type(c_ptr) :: r3_d = C_NULL_PTR

     type(c_ptr) :: p1_d = C_NULL_PTR
     type(c_ptr) :: p2_d = C_NULL_PTR
     type(c_ptr) :: p3_d = C_NULL_PTR

     type(c_ptr) :: z1_d = C_NULL_PTR
     type(c_ptr) :: z2_d = C_NULL_PTR
     type(c_ptr) :: z3_d = C_NULL_PTR

     type(c_ptr) :: tmp_d = C_NULL_PTR

     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     procedure, pass(this) :: init => cg_cpld_device_init
     procedure, pass(this) :: free => cg_cpld_device_free
     procedure, pass(this) :: solve => cg_cpld_device_nop
     procedure, pass(this) :: solve_coupled => cg_cpld_device_solve
  end type cg_cpld_device_t

contains

  !> Initialise a device based PCG solver
  subroutine cg_cpld_device_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(cg_cpld_device_t), target, intent(inout) :: this
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor

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

    call device_map(this%tmp, this%tmp_d, n)
    call device_map(this%z1, this%z1_d, n)
    call device_map(this%z2, this%z2_d, n)
    call device_map(this%z3, this%z3_d, n)
    call device_map(this%p1, this%p1_d, n)
    call device_map(this%p2, this%p2_d, n)
    call device_map(this%p3, this%p3_d, n)
    call device_map(this%r1, this%r1_d, n)
    call device_map(this%r2, this%r2_d, n)
    call device_map(this%r3, this%r3_d, n)
    call device_map(this%w1, this%w1_d, n)
    call device_map(this%w2, this%w2_d, n)
    call device_map(this%w3, this%w3_d, n)

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
  end subroutine cg_cpld_device_init

  !> Deallocate a device based PCG solver
  subroutine cg_cpld_device_free(this)
    class(cg_cpld_device_t), intent(inout) :: this

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

    if (c_associated(this%w1_d)) then
       call device_free(this%w1_d)
    end if

    if (c_associated(this%w2_d)) then
       call device_free(this%w2_d)
    end if

    if (c_associated(this%w3_d)) then
       call device_free(this%w3_d)
    end if

    if (c_associated(this%r1_d)) then
       call device_free(this%r1_d)
    end if

    if (c_associated(this%r2_d)) then
       call device_free(this%r2_d)
    end if

    if (c_associated(this%r3_d)) then
       call device_free(this%r3_d)
    end if

    if (c_associated(this%p1_d)) then
       call device_free(this%p1_d)
    end if

    if (c_associated(this%p2_d)) then
       call device_free(this%p2_d)
    end if

    if (c_associated(this%p3_d)) then
       call device_free(this%p3_d)
    end if

    if (c_associated(this%z1_d)) then
       call device_free(this%z1_d)
    end if

    if (c_associated(this%z2_d)) then
       call device_free(this%z2_d)
    end if

    if (c_associated(this%z3_d)) then
       call device_free(this%z3_d)
    end if

    if (c_associated(this%tmp_d)) then
       call device_free(this%tmp_d)
    end if

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if

  end subroutine cg_cpld_device_free

  function cg_cpld_device_nop(this, Ax, x, f, n, coef, blst, gs_h, niter) &
       result(ksp_results)
    class(cg_cpld_device_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter

    ! Throw and error
    call neko_error('The cpldcg solver is only defined for coupled solves')

    ksp_results%res_final = 0.0
    ksp_results%iter = 0
  end function cg_cpld_device_nop

  !> Standard PCG solve
  function cg_cpld_device_solve(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(cg_cpld_device_t), intent(inout) :: this
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
    integer :: i, iter, max_iter
    real(kind=rp) :: rnorm, rtr, rtr0, rtz2, rtz1
    real(kind=rp) :: beta, pap, alpha, alphm, norm_fac
    integer, parameter :: gdim = 3
    type(c_ptr) :: fx_d
    type(c_ptr) :: fy_d
    type(c_ptr) :: fz_d

    fx_d = device_get_ptr(fx)
    fy_d = device_get_ptr(fy)
    fz_d = device_get_ptr(fz)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / coef%volume

    associate (p1_d => this%p1_d, p2_d => this%p2_d, p3_d => this%p3_d, &
         z1_d => this%z1_d, z2_d => this%z2_d, z3_d => this%z3_d, &
         r1_d => this%r1_d, r2_d => this%r2_d, r3_d => this%r3_d, &
         w1_d => this%w1_d, w2_d => this%w2_d, w3_d => this%w3_d, &
         tmp_d => this%tmp_d)

      rtz1 = 1.0_rp
      call device_rzero(x%x_d, n)
      call device_rzero(y%x_d, n)
      call device_rzero(z%x_d, n)
      call device_rzero(p1_d, n)
      call device_rzero(p2_d, n)
      call device_rzero(p3_d, n)
      call device_rzero(z1_d, n)
      call device_rzero(z2_d, n)
      call device_rzero(z3_d, n)
      call device_copy(r1_d, fx_d, n)
      call device_copy(r2_d, fy_d, n)
      call device_copy(r3_d, fz_d, n)
      call device_vdot3(tmp_d, r1_d, r2_d, r3_d, r1_d, r2_d, r3_d, n)


      rtr = device_glsc3(tmp_d, coef%mult_d, coef%binv_d, n)
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(abscmp(rnorm, 0.0_rp)) then
         ksp_results%converged = .true.
         return
      end if

      call this%monitor_start('device_cpldCG')
      do iter = 1, max_iter
         call this%M%solve(this%z1, this%r1, n)
         call this%M%solve(this%z2, this%r2, n)
         call this%M%solve(this%z3, this%r3, n)
         rtz2 = rtz1

         call device_vdot3(tmp_d, z1_d, z2_d, z3_d, r1_d, r2_d, r3_d, n)

         rtz1 = device_glsc2(tmp_d, coef%mult_d, n)

         beta = rtz1 / rtz2
         if (iter .eq. 1) beta = 0.0_rp
         call device_add2s1(p1_d, z1_d, beta, n)
         call device_add2s1(p2_d, z2_d, beta, n)
         call device_add2s1(p3_d, z3_d, beta, n)

         call Ax%compute_vector(this%w1, this%w2, this%w3, &
              this%p1, this%p2, this%p3, coef, x%msh, x%Xh)
         call gs_h%op(this%w1, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
         call gs_h%op(this%w2, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)
         call gs_h%op(this%w3, n, GS_OP_ADD, this%gs_event)
         call device_event_sync(this%gs_event)

         call blstx%apply(this%w1, n)
         call blsty%apply(this%w2, n)
         call blstz%apply(this%w3, n)

         call device_vdot3(tmp_d, w1_d, w2_d, w3_d, p1_d, p2_d, p3_d, n)

         pap = device_glsc2(tmp_d, coef%mult_d, n)

         alpha = rtz1 / pap
         alphm = -alpha
         call device_opadd2cm(x%x_d, y%x_d, z%x_d, &
              p1_d, p2_d, p3_d, alpha, n, gdim)
         call device_opadd2cm(r1_d, r2_d, r3_d, &
              w1_d, w2_d, w3_d, alphm, n, gdim)
         call device_vdot3(tmp_d, r1_d, r2_d, r3_d, r1_d, r2_d, r3_d, n)

         rtr = device_glsc3(tmp_d, coef%mult_d, coef%binv_d, n)
         if (iter .eq. 1) rtr0 = rtr
         rnorm = sqrt(rtr) * norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol) then
            exit
         end if
      end do
    end associate
    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)

  end function cg_cpld_device_solve

end module cg_cpld_device
