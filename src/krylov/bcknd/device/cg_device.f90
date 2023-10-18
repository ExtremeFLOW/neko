! Copyright (c) 2021-2022, The Neko Authors
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
!> Defines various Conjugate Gradient methods for accelerators
module cg_device
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use device  
  use device_math, only : device_rzero, device_copy, device_glsc3, &
                          device_add2s2, device_add2s1
  implicit none

  !> Device based preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_device_t
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: z(:)
     type(c_ptr) :: w_d = C_NULL_PTR
     type(c_ptr) :: r_d = C_NULL_PTR
     type(c_ptr) :: p_d = C_NULL_PTR
     type(c_ptr) :: z_d = C_NULL_PTR
     type(c_ptr) :: gs_event = C_NULL_PTR
   contains
     procedure, pass(this) :: init => cg_device_init
     procedure, pass(this) :: free => cg_device_free
     procedure, pass(this) :: solve => cg_device_solve
  end type cg_device_t

contains

  !> Initialise a device based PCG solver
  subroutine cg_device_init(this, n, M, rel_tol, abs_tol)
    class(cg_device_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol
    
    call this%free()
    
    allocate(this%w(n))
    allocate(this%r(n))
    allocate(this%p(n))
    allocate(this%z(n))

    call device_map(this%z, this%z_d, n)
    call device_map(this%p, this%p_d, n)
    call device_map(this%r, this%r_d, n)
    call device_map(this%w, this%w_d, n)
    
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

    call device_event_create(this%gs_event, 2)
  end subroutine cg_device_init

  !> Deallocate a device based PCG solver
  subroutine cg_device_free(this)
    class(cg_device_t), intent(inout) :: this

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

    if (c_associated(this%w_d)) then
       call device_free(this%w_d)
    end if

    if (c_associated(this%r_d)) then
       call device_free(this%r_d)
    end if

    if (c_associated(this%p_d)) then
       call device_free(this%p_d)
    end if

    if (c_associated(this%z_d)) then
       call device_free(this%z_d)
    end if

    if (c_associated(this%gs_event)) then
       call device_event_destroy(this%gs_event)
    end if
    
  end subroutine cg_device_free
  
  !> Standard PCG solve
  function cg_device_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(cg_device_t), intent(inout) :: this
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
    integer :: iter, max_iter
    real(kind=rp) :: rnorm, rtr, rtr0, rtz2, rtz1
    real(kind=rp) :: beta, pap, alpha, alphm, norm_fac
    type(c_ptr) :: f_d
    
    f_d = device_get_ptr(f)

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = one/sqrt(coef%volume)

    rtz1 = one
    call device_rzero(x%x_d, n)
    call device_rzero(this%p_d, n)
    call device_copy(this%r_d, f_d, n)

    rtr = device_glsc3(this%r_d, coef%mult_d, this%r_d, n)
    rnorm = sqrt(rtr)*norm_fac
    ksp_results%res_start = rnorm
    ksp_results%res_final = rnorm
    ksp_results%iter = 0
    if(rnorm .eq. zero) return
    do iter = 1, max_iter
       call this%M%solve(this%z, this%r, n)
       rtz2 = rtz1
       rtz1 = device_glsc3(this%r_d, coef%mult_d, this%z_d, n)
       beta = rtz1 / rtz2
       if (iter .eq. 1) beta = zero
       call device_add2s1(this%p_d, this%z_d, beta, n)

       call Ax%compute(this%w, this%p, coef, x%msh, x%Xh)       
       call gs_h%op(this%w, n, GS_OP_ADD, this%gs_event)
       call device_event_sync(this%gs_event)
       call bc_list_apply(blst, this%w, n)       

       pap = device_glsc3(this%w_d, coef%mult_d, this%p_d, n)

       alpha = rtz1 / pap
       alphm = -alpha
       call device_add2s2(x%x_d, this%p_d, alpha, n)
       call device_add2s2(this%r_d, this%w_d, alphm, n)

       rtr = device_glsc3(this%r_d, coef%mult_d, this%r_d, n)
       if (iter .eq. 1) rtr0 = rtr
       rnorm = sqrt(rtr)*norm_fac
       if (rnorm .lt. this%abs_tol) then
          exit
       end if
    end do
    ksp_results%res_final = rnorm
    ksp_results%iter = iter

  end function cg_device_solve

end module cg_device
  

