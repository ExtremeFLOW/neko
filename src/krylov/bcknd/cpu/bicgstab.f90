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
!> Defines various Bi-Conjugate Gradient Stabilized methods
module bicgstab
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, rzero, copy, NEKO_EPS, add2s2, x_update, &
                   p_update
  use utils, only : neko_error
  implicit none
  private

  !> Standard preconditioned Bi-Conjugate Gradient Stabilized method
  type, public, extends(ksp_t) :: bicgstab_t
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: p_hat(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: s_hat(:)
     real(kind=rp), allocatable :: t(:)
     real(kind=rp), allocatable :: v(:)
   contains
     procedure, pass(this) :: init => bicgstab_init
     procedure, pass(this) :: free => bicgstab_free
     procedure, pass(this) :: solve => bicgstab_solve
  end type bicgstab_t

contains

  !> Initialise a standard BiCGSTAB solver
  subroutine bicgstab_init(this, n, M, rel_tol, abs_tol)
    class(bicgstab_t), intent(inout) :: this
    class(pc_t), optional, intent(inout), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol

        
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

    if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(rel_tol, abs_tol)
    else if (present(rel_tol)) then
       call this%ksp_init(rel_tol=rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(abs_tol=abs_tol)
    else
       call this%ksp_init()
    end if
          
  end subroutine bicgstab_init

  !> Deallocate a standard BiCGSTAB solver
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
  
  !> Bi-Conjugate Gradient Stabilized method solve
  function bicgstab_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(bicgstab_t), intent(inout) :: this
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
    real(kind=rp) :: rnorm, rtr, norm_fac, gamma
    real(kind=rp) :: beta, alpha, omega, rho_1, rho_2
    
    if (present(niter)) then
       max_iter = niter
    else
       max_iter = KSP_MAX_ITER
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(r => this%r, t => this%t, s => this%s, v => this%v, p => this%p, &
         s_hat => this%s_hat, p_hat => this%p_hat)
    
      call rzero(x%x, n)
      call copy(r, f, n)

      rtr = sqrt(glsc3(r, coef%mult, r, n))
      rnorm = rtr * norm_fac
      gamma = rnorm * this%rel_tol
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(rnorm .eq. 0.0_rp) return
      do iter = 1, max_iter
         
         rho_1 = glsc3(r, coef%mult, f ,n)
         
         if (abs(rho_1) .lt. NEKO_EPS) then
            call neko_error('Bi-CGStab rho failure')
         end if
   
         if (iter .eq. 1) then
            call copy(p, r, n) 
         else
            beta = (rho_1 / rho_2) * (alpha / omega)
            call p_update(p, r, v, beta, omega, n)
         end if
       
         call this%M%solve(p_hat, p, n)
         call Ax%compute(v, p_hat, coef, x%msh, x%Xh)
         call gs_h%op(v, n, GS_OP_ADD)
         call bc_list_apply(blst, v, n)
         alpha = rho_1 / glsc3(f, coef%mult, v, n)
         call copy(s, r, n)
         call add2s2(s, v, -alpha, n)
         rtr = glsc3(s, coef%mult, s, n)
         rnorm = sqrt(rtr) * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            call add2s2(x%x, p_hat, alpha,n)
            exit
         end if
       
         call this%M%solve(s_hat, s, n)
         call Ax%compute(t, s_hat, coef, x%msh, x%Xh)
         call gs_h%op(t, n, GS_OP_ADD)
         call bc_list_apply(blst, t, n)
         omega = glsc3(t, coef%mult, s, n) &
              / glsc3(t, coef%mult, t, n)
         call x_update(x%x, p_hat, s_hat, alpha, omega, n)
         call copy(r, s, n)
         call add2s2(r, t, -omega, n)
      
         rtr = glsc3(r, coef%mult, r, n)
         rnorm = sqrt(rtr) * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            exit
         end if
    
         if (omega .lt. NEKO_EPS) then
            call neko_error('Bi-CGstab omega failure')
         end if
         rho_2 = rho_1
    
      end do
    end associate
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function bicgstab_solve

end module bicgstab
  

