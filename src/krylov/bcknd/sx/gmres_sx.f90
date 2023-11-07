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
!> Defines various GMRES methods
module gmres_sx
  use krylov, only : ksp_t, ksp_monitor_t
  use precon,  only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc, only : bc_list_t, bc_list_apply
  use math, only : glsc3, rzero, rone, copy, cmult2, col2, col3, add2s2
  use comm
  implicit none
  private

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: sx_gmres_t
     integer :: lgmres
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: c(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: z(:,:)
     real(kind=rp), allocatable :: h(:,:)
     real(kind=rp), allocatable :: ml(:)
     real(kind=rp), allocatable :: v(:,:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: mu(:)
     real(kind=rp), allocatable :: gam(:)
     real(kind=rp), allocatable :: wk1(:)
     real(kind=rp) :: rnorm
   contains
     procedure, pass(this) :: init => sx_gmres_init
     procedure, pass(this) :: free => sx_gmres_free
     procedure, pass(this) :: solve => sx_gmres_solve
  end type sx_gmres_t

contains

  !> Initialise a standard GMRES solver
  subroutine sx_gmres_init(this, n, M, lgmres, rel_tol, abs_tol)
    class(sx_gmres_t), intent(inout) :: this
    integer, intent(in) :: n
    class(pc_t), optional, intent(inout), target :: M
    integer, optional, intent(inout) :: lgmres
    real(kind=rp), optional, intent(inout) :: rel_tol
    real(kind=rp), optional, intent(inout) :: abs_tol

    if (present(lgmres)) then
       this%lgmres = lgmres
    else
       this%lgmres = 30
    end if
    

    call this%free()
    
    if (present(M)) then 
       this%M => M
    end if

    allocate(this%w(n))
    allocate(this%r(n))
    allocate(this%ml(n))
    allocate(this%mu(n))
    allocate(this%wk1(n))
    
    allocate(this%c(this%lgmres))
    allocate(this%s(this%lgmres))
    allocate(this%gam(this%lgmres + 1))
    
    allocate(this%z(n,this%lgmres))
    allocate(this%v(n,this%lgmres))
    
    allocate(this%h(this%lgmres,this%lgmres))
    
       
    if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(rel_tol, abs_tol)
    else if (present(rel_tol)) then
       call this%ksp_init(rel_tol=rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(abs_tol=abs_tol)
    else
       call this%ksp_init(abs_tol)
    end if
          
  end subroutine sx_gmres_init

  !> Deallocate a standard GMRES solver
  subroutine sx_gmres_free(this)
    class(sx_gmres_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%c)) then
       deallocate(this%c)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if
 
    if (allocated(this%z)) then
       deallocate(this%z)
    end if

    if (allocated(this%h)) then
       deallocate(this%h)
    end if
    
    if (allocated(this%ml)) then
       deallocate(this%ml)
    end if
    
    if (allocated(this%v)) then
       deallocate(this%v)
    end if
    
    if (allocated(this%s)) then
       deallocate(this%s)
    end if
    
    if (allocated(this%mu)) then
       deallocate(this%mu)
    end if
    
    if (allocated(this%gam)) then
       deallocate(this%gam)
    end if
    
    if (allocated(this%wk1)) then
       deallocate(this%wk1)
    end if
    
    nullify(this%M)
    
  end subroutine sx_gmres_free
 
  !> Standard PCG solve
  function sx_gmres_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(sx_gmres_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, glb_n
    integer :: i, j, k, ierr 
    real(kind=rp), parameter :: one = 1.0
    real(kind=rp) :: rnorm 
    real(kind=rp) ::  alpha, temp, l
    real(kind=rp) :: ratio, div0, norm_fac
    logical :: conv
    integer outer

    conv = .false.
    iter = 0
    glb_n = n / x%msh%nelv * x%msh%glb_nelv

    call rone(this%ml, n)
    call rone(this%mu, n)
    norm_fac = one / sqrt(coef%volume)
    call rzero(x%x, n)
    call rzero(this%gam, this%lgmres + 1)
    call rone(this%s, this%lgmres)
    call rone(this%c, this%lgmres)
    call rzero(this%h, this%lgmres * this%lgmres)
    outer = 0
    do while (.not. conv .and. iter .lt. niter)
       outer = outer + 1

       if(iter.eq.0) then               
          call col3(this%r,this%ml,f,n) 
       else
          !update residual
          call copy  (this%r,f,n)      
          call Ax%compute(this%w, x%x, coef, x%msh, x%Xh)
          call gs_h%op(this%w, n, GS_OP_ADD)
          call bc_list_apply(blst, this%w, n)
          call add2s2(this%r,this%w,-one,n) 
          call col2(this%r,this%ml,n)       
       endif
       this%gam(1) = sqrt(glsc3(this%r, this%r, coef%mult, n))
       if(iter.eq.0) then
          div0 = this%gam(1) * norm_fac
          ksp_results%res_start = div0
       endif

       if ( this%gam(1) .eq. 0) return

       rnorm = 0.0_rp
       temp = one / this%gam(1)
       call cmult2(this%v(1,1), this%r, temp, n) 
       do j = 1, this%lgmres
          iter = iter+1
          call col3(this%w, this%mu, this%v(1,j), n)

          !Apply precond
          call this%M%solve(this%z(1,j), this%w, n)

          call Ax%compute(this%w, this%z(1,j), coef, x%msh, x%Xh)
          call gs_h%op(this%w, n, GS_OP_ADD)
          call bc_list_apply(blst, this%w, n)
          call col2(this%w, this%ml, n)       

          do i = 1, j
             this%h(i,j) = 0.0_rp
             do k = 1, n
                this%h(i,j) = this%h(i,j) + &
                     this%w(k) * this%v(k,i) * coef%mult(k,1,1,1)
             end do
          end do

          !Could probably be done inplace...
          call MPI_Allreduce(this%h(1,j), this%wk1, j, &
               MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
          call copy(this%h(1,j), this%wk1, j) 

          do i = 1, j
             do k = 1, n
                this%w(k) = this%w(k) - this%h(i,j) * this%v(k,i)
             end do
          end do

          !apply Givens rotations to new column
          do i=1,j-1
             temp = this%h(i,j)                   
             this%h(i  ,j) =  this%c(i)*temp + this%s(i)*this%h(i+1,j)  
             this%h(i+1,j) = -this%s(i)*temp + this%c(i)*this%h(i+1,j)
          end do

          alpha = sqrt(glsc3(this%w, this%w, coef%mult, n))   
          rnorm = 0.0_rp
          if(alpha .eq. 0.0_rp) then 
            conv = .true.
            exit
          end if
          l = sqrt(this%h(j,j) * this%h(j,j) + alpha**2)
          temp = one / l
          this%c(j) = this%h(j,j) * temp
          this%s(j) = alpha  * temp
          this%h(j,j) = l
          this%gam(j+1) = -this%s(j) * this%gam(j)
          this%gam(j)   =  this%c(j) * this%gam(j)

          rnorm = abs(this%gam(j+1)) * norm_fac
          ratio = rnorm / div0
          if (rnorm .lt. this%abs_tol) then 
             conv = .true.
             exit
          end if
         
          if (iter + 1 .gt. niter) exit
          
          if( j .lt. this%lgmres) then
            temp = one / alpha
            call cmult2(this%v(1,j+1), this%w, temp, n)
          endif
       end do
       j = min(j, this%lgmres)
       !back substitution
       do k = j, 1, -1
          temp = this%gam(k)
          do i = j, k+1, -1
             temp = temp - this%h(k,i) * this%c(i)
          enddo
          this%c(k) = temp / this%h(k,k)
       enddo
       !sum up Arnoldi vectors
       do i = 1, j
          do k = 1, n
             x%x(k,1,1,1) = x%x(k,1,1,1) + this%c(i) * this%z(k,i)
          end do
       end do 
    end do

    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function sx_gmres_solve

end module gmres_sx
  

