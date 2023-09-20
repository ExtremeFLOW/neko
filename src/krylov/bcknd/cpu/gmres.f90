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
!> Defines various GMRES methods
module gmres
  use krylov
  use math
  use comm
  implicit none
  private

  !> Base type for standard preconditioned GMRES
  type, private, abstract, extends(ksp_t) :: gmres_base_t
     integer :: lgmres
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: c(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: z(:,:)
     real(kind=rp), allocatable :: h(:,:)
     real(kind=rp), allocatable :: v(:,:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: gam(:)
     real(kind=rp), allocatable :: wk1(:)
  end type gmres_base_t
  
  !> Standard preconditioned generalized minimal residual method 
  type, public, extends(gmres_base_t) :: gmres_t
   contains
     procedure, pass(this) :: init => gmres_init
     procedure, pass(this) :: free => gmres_free
     procedure, pass(this) :: solve => gmres_solve
  end type gmres_t

  !> Standard preconditioned generalized minimal residual method
  !! OpenMP version
  type, public, extends(gmres_base_t) :: gmres_omp_t
   contains
     procedure, pass(this) :: init => gmres_omp_init
     procedure, pass(this) :: free => gmres_omp_free
     procedure, pass(this) :: solve => gmres_omp_solve
  end type gmres_omp_t

contains

  !> Initialise a standard GMRES solver
  subroutine gmres_init(this, n, M, lgmres, rel_tol, abs_tol)
    class(gmres_t), intent(inout) :: this
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
          
  end subroutine gmres_init

  !> Deallocate a standard GMRES solver
  subroutine gmres_free(this)
    class(gmres_t), intent(inout) :: this

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
    
    if (allocated(this%v)) then
       deallocate(this%v)
    end if
    
    if (allocated(this%s)) then
       deallocate(this%s)
    end if
    
    
    if (allocated(this%gam)) then
       deallocate(this%gam)
    end if
    
    if (allocated(this%wk1)) then
       deallocate(this%wk1)
    end if
    
    nullify(this%M)
    
  end subroutine gmres_free
 
  !> Standard GMRES solve
  function gmres_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(gmres_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter 
    integer :: i, j, k, l, ierr 
    real(kind=rp) :: w_plus(NEKO_BLK_SIZE), x_plus(NEKO_BLK_SIZE)
    real(kind=rp) :: rnorm, alpha, temp, lr, alpha2, norm_fac
    logical :: conv

    conv = .false.
    iter = 0

     associate(w => this%w, c => this%c, r => this%r, z => this%z, h => this%h, &
          v => this%v, s => this%s, gam => this%gam, wk1 =>this%wk1)

       norm_fac = 1.0_rp / sqrt(coef%volume)
       call rzero(x%x, n)
       call rzero(gam, this%lgmres + 1)
       call rone(s, this%lgmres)
       call rone(c, this%lgmres)
       call rzero(h, this%lgmres * this%lgmres)
       do while (.not. conv .and. iter .lt. niter)

          if(iter.eq.0) then               
             call copy(r, f, n) 
          else
             call copy(r, f, n)      
             call Ax%compute(w, x%x, coef, x%msh, x%Xh)
             call gs_h%op(w, n, GS_OP_ADD)
             call bc_list_apply(blst, w, n)
             call sub2(r, w, n) 
          end if
          
          gam(1) = sqrt(glsc3(r, r, coef%mult, n))
          if(iter.eq.0) then
             ksp_results%res_start = gam(1) * norm_fac
          end if

          if (gam(1) .eq. 0) return

          rnorm = 0.0_rp
          temp = 1.0_rp / gam(1)
          call cmult2(v(1,1), r, temp, n) 
          do j = 1, this%lgmres
             iter = iter+1
          
             call this%M%solve(z(1,j), v(1,j), n)

             call Ax%compute(w, z(1,j), coef, x%msh, x%Xh)
             call gs_h%op(w, n, GS_OP_ADD)
             call bc_list_apply(blst, w, n)
             
             do l = 1, j
                h(l,j) = 0.0_rp
             enddo

             do i = 0, n, NEKO_BLK_SIZE
                if (i + NEKO_BLK_SIZE .le. n) then
                   do l = 1, j
                      do k = 1, NEKO_BLK_SIZE
                         h(l,j) = h(l,j) + &
                              w(i+k) * v(i+k,l) * coef%mult(i+k,1,1,1)
                      end do
                   end do
                else 
                   do k = 1, n-i
                      do l = 1, j
                         h(l,j) = h(l,j) + &
                              w(i+k) * v(i+k,l) * coef%mult(i+k,1,1,1)
                      end do
                   end do
                end if
             end do
             
             call MPI_Allreduce(h(1,j), wk1, j, &
                  MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
             call copy(h(1,j), wk1, j) 

             alpha2 = 0.0_rp
             do i = 0,n,NEKO_BLK_SIZE
                if (i + NEKO_BLK_SIZE .le. n) then
                   do k = 1, NEKO_BLK_SIZE
                      w_plus(k) = 0.0_rp
                   end do
                   do l = 1,j
                      do k = 1, NEKO_BLK_SIZE
                         w_plus(k) = w_plus(k) - h(l,j) * v(i+k,l)
                      end do
                   end do
                   do k = 1, NEKO_BLK_SIZE
                      w(i+k) = w(i+k) + w_plus(k)
                      alpha2 = alpha2 + w(i+k)**2 * coef%mult(i+k,1,1,1)
                   end do
                else 
                   do k = 1, n-i
                      w_plus(1) = 0.0_rp
                      do l = 1, j
                         w_plus(1) = w_plus(1) - h(l,j) * v(i+k,l)
                      end do
                      w(i+k) = w(i+k) + w_plus(1)
                      alpha2 = alpha2 + (w(i+k)**2) * coef%mult(i+k,1,1,1)
                   end do
                end if
             end do
             
             call MPI_Allreduce(alpha2, temp, 1, &
                  MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
             alpha2 = temp
             alpha = sqrt(alpha2)
             do i=1,j-1
                temp = h(i,j)                   
                h(i  ,j) =  c(i)*temp + s(i) * h(i+1,j)  
                h(i+1,j) = -s(i)*temp + c(i) * h(i+1,j)
             end do
             
             rnorm = 0.0_rp
             if(alpha .eq. 0.0_rp) then 
                conv = .true.
                exit
             end if
             
             lr = sqrt(h(j,j) * h(j,j) + alpha**2)
             temp = 1.0_rp / lr
             c(j) = h(j,j) * temp
             s(j) = alpha  * temp
             h(j,j) = lr
             gam(j+1) = -s(j) * gam(j)
             gam(j)   =  c(j) * gam(j)
             
             rnorm = abs(gam(j+1)) * norm_fac
             if (rnorm .lt. this%abs_tol) then 
                conv = .true.
                exit
             end if
             
             if (iter + 1 .gt. niter) exit
             
             if( j .lt. this%lgmres) then
                temp = 1.0_rp / alpha
                call cmult2(v(1,j+1), w, temp, n)
             end if
             
          end do

          j = min(j, this%lgmres)
          do k = j, 1, -1
             temp = gam(k)
             do i = j, k+1, -1
                temp = temp - h(k,i) * c(i)
             end do
             c(k) = temp / h(k,k)
          end do

          do i = 0, n, NEKO_BLK_SIZE
             if (i + NEKO_BLK_SIZE .le. n) then
                do k = 1, NEKO_BLK_SIZE
                   x_plus(k) = 0.0_rp
                end do
                do l = 1,j
                   do k = 1, NEKO_BLK_SIZE
                      x_plus(k) = x_plus(k) + c(l) * z(i+k,l)
                   end do
                end do
                do k = 1, NEKO_BLK_SIZE
                   x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                end do
             else 
                do k = 1, n-i
                   x_plus(1) = 0.0_rp
                   do l = 1, j
                      x_plus(1) = x_plus(1) + c(l) * z(i+k,l)
                   end do
                   x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                end do
             end if
          end do
       end do

     end associate

     ksp_results%res_final = rnorm
     ksp_results%iter = iter
     
   end function gmres_solve

   !> Initialise a standard GMRES solver (OpenMP version)
   subroutine gmres_omp_init(this, n, M, lgmres, rel_tol, abs_tol)
    class(gmres_omp_t), intent(inout) :: this
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
          
  end subroutine gmres_omp_init

  !> Deallocate a standard GMRES solver (OpenMP version)
  subroutine gmres_omp_free(this)
    class(gmres_omp_t), intent(inout) :: this

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
    
    if (allocated(this%v)) then
       deallocate(this%v)
    end if
    
    if (allocated(this%s)) then
       deallocate(this%s)
    end if
    
    
    if (allocated(this%gam)) then
       deallocate(this%gam)
    end if
    
    if (allocated(this%wk1)) then
       deallocate(this%wk1)
    end if
    
    nullify(this%M)
    
  end subroutine gmres_omp_free

  !> Standard GMRES solve (OpenMP version)
  function gmres_omp_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(gmres_omp_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: i, j, k, iter, ierr 
    real(kind=rp) ::  rnorm, alpha, temp, l, tmp_gam
    real(kind=rp) :: ratio, div0, norm_fac
    logical :: conv
    integer outer

    conv = .false.
    iter = 0
    norm_fac = 1.0_rp / sqrt(coef%volume)
    this%gam(this%lgmres+1) = 0.0_rp
    !$omp parallel
    !$omp do
    do i = 1, n
       x%x(i,1,1,1) = 0.0_rp
    end do
    !$omp end do nowait
    !$omp do
    do i = 1, this%lgmres
       this%s(i) = 1.0_rp
       this%c(i) = 1.0_rp
       this%gam(i) = 0.0_rp
    end do
    !$omp end do

    !$omp do
    do i = 1, this%lgmres**2
       this%h(i, 1) = 0.0_rp
    end do
    !$omp end do
    !$omp end parallel

    outer = 0
    do while (.not. conv .and. iter .lt. niter)
       outer = outer + 1       
       tmp_gam = 0.0_rp ! r^T mult r reduction
       
       !$omp parallel
       !$omp do
       do i = 1, n
          this%r(i) = f(i)
       end do
       !$omp end do          

       if(iter .gt. 0) then          
          call Ax%compute(this%w, x%x, coef, x%msh, x%Xh)
          call gs_h%op(this%w, n, GS_OP_ADD)
          call bc_list_apply(blst, this%w, n)
          !$omp do
          do i = 1, n
             this%r(i) = this%r(i) - this%w(i)
          end do          
          !$omp end do
       endif
       
       !$omp do reduction(+:tmp_gam)
       do i = 1, n
          tmp_gam = tmp_gam + (this%r(i) * coef%mult(i,1,1,1) * this%r(i))
       end do
       !$omp end do
       !$omp single
       call MPI_Allreduce(tmp_gam, this%gam(1), 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
       this%gam(1) = sqrt(this%gam(1))
       if(iter.eq.0) then
          div0 = this%gam(1) * norm_fac
          ksp_results%res_start = div0
       endif
       !$omp end single
       !$omp end parallel

       if ( this%gam(1) .eq. 0) return

       rnorm = 0.0_rp
       temp = 1.0_rp / this%gam(1)
       !$omp parallel do
       do i = 1, n
          this%v(i,1) = temp * this%r(i)
       end do
       !$omp end parallel do
       do j = 1, this%lgmres
          iter = iter+1
          !$omp parallel do
          do i = 1, n
             this%w(i) = this%v(i,j)
          end do
          !$omp end parallel do

          !Apply precond
          call this%M%solve(this%z(1,j), this%w, n)

          !$omp parallel private(i,k)
          call Ax%compute(this%w, this%z(1,j), coef, x%msh, x%Xh)
          call gs_h%op(this%w, n, GS_OP_ADD)
          call bc_list_apply(blst, this%w, n)

          do i = 1, j
             !$omp single
             tmp_gam = 0.0_rp
             !$omp end single
             !$omp do reduction(+:tmp_gam)
             do k = 1, n
                tmp_gam = tmp_gam + this%w(k) * this%v(k,i) * coef%mult(k,1,1,1)
             end do
             !$omp end do
             !$omp single
             this%h(i,j) = tmp_gam
             !$omp end single
          end do
          !$omp single
          
          !Could probably be done inplace...
          call MPI_Allreduce(this%h(1,j), this%wk1, j, &
               MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
          call copy(this%h(1,j), this%wk1, j) 
          !$omp end single
          
          do i = 1, j
             !$omp do
             do k = 1, n
                this%w(k) = this%w(k) - this%h(i,j) * this%v(k,i)
             end do
             !$omp end do
          end do
          
          !apply Givens rotations to new column
          !$omp single
          do i=1,j-1
             temp = this%h(i,j)                   
             this%h(i  ,j) =  this%c(i)*temp + this%s(i)*this%h(i+1,j)  
             this%h(i+1,j) = -this%s(i)*temp + this%c(i)*this%h(i+1,j)
          end do

          alpha = 0.0_rp
          !$omp end single
          !$omp do reduction(+:alpha)
          do i = 1, n
             alpha = alpha + (this%w(i) * coef%mult(i,1,1,1) * this%w(i))
          end do
          !$omp end do
          !$omp single
          call MPI_Allreduce(MPI_IN_PLACE, alpha, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
          alpha = sqrt(alpha)          
          !$omp end single
          !$omp end parallel          

          rnorm = 0.0_rp

          if(alpha .eq. 0.0_rp) then 
            conv = .true.
            exit
          end if
          l = sqrt(this%h(j,j) * this%h(j,j) + alpha**2)
          temp = 1.0_rp / l
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
             temp = 1.0_rp / alpha
             !$omp parallel do
             do i = 1, n
                this%v(i, j+1) = temp * this%w(i)
             end do
             !$omp end parallel do
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
       !$omp parallel private(i)
       do i = 1, j
          !$omp do
          do k = 1, n
             x%x(k,1,1,1) = x%x(k,1,1,1) + this%c(i) * this%z(k,i)
          end do
          !$omp end do
       end do
       !$omp end parallel       
    end do

    ksp_results%res_final = rnorm
    ksp_results%iter = iter
  end function gmres_omp_solve
  
end module gmres
  

