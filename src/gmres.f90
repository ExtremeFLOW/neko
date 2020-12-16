!> Defines various GMRES methods
module gmres
  use krylov
  use comm
  use math
  use operators
  use num_types
  implicit none

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: gmres_t
     integer :: lgmres
     real(kind=dp), allocatable :: w(:)
     real(kind=dp), allocatable :: c(:)
     real(kind=dp), allocatable :: r(:)
     real(kind=dp), allocatable :: z(:,:)
     real(kind=dp), allocatable :: h(:,:)
     real(kind=dp), allocatable :: ml(:)
     real(kind=dp), allocatable :: v(:,:)
     real(kind=dp), allocatable :: s(:)
     real(kind=dp), allocatable :: mu(:)
     real(kind=dp), allocatable :: gam(:)
     real(kind=dp), allocatable :: wk1(:)
   contains
     procedure, pass(this) :: init => gmres_init
     procedure, pass(this) :: free => gmres_free
     procedure, pass(this) :: solve => gmres_solve
  end type gmres_t

contains

  !> Initialise a standard GMRES solver
  subroutine gmres_init(this, n, M, lgmres, rel_tol, abs_tol)
    class(gmres_t), intent(inout) :: this
    integer, intent(in) :: n
    class(pc_t), optional, intent(inout), target :: M
    integer, optional, intent(inout) :: lgmres
    real(kind=dp), optional, intent(inout) :: rel_tol
    real(kind=dp), optional, intent(inout) :: abs_tol

    if (present(lgmres)) then
       this%lgmres = lgmres
    else
       this%lgmres = 50
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
       call this%ksp_init()
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
    
  end subroutine gmres_free
 
  !> Standard PCG solve
  function gmres_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(iter)
    class(gmres_t), intent(inout) :: this
    class(ax_t), intent(inout) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(inout) :: n
    real(kind=dp), dimension(n), intent(inout) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, glb_n
    integer :: i, j, k, ierr 
    real(kind=dp) :: rnorm 
    real(kind=dp) ::  alpha, temp, l
    real(kind=dp) :: ratio, div0, norm_fac, tolpss
    logical :: conv
    integer outer

    conv = .false.
    iter  = 0
    glb_n = n / x%msh%nelv * x%msh%glb_nelv

    call rone(this%ml,n)
    call rone(this%mu ,n)
    norm_fac = 1./sqrt(coef%volume)
    ! Should change when doing real problem
    tolpss = 1d-8
    call rzero(x%x,n)
    call rzero(this%gam,this%lgmres+1)
    call rone(this%s,this%lgmres)
    call rone(this%c,this%lgmres)
    call rzero(this%h,this%lgmres*this%lgmres)
    outer = 0
    do while (.not. conv .and. iter .lt. niter)
       outer = outer+1

       if(iter.eq.0) then               !      -1
          call col3(this%r,this%ml,f,n) ! r = L  res
       else
          !update residual
          call copy  (this%r,f,n)           ! r = f
          call Ax%compute(this%w, x%x, coef, x%msh, x%Xh)
          call gs_op(gs_h, this%w, n, GS_OP_ADD)
          call bc_list_apply(blst, this%w, n)
          call add2s2(this%r,this%w,-1d0,n)  ! r = r - w
          call col2(this%r,this%ml,n)        ! r = L   r
       endif
                                                            !            ______
       this%gam(1) = sqrt(glsc3(this%r,this%r,coef%mult,n)) ! gamma  = \/ (r,r) 
                                                            !      1
       if(iter.eq.0) then
          div0 = this%gam(1)*norm_fac
       endif

       if ( this%gam(1) .eq. 0) return

       rnorm = 0.
       temp = 1d0 / this%gam(1)
       call cmult2(this%v(1,1),this%r,temp,n) ! v  = r / gamma
       do j=1,this%lgmres
          iter = iter+1
          call col3(this%w,this%mu,this%v(1,j),n) ! w  = U   v

          !Apply precond
          call this%M%solve(this%z(1,j), this%w, n)

!          call ortho(this%z(1,j),n,glb_n) ! Orthogonalize wrt null space, if present
          call Ax%compute(this%w, this%z(1,j), coef, x%msh, x%Xh)
          call gs_op(gs_h, this%w, n, GS_OP_ADD)
          call bc_list_apply(blst, this%w, n)
          call col2(this%w,this%ml,n)           ! w = L   w

          do i=1,j
             this%h(i,j)=vlsc3(this%w,this%v(1,i),coef%mult,n) ! h    = (w,v )
          enddo                                                !  i,j       i
         
          !Could prorbably be done inplace...
          call MPI_Allreduce(this%h(1,j), this%wk1, j, &
               MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM, ierr)
          call copy(this%h(1,j), this%wk1, j) 

          do i=1,j
             call add2s2(this%w,this%v(1,i),-this%h(i,j),n) ! w = w - h    v
          enddo                                             !          i,j  i

          !apply Givens rotations to new column
          do i=1,j-1
             temp = this%h(i,j)                   
             this%h(i  ,j)=  this%c(i)*temp + this%s(i)*this%h(i+1,j)  
             this%h(i+1,j)= -this%s(i)*temp + this%c(i)*this%h(i+1,j)
          enddo
                                                         !            ______
          alpha = sqrt(glsc3(this%w,this%w,coef%mult,n)) ! alpha =  \/ (w,w)
          rnorm = 0.
          if(alpha.eq.0.) then 
            conv = .true.
            exit
          end if
          l = sqrt(this%h(j,j)*this%h(j,j)+alpha*alpha)
          temp = 1d0 / l
          this%c(j) = this%h(j,j) * temp
          this%s(j) = alpha  * temp
          this%h(j,j) = l
          this%gam(j+1) = -this%s(j) * this%gam(j)
          this%gam(j)   =  this%c(j) * this%gam(j)

          rnorm = abs(this%gam(j+1))*norm_fac
          ratio = rnorm / div0
          !Should maybe change so that we return that we havent converged if iter > niter
          if (rnorm .lt. tolpss) then 
             conv = .true.
             exit
          end if
         
          if (iter+1.gt.niter) exit
          
          if( j .lt. this%lgmres) then
            temp = 1d0 / alpha
            call cmult2(this%v(1,j+1),this%w,temp,n) ! v    = w / alpha
                                                     !  j+1            
          endif
       enddo
       j = min(j, this%lgmres)
       !back substitution
       !     -1
       !c = H   gamma
       do k=j,1,-1
          temp = this%gam(k)
          do i=j,k+1,-1
             temp = temp - this%h(k,i)*this%c(i)
          enddo
          this%c(k) = temp / this%h(k,k)
       enddo
       !sum up Arnoldi vectors
       do i=1,j
          call add2s2(x%x,this%z(1,i),this%c(i),n) ! x = x + c  z
       enddo                                       !          i  i
       if (pe_rank .eq. 0) write(*,*) "current res", rnorm, iter
    enddo
!    call ortho   (x%x, n, glb_n) ! Orthogonalize wrt null space, if present
    if (pe_rank .eq. 0) write(*,*) "Residual:", rnorm, iter
  end function gmres_solve

end module gmres
  

