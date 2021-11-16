
!> Project x onto X, the space of old solutions and back again
!! @note In this code we assume that the matrix project for the
!! pressure Ax does not vary in time.

module projection
  use num_types
  use math
  use coefs
  use ax_product
  use bc
  use comm
  use gather_scatter
  implicit none
  private

  type, public ::  projection_t
     real(kind=rp), allocatable :: xx(:,:)
     real(kind=rp), allocatable :: bb(:,:)
     real(kind=rp), allocatable :: xbar(:)
     integer :: m, L
     real(kind=rp) :: tol = 1d-7
   contains
     procedure, pass(this) :: project_on => project1
     procedure, pass(this) :: project_back => project2
     procedure, pass(this) :: init => projection_init
     procedure, pass(this) :: free => projection_free
  end type projection_t

contains

  subroutine projection_init(this, n, L)
    class(projection_t), intent(inout) :: this
    integer, intent(in) :: n
    integer, optional, intent(in) :: L
    
    call this%free()
    
    if (present(L)) then
       this%L = L
    else
       this%L = 20
    end if

    this%m = 0 

    allocate(this%xx(n,this%L))
    allocate(this%bb(n,this%L))
    allocate(this%xbar(n))


  end subroutine projection_init

  subroutine projection_free(this)
    class(projection_t), intent(inout) :: this
    if (allocated(this%xx)) then
       deallocate(this%xx)
    end if
    if (allocated(this%bb)) then
       deallocate(this%bb)
    end if
    if (allocated(this%xbar)) then
       deallocate(this%xbar)
    end if

  end subroutine projection_free

  subroutine project1(this, b, coef, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(coef_t), intent(inout) :: coef   
    real(kind=rp), intent(inout), dimension(n) :: b 
    integer :: i, j, k, ierr
    real(kind=rp) :: work(this%L),alpha(this%L)
    associate(xbar => this%xbar, xx => this%xx, &
              bb => this%bb)
    
      if (this%m.le.0) return

      !First round of CGS
      call rzero(alpha,this%m)
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE,n-i+1)
         do k = 1, this%m 
            alpha(k) = alpha(k) + vlsc3(xx(i,k),coef%mult(i,1,1,1),b(i),j)
         end do
      end do
      !First one outside loop to avoid zeroing xbar and bbar
      call MPI_Allreduce(MPI_IN_PLACE, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call rzero(work,this%m)
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE,n-i+1)
         call cmult2(xbar(i),xx(i,1),alpha(1),j)
         call add2s2(b(i),bb(i,1),-alpha(1),j)
         do k = 2,this%m
            call add2s2(xbar(i),xx(i,k),alpha(k),j)
            call add2s2(b(i),bb(i,k),-alpha(k),j)
         end do
         !Second round of CGS
         do k = 1, this%m
            work(k) = work(k) + vlsc3(xx(i,k),coef%mult(i,1,1,1),b(i),j)
         end do
      end do
      call MPI_Allreduce(work, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE,n-i+1)
         do k = 1,this%m
            call add2s2(xbar(i),xx(i,k),alpha(k),j)
            call add2s2(b(i),bb(i,k),-alpha(k),j)
         end do 
      end do
    end associate
  end subroutine project1

  subroutine project2(this,x,Ax,coef, bclst, gs_h, n)
    class(projection_t) :: this
    integer, intent(inout) :: n
    class(Ax_t), intent(inout) :: Ax    
    class(coef_t), intent(inout) :: coef   
    class(bc_list_t), intent(inout) :: bclst
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp), intent(inout), dimension(n) :: x 
    
    if (this%m.gt.0) call add2(x,this%xbar,n)      ! Restore desired solution
    this%m = min(this%m+1,this%L)
    call copy        (this%xx(1,this%m),x,n)   ! Update (X,B)
    call Ax%compute(this%bb(1,this%m), x, coef, coef%msh, coef%Xh)
    call gs_op_vector(gs_h, this%bb(1,this%m), n, GS_OP_ADD)
    call bc_list_apply_scalar(bclst, this%bb(1,this%m), n)
    call proj_ortho  (this,this%xx,this%bb,coef%mult,n) 
  end subroutine project2
  
  subroutine proj_ortho(this, xx, bb, w, n)
    type(projection_t)  :: this
    integer, intent(inout) :: n
    real(kind=rp), dimension(n,this%L), intent(inout) :: xx, bb
    real(kind=rp), dimension(n), intent(inout) :: w
    real(kind=rp) :: nrm, scl1, scl2, c, s
    real(kind=rp) :: work(this%L), alpha(this%L), beta(this%L)
    integer :: i, j, k, h, ierr
    associate(m => this%m)
    if(m.le.0) return !No vectors to ortho-normalize 

    ! AX = B
    ! Calculate dx, db: dx = x-XX^Tb, db=b-BX^Tb     
    call rzero(alpha,m)
    do i = 1, n, NEKO_BLK_SIZE
       j = min(NEKO_BLK_SIZE,n-i+1)
       do k = 1, m !First round CGS
          alpha(k) = alpha(k) + 0.5_rp * (vlsc3(xx(i,k),w(i),bb(i,m),j) &
                   + vlsc3(bb(i,k),w(i),xx(i,m),j))
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, alpha, this%m, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    nrm = sqrt(alpha(m)) !Calculate A-norm of new vector
    call rzero(beta,m)
    do i = 1, n, NEKO_BLK_SIZE
       j = min(NEKO_BLK_SIZE,n-i+1)
       do k = 1,m-1
          call add2s2(xx(i,m),xx(i,k),-alpha(k),j)
          call add2s2(bb(i,m),bb(i,k),-alpha(k),j)
          beta(k) = beta(k) + 0.5_rp * (vlsc3(xx(i,k),w(i),bb(i,m),j) &
                  + vlsc3(bb(i,k),w(i),xx(i,m),j))
       enddo
    end do
    call MPI_Allreduce(MPI_IN_PLACE, beta, this%m-1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    alpha(m) = 0.0_rp
    do i = 1, n, NEKO_BLK_SIZE
       j = min(NEKO_BLK_SIZE,n-i+1)
       do k = 1,m-1
          call add2s2(xx(i,m),xx(i,k),-beta(k),j)
          call add2s2(bb(i,m),bb(i,k),-beta(k),j)
          alpha(k) = alpha(k) + beta(k)
       end do
       alpha(m) = alpha(m) + vlsc3(xx(i,m),w(i),bb(i,m), j)
    end do

    !alpha(m) = glsc3(xx(1,m), w, bb(1,m), n) 
    call MPI_Allreduce(MPI_IN_PLACE, alpha(m), 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    alpha(m) = sqrt(alpha(m))
    !dx and db now stored in last column of xx and bb

    if(alpha(m).gt.this%tol*nrm) then !New vector is linearly independent    
       !Normalize dx and db
       scl1 = 1.0_rp / alpha(m) 
       call cmult(xx(1,m), scl1, n)   
       call cmult(bb(1,m), scl1, n)   

       !We want to throw away the oldest information
       !The below propagates newest information to first vector.
       !This will make the first vector a scalar 
       !multiple of x.
       do k = m, 2, -1
          h = k - 1   
          call givens_rotation(alpha(h),alpha(k),c,s,nrm)
          alpha(h) = nrm     
          do i = 1, n !Apply rotation to xx and bb
             scl1 = c*xx(i,h) + s*xx(i,k)
             xx(i,k) = -s*xx(i,h) + c*xx(i,k)
             xx(i,h) = scl1       
             scl2 = c*bb(i,h) + s*bb(i,k)
             bb(i,k) = -s*bb(i,h) + c*bb(i,k)    
             bb(i,h) = scl2        
          enddo
       enddo

    else !New vector is not linearly independent, forget about it
       k = m !location of rank deficient column
       if(pe_rank .eq. 0) call neko_warning('New vector not linearly indepependent!')
       m = m - 1 !Remove column
    endif   
    end associate
  end subroutine proj_ortho
  subroutine givens_rotation(a, b, c, s, r)
    real(kind=rp), intent(inout) :: a, b, c, s, r
    real(kind=rp) ::  h, d

    if(b .ne. 0.0_rp) then
       h = hypot(a,b) 
       d = 1.0_rp / h
       c = abs(a) * d
       s = sign(d,a) * b
       r = sign(1.0_rp,a) * h
    else
       c = 1.0_rp
       s = 0.0_rp
       r = a
    endif
      
    return
  end subroutine givens_rotation

end module projection
