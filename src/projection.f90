
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

  type projection_t
     real(kind=rp), allocatable :: xx(:,:)
     real(kind=rp), allocatable :: bb(:,:)
     real(kind=rp), allocatable :: xbar(:), bbar(:)
     integer :: m, L
     real(kind=rp) :: tol = real(1d-7,rp)
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
    allocate(this%bbar(n))


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
    if (allocated(this%bbar)) then
       deallocate(this%bbar)
    end if

  end subroutine projection_free

  subroutine project1(this, b, Ax, coef, bclst, gs_h, n)
    class(projection_t) :: this
    integer, intent(inout) :: n
    class(Ax_t), intent(inout) :: Ax    
    class(coef_t), intent(inout) :: coef   
    class(bc_list_t), intent(inout) :: bclst
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp), intent(inout), dimension(n) :: b 
    integer :: i, j, k, ierr
    real(kind=rp) :: work(this%L),alpha(this%L)
    associate(xbar => this%xbar, xx => this%xx, &
              bbar => this%bbar, bb => this%bb)
    
      if (this%m.le.0) return

      !First round of CGS
      do k = 1, this%m 
         alpha(k) = vlsc3(xx(1,k),coef%mult,b,n)
      enddo
      !First one outside loop to avoid zeroing xbar and bbar
      !Could prorbably be done inplace...
      call MPI_Allreduce(alpha, work, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call copy(alpha, work, this%m) 

      call cmult2(xbar,xx(1,1),alpha(1),n)
      call cmult2(bbar,bb(1,1),alpha(1),n)
      call add2s2(b,bb(1,1),-alpha(1),n)
      do k = 2,this%m
         call add2s2(xbar,xx(1,k),alpha(k),n)
         call add2s2(bbar,bb(1,k),alpha(k),n)
         call add2s2(b,bb(1,k),-alpha(k),n)
      enddo
      !Second round of CGS
      do k = 1, this%m
         alpha(k) = vlsc3(xx(1,k),coef%mult,b,n)
      enddo
      call MPI_Allreduce(alpha, work, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call copy(alpha, work, this%m) 
      do k = 1,this%m
         call add2s2(xbar,xx(1,k),alpha(k),n)
         call add2s2(bbar,bb(1,k),alpha(k),n)
         call add2s2(b,bb(1,k),-alpha(k),n)
      enddo 
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
       
    do k = 1, m !First round CGS
       alpha(k) = 0.5*(vlsc3(xx(1,k),w,bb(1,m),n) &
                + vlsc3(bb(1,k),w,xx(1,m),n))
    enddo
    call MPI_Allreduce(alpha, work, this%m, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    call copy(alpha, work, this%m) 
    nrm = sqrt(alpha(m)) !Calculate A-norm of new vector
    do k = 1,m-1
       call add2s2(xx(1,m),xx(1,k),-alpha(k),n)
       call add2s2(bb(1,m),bb(1,k),-alpha(k),n)
    enddo
    
    do k = 1, m-1 !Second round CGS
       beta(k) = 0.5*(vlsc3(xx(1,k),w,bb(1,m),n) &
               + vlsc3(bb(1,k),w,xx(1,m),n))
    enddo
    call MPI_Allreduce(beta, work, this%m-1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    call copy(beta, work, this%m) 
    do k = 1,m-1
       call add2s2(xx(1,m),xx(1,k),-beta(k),n)
       call add2s2(bb(1,m),bb(1,k),-beta(k),n)
       !While we're at it,
       !Sum weights from each round to get the total alpha
       alpha(k) = alpha(k) + beta(k)
    enddo

    !Calculate A-norm of newest solution
    alpha(m) = glsc3(xx(1,m), w, bb(1,m), n) 
    alpha(m) = sqrt(alpha(m))
    !dx and db now stored in last column of xx and bb


    if(alpha(m).gt.this%tol*nrm) then !New vector is linearly independent    
       !Normalize dx and db
       scl1 = 1.0/alpha(m) 
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

    if(b.ne.0d0) then
       h = hypot(a,b) 
       d = 1d0/h
       c = abs(a)*d
       s = sign(d,a)*b
       r = sign(real(1d0,rp),a)*h
    else
       c = 1d0
       s = 0d0
       r = a
    endif
      
    return
  end subroutine givens_rotation

end module projection
