!> Overlapping schwarz solves
module schwarz
  use num_types
  use speclib
  use math
  use space
  use dofmap
  use bc
  use dirichlet
  use gather_scatter
  use fast3d
  use fdm
  implicit none  
  type, public :: schwarz_t
    real(kind=rp), allocatable :: work1(:)
    real(kind=rp), allocatable :: work2(:)
    real(kind=rp), allocatable :: wt(:,:,:,:,:)
    type(space_t) :: Xh_schwarz !< needed to init gs
    type(gs_t) :: gs_schwarz !< We are only interested in the gather-scatter!
    type(dofmap_t) :: dm_schwarz !< needed to init gs
    type(fdm_t) :: fdm
    type(space_t), pointer :: Xh
    type(bc_list_t), pointer :: bclst
    type(dofmap_t), pointer :: dm
    type(gs_t), pointer :: gs_h
    type(mesh_t), pointer :: msh
  contains 
    procedure, pass(this) :: init => schwarz_init
    procedure, pass(this) :: free => schwarz_free
    procedure, pass(this) :: compute => schwarz_compute
  end type schwarz_t
contains
  subroutine schwarz_init(this, Xh, dm, gs_h, bclst, msh)
    class(schwarz_t), intent(inout) :: this
    type(space_t), target, intent(inout) :: Xh
    type(dofmap_t), target, intent(inout) :: dm
    type(gs_t), target, intent(inout) :: gs_h
    type(mesh_t), target, intent(inout) :: msh
    type(bc_list_t), target, intent(inout):: bclst
    integer :: nl, n
    call this%free()
    
    call space_init(this%Xh_schwarz, GLL, Xh%lx+2, Xh%lx+2, Xh%lx+2)
    this%dm_schwarz = dofmap_t(msh, this%Xh_schwarz) 
    call gs_init(this%gs_schwarz, this%dm_schwarz)

    allocate(this%work1(this%dm_schwarz%n_dofs))
    allocate(this%work2(this%dm_schwarz%n_dofs))
    allocate(this%wt(Xh%lx, Xh%lx, 4, msh%gdim, msh%nelv))
    
    call fdm_init(this%fdm,Xh, dm, gs_h, bclst)


    this%msh => msh
    this%Xh => Xh
    this%bclst => bclst
    this%dm => dm
    this%gs_h => gs_h


    call schwarz_setup_wt(this)

  end subroutine schwarz_init
 
  subroutine schwarz_free(this)
    class(schwarz_t), intent(inout) :: this
    
    if(allocated(this%work1)) deallocate(this%work1)
    if(allocated(this%work2)) deallocate(this%work2)
    if(allocated(this%wt)) deallocate(this%wt)
    
    call space_free(this%Xh_schwarz)
    call gs_free(this%gs_schwarz)
    !why cant I do this?
    !call dofmap_free(this%dm_schwarz)
    call this%fdm%free()

    nullify(this%Xh)
    nullify(this%bclst)
    nullify(this%dm)
    nullify(this%gs_h)
    nullify(this%msh)
  end subroutine schwarz_free
  !> setup weights 
  subroutine schwarz_setup_wt(this)
    class(schwarz_t), intent(inout) :: this
    integer :: enx,eny,enz, n, ie, k, ns
    real(kind=rp), parameter :: zero = 0.0
    real(kind=rp), parameter :: one = 1.0
    associate(work1 => this%work1, work2 => this%work2)
    n  = this%dm%n_dofs

    enx=this%Xh_schwarz%lx
    eny=this%Xh_schwarz%ly
    enz=this%Xh_schwarz%lz
    if(.not. this%msh%gdim .eq. 3) enz=1
    ns = enx*eny*enz*this%msh%nelv

    call rone(this%work2,ns)
 
!   Sum overlap region (border excluded)
!   Cred to PFF for this, very clever
    call schwarz_extrude(work1,0,zero,work2,0,one ,enx,eny,enz, this%msh%nelv)
    call gs_op_vector(this%gs_schwarz, work2, ns, GS_OP_ADD) 
    call schwarz_extrude(work2,0,one ,work1,0,-one,enx,eny,enz, this%msh%nelv)
    call schwarz_extrude(work2,2,one,work2,0,one,enx,eny,enz, this%msh%nelv)

   ! if(.not.if3d) then ! Go back to regular size array
   !    call hsmg_schwarz_toreg2d(mg_work,mg_work(i),mg_nh(l))
   ! else
       call schwarz_toreg3d(work1,work2,this%Xh%lx, this%msh%nelv)
   ! endif

   call gs_op_vector(this%gs_h, work1, n, GS_OP_ADD) 

    k = 1
    do ie=1,this%msh%nelv
       if (this%msh%gdim .eq. 2) call schwarz_setup_schwarz_wt2d_2(this%wt,ie,this%Xh%lx,work1(k), this%msh%nelv)
       if (this%msh%gdim.eq. 3) call schwarz_setup_schwarz_wt3d_2(this%wt,ie,this%Xh%lx,work1(k), this%msh%nelv)
       k = k+this%Xh%lxyz
    enddo
    end associate
  end subroutine schwarz_setup_wt

  !>Setup schwarz weights, 2d, second step
  subroutine schwarz_setup_schwarz_wt2d_2(wt,ie,n,work, nelv)
    integer, intent(in) :: n, nelv
    real(kind=rp), intent(inout) :: wt(n,4,2,nelv)
    real(kind=rp), intent(inout) :: work(n,n)
    integer :: ie,i,j
    do j=1,n
       wt(j,1,1,ie)=1d0/work(1,j)
       wt(j,2,1,ie)=1d0/work(2,j)
       wt(j,3,1,ie)=1d0/work(n-1,j)
       wt(j,4,1,ie)=1d0/work(n,j)
    enddo
    do i=1,n
       wt(i,1,2,ie)=1d0/work(i,1)
       wt(i,2,2,ie)=1d0/work(i,2)
       wt(i,3,2,ie)=1d0/work(i,n-1)
       wt(i,4,2,ie)=1d0/work(i,n)
    enddo

    return
  end subroutine schwarz_setup_schwarz_wt2d_2

  !>Setup schwarz weights, 3d, second step
  subroutine schwarz_setup_schwarz_wt3d_2(wt,ie,n,work, nelv)
      integer, intent(in) ::n, nelv, ie
      real(kind=rp), intent(inout) :: wt(n,n,4,3,nelv)
      real(kind=rp), intent(inout) :: work(n,n,n)
      
      integer :: i,j,k
      integer :: lbr,rbr,lbs,rbs,lbt,rbt

      do k=1,n
      do j=1,n
         wt(j,k,1,1,ie)=1d0/work(1,j,k)
         wt(j,k,2,1,ie)=1d0/work(2,j,k)
         wt(j,k,3,1,ie)=1d0/work(n-1,j,k)
         wt(j,k,4,1,ie)=1d0/work(n,j,k)
      enddo
      enddo
      do k=1,n
      do i=1,n
         wt(i,k,1,2,ie)=1d0/work(i,1,k)
         wt(i,k,2,2,ie)=1d0/work(i,2,k)
         wt(i,k,3,2,ie)=1d0/work(i,n-1,k)
         wt(i,k,4,2,ie)=1d0/work(i,n,k)
      enddo
      enddo
      do j=1,n
      do i=1,n
         wt(i,j,1,3,ie)=1d0/work(i,j,1)
         wt(i,j,2,3,ie)=1d0/work(i,j,2)
         wt(i,j,3,3,ie)=1d0/work(i,j,n-1)
         wt(i,j,4,3,ie)=1d0/work(i,j,n)
      enddo
      enddo
  end subroutine schwarz_setup_schwarz_wt3d_2

  !> convert array a from extended size to regular
  subroutine schwarz_toreg3d(b,a,n, nelv)
    integer, intent(in) :: n, nelv
    real (kind=rp), intent(inout) :: a(0:n+1,0:n+1,0:n+1,nelv),b(n,n,n,nelv)
    integer :: i,j,k,ie
    do ie=1,nelv
    do k=1,n
    do j=1,n
    do i=1,n
       b(i,j,k,ie)=a(i,j,k,ie)
    enddo
    enddo
    enddo
    enddo
    return
  end subroutine schwarz_toreg3d

  !> convert array a from original size to size extended array with border
  subroutine schwarz_toext3d(a,b,n, nelv)
    integer, intent(in) :: n, nelv
    real (kind=rp), intent(inout) :: a(0:n+1,0:n+1,0:n+1,nelv),b(n,n,n,nelv)
    integer :: i,j,k,ie

    call rzero(a,(n+2)*(n+2)*(n+2)*nelv)
    do ie=1,nelv
    do k=1,n
    do j=1,n
    do i=1,n
       a(i,j,k,ie)=b(i,j,k,ie)
    enddo
    enddo
    enddo
    enddo
  end subroutine schwarz_toext3d

  !> Sum values along rows l1, l2 with weights f1, f2 and store along row l1. 
  !! Helps us avoid complicated communcation to get neighbor values.
  !! Simply copy interesting values to the boundary and then do gs_op on extended array.
  subroutine schwarz_extrude(arr1,l1,f1,arr2,l2,f2,nx,ny,nz, nelv)
    integer, intent(in) :: l1,l2,nx,ny,nz, nelv
    real(kind=rp), intent(inout) :: arr1(nx,ny,nz,nelv),arr2(nx,ny,nz,nelv)
    real(kind=rp), intent(in) :: f1,f2
    integer :: i,j,k,ie,i0,i1
    i0=2
    i1=nx-1
    
    if(nelv .ne. 3) then
       do ie=1,nelv
          do j=i0,i1
             arr1(l1+1 ,j,1,ie) = f1*arr1(l1+1 ,j,1,ie) &
                                 +f2*arr2(l2+1 ,j,1,ie)
             arr1(nx-l1,j,1,ie) = f1*arr1(nx-l1,j,1,ie) &
                                 +f2*arr2(nx-l2,j,1,ie)
          enddo
          do i=i0,i1
             arr1(i,l1+1 ,1,ie) = f1*arr1(i,l1+1 ,1,ie) &
                                 +f2*arr2(i,l2+1 ,1,ie)
             arr1(i,ny-l1,1,ie) = f1*arr1(i,ny-l1,1,ie) &
                                 +f2*arr2(i,nx-l2,1,ie)
          enddo
       enddo
    else
       do ie=1,nelv
          do k=i0,i1
          do j=i0,i1
             arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie) &
                                 +f2*arr2(l2+1 ,j,k,ie)
             arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie) &
                                 +f2*arr2(nx-l2,j,k,ie)
          enddo
          enddo
          do k=i0,i1
          do i=i0,i1
             arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie) &
                                 +f2*arr2(i,l2+1 ,k,ie)
             arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie) &
                                 +f2*arr2(i,nx-l2,k,ie)
          enddo
          enddo
          do j=i0,i1
          do i=i0,i1
             arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie) &
                                 +f2*arr2(i,j,l2+1 ,ie)
             arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie) &
                                 +f2*arr2(i,j,nx-l2,ie)
          enddo
          enddo
       enddo
    endif
  end subroutine schwarz_extrude
  
  subroutine schwarz_compute(this, e, r)
    class(schwarz_t), intent(inout) :: this
    real(kind=rp), dimension(this%dm%n_dofs), intent(inout) :: e, r
    integer :: n, enx, eny, enz, ns
    real(kind=rp) :: zero, one
    associate(work1 => this%work1, work2 => this%work2)

    n  = this%dm%n_dofs
    enx=this%Xh_schwarz%lx
    eny=this%Xh_schwarz%ly
    enz=this%Xh_schwarz%lz
    if(.not. this%msh%gdim .eq. 3) enz=1
    ns = enx*eny*enz*this%msh%nelv
    zero = real(0,rp)
    one = real(1, rp)

    call bc_list_apply_scalar(this%bclst, r, n)
    !if (if3d) then ! extended array 
      call schwarz_toext3d(work1,r,this%Xh%lx, this%msh%nelv)
    !  else
    !     call hsmg_schwarz_toext2d(mg_work,r,mg_nh(l))
    !  endif

 

!  exchange interior nodes
   call schwarz_extrude(work1,0,zero,work1,2,one ,enx,eny,enz, this%msh%nelv)
   call gs_op_vector(this%gs_schwarz, work1, ns, GS_OP_ADD) 
   call schwarz_extrude(work1,0,one,work1,2,-one,enx,eny,enz, this%msh%nelv)
   
   call this%fdm%compute(work2, work1) ! do local solves

   !   Sum overlap region (border excluded)
    call schwarz_extrude(work1,0,zero,work2,0,one,enx,eny,enz, this%msh%nelv)
    call gs_op_vector(this%gs_schwarz, work2, ns, GS_OP_ADD) 
    call schwarz_extrude(work2,0,one,work1,0,-one,enx,eny,enz, this%msh%nelv)
    call schwarz_extrude(work2,2,one,work2,0,one,enx,eny,enz, this%msh%nelv)

   ! if(.not.if3d) then ! Go back to regular size array
   !    call hsmg_schwarz_toreg2d(mg_work,mg_work(i),mg_nh(l))
   ! else
       call schwarz_toreg3d(e,work2,this%Xh%lx, this%msh%nelv)
   ! endif


   ! sum border nodes
   call gs_op_vector(this%gs_h, e, n, GS_OP_ADD) 
   call bc_list_apply_scalar(this%bclst, e, n)

  ! if(.not.if3d) call hsmg_schwarz_wt2d(e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
   !if(if3d) 
   call schwarz_wt3d(e,this%wt,this%Xh%lx,this%msh%nelv)
  end associate
  end subroutine schwarz_compute

  !Apply schwarz weights along the boundary of each element.
  subroutine schwarz_wt3d(e,wt,n, nelv)
    integer, intent(in) :: n, nelv
    real(kind=rp), intent(inout) :: e(n,n,n,nelv)
    real(kind=rp), intent(inout) ::  wt(n,n,4,3,nelv)
    integer :: ie,i,j,k

    do ie=1,nelv
       do k=1,n
       do j=1,n
          e(1  ,j,k,ie)=e(1  ,j,k,ie)*wt(j,k,1,1,ie)
          e(2  ,j,k,ie)=e(2  ,j,k,ie)*wt(j,k,2,1,ie)
          e(n-1,j,k,ie)=e(n-1,j,k,ie)*wt(j,k,3,1,ie)
          e(n  ,j,k,ie)=e(n  ,j,k,ie)*wt(j,k,4,1,ie)
       enddo
       enddo
       do k=1,n
       do i=3,n-2
          e(i,1  ,k,ie)=e(i,1  ,k,ie)*wt(i,k,1,2,ie)
          e(i,2  ,k,ie)=e(i,2  ,k,ie)*wt(i,k,2,2,ie)
          e(i,n-1,k,ie)=e(i,n-1,k,ie)*wt(i,k,3,2,ie)
          e(i,n  ,k,ie)=e(i,n  ,k,ie)*wt(i,k,4,2,ie)
       enddo
       enddo
       do j=3,n-2
       do i=3,n-2
          e(i,j,1  ,ie)=e(i,j,1  ,ie)*wt(i,j,1,3,ie)
          e(i,j,2  ,ie)=e(i,j,2  ,ie)*wt(i,j,2,3,ie)
          e(i,j,n-1,ie)=e(i,j,n-1,ie)*wt(i,j,3,3,ie)
          e(i,j,n  ,ie)=e(i,j,n  ,ie)*wt(i,j,4,3,ie)
       enddo
       enddo
    enddo
  end subroutine schwarz_wt3d
end module schwarz
