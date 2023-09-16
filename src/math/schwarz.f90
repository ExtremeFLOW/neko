! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC. 
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF 
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract 
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE, 
! LLC nor any of their employees, makes any warranty, 
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process, 
! or services by trade name, trademark, manufacturer or otherwise does 
! not necessarily constitute or imply its endorsement, recommendation, 
! or favoring by the United States Government or UCHICAGO ARGONNE LLC. 
! The views and opinions of authors expressed 
! herein do not necessarily state or reflect those of the United States 
! Government or UCHICAGO ARGONNE, LLC, and shall 
! not be used for advertising or product endorsement purposes.
!
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
  use device_schwarz
  use fdm
  use device_math
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private
  
  type, public :: schwarz_t
    real(kind=rp), allocatable :: work1(:)
    real(kind=rp), allocatable :: work2(:)
    real(kind=rp), allocatable :: wt(:,:,:,:,:)
    type(c_ptr) :: work1_d = C_NULL_PTR
    type(c_ptr) :: work2_d = C_NULL_PTR
    type(c_ptr) :: wt_d = C_NULL_PTR
    type(space_t) :: Xh_schwarz !< needed to init gs
    type(gs_t) :: gs_schwarz !< We are only interested in the gather-scatter!
    type(dofmap_t) :: dm_schwarz !< needed to init gs
    type(fdm_t) :: fdm
    type(space_t), pointer :: Xh
    type(bc_list_t), pointer :: bclst
    type(dofmap_t), pointer :: dm
    type(gs_t), pointer :: gs_h
    type(mesh_t), pointer :: msh
    type(c_ptr) :: event
  contains 
    procedure, pass(this) :: init => schwarz_init
    procedure, pass(this) :: free => schwarz_free
    procedure, pass(this) :: compute => schwarz_compute
  end type schwarz_t
 
contains
  
  subroutine schwarz_init(this, Xh, dm, gs_h, bclst, msh)
    class(schwarz_t), target, intent(inout) :: this
    type(space_t), target, intent(inout) :: Xh
    type(dofmap_t), target, intent(inout) :: dm
    type(gs_t), target, intent(inout) :: gs_h
    type(mesh_t), target, intent(inout) :: msh
    type(bc_list_t), target, intent(inout):: bclst

    call this%free()
    
    call space_init(this%Xh_schwarz, GLL, Xh%lx+2, Xh%lx+2, Xh%lx+2)
    this%dm_schwarz = dofmap_t(msh, this%Xh_schwarz) 
    call gs_init(this%gs_schwarz, this%dm_schwarz)

    allocate(this%work1(this%dm_schwarz%size()))
    allocate(this%work2(this%dm_schwarz%size()))
    allocate(this%wt(Xh%lx, Xh%lx, 4, msh%gdim, msh%nelv))
    
    call this%fdm%init(Xh, dm, gs_h)


    this%msh => msh
    this%Xh => Xh
    this%bclst => bclst
    this%dm => dm
    this%gs_h => gs_h
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%work1, this%work1_d,this%dm_schwarz%size()) 
       call device_map(this%work2, this%work2_d,this%dm_schwarz%size()) 
    end if


    call schwarz_setup_wt(this)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_alloc(this%wt_d,int(this%dm%size()*rp, i8)) 
       call rone(this%work1, this%dm%size())
       call schwarz_wt3d(this%work1, this%wt, Xh%lx, msh%nelv)
       call device_memcpy(this%work1, this%wt_d, this%dm%size(), HOST_TO_DEVICE)
       call device_event_create(this%event, 2)
    end if
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
    associate(work1 => this%work1, work2 => this%work2, msh => this%msh, &
         Xh => this%Xh, Xh_schwarz => this%Xh_schwarz)

      n  = this%dm%size()

      enx = Xh_schwarz%lx
      eny = Xh_schwarz%ly
      enz = Xh_schwarz%lz
      if(.not. msh%gdim .eq. 3) enz=1
      ns = enx*eny*enz*msh%nelv
      
      call rone(work2, ns)
      call rzero(work1, ns)
 
      !   Sum overlap region (border excluded)
      !   Cred to PFF for this, very clever
      call schwarz_extrude(work1, 0, zero, work2, 0, one , enx, eny, enz, msh%nelv)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_memcpy(work2, this%work2_d, ns, HOST_TO_DEVICE)
         call gs_op(this%gs_schwarz, work2, ns, GS_OP_ADD) 
         call device_memcpy(work2, this%work2_d, ns, DEVICE_TO_HOST)
      else
         call gs_op(this%gs_schwarz, work2, ns, GS_OP_ADD) 
      end if
      call schwarz_extrude(work2, 0, one, work1, 0, -one, enx, eny, enz, msh%nelv)
      call schwarz_extrude(work2, 2, one, work2, 0, one, enx, eny, enz, msh%nelv)

      ! if(.not.if3d) then ! Go back to regular size array
      !    call hsmg_schwarz_toreg2d(mg_work,mg_work(i),mg_nh(l))
      ! else
      call schwarz_toreg3d(work1, work2, Xh%lx, msh%nelv)
      ! endif
      
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_memcpy(work1, this%work1_d, n, HOST_TO_DEVICE)
         call gs_op(this%gs_h, work1, n, GS_OP_ADD) 
         call device_memcpy(work1, this%work1_d, n, DEVICE_TO_HOST)
      else
          call gs_op(this%gs_h, work1, n, GS_OP_ADD) 
      end if

      k = 1
      do ie = 1,msh%nelv
         if (msh%gdim .eq. 2) then
            call schwarz_setup_schwarz_wt2d_2(this%wt,ie,Xh%lx, work1(k), msh%nelv)
         end if
         if (this%msh%gdim.eq. 3) then
            call schwarz_setup_schwarz_wt3d_2(this%wt,ie,Xh%lx, work1(k), msh%nelv)
            k = k + Xh%lxyz
         end if
      end do
    end associate
  end subroutine schwarz_setup_wt

  !>Setup schwarz weights, 2d, second step
  subroutine schwarz_setup_schwarz_wt2d_2(wt,ie,n,work, nelv)
    integer, intent(in) :: n, nelv
    real(kind=rp), intent(inout) :: wt(n,4,2,nelv)
    real(kind=rp), intent(inout) :: work(n,n)
    integer :: ie,i,j
    do j = 1, n
       wt(j,1,1,ie) = 1.0_rp / work(1,j)
       wt(j,2,1,ie) = 1.0_rp / work(2,j)
       wt(j,3,1,ie) = 1.0_rp / work(n-1,j)
       wt(j,4,1,ie) = 1.0_rp / work(n,j)
    end do
    do i = 1, n
       wt(i,1,2,ie) = 1.0_rp / work(i,1)
       wt(i,2,2,ie) = 1.0_rp / work(i,2)
       wt(i,3,2,ie) = 1.0_rp / work(i,n-1)
       wt(i,4,2,ie) = 1.0_rp / work(i,n)
    end do

    return
  end subroutine schwarz_setup_schwarz_wt2d_2

  !>Setup schwarz weights, 3d, second step
  subroutine schwarz_setup_schwarz_wt3d_2(wt, ie, n, work, nelv)
    integer, intent(in) ::n, nelv, ie
    real(kind=rp), intent(inout) :: wt(n,n,4,3,nelv)
    real(kind=rp), intent(inout) :: work(n,n,n)      
    integer :: i,j,k
    
    do k = 1, n
       do j = 1, n
          wt(j,k,1,1,ie) = 1.0_rp / work(1,j,k)
          wt(j,k,2,1,ie) = 1.0_rp / work(2,j,k)
          wt(j,k,3,1,ie) = 1.0_rp / work(n-1,j,k)
          wt(j,k,4,1,ie) = 1.0_rp / work(n,j,k)
       end do
    end do
    
    do k = 1, n
       do i = 1, n
          wt(i,k,1,2,ie) = 1.0_rp / work(i,1,k)
          wt(i,k,2,2,ie) = 1.0_rp / work(i,2,k)
          wt(i,k,3,2,ie) = 1.0_rp / work(i,n-1,k)
          wt(i,k,4,2,ie) = 1.0_rp / work(i,n,k)
       end do
    end do
    
    do j = 1, n
       do i = 1, n
          wt(i,j,1,3,ie) = 1.0_rp / work(i,j,1)
          wt(i,j,2,3,ie) = 1.0_rp / work(i,j,2)
          wt(i,j,3,3,ie) = 1.0_rp / work(i,j,n-1)
          wt(i,j,4,3,ie) = 1.0_rp / work(i,j,n)
       end do
    end do
    
  end subroutine schwarz_setup_schwarz_wt3d_2

  !> convert array a from extended size to regular
  subroutine schwarz_toreg3d(b, a, n, nelv)
    integer, intent(in) :: n, nelv
    real(kind=rp), intent(inout) :: a(0:n+1, 0:n+1, 0:n+1, nelv)
    real(kind=rp), intent(inout) :: b(n,n,n,nelv)
    integer :: i, j, k, ie
    do ie = 1, nelv
       do k = 1, n
          do j = 1, n
             do i = 1, n
                b(i,j,k,ie) = a(i,j,k,ie)
             end do
          end do
       end do
    end do
  end subroutine schwarz_toreg3d

  !> convert array a from original size to size extended array with border
  subroutine schwarz_toext3d(a, b, n, nelv)
    integer, intent(in) :: n, nelv
    real (kind=rp), intent(inout) :: a(0:n+1,0:n+1,0:n+1,nelv), b(n,n,n,nelv)
    integer :: i,j,k,ie

    call rzero(a, (n+2)*(n+2)*(n+2)*nelv)
    do ie = 1, nelv
       do k = 1, n
          do j = 1, n
             do i = 1, n
                a(i,j,k,ie) = b(i,j,k,ie)
             end do
          end do
       end do
    end do
  end subroutine schwarz_toext3d

  !> Sum values along rows l1, l2 with weights f1, f2 and store along row l1. 
  !! Helps us avoid complicated communcation to get neighbor values.
  !! Simply copy interesting values to the boundary and then do gs_op on extended array.
  subroutine schwarz_extrude(arr1, l1, f1, arr2, l2, f2, nx, ny, nz, nelv)
    integer, intent(in) :: l1, l2, nx, ny, nz, nelv
    real(kind=rp), intent(inout) :: arr1(nx,ny,nz,nelv), arr2(nx,ny,nz,nelv)
    real(kind=rp), intent(in) :: f1, f2
    integer :: i, j, k, ie, i0, i1
    i0=2
    i1=nx-1
    
    if(nz .eq. 1) then
       do ie = 1, nelv
          do j = i0, i1
             arr1(l1+1 ,j,1,ie) = f1*arr1(l1+1 ,j,1,ie) &
                                 +f2*arr2(l2+1 ,j,1,ie)
             arr1(nx-l1,j,1,ie) = f1*arr1(nx-l1,j,1,ie) &
                                 +f2*arr2(nx-l2,j,1,ie)
          end do
          do i = i0, i1
             arr1(i,l1+1 ,1,ie) = f1*arr1(i,l1+1 ,1,ie) &
                                 +f2*arr2(i,l2+1 ,1,ie)
             arr1(i,ny-l1,1,ie) = f1*arr1(i,ny-l1,1,ie) &
                                 +f2*arr2(i,nx-l2,1,ie)
          end do
       end do
    else
       do ie = 1, nelv
          do k = i0, i1
             do j = i0, i1
                arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie) &
                                    +f2*arr2(l2+1 ,j,k,ie)
                arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie) &
                                    +f2*arr2(nx-l2,j,k,ie)
             end do
          end do
          do k = i0, i1
             do i = i0, i1
                arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie) &
                                    +f2*arr2(i,l2+1 ,k,ie)
                arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie) &
                                    +f2*arr2(i,nx-l2,k,ie)
             end do
          end do
          do j = i0, i1
             do i = i0, i1
                arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie) &
                                    +f2*arr2(i,j,l2+1 ,ie)
                arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie) &
                                    +f2*arr2(i,j,nx-l2,ie)
             end do
          end do
       end do
    endif
  end subroutine schwarz_extrude
  
  subroutine schwarz_compute(this, e, r)
    class(schwarz_t), intent(inout) :: this
    real(kind=rp), dimension(this%dm%size()), intent(inout) :: e, r
    integer :: n, enx, eny, enz, ns
    real(kind=rp), parameter :: zero = 0.0_rp
    real(kind=rp), parameter :: one = 1.0_rp
    type(c_ptr) :: e_d, r_d
    associate(work1 => this%work1, work1_d => this%work1_d,&
              work2 => this%work2, work2_d => this%work2_d)

    n  = this%dm%size()
    enx=this%Xh_schwarz%lx
    eny=this%Xh_schwarz%ly
    enz=this%Xh_schwarz%lz
    if(.not. this%msh%gdim .eq. 3) enz=1
    ns = enx*eny*enz*this%msh%nelv
    if (NEKO_BCKND_DEVICE .eq. 1) then
       r_d = device_get_ptr(r)
       e_d = device_get_ptr(e)
       call device_event_record(this%event, glb_cmd_queue)
       call device_stream_wait_event(aux_cmd_queue, this%event, 0)
       call device_schwarz_toext3d(work1_d, r_d, this%Xh%lx, &
                                   this%msh%nelv, aux_cmd_queue)
       call device_schwarz_extrude(work1_d, 0, zero, work1_d, 2, one, &
                                   enx,eny,enz, this%msh%nelv,aux_cmd_queue)

       this%gs_schwarz%bcknd%gs_stream = aux_cmd_queue
       call gs_op(this%gs_schwarz, work1, ns, GS_OP_ADD,this%event) 
       call device_event_sync(this%event)
       call device_schwarz_extrude(work1_d, 0, one, work1_d, 2, -one, &
                                   enx, eny, enz, this%msh%nelv, aux_cmd_queue)
       
       call this%fdm%compute(work2, work1,aux_cmd_queue) ! do local solves

       call device_schwarz_extrude(work1_d, 0, zero, work2_d, 0, one, &
                                   enx, eny, enz, this%msh%nelv, aux_cmd_queue)
       call gs_op(this%gs_schwarz, work2, ns, GS_OP_ADD,this%event)
       call device_event_sync(this%event)

       call device_schwarz_extrude(work2_d, 0, one, work1_d, 0, -one, &
                                   enx, eny, enz, this%msh%nelv, aux_cmd_queue)
       call device_schwarz_extrude(work2_d, 2, one, work2_d, 0, one, &
                                   enx, eny, enz, this%msh%nelv, aux_cmd_queue)
       call device_schwarz_toreg3d(e_d, work2_d, this%Xh%lx, &
                                   this%msh%nelv, aux_cmd_queue)

       call device_event_record(this%event,aux_cmd_queue)
       call device_event_sync(this%event)

       call gs_op(this%gs_h, e, n, GS_OP_ADD, this%event)
       call bc_list_apply_scalar(this%bclst, e, n)
       call device_col2(e_d,this%wt_d, n)
       call device_stream_wait_event(aux_cmd_queue, this%event, 0)
    else
       call bc_list_apply_scalar(this%bclst, r, n)
       call schwarz_toext3d(work1, r, this%Xh%lx, this%msh%nelv)

       !  exchange interior nodes
       call schwarz_extrude(work1, 0, zero, work1, 2, one, &
                            enx, eny, enz, this%msh%nelv)
       call gs_op(this%gs_schwarz, work1, ns, GS_OP_ADD) 
       call schwarz_extrude(work1, 0, one, work1, 2, -one, &
                            enx, eny, enz, this%msh%nelv)
       
       call this%fdm%compute(work2, work1) ! do local solves

       !   Sum overlap region (border excluded)
       call schwarz_extrude(work1, 0, zero, work2, 0, one, &
                            enx, eny, enz, this%msh%nelv)
       call gs_op(this%gs_schwarz, work2, ns, GS_OP_ADD) 
       call schwarz_extrude(work2, 0, one, work1, 0, -one, &
                            enx, eny, enz, this%msh%nelv)
       call schwarz_extrude(work2, 2, one, work2, 0, one, &
                            enx, eny, enz, this%msh%nelv)

       call schwarz_toreg3d(e, work2, this%Xh%lx, this%msh%nelv)

       ! sum border nodes
       call gs_op(this%gs_h, e, n, GS_OP_ADD) 
       call bc_list_apply_scalar(this%bclst, e, n)

       call schwarz_wt3d(e, this%wt, this%Xh%lx, this%msh%nelv)
    end if
  end associate
  end subroutine schwarz_compute

  !Apply schwarz weights along the boundary of each element.
  subroutine schwarz_wt3d(e,wt,n, nelv)
    integer, intent(in) :: n, nelv
    real(kind=rp), intent(inout) :: e(n,n,n,nelv)
    real(kind=rp), intent(inout) ::  wt(n,n,4,3,nelv)
    integer :: ie, i, j, k

    do ie = 1, nelv
       do k = 1, n
          do j = 1, n
             e(1  ,j,k,ie) = e(1  ,j,k,ie) * wt(j,k,1,1,ie)
             e(2  ,j,k,ie) = e(2  ,j,k,ie) * wt(j,k,2,1,ie)
             e(n-1,j,k,ie) = e(n-1,j,k,ie) * wt(j,k,3,1,ie)
             e(n  ,j,k,ie) = e(n  ,j,k,ie) * wt(j,k,4,1,ie)
          end do
       end do
       do k = 1, n
          do i = 3, n-2
             e(i,1  ,k,ie) = e(i,1  ,k,ie) * wt(i,k,1,2,ie)
             e(i,2  ,k,ie) = e(i,2  ,k,ie) * wt(i,k,2,2,ie)
             e(i,n-1,k,ie) = e(i,n-1,k,ie) * wt(i,k,3,2,ie)
             e(i,n  ,k,ie) = e(i,n  ,k,ie) * wt(i,k,4,2,ie)
          end do
       end do
       do j = 3, n-2
          do i = 3, n-2
             e(i,j,1  ,ie) = e(i,j,1  ,ie) * wt(i,j,1,3,ie)
             e(i,j,2  ,ie) = e(i,j,2  ,ie) * wt(i,j,2,3,ie)
             e(i,j,n-1,ie) = e(i,j,n-1,ie) * wt(i,j,3,3,ie)
             e(i,j,n  ,ie) = e(i,j,n  ,ie) * wt(i,j,4,3,ie)
          end do
       end do
    end do
  end subroutine schwarz_wt3d
end module schwarz
