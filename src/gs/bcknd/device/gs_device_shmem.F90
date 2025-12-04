! Copyright (c) 2020-2024, The Neko Authors
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
!> Defines GPU aware MPI gather-scatter communication
module gs_device_shmem
  use num_types, only : rp, c_rp
  use gs_comm, only : gs_comm_t
  use gs_ops
  use stack, only : stack_i4_t
  use htable, only : htable_i4_t
  use device
  use comm, only : pe_size, pe_rank, NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, &
       MPI_MAX, MPI_Sendrecv, MPI_STATUS_IGNORE
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_sizeof, c_int32_t, &
       c_ptr, C_NULL_PTR, c_size_t, c_associated
  implicit none
  private

  !> Buffers for non-blocking communication and packing/unpacking
  type, private :: gs_device_shmem_buf_t
     integer, allocatable :: ndofs(:) !< Number of dofs
     integer, allocatable :: offset(:) !< Offset into buf
     integer, allocatable :: remote_offset(:) !< Offset into buf for remote rank
     integer :: total !< Total number of dofs
     type(c_ptr) :: buf_d = C_NULL_PTR !< Device buffer
     type(c_ptr) :: dof_d = C_NULL_PTR !< Dof mapping for pack/unpack
   contains
     procedure, pass(this) :: init => gs_device_shmem_buf_init
     procedure, pass(this) :: free => gs_device_shmem_buf_free
  end type gs_device_shmem_buf_t

  !> Gather-scatter communication using device SHMEM.
  !! The arrays are indexed per PE like @a send_pe and @ recv_pe.
  type, public, extends(gs_comm_t) :: gs_device_shmem_t
     type(gs_device_shmem_buf_t) :: send_buf
     type(gs_device_shmem_buf_t) :: recv_buf
     type(c_ptr), allocatable :: stream(:)
     type(c_ptr), allocatable :: event(:)
     integer :: nvshmem_counter = 1
     type(c_ptr), allocatable :: notifyDone(:)
     type(c_ptr), allocatable :: notifyReady(:)
   contains
     procedure, pass(this) :: init => gs_device_shmem_init
     procedure, pass(this) :: free => gs_device_shmem_free
     procedure, pass(this) :: nbsend => gs_device_shmem_nbsend
     procedure, pass(this) :: nbrecv => gs_device_shmem_nbrecv
     procedure, pass(this) :: nbwait => gs_device_shmem_nbwait
  end type gs_device_shmem_t


#if defined (HAVE_CUDA) && defined(HAVE_NVSHMEM)

  interface
     subroutine cudamalloc_nvshmem(ptr, size) &
          bind(c, name='cudamalloc_nvshmem')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr
       integer(c_size_t), value :: size
     end subroutine cudamalloc_nvshmem
  end interface

  interface
     subroutine cudafree_nvshmem(ptr) &
          bind(c, name='cudafree_nvshmem')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr
     end subroutine cudafree_nvshmem
  end interface

  interface
     subroutine cuda_gs_pack_and_push(u_d, buf_d, dof_d, offset, n, stream, &
          srank, rbuf_d, roffset, remote_offset, &
          rrank, nvshmem_counter, notifyDone, &
          notifyReady, iter) &
          bind(c, name='cuda_gs_pack_and_push')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset, srank, roffset, rrank, iter
       integer(c_int), value :: nvshmem_counter
       type(c_ptr), value :: u_d, buf_d, dof_d, stream, rbuf_d, notifyDone, notifyReady
       integer(c_int),dimension(*) :: remote_offset
     end subroutine cuda_gs_pack_and_push
  end interface

  interface
     subroutine cuda_gs_pack_and_push_wait(stream, nvshmem_counter, notifyDone) &
          bind(c, name='cuda_gs_pack_and_push_wait')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: nvshmem_counter
       type(c_ptr), value :: stream, notifyDone
     end subroutine cuda_gs_pack_and_push_wait
  end interface

  interface
     subroutine cuda_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name='cuda_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack
  end interface
#endif

contains

  subroutine gs_device_shmem_buf_init(this, pe_order, dof_stack, mark_dupes)
    class(gs_device_shmem_buf_t), intent(inout) :: this
    integer, allocatable, intent(inout) :: pe_order(:)
    type(stack_i4_t), allocatable, intent(inout) :: dof_stack(:)
    logical, intent(in) :: mark_dupes
    integer, allocatable :: dofs(:)
    integer :: i, j, total, max_total
    integer(c_size_t) :: sz
    type(htable_i4_t) :: doftable
    integer :: dupe, marked, k
    real(c_rp) :: rp_dummy
    integer(c_int32_t) :: i4_dummy

    allocate(this%ndofs(size(pe_order)))
    allocate(this%offset(size(pe_order)))
    allocate(this%remote_offset(size(pe_order)))

    do i = 1, size(pe_order)
       this%remote_offset(i)=-1
    end do

    total = 0
    do i = 1, size(pe_order)
       this%ndofs(i) = dof_stack(pe_order(i))%size()
       this%offset(i) = total
       total = total + this%ndofs(i)
    end do

    call MPI_Allreduce(total, max_total, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)

    this%total = total

    sz = c_sizeof(rp_dummy) * max_total
#ifdef HAVE_NVSHMEM
    call cudamalloc_nvshmem(this%buf_d, sz)
#endif

    sz = c_sizeof(i4_dummy) * total
    call device_alloc(this%dof_d, sz)

    if (mark_dupes) call doftable%init(2*total)
    allocate(dofs(total))

    ! Copy from dof_stack into dofs, optionally marking duplicates with doftable
    marked = 0
    do i = 1, size(pe_order)
       ! %array() breaks on cray
       select type (arr => dof_stack(pe_order(i))%data)
       type is (integer)
          do j = 1, this%ndofs(i)
             k = this%offset(i) + j
             if (mark_dupes) then
                if (doftable%get(arr(j), dupe) .eq. 0) then
                   if (dofs(dupe) .gt. 0) then
                      dofs(dupe) = -dofs(dupe)
                      marked = marked + 1
                   end if
                   dofs(k) = -arr(j)
                   marked = marked + 1
                else
                   call doftable%set(arr(j), k)
                   dofs(k) = arr(j)
                end if
             else
                dofs(k) = arr(j)
             end if
          end do
       end select
    end do

    call device_memcpy(dofs, this%dof_d, total, HOST_TO_DEVICE, sync=.true.)

    deallocate(dofs)
    call doftable%free()

  end subroutine gs_device_shmem_buf_init

  subroutine gs_device_shmem_buf_free(this)
    class(gs_device_shmem_buf_t), intent(inout) :: this


    if (allocated(this%ndofs)) deallocate(this%ndofs)
    if (allocated(this%offset)) deallocate(this%offset)

#ifdef HAVE_NVSHMEM
    if (c_associated(this%buf_d)) call cudafree_nvshmem(this%buf_d)
#endif
    if (c_associated(this%dof_d)) call device_free(this%dof_d)

  end subroutine gs_device_shmem_buf_free

  !> Initialise MPI based communication method
  subroutine gs_device_shmem_init(this, send_pe, recv_pe)
    class(gs_device_shmem_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer :: i

    call this%init_order(send_pe, recv_pe)

    call this%send_buf%init(this%send_pe, this%send_dof, .false.)
    call this%recv_buf%init(this%recv_pe, this%recv_dof, .true.)

#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    ! Create a set of non-blocking streams
    allocate(this%stream(size(this%recv_pe)))
    do i = 1, size(this%recv_pe)
       call device_stream_create_with_priority(this%stream(i), 1, STRM_HIGH_PRIO)
    end do

    allocate(this%event(size(this%recv_pe)))
    do i = 1, size(this%recv_pe)
       call device_event_create(this%event(i), 2)
    end do

#ifdef HAVE_NVSHMEM
    allocate(this%notifyDone(size(this%recv_pe)))
    allocate(this%notifyReady(size(this%recv_pe)))
    do i = 1, size(this%recv_pe)
       call cudamalloc_nvshmem(this%notifyDone(i), 8_8)
       call cudamalloc_nvshmem(this%notifyReady(i), 8_8)
    end do
#endif
#endif

  end subroutine gs_device_shmem_init

  !> Deallocate MPI based communication method
  subroutine gs_device_shmem_free(this)
    class(gs_device_shmem_t), intent(inout) :: this
    integer :: i

    call this%send_buf%free()
    call this%recv_buf%free()

    call this%free_order()
    call this%free_dofs()

#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    if (allocated(this%stream)) then
       do i = 1, size(this%stream)
          call device_stream_destroy(this%stream(i))
       end do
       deallocate(this%stream)
    end if
#endif

  end subroutine gs_device_shmem_free

  !> Post non-blocking send operations
  subroutine gs_device_shmem_nbsend(this, u, n, deps, strm)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u)

    do i = 1, size(this%send_pe)
       call device_stream_wait_event(this%stream(i), deps, 0)
       ! Not clear why this sync is required, but there seems to be a race condition
       ! without it for certain run configs
       call device_sync(this%stream(i))
    end do

    ! We do the rest in the "wait" routine below

  end subroutine gs_device_shmem_nbsend

  !> Post non-blocking receive operations
  subroutine gs_device_shmem_nbrecv(this)
    class(gs_device_shmem_t), intent(inout) :: this
    integer :: i

    ! We do everything in the "wait" routine below

  end subroutine gs_device_shmem_nbrecv

  !> Wait for non-blocking operations
  subroutine gs_device_shmem_nbwait(this, u, n, op, strm)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op, done_req, i
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u)
#ifdef HAVE_NVSHMEM
    do i = 1, size(this%send_pe)
       if (this%recv_buf%remote_offset(i) .eq. -1) then
          call MPI_Sendrecv(this%recv_buf%offset(i), 1, MPI_INTEGER, &
               this%recv_pe(i), 0, &
               this%recv_buf%remote_offset(i), 1, MPI_INTEGER, &
               this%send_pe(i), 0, NEKO_COMM, MPI_STATUS_IGNORE)
       end if

       call cuda_gs_pack_and_push(u_d, &
            this%send_buf%buf_d, &
            this%send_buf%dof_d, &
            this%send_buf%offset(i), &
            this%send_buf%ndofs(i), &
            this%stream(i), &
            this%send_pe(i), &
            this%recv_buf%buf_d, &
            this%recv_buf%offset(i), &
            this%recv_buf%remote_offset, &
            this%recv_pe(i), &
            this%nvshmem_counter, &
            this%notifyDone(i), &
            this%notifyReady(i), &
            i)
       this%nvshmem_counter = this%nvshmem_counter + 1
    end do

    do i = 1, size(this%send_pe)
       call cuda_gs_pack_and_push_wait(this%stream(i), &
            this%nvshmem_counter - size(this%send_pe) + i - 1, &
            this%notifyDone(i))
    end do

    do done_req = 1, size(this%recv_pe)
       call cuda_gs_unpack(u_d, op, &
            this%recv_buf%buf_d, &
            this%recv_buf%dof_d, &
            this%recv_buf%offset(done_req), &
            this%recv_buf%ndofs(done_req), &
            this%stream(done_req))
       call device_event_record(this%event(done_req), this%stream(done_req))
    end do

    ! Sync non-blocking streams
    do done_req = 1, size(this%recv_pe)
       call device_stream_wait_event(strm, &
            this%event(done_req), 0)
    end do
#endif
  end subroutine gs_device_shmem_nbwait

end module gs_device_shmem
