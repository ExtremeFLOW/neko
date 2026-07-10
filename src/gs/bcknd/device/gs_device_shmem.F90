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
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use stack, only : stack_i4_t
  use htable, only : htable_i4_t
  use device
  use comm, only : pe_size, pe_rank, NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_Alltoall, MPI_INTEGER, MPI_MAX
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_sizeof, c_int32_t, c_int64_t, &
       c_ptr, C_NULL_PTR, c_size_t, c_associated
  implicit none
  private

  !> Buffers for non-blocking communication and packing/unpacking
  type, private :: gs_device_shmem_buf_t
     integer, allocatable :: ndofs(:) !< Number of dofs
     integer, allocatable :: offset(:) !< Offset into buf
     integer, allocatable :: remote_offset(:) !< Offset into buf for remote rank
     integer :: total !< Total number of dofs
     !> Elements per parity slab of buf_d (the global max dof count); the
     !! stride of a buf_v_d parity slab is GS_VEC_NC times this. Send buffers
     !! hold two parity slabs, recv buffers one (see nslabs in init).
     integer :: slab_stride
     type(c_ptr) :: buf_d = C_NULL_PTR !< Device buffer
     type(c_ptr) :: buf_v_d = C_NULL_PTR !< Symmetric buffer, fused vector (x nc)
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
     !> Round counter for the rank-indexed signal protocol. Advances once per
     !! gs op (lockstep across ranks, SPMD); all waits use CMP_GE, so no
     !! cross-rank counter matching is required and peer counts may differ
     !! between ranks.
     integer :: iter = 0
     !> Symmetric rank-indexed signal arrays, pe_size slots each (a single
     !! collective allocation per array, so nvshmem_malloc collectivity is
     !! independent of the local peer count). done_sig(r) on our PE is set to
     !! iter by rank r's put_signal when its slab has landed in our recv
     !! buffer; ready_sig(r) on our PE is set to iter by rank r once it has
     !! consumed our round-iter slab (so we may put into its recv slab again).
     type(c_ptr) :: done_sig_d = C_NULL_PTR
     type(c_ptr) :: ready_sig_d = C_NULL_PTR
     !> Records the bulk pack on the main stream in nbsend; every per-peer
     !! stream waits on it, ordering all reads of the shared buffer u before
     !! any unpack write to u.
     type(c_ptr) :: pack_event = C_NULL_PTR
   contains
     procedure, pass(this) :: init => gs_device_shmem_init
     procedure, pass(this) :: free => gs_device_shmem_free
     procedure, pass(this) :: nbsend => gs_device_shmem_nbsend
     procedure, pass(this) :: nbrecv => gs_device_shmem_nbrecv
     procedure, pass(this) :: nbwait => gs_device_shmem_nbwait
     procedure, pass(this) :: nbsend_vec => gs_device_shmem_nbsend_vec
     procedure, pass(this) :: nbrecv_vec => gs_device_shmem_nbrecv_vec
     procedure, pass(this) :: nbwait_vec => gs_device_shmem_nbwait_vec
  end type gs_device_shmem_t


#if defined (HAVE_CUDA) && defined(HAVE_NVSHMEM)

  interface
     subroutine cudamalloc_nvshmem(ptr, size) &
          bind(c, name = 'cudamalloc_nvshmem')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr
       integer(c_size_t), value :: size
     end subroutine cudamalloc_nvshmem
  end interface

  interface
     subroutine cudafree_nvshmem(ptr) &
          bind(c, name = 'cudafree_nvshmem')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr
     end subroutine cudafree_nvshmem
  end interface

  interface
     subroutine cuda_gs_nvshmem_pack(u_d, buf_d, dof_d, boffset, n, stream) &
          bind(c, name = 'cuda_gs_nvshmem_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: boffset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_nvshmem_pack
  end interface

  interface
     subroutine cuda_gs_nvshmem_pack_vec(u_d, buf_d, dof_d, boffset, n, nc, &
          ns, stream) bind(c, name = 'cuda_gs_nvshmem_pack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: boffset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_nvshmem_pack_vec
  end interface

  interface
     subroutine cuda_gs_push(buf_d, offset, n, stream, dest_rank, rbuf_d, &
          roffset, iter, done_d, ready_d, mype) &
          bind(c, name = 'cuda_gs_push')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset, dest_rank, roffset, iter, mype
       type(c_ptr), value :: buf_d, stream, rbuf_d, done_d, ready_d
     end subroutine cuda_gs_push
  end interface

  interface
     subroutine cuda_gs_push_wait(stream, iter, done_d, src_rank) &
          bind(c, name = 'cuda_gs_push_wait')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: iter, src_rank
       type(c_ptr), value :: stream, done_d
     end subroutine cuda_gs_push_wait
  end interface

  interface
     subroutine cuda_gs_post_ready(stream, iter, ready_d, mype, src_rank) &
          bind(c, name = 'cuda_gs_post_ready')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: iter, mype, src_rank
       type(c_ptr), value :: stream, ready_d
     end subroutine cuda_gs_post_ready
  end interface

  interface
     subroutine cuda_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'cuda_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack
  end interface

  interface
     subroutine cuda_gs_unpack_vec(u_d, op, buf_d, dof_d, offset, n, nc, ns, &
          stream) bind(c, name = 'cuda_gs_unpack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack_vec
  end interface
#endif

contains

  !> @param nslabs number of parity slabs (2 for send buffers, which are
  !! round-parity double-buffered so a slab is never repacked while a
  !! previous non-blocking put may still be draining it; 1 for recv buffers,
  !! whose overwrite is gated by the ready signal).
  subroutine gs_device_shmem_buf_init(this, pe_order, dof_stack, mark_dupes, &
       nslabs)
    class(gs_device_shmem_buf_t), intent(inout) :: this
    integer, allocatable, intent(inout) :: pe_order(:)
    type(stack_i4_t), allocatable, intent(inout) :: dof_stack(:)
    logical, intent(in) :: mark_dupes
    integer, intent(in) :: nslabs
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
       this%remote_offset(i) = -1
    end do

    total = 0
    do i = 1, size(pe_order)
       this%ndofs(i) = dof_stack(pe_order(i))%size()
       this%offset(i) = total
       total = total + this%ndofs(i)
    end do

    call MPI_Allreduce(total, max_total, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)

    this%total = total
    this%slab_stride = max_total

    sz = c_sizeof(rp_dummy) * nslabs * max_total
#ifdef HAVE_NVSHMEM
    call cudamalloc_nvshmem(this%buf_d, sz)
    ! Fused vector symmetric buffer, sized for up to GS_VEC_NC components.
    sz = c_sizeof(rp_dummy) * nslabs * GS_VEC_NC * max_total
    call cudamalloc_nvshmem(this%buf_v_d, sz)
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

    call device_memcpy(dofs, this%dof_d, total, HOST_TO_DEVICE, sync = .true.)

    deallocate(dofs)
    call doftable%free()

  end subroutine gs_device_shmem_buf_init

  subroutine gs_device_shmem_buf_free(this)
    class(gs_device_shmem_buf_t), intent(inout) :: this


    if (allocated(this%ndofs)) deallocate(this%ndofs)
    if (allocated(this%offset)) deallocate(this%offset)
    if (allocated(this%remote_offset)) deallocate(this%remote_offset)

#ifdef HAVE_NVSHMEM
    if (c_associated(this%buf_d)) call cudafree_nvshmem(this%buf_d)
    if (c_associated(this%buf_v_d)) call cudafree_nvshmem(this%buf_v_d)
#endif
    if (c_associated(this%dof_d)) call device_free(this%dof_d)

  end subroutine gs_device_shmem_buf_free

  !> Initialise MPI based communication method
  subroutine gs_device_shmem_init(this, send_pe, recv_pe)
    class(gs_device_shmem_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer :: i, nstrm, ierr
    integer, allocatable :: local_offsets(:), remote_offsets(:)
    integer(c_size_t) :: sz
    integer(c_int64_t) :: i64_dummy

    call this%init_order(send_pe, recv_pe)

    ! The gs schedule registers every sharing peer for both directions, so
    ! the two lists are identical and index-aligned. Both the per-peer
    ! stream reuse below (stream(i) serves send_pe(i) and recv_pe(i)) and
    ! the ungated repack of the parity send slab in nbsend (see there)
    ! depend on this; fail loudly if a future schedule breaks it.
    if (size(this%send_pe) .ne. size(this%recv_pe)) then
       call neko_error('gs_device_shmem requires symmetric peer lists')
    end if
    do i = 1, size(this%send_pe)
       if (this%send_pe(i) .ne. this%recv_pe(i)) then
          call neko_error('gs_device_shmem requires aligned peer lists')
       end if
    end do

    ! Send buffers are parity double-buffered (nslabs = 2), recv buffers
    ! single (nslabs = 1); see gs_device_shmem_buf_init.
    call this%send_buf%init(this%send_pe, this%send_dof, .false., 2)
    call this%recv_buf%init(this%recv_pe, this%recv_dof, .true., 1)

    ! Exchange, for every send peer, the offset in the receiver's recv buffer
    ! where our slab must land. A single Alltoall keeps the exchange
    ! deadlock-free for arbitrary, non-uniform peer sets (same rationale as
    ! the host OpenSHMEM backend); the lazy pairwise exchange this replaces
    ! indexed send_pe/recv_pe with the same loop index and required uniform,
    ! index-aligned peer lists on every rank.
    allocate(local_offsets(0:pe_size - 1))
    allocate(remote_offsets(0:pe_size - 1))
    local_offsets = -1
    do i = 1, size(this%recv_pe)
       local_offsets(this%recv_pe(i)) = this%recv_buf%offset(i)
    end do
    call MPI_Alltoall(local_offsets, 1, MPI_INTEGER, &
         remote_offsets, 1, MPI_INTEGER, NEKO_COMM, ierr)
    do i = 1, size(this%send_pe)
       this%send_buf%remote_offset(i) = remote_offsets(this%send_pe(i))
    end do
    deallocate(local_offsets)
    deallocate(remote_offsets)

#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    ! Create a set of non-blocking streams. The per-peer streams, events and
    ! notify signals are indexed over both send_pe (pack-and-push) and
    ! recv_pe (unpack, sync), so size them for the larger of the two peer
    ! lists.
    nstrm = max(size(this%send_pe), size(this%recv_pe))
    allocate(this%stream(nstrm))
    do i = 1, nstrm
       call device_stream_create_with_priority(this%stream(i), 1, &
            STRM_HIGH_PRIO)
    end do

    allocate(this%event(nstrm))
    do i = 1, nstrm
       call device_event_create(this%event(i), 2)
    end do

    call device_event_create(this%pack_event, 2)

#ifdef HAVE_NVSHMEM
    ! Rank-indexed symmetric signal arrays; zero-initialised by
    ! cudamalloc_nvshmem, so the first round's CMP_GE waits on iter-1 = 0
    ! pass immediately.
    sz = c_sizeof(i64_dummy) * pe_size
    call cudamalloc_nvshmem(this%done_sig_d, sz)
    call cudamalloc_nvshmem(this%ready_sig_d, sz)
#endif
#endif

    this%iter = 0
    this%vec_supported = .true.

  end subroutine gs_device_shmem_init

  !> Deallocate MPI based communication method
  subroutine gs_device_shmem_free(this)
    class(gs_device_shmem_t), intent(inout) :: this
    integer :: i

    call this%send_buf%free()
    call this%recv_buf%free()

#ifdef HAVE_NVSHMEM
    ! Collective frees; every rank allocated both arrays exactly once.
    if (c_associated(this%done_sig_d)) then
       call cudafree_nvshmem(this%done_sig_d)
       this%done_sig_d = C_NULL_PTR
    end if
    if (c_associated(this%ready_sig_d)) then
       call cudafree_nvshmem(this%ready_sig_d)
       this%ready_sig_d = C_NULL_PTR
    end if
#endif

    call this%free_order()
    call this%free_dofs()

#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    if (allocated(this%stream)) then
       do i = 1, size(this%stream)
          call device_stream_destroy(this%stream(i))
       end do
       deallocate(this%stream)
    end if

    if (allocated(this%event)) then
       do i = 1, size(this%event)
          call device_event_destroy(this%event(i))
       end do
       deallocate(this%event)
    end if

    if (c_associated(this%pack_event)) then
       call device_event_destroy(this%pack_event)
       this%pack_event = C_NULL_PTR
    end if
#endif

  end subroutine gs_device_shmem_free

  !> Post non-blocking send operations
  subroutine gs_device_shmem_nbsend(this, u, n, tag, deps, strm)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, parity
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u)

#ifdef HAVE_NVSHMEM
    ! One round counter per gs op; advanced here so nbsend and nbwait agree
    ! on the round's parity slab. The rank-indexed signal slots make the
    ! exchange independent of peer counts and peer-list ordering.
    this%iter = this%iter + 1
    parity = mod(this%iter, 2)

    ! Bulk-pack every peer's slab on the main stream, ordered after the
    ! gather (deps), BEFORE any per-peer work. This orders all reads of the
    ! shared buffer u before any unpack write to u. Packing inside the push
    ! kernel, gated on the peer's ready signal, let an unpack from a fast
    ! peer modify u before the pack for a slow peer had read it, so a dof
    ! shared with both peers reached the slow peer already partially
    ! reduced (divergence at large rank counts, where multi-peer dofs and
    ! round-level skew are common).
    !
    ! The pack needs no remote gate: the parity slab was last put in round
    ! iter-2, and this pack is stream-ordered (via strm's event joins at
    ! the end of the previous nbwait) after the unpack of every peer's
    ! round iter-1 slab. Peer i putting its round iter-1 slab here implies,
    ! by stream order on peer i -- its push to us and its unpack of our
    ! data run on one stream, since peer lists are index-aligned (checked
    ! in init) -- that peer i consumed our round iter-2 slab, which in turn
    ! implies that put has drained the parity slab.
    call device_stream_wait_event(strm, deps, 0)
    call cuda_gs_nvshmem_pack(u_d, this%send_buf%buf_d, &
         this%send_buf%dof_d, parity * this%send_buf%slab_stride, &
         this%send_buf%total, strm)
    call device_event_record(this%pack_event, strm)

    ! Order each per-peer stream after the pack; the push kernels in nbwait
    ! read the packed parity slab.
    do i = 1, size(this%send_pe)
       call device_stream_wait_event(this%stream(i), this%pack_event, 0)
    end do
#endif

    ! We do the rest in the "wait" routine below

  end subroutine gs_device_shmem_nbsend

  !> Post non-blocking receive operations
  subroutine gs_device_shmem_nbrecv(this, tag)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: tag
    integer :: i

    ! We do everything in the "wait" routine below

  end subroutine gs_device_shmem_nbrecv

  !> Wait for non-blocking operations
  subroutine gs_device_shmem_nbwait(this, u, n, op, strm)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op, done_req, i, parity
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u)
#ifdef HAVE_NVSHMEM
    parity = mod(this%iter, 2)

    ! Push the packed parity slab to every send peer. The kernel waits for
    ! the peer's ready signal (previous round consumed, so its recv slab is
    ! free) before the put; the pack already happened in nbsend.
    do i = 1, size(this%send_pe)
       call cuda_gs_push(this%send_buf%buf_d, &
            parity * this%send_buf%slab_stride + this%send_buf%offset(i), &
            this%send_buf%ndofs(i), &
            this%stream(i), &
            this%send_pe(i), &
            this%recv_buf%buf_d, &
            this%send_buf%remote_offset(i), &
            this%iter, this%done_sig_d, this%ready_sig_d, pe_rank)
    end do

    ! For every recv peer: wait until its slab has landed, reduce it into u,
    ! then post our ready signal so the peer may start its next round.
    do done_req = 1, size(this%recv_pe)
       call cuda_gs_push_wait(this%stream(done_req), this%iter, &
            this%done_sig_d, this%recv_pe(done_req))
       call cuda_gs_unpack(u_d, op, &
            this%recv_buf%buf_d, &
            this%recv_buf%dof_d, &
            this%recv_buf%offset(done_req), &
            this%recv_buf%ndofs(done_req), &
            this%stream(done_req))
       call cuda_gs_post_ready(this%stream(done_req), this%iter, &
            this%ready_sig_d, pe_rank, this%recv_pe(done_req))
       call device_event_record(this%event(done_req), this%stream(done_req))
    end do

    ! Sync non-blocking streams
    do done_req = 1, size(this%recv_pe)
       call device_stream_wait_event(strm, &
            this%event(done_req), 0)
    end do
#endif
  end subroutine gs_device_shmem_nbwait

  !> Fused nc-component send. @a u is the compact shared device buffer
  !! (component-outer, per-component stride n = nshared). Bulk-packs the
  !! parity slab of the fused vector send buffer; see gs_device_shmem_nbsend
  !! for the ordering and parity-slab rationale.
  subroutine gs_device_shmem_nbsend_vec(this, u, n, nc, tag, deps, strm)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, parity
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u)

#ifdef HAVE_NVSHMEM
    ! Shared round counter with the scalar path; all ranks execute the same
    ! op sequence, so counters stay in lockstep.
    this%iter = this%iter + 1
    parity = mod(this%iter, 2)

    call device_stream_wait_event(strm, deps, 0)
    call cuda_gs_nvshmem_pack_vec(u_d, this%send_buf%buf_v_d, &
         this%send_buf%dof_d, &
         parity * GS_VEC_NC * this%send_buf%slab_stride, &
         this%send_buf%total, nc, n, strm)
    call device_event_record(this%pack_event, strm)

    ! Order each per-peer stream after the pack; the push kernels in
    ! nbwait_vec read the packed parity slab.
    do i = 1, size(this%send_pe)
       call device_stream_wait_event(this%stream(i), this%pack_event, 0)
    end do
#endif

  end subroutine gs_device_shmem_nbsend_vec

  !> No-op: everything happens in nbwait_vec.
  subroutine gs_device_shmem_nbrecv_vec(this, tag, nc)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: tag, nc
  end subroutine gs_device_shmem_nbrecv_vec

  !> Fused nc-component push + unpack (the pack happens in nbsend_vec).
  subroutine gs_device_shmem_nbwait_vec(this, u, n, nc, op, strm)
    class(gs_device_shmem_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op, done_req, i, parity
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u)
#ifdef HAVE_NVSHMEM
    parity = mod(this%iter, 2)

    ! Push the packed nc-component parity slab to every send peer. Offsets
    ! and counts are in elements (nc values per packed position); the recv
    ! buffer holds a single slab, so only the local parity offset applies.
    do i = 1, size(this%send_pe)
       call cuda_gs_push(this%send_buf%buf_v_d, &
            parity * GS_VEC_NC * this%send_buf%slab_stride &
            + nc * this%send_buf%offset(i), &
            nc * this%send_buf%ndofs(i), &
            this%stream(i), &
            this%send_pe(i), &
            this%recv_buf%buf_v_d, &
            nc * this%send_buf%remote_offset(i), &
            this%iter, this%done_sig_d, this%ready_sig_d, pe_rank)
    end do

    ! For every recv peer: wait until its slab has landed, reduce it into u,
    ! then post our ready signal so the peer may start its next round.
    do done_req = 1, size(this%recv_pe)
       call cuda_gs_push_wait(this%stream(done_req), this%iter, &
            this%done_sig_d, this%recv_pe(done_req))
       call cuda_gs_unpack_vec(u_d, op, &
            this%recv_buf%buf_v_d, &
            this%recv_buf%dof_d, &
            this%recv_buf%offset(done_req), &
            this%recv_buf%ndofs(done_req), &
            nc, n, &
            this%stream(done_req))
       call cuda_gs_post_ready(this%stream(done_req), this%iter, &
            this%ready_sig_d, pe_rank, this%recv_pe(done_req))
       call device_event_record(this%event(done_req), this%stream(done_req))
    end do

    ! Sync non-blocking streams
    do done_req = 1, size(this%recv_pe)
       call device_stream_wait_event(strm, &
            this%event(done_req), 0)
    end do
#endif
  end subroutine gs_device_shmem_nbwait_vec

end module gs_device_shmem
