! Copyright (c) 2026, The Neko Authors
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
!> Defines a gather-scatter backend using the native Tofu interconnect (uTofu).
!! Each rank registers its receive buffer once and a send is a single one-sided
!! RDMA put straight into the neighbour's buffer, with the arrival signal fused
!! into the descriptor (`edata`) instead of a separate atomic/credit message.
!! Injection is non-blocking and thread-parallel: each OpenMP thread owns one
!! injection VCQ (created thread-safe, dealt over the TNIs) and fires the puts
!! for its share of the peers itself. The fused vector (multi-component)
!! exchange of gs_op_r3 is supported through a second, independent context
!! with its own buffers, receive VCQ, and parity. See gs_utofu_aux.c for the
!! transport.
module gs_utofu
  use num_types, only : rp, i8
  use gs_comm, only : gs_comm_t, GS_COMM_UTOFU, GS_VEC_NC
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use comm, only : NEKO_COMM, pe_rank
  use logger, only : neko_log, LOG_SIZE
#ifdef HAVE_UTOFU
  use mpi_f08, only : MPI_Request, MPI_Isend, MPI_Irecv, MPI_Waitall, &
       MPI_INTEGER8, MPI_STATUSES_IGNORE
#endif
  use utils, only : neko_error
  !$ use omp_lib, only : omp_get_max_threads, omp_get_thread_num, &
  !$      omp_get_num_threads
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Whether uTofu support was built into this Neko (--with-utofu). Lets
  !! callers (e.g. the gs comm. autotuner) skip the backend rather than
  !! aborting in init on builds without it.
#ifdef HAVE_UTOFU
  logical, parameter, public :: GS_UTOFU_AVAIL = .true.
#else
  logical, parameter, public :: GS_UTOFU_AVAIL = .false.
#endif

  !> MPI tag used for the one-off neighbour metadata exchange at init.
  integer, parameter :: GS_UTOFU_XCHG_TAG = 8123
  !> Ditto for the fused vector path's metadata exchange.
  integer, parameter :: GS_UTOFU_XCHG_TAG_V = 8124

  !> Gather-scatter communication using one-sided uTofu puts.
  type, public, extends(gs_comm_t) :: gs_utofu_t
     !> Concatenated send slabs, one per peer (the put sources). Registered
     !! with uTofu; a pointer (not allocatable) so c_loc can take its address
     !! -- the standard forbids the target attribute on a type component.
     real(kind=rp), pointer :: send_buf(:) => null()
     !> Double-buffered receive slabs: two halves of buf_size elements each,
     !! selected by parity. Registered with uTofu; pointer for the same reason
     !! as send_buf.
     real(kind=rp), pointer :: recv_buf(:) => null()
     !> Per-peer slab lengths and 0-based offsets into send_buf / recv_buf.
     integer, allocatable :: send_len(:), recv_len(:)
     integer, allocatable :: send_offset(:), recv_offset(:)
     !> Int64 views of the send layout, passed to the C transport.
     integer(c_int64_t), allocatable :: send_off64(:), send_len64(:)
     !> Per send-peer remote put targets, learned from each receiver at init:
     !! its receive VCQ id, buffer STADD, the byte-offsets of our slab in its
     !! two parity halves, and the edata tag it wants us to stamp.
     integer(c_int64_t), allocatable :: rmt_vcq_id(:), rmt_stadd(:)
     integer(c_int64_t), allocatable :: dst_off0(:), dst_off1(:)
     integer(c_int), allocatable :: dst_tag(:)
     !> Scratch for the unpack-as-ready loop: receive-peer indices completed
     !! by the latest poll and how many (c_int to match the C transport).
     integer(c_int), allocatable :: recv_indices(:)
     integer(c_int) :: ncompleted = 0
     !> Size in elements of one receive half (= total receive count).
     integer :: buf_size = 0
     !> Size in elements of one send half (= total send count); send_buf is
     !! double-buffered (2 * send_size) so a round can pack while the previous
     !! round's puts are still draining.
     integer :: send_size = 0
     !> Storage size of one element in bytes.
     integer :: elem_size = 0
     !> Active double-buffer half for the next round.
     integer :: parity = 0
     !> Opaque per-instance uTofu context (registrations + notice stash).
     type(c_ptr) :: ctx = c_null_ptr
     !> Fused vector (multi-component) path: double-buffered send/recv slabs
     !! sized for GS_VEC_NC components. Each peer slab holds nc consecutive
     !! component blocks, so the scalar layout is reused scaled by nc. The
     !! path has its OWN uTofu context (registrations, receive VCQ and thus
     !! private MRQ, completion counters) and its own parity: a neighbour
     !! racing from a vector round into a scalar round (or vice versa) must
     !! not land notices in the other path's queue.
     real(kind=rp), pointer :: send_buf_v(:) => null()
     real(kind=rp), pointer :: recv_buf_v(:) => null()
     !> Per send-peer vector-path put targets: the receiver's vector recv
     !! VCQ id and buffer STADD, our slab's scalar element offset within one
     !! of its halves, and its scalar half size (both unscaled; the sender
     !! scales by nc and folds the parity in at round time).
     integer(c_int64_t), allocatable :: rmt_vcq_id_v(:), rmt_stadd_v(:)
     integer(c_int64_t), allocatable :: rmt_off_v(:), rmt_bufsz_v(:)
     !> Per-round scratch, computed by the master each nbsend_vec: the send
     !! layout scaled by nc and the absolute destination element offsets
     !! (vector parity already folded in).
     integer(c_int64_t), allocatable :: send_off64_v(:), send_len64_v(:)
     integer(c_int64_t), allocatable :: dst_off_v(:)
     !> Active vector double-buffer half for the next vector round.
     integer :: parity_v = 0
     !> Opaque uTofu context of the vector path.
     type(c_ptr) :: ctx_v = c_null_ptr
   contains
     procedure, pass(this) :: init => gs_utofu_init
     procedure, pass(this) :: free => gs_utofu_free
     procedure, pass(this) :: nbsend => gs_nbsend_utofu
     procedure, pass(this) :: nbrecv => gs_nbrecv_utofu
     procedure, pass(this) :: nbwait => gs_nbwait_utofu
     procedure, pass(this) :: nbsend_vec => gs_nbsend_vec_utofu
     procedure, pass(this) :: nbrecv_vec => gs_nbrecv_vec_utofu
     procedure, pass(this) :: nbwait_vec => gs_nbwait_vec_utofu
  end type gs_utofu_t

#ifdef HAVE_UTOFU
  !> Largest edata value the hardware will carry (set once at first init).
  integer(c_int64_t) :: gs_utofu_edata_max = 255
  !> Receive VCQs requested per instance path (set once at first init;
  !! instances may be granted fewer when the VCQ budget runs dry).
  integer(c_int) :: gs_utofu_nrvcq = 1
  !> Log the granted VCQ/TNI counts only once, on the first instance.
  logical, save :: gs_utofu_logged = .false.

  interface
     subroutine gs_utofu_init_ep(ntni_req, nthreads, ntni_out, nvcq_out, &
          nrvcq_out, edata_max, ierr) bind(c, name = 'gs_utofu_init')
       import :: c_int, c_int64_t
       integer(c_int), value :: ntni_req, nthreads
       integer(c_int) :: ntni_out, nvcq_out, nrvcq_out
       integer(c_int64_t) :: edata_max
       integer(c_int) :: ierr
     end subroutine gs_utofu_init_ep

     subroutine gs_utofu_ctx_create(send_buf, send_bytes, recv_buf, &
          recv_bytes, ctx, nrvcq, recv_vcq_id, recv_stadd, ierr) &
          bind(c, name = 'gs_utofu_ctx_create')
       import :: c_ptr, c_size_t, c_int64_t, c_int
       type(c_ptr), value :: send_buf, recv_buf
       integer(c_size_t), value :: send_bytes, recv_bytes
       type(c_ptr) :: ctx
       integer(c_int) :: nrvcq
       integer(c_int64_t) :: recv_vcq_id(*)
       integer(c_int64_t) :: recv_stadd(*)
       integer(c_int) :: ierr
     end subroutine gs_utofu_ctx_create

     subroutine gs_utofu_ctx_free(ctx, ierr) bind(c, name = 'gs_utofu_ctx_free')
       import :: c_ptr, c_int
       type(c_ptr), value :: ctx
       integer(c_int) :: ierr
     end subroutine gs_utofu_ctx_free

     subroutine gs_utofu_post_puts(ctx, tid, nteam, nsend, send_off, &
          send_len, rmt_vcq_id, rmt_stadd, dst_off0, dst_off1, dst_tag, &
          parity, send_base, elem_size, ierr) &
          bind(c, name = 'gs_utofu_post_puts')
       import :: c_ptr, c_int, c_int64_t
       type(c_ptr), value :: ctx
       integer(c_int), value :: tid, nteam, nsend, parity, elem_size
       integer(c_int64_t) :: send_off(*), send_len(*)
       integer(c_int64_t) :: rmt_vcq_id(*), rmt_stadd(*)
       integer(c_int64_t) :: dst_off0(*), dst_off1(*)
       integer(c_int) :: dst_tag(*)
       integer(c_int64_t), value :: send_base
       integer(c_int) :: ierr
     end subroutine gs_utofu_post_puts

     subroutine gs_utofu_poll_recv(ctx, parity, cap, ncompleted, out_idx, &
          ierr) bind(c, name = 'gs_utofu_poll_recv')
       import :: c_ptr, c_int
       type(c_ptr), value :: ctx
       integer(c_int), value :: parity, cap
       integer(c_int) :: ncompleted
       integer(c_int) :: out_idx(*)
       integer(c_int) :: ierr
     end subroutine gs_utofu_poll_recv

     subroutine gs_utofu_drain_half(ctx, half, tid, nteam, ierr) &
          bind(c, name = 'gs_utofu_drain_half')
       import :: c_ptr, c_int
       type(c_ptr), value :: ctx
       integer(c_int), value :: half, tid, nteam
       integer(c_int) :: ierr
     end subroutine gs_utofu_drain_half
  end interface
#endif

contains

  !> Initialise the uTofu communication method. See gs_comm.f90 for details.
  subroutine gs_utofu_init(this, send_pe, recv_pe)
    class(gs_utofu_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
#ifdef HAVE_UTOFU
    integer :: i, nsend, nrecv, send_total, recv_total
    integer :: ntni_req, env_len, max_tag, nthreads
    integer(c_int) :: ntni, nvcq, nrv, ierr
    integer(c_int64_t), allocatable :: recv_vcq_id(:), recv_stadd(:)
    integer(c_size_t) :: send_bytes, recv_bytes
    character(len=32) :: env_val
    character(len=LOG_SIZE) :: log_buf
    integer(c_int64_t), allocatable :: msg_out(:,:), msg_in(:,:)
    type(MPI_Request), allocatable :: sreq(:), rreq(:)

    call this%init_order(send_pe, recv_pe)

    nsend = size(this%send_pe)
    nrecv = size(this%recv_pe)
    this%elem_size = int(storage_size(1.0_rp) / 8)

    allocate(this%send_len(nsend), this%send_offset(nsend))
    allocate(this%send_off64(nsend), this%send_len64(nsend))
    allocate(this%rmt_vcq_id(nsend), this%rmt_stadd(nsend))
    allocate(this%dst_off0(nsend), this%dst_off1(nsend), this%dst_tag(nsend))
    allocate(this%recv_len(nrecv), this%recv_offset(nrecv))
    allocate(this%recv_indices(nrecv))

    ! Local receive layout (one concatenated half; the other half mirrors it).
    recv_total = 0
    do i = 1, nrecv
       this%recv_len(i) = this%recv_dof(this%recv_pe(i))%size()
       this%recv_offset(i) = recv_total
       recv_total = recv_total + this%recv_len(i)
    end do
    this%buf_size = recv_total

    ! Local send layout (concatenated per-peer slabs).
    send_total = 0
    do i = 1, nsend
       this%send_len(i) = this%send_dof(this%send_pe(i))%size()
       this%send_offset(i) = send_total
       this%send_off64(i) = int(send_total, c_int64_t)
       this%send_len64(i) = int(this%send_len(i), c_int64_t)
       send_total = send_total + this%send_len(i)
    end do
    this%send_size = send_total

    allocate(this%send_buf(max(1, 2 * send_total)))
    allocate(this%recv_buf(max(1, 2 * recv_total)))

    ! Create the (process-global) uTofu endpoint and register our buffers.
    ntni_req = 1
    call get_environment_variable("NEKO_GS_UTOFU_NTNI", env_val, env_len)
    if (env_len .gt. 0) read(env_val(1:env_len), *) ntni_req

    ! One injection VCQ is requested per OpenMP thread, dealt over the TNIs.
    nthreads = 1
    !$ nthreads = omp_get_max_threads()

    call gs_utofu_init_ep(int(ntni_req, c_int), int(nthreads, c_int), &
         ntni, nvcq, gs_utofu_nrvcq, gs_utofu_edata_max, ierr)
    if (ierr .ne. 0) call neko_error("gs_utofu: endpoint init failed")

    send_bytes = int(size(this%send_buf), c_size_t) * this%elem_size
    recv_bytes = int(size(this%recv_buf), c_size_t) * this%elem_size
    ! Filled by gs_utofu_ctx_create: one (vcq_id, stadd) pair per receive
    ! VCQ granted to this instance; each receive peer is bound to one pair
    ! below. Initialised to keep the compiler from warning (it cannot see
    ! the C side write).
    allocate(recv_vcq_id(gs_utofu_nrvcq), recv_stadd(gs_utofu_nrvcq))
    recv_vcq_id = 0_c_int64_t
    recv_stadd = 0_c_int64_t
    nrv = 1
    call gs_utofu_ctx_create(c_loc(this%send_buf), send_bytes, &
         c_loc(this%recv_buf), recv_bytes, this%ctx, nrv, recv_vcq_id, &
         recv_stadd, ierr)
    if (ierr .eq. 1) then
       call neko_error("gs_utofu: out of VCQs for the receive queues " // &
            "(lower OMP_NUM_THREADS or NEKO_GS_UTOFU_NRVCQ, or set " // &
            "NEKO_GS_UTOFU_NVCQ)")
    else if (ierr .ne. 0) then
       call neko_error("gs_utofu: buffer registration failed")
    end if

    ! Report what the hardware actually granted: requested vs granted can
    ! differ when ranks on a node share its TNIs, or when the VCQ budget
    ! runs dry. Once, on rank 0.
    if (.not. gs_utofu_logged .and. pe_rank .eq. 0) then
       write(log_buf, '(A,I0,A,I0,A,I0,A,I0,A)') 'uTofu inj.   : ', &
            int(nvcq), ' VCQs (', nthreads, ' threads) over ', int(ntni), &
            ' TNIs (', ntni_req, ' requested)'
       call neko_log%message(log_buf)
       write(log_buf, '(A,I0,A,I0,A)') 'uTofu recv   : ', int(nrv), &
            ' of ', int(gs_utofu_nrvcq), ' VCQs (first instance)'
       call neko_log%message(log_buf)
       gs_utofu_logged = .true.
    end if

    ! Tell each sender where, and with which tag, to put our slab; learn the
    ! same for each receiver we send to. Message layout per peer:
    !   [recv_stadd, recv_vcq_id, byte_off_parity0, byte_off_parity1, tag].
    allocate(msg_out(5, max(1, nrecv)), msg_in(5, max(1, nsend)))
    allocate(sreq(max(1, nrecv)), rreq(max(1, nsend)))

    ! Receive peer i is bound to receive VCQ mod(i-1, nrv) + 1; advertising
    ! that VCQ's STADD/id spreads the incoming slabs over the TNIs.
    do i = 1, nrecv
       msg_out(1, i) = recv_stadd(mod(i - 1, int(nrv)) + 1)
       msg_out(2, i) = recv_vcq_id(mod(i - 1, int(nrv)) + 1)
       msg_out(3, i) = int(this%recv_offset(i), c_int64_t)
       msg_out(4, i) = int(this%buf_size + this%recv_offset(i), c_int64_t)
       msg_out(5, i) = int(i - 1, c_int64_t)
    end do

    do i = 1, nsend
       call MPI_Irecv(msg_in(:, i), 5, MPI_INTEGER8, this%send_pe(i), &
            GS_UTOFU_XCHG_TAG, NEKO_COMM, rreq(i))
    end do
    do i = 1, nrecv
       call MPI_Isend(msg_out(:, i), 5, MPI_INTEGER8, this%recv_pe(i), &
            GS_UTOFU_XCHG_TAG, NEKO_COMM, sreq(i))
    end do
    if (nsend .gt. 0) call MPI_Waitall(nsend, rreq, MPI_STATUSES_IGNORE)
    if (nrecv .gt. 0) call MPI_Waitall(nrecv, sreq, MPI_STATUSES_IGNORE)

    max_tag = 0
    do i = 1, nsend
       this%rmt_stadd(i) = msg_in(1, i)
       this%rmt_vcq_id(i) = msg_in(2, i)
       this%dst_off0(i) = msg_in(3, i)
       this%dst_off1(i) = msg_in(4, i)
       this%dst_tag(i) = int(msg_in(5, i), c_int)
       max_tag = max(max_tag, this%dst_tag(i))
    end do

    deallocate(msg_out, msg_in, sreq, rreq)

    ! The edata carries the slab tag in its high bits and the parity bit in
    ! bit 0, so a neighbour's count must fit (2*tag + 1) in the edata range.
    ! The vector path reuses the same tags, so this check covers it too.
    if (int(2 * max_tag + 1, i8) .gt. int(gs_utofu_edata_max, i8)) &
         call neko_error("gs_utofu: neighbour count exceeds edata capacity")

    ! Fused vector path: double-buffered slabs sized for GS_VEC_NC
    ! components, registered under their own context (own receive VCQs).
    allocate(this%send_buf_v(max(1, 2 * GS_VEC_NC * send_total)))
    allocate(this%recv_buf_v(max(1, 2 * GS_VEC_NC * recv_total)))
    allocate(this%rmt_vcq_id_v(nsend), this%rmt_stadd_v(nsend))
    allocate(this%rmt_off_v(nsend), this%rmt_bufsz_v(nsend))
    allocate(this%send_off64_v(nsend), this%send_len64_v(nsend))
    allocate(this%dst_off_v(nsend))

    send_bytes = int(size(this%send_buf_v), c_size_t) * this%elem_size
    recv_bytes = int(size(this%recv_buf_v), c_size_t) * this%elem_size
    recv_vcq_id = 0_c_int64_t
    recv_stadd = 0_c_int64_t
    nrv = 1
    call gs_utofu_ctx_create(c_loc(this%send_buf_v), send_bytes, &
         c_loc(this%recv_buf_v), recv_bytes, this%ctx_v, nrv, recv_vcq_id, &
         recv_stadd, ierr)
    if (ierr .eq. 1) then
       call neko_error("gs_utofu: out of VCQs for the vector receive " // &
            "queues (lower OMP_NUM_THREADS or NEKO_GS_UTOFU_NRVCQ, or " // &
            "set NEKO_GS_UTOFU_NVCQ)")
    else if (ierr .ne. 0) then
       call neko_error("gs_utofu: vector buffer registration failed")
    end if

    ! Second metadata exchange for the vector path. Because the component
    ! count nc varies per round, the receiver advertises its UNSCALED layout
    ! -- base STADD, VCQ id, our slab's scalar offset, and its scalar half
    ! size -- and the sender scales by nc and folds the parity in at round
    ! time. Tags are the scalar ones (same peer ordering).
    allocate(msg_out(4, max(1, nrecv)), msg_in(4, max(1, nsend)))
    allocate(sreq(max(1, nrecv)), rreq(max(1, nsend)))

    ! Same peer-to-receive-VCQ binding as the scalar path (its own VCQ set).
    do i = 1, nrecv
       msg_out(1, i) = recv_stadd(mod(i - 1, int(nrv)) + 1)
       msg_out(2, i) = recv_vcq_id(mod(i - 1, int(nrv)) + 1)
       msg_out(3, i) = int(this%recv_offset(i), c_int64_t)
       msg_out(4, i) = int(this%buf_size, c_int64_t)
    end do

    do i = 1, nsend
       call MPI_Irecv(msg_in(:, i), 4, MPI_INTEGER8, this%send_pe(i), &
            GS_UTOFU_XCHG_TAG_V, NEKO_COMM, rreq(i))
    end do
    do i = 1, nrecv
       call MPI_Isend(msg_out(:, i), 4, MPI_INTEGER8, this%recv_pe(i), &
            GS_UTOFU_XCHG_TAG_V, NEKO_COMM, sreq(i))
    end do
    if (nsend .gt. 0) call MPI_Waitall(nsend, rreq, MPI_STATUSES_IGNORE)
    if (nrecv .gt. 0) call MPI_Waitall(nrecv, sreq, MPI_STATUSES_IGNORE)

    do i = 1, nsend
       this%rmt_stadd_v(i) = msg_in(1, i)
       this%rmt_vcq_id_v(i) = msg_in(2, i)
       this%rmt_off_v(i) = msg_in(3, i)
       this%rmt_bufsz_v(i) = msg_in(4, i)
    end do

    deallocate(msg_out, msg_in, sreq, rreq)

    ! Runtime kill-switch for the fused vector path (validation aid):
    ! NEKO_GS_UTOFU_VEC=0 makes gs_op_r3 fall back to three scalar rounds.
    ! The vector context stays registered either way; only its use is gated.
    call get_environment_variable("NEKO_GS_UTOFU_VEC", env_val, env_len)
    this%vec_supported = .not. (env_len .gt. 0 .and. env_val(1:1) .eq. '0')
#else
    call neko_error("uTofu support not built; reconfigure with --with-utofu")
#endif
  end subroutine gs_utofu_init

  !> Deallocate the uTofu communication method.
  subroutine gs_utofu_free(this)
    class(gs_utofu_t), intent(inout) :: this
#ifdef HAVE_UTOFU
    integer(c_int) :: ierr

    if (c_associated(this%ctx)) then
       call gs_utofu_ctx_free(this%ctx, ierr)
       this%ctx = c_null_ptr
    end if

    if (c_associated(this%ctx_v)) then
       call gs_utofu_ctx_free(this%ctx_v, ierr)
       this%ctx_v = c_null_ptr
    end if

    if (associated(this%send_buf)) then
       deallocate(this%send_buf)
    end if

    if (associated(this%recv_buf)) then
       deallocate(this%recv_buf)
    end if

    if (associated(this%send_buf_v)) then
       deallocate(this%send_buf_v)
    end if

    if (associated(this%recv_buf_v)) then
       deallocate(this%recv_buf_v)
    end if

    if (allocated(this%send_len)) then
       deallocate(this%send_len)
    end if

    if (allocated(this%recv_len)) then
       deallocate(this%recv_len)
    end if

    if (allocated(this%send_offset)) then
       deallocate(this%send_offset)
    end if

    if (allocated(this%recv_offset)) then
       deallocate(this%recv_offset)
    end if

    if (allocated(this%send_off64)) then
       deallocate(this%send_off64)
    end if

    if (allocated(this%send_len64)) then
       deallocate(this%send_len64)
    end if

    if (allocated(this%rmt_vcq_id)) then
       deallocate(this%rmt_vcq_id)
    end if

    if (allocated(this%rmt_stadd)) then
       deallocate(this%rmt_stadd)
    end if

    if (allocated(this%dst_off0)) then
       deallocate(this%dst_off0)
    end if

    if (allocated(this%dst_off1)) then
       deallocate(this%dst_off1)
    end if

    if (allocated(this%dst_tag)) then
       deallocate(this%dst_tag)
    end if

    if (allocated(this%recv_indices)) then
       deallocate(this%recv_indices)
    end if

    if (allocated(this%rmt_vcq_id_v)) then
       deallocate(this%rmt_vcq_id_v)
    end if

    if (allocated(this%rmt_stadd_v)) then
       deallocate(this%rmt_stadd_v)
    end if

    if (allocated(this%rmt_off_v)) then
       deallocate(this%rmt_off_v)
    end if

    if (allocated(this%rmt_bufsz_v)) then
       deallocate(this%rmt_bufsz_v)
    end if

    if (allocated(this%send_off64_v)) then
       deallocate(this%send_off64_v)
    end if

    if (allocated(this%send_len64_v)) then
       deallocate(this%send_len64_v)
    end if

    if (allocated(this%dst_off_v)) then
       deallocate(this%dst_off_v)
    end if

    call this%free_order()
    call this%free_dofs()
#endif
  end subroutine gs_utofu_free

  !> Pack the send slabs and fire all one-sided puts (non-blocking).
  subroutine gs_nbsend_utofu(this, u, n, tag, deps, strm)
    class(gs_utofu_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
#ifdef HAVE_UTOFU
    integer :: i, j, dst, off, ndst, send_base, tid, nteam
    integer(c_int) :: ierr
    integer, pointer :: sp(:)

    ! Entered from inside the gs_op_vector OpenMP parallel region, and every
    ! phase is work-shared: each thread owns one injection VCQ (created
    ! THREAD_SAFE) and fires the puts for its share of the peers itself --
    ! the uTofu substitute for the per-thread MPI_Isend that Fujitsu MPI's
    ! missing MPI_THREAD_MULTIPLE support rules out.
    tid = 0
    nteam = 1
    !$ tid = omp_get_thread_num()
    !$ nteam = omp_get_num_threads()

    ! Reclaim the send-buffer half we are about to repack: the team strides
    ! over the injection VCQs so each one is drained by exactly one thread
    ! (no atomics on the completion counters). The half was last used two
    ! rounds ago, so its puts have long completed -- this just reaps their
    ! TCQ notices. The barrier then guarantees the whole half is free before
    ! anyone packs into it.
    call gs_utofu_drain_half(this%ctx, int(this%parity, c_int), &
         int(tid, c_int), int(nteam, c_int), ierr)
    if (ierr .ne. 0) call neko_error("gs_utofu: send drain failed")
    !$omp barrier

    ! Pack each peer's slab in one work-shared loop. Subroutine locals are
    ! per-thread private (orphaned region); the implicit barrier at end do
    ! guarantees every slab is packed before any thread fires its puts (a
    ! peer packed by one thread may be sent by another -- the pack and post
    ! partitions differ).
    send_base = this%parity * this%send_size
    !$omp do
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       off = this%send_offset(i)
       ndst = this%send_len(i)
       sp => this%send_dof(dst)%array()
       !OCL NORECURRENCE, NOVREC, NOALIAS
       !DIR$ CONCURRENT
       !DIR$ IVDEP
       !GCC$ ivdep
       !NEC$ IVDEP
       !$omp simd
       do j = 1, ndst
          this%send_buf(send_base + off + j) = u(sp(j))
       end do
    end do
    !$omp end do

    ! Fire the one-sided puts: each thread injects its modulo share of the
    ! peers on its own VCQ (or thread 0 injects everything when
    ! NEKO_GS_UTOFU_MASTER_INJECT is set, as an A/B reference). The trailing
    ! barrier holds the team until every put has been issued.
    call gs_utofu_post_puts(this%ctx, int(tid, c_int), int(nteam, c_int), &
         int(size(this%send_pe), c_int), &
         this%send_off64, this%send_len64, this%rmt_vcq_id, this%rmt_stadd, &
         this%dst_off0, this%dst_off1, this%dst_tag, &
         int(this%parity, c_int), int(send_base, c_int64_t), &
         int(this%elem_size, c_int), ierr)
    if (ierr .ne. 0) call neko_error("gs_utofu: put failed")
    !$omp barrier
#else
    call neko_error("uTofu support not built")
#endif
  end subroutine gs_nbsend_utofu

  !> No-op: the one-sided receive buffer is registered once at init, so there
  !! is nothing to post per round.
  subroutine gs_nbrecv_utofu(this, tag)
    class(gs_utofu_t), intent(inout) :: this
    integer, intent(in) :: tag
  end subroutine gs_nbrecv_utofu

  !> Poll for incoming puts, reducing each slab into u as it lands, then wait
  !! for our own sends to drain and flip the double-buffer parity.
  subroutine gs_nbwait_utofu(this, u, n, op, strm)
    class(gs_utofu_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
#ifdef HAVE_UTOFU
    integer :: i, j, k, src, off, nsrc, nreqs, base
    integer(c_int) :: ierr
    integer, pointer :: sp(:)
    ! Bind a local pointer to the receive buffer and index it in the unpack
    ! loops, mirroring how the CAF backend uses its module-level buffer.
    real(kind=rp), pointer :: rbuf(:)

    base = this%parity * this%buf_size
    nreqs = size(this%recv_pe)
    rbuf => this%recv_buf

    do while (nreqs .gt. 0)
       !$omp master
       ! Spin on the MRQ from the master alone until at least one slab lands.
       ! Polling here (rather than around the whole team) keeps the idle
       ! threads from paying a barrier on every empty poll -- the dominant
       ! cost under interconnect latency.
       this%ncompleted = 0
       do while (this%ncompleted .eq. 0)
          call gs_utofu_poll_recv(this%ctx, int(this%parity, c_int), &
               int(size(this%recv_pe), c_int), this%ncompleted, &
               this%recv_indices, ierr)
          if (ierr .ne. 0) call neko_error("gs_utofu: mrq poll failed")
       end do
       !$omp end master
       !$omp barrier

       do k = 1, this%ncompleted
          i = this%recv_indices(k)
          src = this%recv_pe(i)
          off = this%recv_offset(i)
          nsrc = this%recv_len(i)
          sp => this%recv_dof(src)%array()
          select case (op)
          case (GS_OP_ADD)
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ CONCURRENT
             !DIR$ IVDEP
             !GCC$ ivdep
             !NEC$ IVDEP
             !$omp do
             do j = 1, nsrc
                u(sp(j)) = u(sp(j)) + rbuf(base + off + j)
             end do
             !$omp end do
          case (GS_OP_MUL)
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ CONCURRENT
             !DIR$ IVDEP
             !GCC$ ivdep
             !NEC$ IVDEP
             !$omp do
             do j = 1, nsrc
                u(sp(j)) = u(sp(j)) * rbuf(base + off + j)
             end do
             !$omp end do
          case (GS_OP_MIN)
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ CONCURRENT
             !DIR$ IVDEP
             !GCC$ ivdep
             !NEC$ IVDEP
             !$omp do
             do j = 1, nsrc
                u(sp(j)) = min(u(sp(j)), rbuf(base + off + j))
             end do
             !$omp end do
          case (GS_OP_MAX)
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ CONCURRENT
             !DIR$ IVDEP
             !GCC$ ivdep
             !NEC$ IVDEP
             !$omp do
             do j = 1, nsrc
                u(sp(j)) = max(u(sp(j)), rbuf(base + off + j))
             end do
             !$omp end do
          case default
             call neko_error("Unknown operation in gs_nbwait_utofu")
          end select
       end do

       nreqs = nreqs - this%ncompleted
       ! Hold the team until every thread has read ncompleted (above) before
       ! the master re-enters the poll and overwrites recv_indices/ncompleted.
       !$omp barrier
    end do

    !$omp master
    ! Flip to the other double-buffer half for the next round. Send completions
    ! are NOT waited on here; nbsend drains the half it reuses, keeping the
    ! send-completion wait off the critical path.
    this%parity = 1 - this%parity
    !$omp end master
    !$omp barrier
#else
    call neko_error("uTofu support not built")
#endif
  end subroutine gs_nbwait_utofu

  !> Pack the fused nc-component send slabs and fire the puts (non-blocking).
  !! @param u compact shared buffer, component-outer: u((c-1)*n + idx).
  !! Peer i's slab holds nc consecutive blocks of send_len(i) values, so the
  !! scalar layout is reused scaled by nc. Same protocol as gs_nbsend_utofu,
  !! on the vector context and parity.
  subroutine gs_nbsend_vec_utofu(this, u, n, nc, tag, deps, strm)
    class(gs_utofu_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
#ifdef HAVE_UTOFU
    integer :: i, j, c, dst, off, ndst, send_base, tid, nteam
    integer(c_int) :: ierr
    integer, pointer :: sp(:)

    tid = 0
    nteam = 1
    !$ tid = omp_get_thread_num()
    !$ nteam = omp_get_num_threads()

    ! Reclaim the vector send-buffer half we are about to repack (see
    ! gs_nbsend_utofu for the ownership rules).
    call gs_utofu_drain_half(this%ctx_v, int(this%parity_v, c_int), &
         int(tid, c_int), int(nteam, c_int), ierr)
    if (ierr .ne. 0) call neko_error("gs_utofu: vector send drain failed")
    !$omp barrier

    ! Scale the scalar per-peer layout by this round's nc and fold the
    ! vector parity into the destination offsets. Master-only; the pack
    ! loop's implicit end-do barrier orders these writes before any thread
    ! posts.
    send_base = this%parity_v * (GS_VEC_NC * this%send_size)
    !$omp master
    do i = 1, size(this%send_pe)
       this%send_off64_v(i) = int(nc, c_int64_t) * this%send_off64(i)
       this%send_len64_v(i) = int(nc, c_int64_t) * this%send_len64(i)
       this%dst_off_v(i) = int(this%parity_v, c_int64_t) &
            * (GS_VEC_NC * this%rmt_bufsz_v(i)) &
            + int(nc, c_int64_t) * this%rmt_off_v(i)
    end do
    !$omp end master

    ! Pack each peer's slab: nc consecutive component blocks.
    !$omp do
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       off = this%send_offset(i)
       ndst = this%send_len(i)
       sp => this%send_dof(dst)%array()
       do c = 1, nc
          !OCL NORECURRENCE, NOVREC, NOALIAS
          !DIR$ CONCURRENT
          !DIR$ IVDEP
          !GCC$ ivdep
          !NEC$ IVDEP
          !$omp simd
          do j = 1, ndst
             this%send_buf_v(send_base + nc*off + (c-1)*ndst + j) = &
                  u((c-1)*n + sp(j))
          end do
       end do
    end do
    !$omp end do

    ! Fire the puts (see gs_nbsend_utofu). The same dst_off_v serves both
    ! parity slots of the C transport -- the vector parity is already folded
    ! into the offsets; the parity argument still stamps the edata bit.
    call gs_utofu_post_puts(this%ctx_v, int(tid, c_int), &
         int(nteam, c_int), int(size(this%send_pe), c_int), &
         this%send_off64_v, this%send_len64_v, this%rmt_vcq_id_v, &
         this%rmt_stadd_v, this%dst_off_v, this%dst_off_v, this%dst_tag, &
         int(this%parity_v, c_int), int(send_base, c_int64_t), &
         int(this%elem_size, c_int), ierr)
    if (ierr .ne. 0) call neko_error("gs_utofu: vector put failed")
    !$omp barrier
#else
    call neko_error("uTofu support not built")
#endif
  end subroutine gs_nbsend_vec_utofu

  !> No-op: the vector receive buffer is registered once at init, so there
  !! is nothing to post per round.
  subroutine gs_nbrecv_vec_utofu(this, tag, nc)
    class(gs_utofu_t), intent(inout) :: this
    integer, intent(in) :: tag, nc
  end subroutine gs_nbrecv_vec_utofu

  !> Poll for incoming vector puts, reducing each nc-component slab into u
  !! as it lands, then flip the vector double-buffer parity. Same protocol
  !! as gs_nbwait_utofu, on the vector context and parity.
  subroutine gs_nbwait_vec_utofu(this, u, n, nc, op, strm)
    class(gs_utofu_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
#ifdef HAVE_UTOFU
    integer :: i, j, c, k, src, off, nsrc, nreqs, base
    integer(c_int) :: ierr
    integer, pointer :: sp(:)
    real(kind=rp), pointer :: rbuf(:)

    base = this%parity_v * (GS_VEC_NC * this%buf_size)
    nreqs = size(this%recv_pe)
    rbuf => this%recv_buf_v

    do while (nreqs .gt. 0)
       !$omp master
       ! Spin on the vector MRQ from the master alone (see gs_nbwait_utofu).
       this%ncompleted = 0
       do while (this%ncompleted .eq. 0)
          call gs_utofu_poll_recv(this%ctx_v, int(this%parity_v, c_int), &
               int(size(this%recv_pe), c_int), this%ncompleted, &
               this%recv_indices, ierr)
          if (ierr .ne. 0) &
               call neko_error("gs_utofu: vector mrq poll failed")
       end do
       !$omp end master
       !$omp barrier

       do k = 1, this%ncompleted
          i = this%recv_indices(k)
          src = this%recv_pe(i)
          off = this%recv_offset(i)
          nsrc = this%recv_len(i)
          sp => this%recv_dof(src)%array()
          select case (op)
          case (GS_OP_ADD)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + sp(j)) = u((c-1)*n + sp(j)) + &
                        rbuf(base + nc*off + (c-1)*nsrc + j)
                end do
             end do
             !$omp end do
          case (GS_OP_MUL)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + sp(j)) = u((c-1)*n + sp(j)) * &
                        rbuf(base + nc*off + (c-1)*nsrc + j)
                end do
             end do
             !$omp end do
          case (GS_OP_MIN)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + sp(j)) = min(u((c-1)*n + sp(j)), &
                        rbuf(base + nc*off + (c-1)*nsrc + j))
                end do
             end do
             !$omp end do
          case (GS_OP_MAX)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + sp(j)) = max(u((c-1)*n + sp(j)), &
                        rbuf(base + nc*off + (c-1)*nsrc + j))
                end do
             end do
             !$omp end do
          case default
             call neko_error("Unknown operation in gs_nbwait_vec_utofu")
          end select
       end do

       nreqs = nreqs - this%ncompleted
       ! Hold the team until every thread has read ncompleted (above) before
       ! the master re-enters the poll and overwrites recv_indices/ncompleted.
       !$omp barrier
    end do

    !$omp master
    ! Flip to the other vector double-buffer half for the next vector round.
    ! Send completions are reaped by the next nbsend_vec's drain.
    this%parity_v = 1 - this%parity_v
    !$omp end master
    !$omp barrier
#else
    call neko_error("uTofu support not built")
#endif
  end subroutine gs_nbwait_vec_utofu

end module gs_utofu
