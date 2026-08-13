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
!> Defines Coarray Fortran gather-scatter communication
module gs_caf
  use num_types, only : rp
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use comm, only : pe_size
  use, intrinsic :: iso_c_binding
#ifdef HAVE_COARRAY_EVENTS
  use, intrinsic :: iso_fortran_env, only : atomic_int_kind, event_type
#else
  use, intrinsic :: iso_fortran_env, only : atomic_int_kind
#endif
  use utils, only : neko_error
  implicit none
  private

  ! Signaling mode constants. Selected at first init via the
  ! NEKO_GS_CAF_SIGNALING environment variable
  ! ("sync", "atomic", or "event").
  integer, parameter, public :: GS_CAF_SIGNAL_SYNC = 1
  integer, parameter, public :: GS_CAF_SIGNAL_ATOMIC = 2
  integer, parameter, public :: GS_CAF_SIGNAL_EVENT = 3

#ifdef HAVE_COARRAY
  ! Module-level receive coarray, shared by all gs_caf_t instances.
  ! F2008 forbids a derived type from adding a coarray ultimate
  ! component when its parent type has none, so the coarray buffer is
  ! held at module scope rather than as a component of gs_caf_t.
  !
  ! The buffer is double-buffered: it is allocated to twice the global
  ! max receive count so that consecutive rounds write to alternating
  ! halves. In sync mode this eliminates back-pressure entirely (the
  ! receiver may still be unpacking the previous round, but from a
  ! different half, so no overwrite hazard exists). In atomic and event
  ! modes the same property is used to relax the back-pressure spin to
  ! a one-round tolerance -- the next overwrite is two rounds away, so
  ! the receiver only needs to be at most one round behind.
  ! gs_caf_buf_size is the size of one half.
  !
  ! Multiple gs_caf_t instances may coexist (each carrying its own
  ! offset bookkeeping) provided they are used strictly sequentially
  ! -- no overlapping nbsend/nbwait rounds across instances. The buffer
  ! is grown on demand to fit the largest gs ever initialised; it is
  ! never shrunk and is retained for the program lifetime.
  real(kind=rp), allocatable :: gs_caf_recv_buf(:)[:]
  integer :: gs_caf_buf_size = 0

  ! Active signaling mode; bound on the first gs_caf_t init from the
  ! NEKO_GS_CAF_SIGNALING environment variable. Subsequent instances
  ! must use the same mode (the env var is read once).
  integer :: gs_caf_mode = 0

  ! Atomic-mode signaling counters, indexed by remote rank.
  ! gs_caf_data_ready(s_rank) on image r counts rounds image s has put
  ! into r so far. gs_caf_buf_ready(r_rank) on image s counts rounds
  ! image r has finished unpacking from s so far. Allocated only in
  ! atomic mode and shared by all instances.
  integer(kind=atomic_int_kind), allocatable :: gs_caf_data_ready(:)[:]
  integer(kind=atomic_int_kind), allocatable :: gs_caf_buf_ready(:)[:]

  ! Local caches of "rounds we have sent to / received from each remote
  ! rank" -- size pe_size per image, indexed by remote rank. Updated
  ! locally on every atomic_define / wait completion in nbsend / nbwait.
  ! Reading these to baseline a new gs_caf_t avoids any remote
  ! atomic_ref during init, which Cray CCE has historically deadlocked
  ! on. The values match the remote atomic counters at quiescent
  ! points (i.e. between gs ops) for symmetric, lockstep gs traffic.
  integer, allocatable :: gs_caf_send_count(:)
  integer, allocatable :: gs_caf_recv_count(:)

#ifdef HAVE_COARRAY_EVENTS
  ! Event-mode signaling. The data_ready event accumulates one post per
  ! sender per round; buf_ready is the back-channel from receiver to
  ! sender. Events are scalar coarrays whose count cannot distinguish
  ! posts coming from different gs_caf_t instances, so event mode is
  ! restricted to a single live instance at a time.
  !
  ! The events are allocatable rather than static module-scope coarrays
  ! because Cray CCE has historically had layout issues with mixing
  ! module-scope static coarrays of derived type with allocatable
  ! coarrays on the symmetric heap; an explicit allocate side-steps
  ! that.
  type(event_type), allocatable :: gs_caf_data_ready_ev[:]
  type(event_type), allocatable :: gs_caf_buf_ready_ev[:]
  logical :: gs_caf_event_in_use = .false.
#endif
#endif

  !> Gather-scatter communication using Coarray Fortran (F2008).
  !! Each image puts directly into the (module-level) receive coarray on
  !! the destination image. The signaling mode (`NEKO_GS_CAF_SIGNALING`)
  !! selects how segment ordering is established: `sync` uses `sync
  !! images` over the union of send and recv peers, `atomic` uses
  !! per-peer `atomic_define`/`atomic_ref` on `data_ready`/`buf_ready`
  !! counters, and `event` uses coarray events.
  !!
  !! Under OpenMP the pack and unpack loops are work-shared across the
  !! calling team, while every coarray operation is funnelled through the
  !! master thread -- see gs_nbsend_caf.
  type, public, extends(gs_comm_t) :: gs_caf_t
     !> Local gathered send buffer (concatenated slabs, one per send peer).
     real(kind=rp), allocatable :: send_buf(:)
     !> Number of dofs to send to / receive from each peer.
     integer, allocatable :: send_len(:), recv_len(:)
     !> 0-based offset into send_buf / recv_buf for each peer.
     integer, allocatable :: send_offset(:), recv_offset(:)
     !> 0-based offset in the remote peer's recv_buf where our slab is placed.
     integer, allocatable :: dest_offset(:)
     !> 1-based image numbers for the send and recv peer arrays.
     integer, allocatable :: send_img(:), recv_img(:)
     !> sync-mode only: image numbers of the union of send and recv peers,
     !! used for the pairwise sync images that bracket the put. Unallocated
     !! in atomic and event modes.
     integer, allocatable :: sync_img(:)
     !> event-mode only: false on the very first nbsend so the buf_ready
     !! wait is skipped (there are no credits posted yet).
     logical :: send_started = .false.
     !> Double-buffer parity (0 or 1), flipped after every nbwait. The
     !! current round writes to and reads from the parity*buf_size half
     !! of gs_caf_recv_buf, the next round uses the other half.
     integer :: parity = 0
   contains
     procedure, pass(this) :: init => gs_caf_init
     procedure, pass(this) :: free => gs_caf_free
     procedure, pass(this) :: nbsend => gs_nbsend_caf
     procedure, pass(this) :: nbrecv => gs_nbrecv_caf
     procedure, pass(this) :: nbwait => gs_nbwait_caf
     procedure, pass(this) :: nbsend_vec => gs_nbsend_vec_caf
     procedure, pass(this) :: nbrecv_vec => gs_nbrecv_vec_caf
     procedure, pass(this) :: nbwait_vec => gs_nbwait_vec_caf
  end type gs_caf_t

contains

  !> Initialise Coarray Fortran based communication method
  subroutine gs_caf_init(this, send_pe, recv_pe)
    class(gs_caf_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
#ifdef HAVE_COARRAY
    integer, allocatable :: dest_xchg(:)[:]
    logical, allocatable :: in_neigh(:)
    integer :: i, nsend, nrecv, send_total, recv_total, max_total, n_neigh
    integer :: me, env_len
    character(len=64) :: env_val

    ! Bind the signaling mode on the first init.
    if (gs_caf_mode .eq. 0) then
       call get_environment_variable("NEKO_GS_CAF_SIGNALING", env_val, env_len)
       if (env_len .gt. 0 .and. env_val(1:env_len) .eq. "atomic") then
          gs_caf_mode = GS_CAF_SIGNAL_ATOMIC
          allocate(gs_caf_data_ready(0:pe_size - 1)[*])
          allocate(gs_caf_buf_ready(0:pe_size - 1)[*])
          allocate(gs_caf_send_count(0:pe_size - 1))
          allocate(gs_caf_recv_count(0:pe_size - 1))
          gs_caf_send_count = 0
          gs_caf_recv_count = 0
          ! F2008 forbids mixing atomic and non-atomic accesses on the
          ! same variable, so initialise via atomic_define rather than
          ! a regular array assignment.
          do i = 0, pe_size - 1
             call atomic_define(gs_caf_data_ready(i), 0_atomic_int_kind)
             call atomic_define(gs_caf_buf_ready(i), 0_atomic_int_kind)
          end do
       else if (env_len .gt. 0 .and. env_val(1:env_len) .eq. "event") then
#ifdef HAVE_COARRAY_EVENTS
          gs_caf_mode = GS_CAF_SIGNAL_EVENT
#else
          call neko_error("NEKO_GS_CAF_SIGNALING=event requires a Fortran " // &
               "compiler with coarray events support")
#endif
       else
          gs_caf_mode = GS_CAF_SIGNAL_SYNC
       end if
    end if

#ifdef HAVE_COARRAY_EVENTS
    ! Allocate the shared event coarrays once, lazily, on the first
    ! event-mode init. The in-use guard against overlapping gs ops is
    ! enforced in nbsend/nbwait, not here -- multiple gs_caf_t instances
    ! may be initialised back-to-back.
    if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       if (.not. allocated(gs_caf_data_ready_ev)) then
          allocate(gs_caf_data_ready_ev[*])
          allocate(gs_caf_buf_ready_ev[*])
       end if
    end if
#endif

    call this%init_order(send_pe, recv_pe)

    nsend = size(this%send_pe)
    nrecv = size(this%recv_pe)

    allocate(this%send_len(nsend), this%send_offset(nsend), &
         this%send_img(nsend), this%dest_offset(nsend))
    allocate(this%recv_len(nrecv), this%recv_offset(nrecv), &
         this%recv_img(nrecv))

    ! Local receive layout
    recv_total = 0
    do i = 1, nrecv
       this%recv_len(i) = this%recv_dof(this%recv_pe(i))%size()
       this%recv_offset(i) = recv_total
       recv_total = recv_total + this%recv_len(i)
       this%recv_img(i) = this%recv_pe(i) + 1
    end do

    ! Local send layout (concatenated per-peer slabs in one buffer)
    send_total = 0
    do i = 1, nsend
       this%send_len(i) = this%send_dof(this%send_pe(i))%size()
       this%send_offset(i) = send_total
       send_total = send_total + this%send_len(i)
       this%send_img(i) = this%send_pe(i) + 1
    end do
    ! Sized for up to GS_VEC_NC components; scalar path uses the first
    ! send_total elements, the fused vector path uses nc*send_total.
    allocate(this%send_buf(max(1, GS_VEC_NC*send_total)))

    ! Symmetric coarray sized to twice the global max total receive
    ! count (double buffering). gs_caf_buf_size tracks the size of one
    ! half. Grow the shared buffer on demand; allocate / deallocate of
    ! an allocatable coarray is implicitly collective and acts as a
    ! global sync.
    max_total = recv_total
    call co_max(max_total)
    max_total = max(1, max_total)
    ! Half size = GS_VEC_NC * max_total, so the single double-buffered coarray
    ! serves both the scalar path (which uses only the low max_total of each
    ! half) and the fused vector path (which uses the full half). gs_caf_buf_size
    ! is the half size; half_off = parity * gs_caf_buf_size in both paths.
    if (GS_VEC_NC * max_total .gt. gs_caf_buf_size) then
       if (allocated(gs_caf_recv_buf)) deallocate(gs_caf_recv_buf)
       allocate(gs_caf_recv_buf(2 * GS_VEC_NC * max_total)[*])
       gs_caf_buf_size = GS_VEC_NC * max_total
    end if
    this%vec_supported = .true.

    ! Tell each sender at what offset in our recv_buf to place their slab,
    ! and learn at what offset in each receiver's recv_buf our slab should go.
    ! Each image puts its own offset for each sender into a slot on the
    ! sender's image indexed by our rank; after sync_all, each image reads
    ! the offsets directly from its local copy.
    me = this_image()
    allocate(dest_xchg(0:pe_size - 1)[*])
    do i = 1, nrecv
       dest_xchg(me - 1)[this%recv_img(i)] = this%recv_offset(i)
    end do
    sync all
    do i = 1, nsend
       this%dest_offset(i) = dest_xchg(this%send_pe(i))
    end do
    deallocate(dest_xchg)

    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       ! Sync image set = union of send and recv peers. Both endpoints of
       ! every neighbour pair must include each other so the pairwise
       ! sync images statements match up.
       allocate(in_neigh(0:pe_size - 1))
       in_neigh = .false.
       do i = 1, nsend
          in_neigh(this%send_pe(i)) = .true.
       end do
       do i = 1, nrecv
          in_neigh(this%recv_pe(i)) = .true.
       end do
       n_neigh = count(in_neigh)
       allocate(this%sync_img(n_neigh))
       n_neigh = 0
       do i = 0, pe_size - 1
          if (in_neigh(i)) then
             n_neigh = n_neigh + 1
             this%sync_img(n_neigh) = i + 1
          end if
       end do
       deallocate(in_neigh)
    end if ! atomic & event modes: no per-instance state to allocate

    ! Ensure recv_buf is allocated and (atomic mode) baselines are stable
    ! on every image before any signalling activity begins.
    sync all
#else
    call neko_error("Coarray Fortran support not built; reconfigure with " // &
         "a coarray-capable Fortran compiler")
#endif
  end subroutine gs_caf_init

  !> Deallocate Coarray Fortran based communication method.
  !! The shared module-level recv coarray is intentionally retained so
  !! it can be reused by subsequent gs_caf_t instances.
  subroutine gs_caf_free(this)
    class(gs_caf_t), intent(inout) :: this
#ifdef HAVE_COARRAY
    if (allocated(this%send_buf)) deallocate(this%send_buf)
    if (allocated(this%send_len)) deallocate(this%send_len)
    if (allocated(this%recv_len)) deallocate(this%recv_len)
    if (allocated(this%send_offset)) deallocate(this%send_offset)
    if (allocated(this%recv_offset)) deallocate(this%recv_offset)
    if (allocated(this%dest_offset)) deallocate(this%dest_offset)
    if (allocated(this%send_img)) deallocate(this%send_img)
    if (allocated(this%recv_img)) deallocate(this%recv_img)
    if (allocated(this%sync_img)) deallocate(this%sync_img)

    call this%free_order()
    call this%free_dofs()
#endif
  end subroutine gs_caf_free

  !> Pack u into per-peer slabs and put each slab into the remote image's
  !! recv_buf. Double buffering means each round writes to a different
  !! half of the recv coarray, so no back-pressure synchronisation is
  !! needed in sync mode -- the visibility synchronisation in nbwait
  !! suffices. Atomic and event modes still use their per-pair signalling.
  subroutine gs_nbsend_caf(this, u, n, tag, deps, strm)
    class(gs_caf_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
#ifdef HAVE_COARRAY
    integer :: i, j, dst, off, dimg, ndst, doff, half_off
    integer, pointer :: sp(:)
    integer(kind=atomic_int_kind) :: flag
    integer :: me_rank

    ! These routines are entered from inside the gs_op_vector OpenMP
    ! parallel region (the CAF backend is a host comm method), so the whole
    ! thread team executes them. Coarray Fortran image control is not
    ! thread-safe: every put, sync, atomic_* and event statement -- and the
    ! shared module/instance signalling state -- is performed by the master
    ! thread alone, inside !$omp master regions. Only the local pack of the
    ! send buffer is work-shared across the team with !$omp do. The implicit
    ! barrier at the end of each !$omp do guarantees the slab is fully packed
    ! before the master thread reads it for the put.
    !
    ! parity is flipped by the master at the end of nbwait and published by
    ! the barrier there, so it is stable for every thread here.
    half_off = this%parity * gs_caf_buf_size

    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          sp => this%send_dof(dst)%array()
          !OCL NORECURRENCE, NOVREC, NOALIAS
          !DIR$ CONCURRENT
          !DIR$ IVDEP
          !GCC$ ivdep
          !NEC$ IVDEP
          !$omp do
          do j = 1, ndst
             this%send_buf(off + j) = u(sp(j))
          end do
          !$omp end do
          ! No barrier needed between the pack and the put: end do carries
          ! one. The team meanwhile packs the next peer's (disjoint) slab
          ! while the master injects this one.
          !$omp master
          gs_caf_recv_buf(half_off + doff + 1 : half_off + doff + ndst)[dimg] &
               = this%send_buf(off + 1 : off + ndst)
          !$omp end master
       end do
       !$omp barrier
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       ! Event-mode signalling is a sequence of image-control statements and
       ! is executed by the master thread alone.
       !$omp master
       ! Event mode shares one set of module-level event coarrays among
       ! all instances and cannot disambiguate posts from concurrent gs
       ! ops, so we must guarantee non-overlapping nbsend/nbwait windows.
       if (gs_caf_event_in_use) then
          call neko_error("Event-mode coarray gather-scatter does not " // &
               "support overlapping gs ops on different instances")
       end if
       gs_caf_event_in_use = .true.

       ! Wait for all receivers to have credited their buffers (skipped
       ! on the first nbsend; there are no credits posted yet).
       if (this%send_started) then
          if (size(this%send_pe) .gt. 0) then
             event wait(gs_caf_buf_ready_ev, until_count=size(this%send_pe))
          end if
       else
          this%send_started = .true.
       end if
       !$omp end master
       ! The back-pressure wait above must complete before any put below
       ! overwrites the receivers' buffers.
       !$omp barrier

       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          sp => this%send_dof(dst)%array()
          !OCL NORECURRENCE, NOVREC, NOALIAS
          !DIR$ CONCURRENT
          !DIR$ IVDEP
          !GCC$ ivdep
          !NEC$ IVDEP
          !$omp do
          do j = 1, ndst
             this%send_buf(off + j) = u(sp(j))
          end do
          !$omp end do
          !$omp master
          gs_caf_recv_buf(half_off + doff + 1 : half_off + doff + ndst)[dimg] &
               = this%send_buf(off + 1 : off + ndst)
          ! event post is meant to act as an image-control statement
          ! that establishes segment ordering with the matching event
          ! wait, but real-world coarray runtimes can let a small event
          ! message race past a still-in-flight RDMA put -- the
          ! receiver's wait then completes before the data has landed.
          ! sync memory forces the put to commit locally before the post.
          sync memory
          event post(gs_caf_data_ready_ev[dimg])
          !$omp end master
       end do
       !$omp barrier
#endif
    else
       ! Pack all peers up front (work-shared) so the subsequent network
       ! waits and puts can overlap with each other rather than serialising
       ! behind per-peer pack work.
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
          !$omp do
          do j = 1, ndst
             this%send_buf(off + j) = u(sp(j))
          end do
          !$omp end do
       end do

       ! Back-pressure, put and signal per peer -- all coarray operations,
       ! so master only. With double-buffering the half we are about to
       ! write last carried round (send_count - 2), so we only need the
       ! receiver to have unpacked through (send_count - 1).
       !$omp master
       me_rank = this_image() - 1
       do i = 1, size(this%send_pe)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)

          do
             call atomic_ref(flag, gs_caf_buf_ready(this%send_pe(i)))
             if (int(flag) .ge. gs_caf_send_count(this%send_pe(i)) - 1) exit
          end do

          gs_caf_recv_buf(half_off + doff + 1 : half_off + doff + ndst)[dimg] &
               = this%send_buf(off + 1 : off + ndst)

          gs_caf_send_count(this%send_pe(i)) = &
               gs_caf_send_count(this%send_pe(i)) + 1
          call atomic_define(gs_caf_data_ready(me_rank)[dimg], &
               int(gs_caf_send_count(this%send_pe(i)), atomic_int_kind))
       end do
       !$omp end master
       !$omp barrier
    end if
#else
    call neko_error("Coarray Fortran support not built")
#endif
  end subroutine gs_nbsend_caf

  !> No-op for coarrays: senders push into the receiver's buffer, so
  !! the receive side does not need to post anything.
  subroutine gs_nbrecv_caf(this, tag)
    class(gs_caf_t), intent(inout) :: this
    integer, intent(in) :: tag
  end subroutine gs_nbrecv_caf

  !> Wait for all incoming puts and reduce them into u. In sync mode a
  !! sync_images bracket pairs with the senders' nbsend; in atomic mode
  !! each sender is awaited via its data_ready counter and credited via
  !! buf_ready after unpack.
  subroutine gs_nbwait_caf(this, u, n, op, strm)
    class(gs_caf_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
#ifdef HAVE_COARRAY
    integer :: i, j, src, off, nsrc, half_off
    integer, pointer :: sp(:)
    integer(kind=atomic_int_kind) :: flag
    integer :: me_rank

    half_off = this%parity * gs_caf_buf_size

    ! All incoming data is awaited by the master (see gs_nbsend_caf on why
    ! coarray operations are funnelled); the barrier then releases the team
    ! into the unpack once every slab has landed.
    !$omp master
    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       if (allocated(this%sync_img)) then
          if (size(this%sync_img) .gt. 0) then
             sync images(this%sync_img)
          end if
       end if
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       if (size(this%recv_pe) .gt. 0) then
          event wait(gs_caf_data_ready_ev, until_count=size(this%recv_pe))
       end if
#endif
    else
       ! Atomic mode: spin per-sender on data_ready until the expected
       ! round count is observed.
       do i = 1, size(this%recv_pe)
          gs_caf_recv_count(this%recv_pe(i)) = &
               gs_caf_recv_count(this%recv_pe(i)) + 1
          do
             call atomic_ref(flag, gs_caf_data_ready(this%recv_pe(i)))
             if (int(flag) .ge. gs_caf_recv_count(this%recv_pe(i))) exit
          end do
       end do
    end if
    !$omp end master
    !$omp barrier

    ! Reduce each received slab into u. The loop over peers stays serial: a
    ! dof shared by 3+ ranks appears in several recv_dof lists, so reducing
    ! two slabs concurrently would race on that dof. Parallelism is taken
    ! within each slab instead.
    do i = 1, size(this%recv_pe)
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
             u(sp(j)) = u(sp(j)) + gs_caf_recv_buf(half_off + off + j)
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
             u(sp(j)) = u(sp(j)) * gs_caf_recv_buf(half_off + off + j)
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
             u(sp(j)) = min(u(sp(j)), gs_caf_recv_buf(half_off + off + j))
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
             u(sp(j)) = max(u(sp(j)), gs_caf_recv_buf(half_off + off + j))
          end do
          !$omp end do
       case default
          call neko_error("Unknown operation in gs_nbwait_caf")
       end select
    end do

    ! The implicit barrier at the last end do guarantees every slab has been
    ! consumed before the master credits the senders and flips the parity;
    ! the trailing barrier publishes the new parity to the whole team.
    !$omp master
    if (gs_caf_mode .eq. GS_CAF_SIGNAL_ATOMIC) then
       ! Credit each sender that we have unpacked their slab so they
       ! may proceed with their next round.
       me_rank = this_image() - 1
       do i = 1, size(this%recv_pe)
          call atomic_define(gs_caf_buf_ready(me_rank)[this%recv_img(i)], &
               int(gs_caf_recv_count(this%recv_pe(i)), atomic_int_kind))
       end do
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       do i = 1, size(this%recv_pe)
          event post(gs_caf_buf_ready_ev[this%recv_img(i)])
       end do
       gs_caf_event_in_use = .false.
#endif
    end if

    ! Flip the double-buffer parity for the next round.
    this%parity = 1 - this%parity
    !$omp end master
    !$omp barrier
#else
    call neko_error("Coarray Fortran support not built")
#endif
  end subroutine gs_nbwait_caf

  !> Fused nc-component put. Each peer slab is nc consecutive component
  !! blocks; the remote placement offset and slab length scale by nc, the
  !! per-peer signalling (sync/atomic/event) is unchanged.
  !! @param u compact shared buffer, component-outer: u((c-1)*n + idx).
  subroutine gs_nbsend_vec_caf(this, u, n, nc, tag, deps, strm)
    class(gs_caf_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
#ifdef HAVE_COARRAY
    integer :: i, j, c, dst, off, dimg, ndst, doff, half_off
    integer, pointer :: sp(:)
    integer(kind=atomic_int_kind) :: flag
    integer :: me_rank

    ! Same threading split as the scalar path: the nc-component pack is
    ! work-shared over the dofs of one peer at a time, every coarray
    ! operation is issued by the master alone.
    half_off = this%parity * gs_caf_buf_size

    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          sp => this%send_dof(dst)%array()
          !$omp do
          do j = 1, ndst
             do c = 1, nc
                this%send_buf(nc*off + (c-1)*ndst + j) = u((c-1)*n + sp(j))
             end do
          end do
          !$omp end do
          !$omp master
          gs_caf_recv_buf(half_off + nc*doff + 1 : half_off + nc*doff + nc*ndst) &
               [dimg] = this%send_buf(nc*off + 1 : nc*off + nc*ndst)
          !$omp end master
       end do
       !$omp barrier
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       !$omp master
       if (gs_caf_event_in_use) then
          call neko_error("Event-mode coarray gather-scatter does not " // &
               "support overlapping gs ops on different instances")
       end if
       gs_caf_event_in_use = .true.

       if (this%send_started) then
          if (size(this%send_pe) .gt. 0) then
             event wait(gs_caf_buf_ready_ev, until_count=size(this%send_pe))
          end if
       else
          this%send_started = .true.
       end if
       !$omp end master
       !$omp barrier

       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          sp => this%send_dof(dst)%array()
          !$omp do
          do j = 1, ndst
             do c = 1, nc
                this%send_buf(nc*off + (c-1)*ndst + j) = u((c-1)*n + sp(j))
             end do
          end do
          !$omp end do
          !$omp master
          gs_caf_recv_buf(half_off + nc*doff + 1 : half_off + nc*doff + nc*ndst) &
               [dimg] = this%send_buf(nc*off + 1 : nc*off + nc*ndst)
          sync memory
          event post(gs_caf_data_ready_ev[dimg])
          !$omp end master
       end do
       !$omp barrier
#endif
    else
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          sp => this%send_dof(dst)%array()
          !$omp do
          do j = 1, ndst
             do c = 1, nc
                this%send_buf(nc*off + (c-1)*ndst + j) = u((c-1)*n + sp(j))
             end do
          end do
          !$omp end do
       end do

       !$omp master
       me_rank = this_image() - 1
       do i = 1, size(this%send_pe)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)

          do
             call atomic_ref(flag, gs_caf_buf_ready(this%send_pe(i)))
             if (int(flag) .ge. gs_caf_send_count(this%send_pe(i)) - 1) exit
          end do

          gs_caf_recv_buf(half_off + nc*doff + 1 : half_off + nc*doff + nc*ndst) &
               [dimg] = this%send_buf(nc*off + 1 : nc*off + nc*ndst)

          gs_caf_send_count(this%send_pe(i)) = &
               gs_caf_send_count(this%send_pe(i)) + 1
          call atomic_define(gs_caf_data_ready(me_rank)[dimg], &
               int(gs_caf_send_count(this%send_pe(i)), atomic_int_kind))
       end do
       !$omp end master
       !$omp barrier
    end if
#else
    call neko_error("Coarray Fortran support not built")
#endif
  end subroutine gs_nbsend_vec_caf

  !> No-op: senders push into the receiver's buffer.
  subroutine gs_nbrecv_vec_caf(this, tag, nc)
    class(gs_caf_t), intent(inout) :: this
    integer, intent(in) :: tag, nc
  end subroutine gs_nbrecv_vec_caf

  !> Fused nc-component wait/reduce for the coarray backend.
  subroutine gs_nbwait_vec_caf(this, u, n, nc, op, strm)
    class(gs_caf_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
#ifdef HAVE_COARRAY
    integer :: i, j, c, src, off, nsrc, half_off
    integer, pointer :: sp(:)
    integer(kind=atomic_int_kind) :: flag
    integer :: me_rank

    half_off = this%parity * gs_caf_buf_size

    !$omp master
    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       if (allocated(this%sync_img)) then
          if (size(this%sync_img) .gt. 0) then
             sync images(this%sync_img)
          end if
       end if
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       if (size(this%recv_pe) .gt. 0) then
          event wait(gs_caf_data_ready_ev, until_count=size(this%recv_pe))
       end if
#endif
    else
       do i = 1, size(this%recv_pe)
          gs_caf_recv_count(this%recv_pe(i)) = &
               gs_caf_recv_count(this%recv_pe(i)) + 1
          do
             call atomic_ref(flag, gs_caf_data_ready(this%recv_pe(i)))
             if (int(flag) .ge. gs_caf_recv_count(this%recv_pe(i))) exit
          end do
       end do
    end if
    !$omp end master
    !$omp barrier

    ! Serial over peers (a dof shared by 3+ ranks appears in several recv
    ! lists); parallelism is taken within each slab.
    do i = 1, size(this%recv_pe)
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
                     gs_caf_recv_buf(half_off + nc*off + (c-1)*nsrc + j)
             end do
          end do
          !$omp end do
       case (GS_OP_MUL)
          !$omp do
          do j = 1, nsrc
             do c = 1, nc
                u((c-1)*n + sp(j)) = u((c-1)*n + sp(j)) * &
                     gs_caf_recv_buf(half_off + nc*off + (c-1)*nsrc + j)
             end do
          end do
          !$omp end do
       case (GS_OP_MIN)
          !$omp do
          do j = 1, nsrc
             do c = 1, nc
                u((c-1)*n + sp(j)) = min(u((c-1)*n + sp(j)), &
                     gs_caf_recv_buf(half_off + nc*off + (c-1)*nsrc + j))
             end do
          end do
          !$omp end do
       case (GS_OP_MAX)
          !$omp do
          do j = 1, nsrc
             do c = 1, nc
                u((c-1)*n + sp(j)) = max(u((c-1)*n + sp(j)), &
                     gs_caf_recv_buf(half_off + nc*off + (c-1)*nsrc + j))
             end do
          end do
          !$omp end do
       case default
          call neko_error("Unknown operation in gs_nbwait_vec_caf")
       end select
    end do

    !$omp master
    if (gs_caf_mode .eq. GS_CAF_SIGNAL_ATOMIC) then
       me_rank = this_image() - 1
       do i = 1, size(this%recv_pe)
          call atomic_define(gs_caf_buf_ready(me_rank)[this%recv_img(i)], &
               int(gs_caf_recv_count(this%recv_pe(i)), atomic_int_kind))
       end do
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       do i = 1, size(this%recv_pe)
          event post(gs_caf_buf_ready_ev[this%recv_img(i)])
       end do
       gs_caf_event_in_use = .false.
#endif
    end if

    this%parity = 1 - this%parity
    !$omp end master
    !$omp barrier
#else
    call neko_error("Coarray Fortran support not built")
#endif
  end subroutine gs_nbwait_vec_caf

end module gs_caf
