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
  !$ use omp_lib
  implicit none
  private

  !> Whether coarray support was built into this Neko. Lets callers (e.g.
  !! the gs comm. autotuner) skip the backend rather than aborting in init
  !! on builds without it.
#ifdef HAVE_COARRAY
  logical, parameter, public :: GS_CAF_AVAIL = .true.
#else
  logical, parameter, public :: GS_CAF_AVAIL = .false.
#endif

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

  ! Double-buffer parity (0 or 1) per remote rank, flipped at the end of
  ! every nbwait for the peers that round involved: the round writes to and
  ! reads from the gs_caf_peer_parity(peer)*gs_caf_buf_size half, the next
  ! round with that same peer uses the other one.
  !
  ! The parity belongs to the buffer, which is module state shared by every
  ! instance, so it cannot live in gs_caf_t. It has to be per *pair* rather
  ! than one counter for the whole buffer, because what makes it safe for a
  ! sender to run ahead of a receiver is "the next round I exchange with this
  ! peer writes the other half", and only traffic between that pair advances
  ! it. A single global parity holds only while every instance shares one
  ! peer set: with instances whose peer sets differ -- hsmg, where the
  ! Schwarz smoother's gs lives on the extended lx+2 dofmap and so does not
  ! have the connectivity of the level's own gs -- two consecutive rounds
  ! between I and J can sit an even number of global rounds apart and land on
  ! the same half. Nothing then orders I's put against J's still-running
  ! unpack: in sync mode the intervening round does not name the pair in its
  ! sync images set, and in atomic mode the back-pressure spin explicitly
  ! tolerates a receiver one round behind. p-multigrid never shows it, since
  ! every level shares the partition's peer set and so every pair
  ! synchronises every round.
  !
  ! Sender and receiver stay in step because gs schedules are symmetric: both
  ! endpoints of a pair take part in the same rounds and flip at the same
  ! point, at the end of nbwait, so a round's put and its unpack agree on the
  ! half. Flipping in nbsend as well would make the two ends disagree.
  integer, allocatable :: gs_caf_peer_parity(:)

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
     !> Ranks of the union of send and recv peers, i.e. every peer this
     !! instance exchanges with in a round. Used to advance the per-peer
     !! double-buffer parity in nbwait.
     integer, allocatable :: peer(:)
     !> sync-mode only: the same set as image numbers, used for the pairwise
     !! sync images that bracket the put. Unallocated in atomic and event
     !! modes.
     integer, allocatable :: sync_img(:)
     !> event-mode only: false on the very first nbsend so the buf_ready
     !! wait is skipped (there are no credits posted yet).
     logical :: send_started = .false.
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

  public :: gs_caf_usable, gs_caf_signal_auto, gs_caf_signal_modes, &
       gs_caf_mode_name, gs_caf_mode_get, gs_caf_set_mode

contains

  !> Whether the coarray backend can actually run in this job. GS_CAF_AVAIL
  !! only says that the compiler accepted coarrays when Neko was configured;
  !! a compiler may well do so and still provide a single image per process
  !! (gfortran's -fcoarray=single, a coarray runtime that was not linked in),
  !! in which case every image sees num_images() == 1 and no halo can be
  !! exchanged. The images must map one-to-one onto the ranks of NEKO_COMM,
  !! as the backend addresses its peers by image number (rank + 1) -- which
  !! also rules out a run split into several communicators (NEKO_COMM_ID).
  !! @return whether a gs_caf_t can be used by this run
  function gs_caf_usable() result(usable)
    logical :: usable

#ifdef HAVE_COARRAY
    usable = (num_images() .eq. pe_size)
#else
    usable = .false.
#endif

  end function gs_caf_usable

  !> Whether the signaling mode should be selected by benchmarking, i.e.
  !! NEKO_GS_CAF_SIGNALING=auto. The mode is a program-wide binding shared
  !! by every gs_caf_t instance, so the caller (the gs comm. autotuner) must
  !! bind it once and keep it: see gs_caf_set_mode.
  !! @return whether the mode should be selected by benchmarking
  function gs_caf_signal_auto() result(auto)
    logical :: auto
    character(len=64) :: env_val
    integer :: env_len

    call get_environment_variable("NEKO_GS_CAF_SIGNALING", env_val, env_len)
    auto = (env_len .gt. 0)
    if (auto) auto = (env_val(1:env_len) .eq. "auto")

  end function gs_caf_signal_auto

  !> The signaling modes this build can run, in the order they should be
  !! benchmarked. Events are only available with a compiler that implements
  !! them.
  !! @return the runnable modes, as GS_CAF_SIGNAL_* constants
  function gs_caf_signal_modes() result(mode)
    integer, allocatable :: mode(:)
    integer :: m(3), n

    n = 2
    m(1) = GS_CAF_SIGNAL_SYNC
    m(2) = GS_CAF_SIGNAL_ATOMIC
#ifdef HAVE_COARRAY_EVENTS
    n = n + 1
    m(n) = GS_CAF_SIGNAL_EVENT
#endif

    allocate(mode(n))
    mode = m(1:n)

  end function gs_caf_signal_modes

  !> The signaling mode currently in force, or 0 if none has been bound yet
  !! (no gs_caf_t has been initialised and the autotuner has not run).
  !! @return the bound mode, a GS_CAF_SIGNAL_* constant, or 0
  function gs_caf_mode_get() result(mode)
    integer :: mode

#ifdef HAVE_COARRAY
    mode = gs_caf_mode
#else
    mode = 0
#endif

  end function gs_caf_mode_get

  !> Name of the signaling mode @a mode, right-adjusted for the log
  !! @param mode the mode to name, one of the GS_CAF_SIGNAL_* constants
  function gs_caf_mode_name(mode) result(name)
    integer, intent(in) :: mode
    character(len=12) :: name

    select case (mode)
    case (GS_CAF_SIGNAL_SYNC)
       name = '        sync'
    case (GS_CAF_SIGNAL_ATOMIC)
       name = '      atomic'
    case (GS_CAF_SIGNAL_EVENT)
       name = '       event'
    case default
       name = '     unknown'
       call neko_error('Unknown coarray gather-scatter signaling mode')
    end select

  end function gs_caf_mode_name

  !> Bind the signaling mode shared by every gs_caf_t instance, allocating
  !! whatever module-level state the mode needs. Idempotent per mode, so it
  !! may be called again to switch modes while no gs op is in flight and no
  !! gs_caf_t instance is live -- which is what the autotuner does to
  !! benchmark the modes against each other.
  !! @param mode the mode to bind, one of the GS_CAF_SIGNAL_* constants
  !! @note Collective: the atomic counters and the event coarrays are
  !! coarrays, so every image must call this in the same order.
  subroutine gs_caf_set_mode(mode)
    integer, intent(in) :: mode
#ifdef HAVE_COARRAY
    integer :: i

    select case (mode)
    case (GS_CAF_SIGNAL_SYNC)
       ! Nothing to set up; the image set is built per instance in init
    case (GS_CAF_SIGNAL_ATOMIC)
       if (.not. allocated(gs_caf_data_ready)) then
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
       end if
    case (GS_CAF_SIGNAL_EVENT)
#ifdef HAVE_COARRAY_EVENTS
       ! One set of event coarrays shared by all instances. The guard
       ! against overlapping gs ops is enforced in nbsend/nbwait, not here.
       if (.not. allocated(gs_caf_data_ready_ev)) then
          allocate(gs_caf_data_ready_ev[*])
          allocate(gs_caf_buf_ready_ev[*])
       end if
#else
       call neko_error("NEKO_GS_CAF_SIGNALING=event requires a Fortran " // &
            "compiler with coarray events support")
#endif
    case default
       call neko_error('Unknown coarray gather-scatter signaling mode')
    end select

    gs_caf_mode = mode
#else
    call neko_error("Coarray Fortran support not built; reconfigure with " // &
         "a coarray-capable Fortran compiler")
#endif
  end subroutine gs_caf_set_mode

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

    ! A build with coarray support is not necessarily a build that can run
    ! them across the job (see gs_caf_usable). Catch it here rather than
    ! silently exchanging nothing, or deadlocking in the first put.
    if (.not. gs_caf_usable()) then
       call neko_error("Coarray gather-scatter needs one image per rank; " // &
            "this build or run has num_images() /= pe_size")
    end if

    ! Bind the signaling mode on the first init. With
    ! NEKO_GS_CAF_SIGNALING=auto the mode is instead selected by
    ! benchmarking, and the gs comm. autotuner has already bound it with
    ! gs_caf_set_mode by the time we get here; falling through to sync
    ! covers the case where nothing tuned it (CAF requested explicitly).
    if (gs_caf_mode .eq. 0) then
       call get_environment_variable("NEKO_GS_CAF_SIGNALING", env_val, env_len)
       if (env_len .gt. 0 .and. env_val(1:env_len) .eq. "atomic") then
          call gs_caf_set_mode(GS_CAF_SIGNAL_ATOMIC)
       else if (env_len .gt. 0 .and. env_val(1:env_len) .eq. "event") then
          call gs_caf_set_mode(GS_CAF_SIGNAL_EVENT)
       else
          call gs_caf_set_mode(GS_CAF_SIGNAL_SYNC)
       end if
    end if

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

    ! Peer set = union of send and recv peers. Both endpoints of every
    ! neighbour pair include each other, so the two sides advance the pair's
    ! double-buffer parity on the same rounds (and, in sync mode, their
    ! pairwise sync images statements match up).
    allocate(in_neigh(0:pe_size - 1))
    in_neigh = .false.
    do i = 1, nsend
       in_neigh(this%send_pe(i)) = .true.
    end do
    do i = 1, nrecv
       in_neigh(this%recv_pe(i)) = .true.
    end do
    n_neigh = count(in_neigh)
    allocate(this%peer(n_neigh))
    n_neigh = 0
    do i = 0, pe_size - 1
       if (in_neigh(i)) then
          n_neigh = n_neigh + 1
          this%peer(n_neigh) = i
       end if
    end do
    deallocate(in_neigh)

    ! Shared by every instance, like the buffer it indexes, so allocate it
    ! once and never reset it: the parity of a pair has to keep advancing
    ! across instances.
    if (.not. allocated(gs_caf_peer_parity)) then
       allocate(gs_caf_peer_parity(0:pe_size - 1))
       gs_caf_peer_parity = 0
    end if

    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       allocate(this%sync_img(size(this%peer)))
       this%sync_img = this%peer + 1
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       ! Start from zero event counts; see gs_caf_event_drain for why an
       ! instance leaves credits behind and what accumulating them costs.
       call gs_caf_event_drain()
       ! No gs op can be in flight here (init is collective), so any in-use
       ! flag left set by a torn-down instance is stale.
       gs_caf_event_in_use = .false.
#endif
    end if ! atomic mode: no per-instance state to allocate

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
    if (allocated(this%peer)) deallocate(this%peer)
    if (allocated(this%sync_img)) deallocate(this%sync_img)

    call this%free_order()
    call this%free_dofs()
#endif
  end subroutine gs_caf_free

#ifdef HAVE_COARRAY

#ifdef HAVE_COARRAY_EVENTS

  !> Consume the event posts left over by a previous gs_caf_t, so that each
  !! instance starts counting from zero.
  !!
  !! The event coarrays are module state that outlives an instance, while
  !! the credit accounting is per-instance. Over a run of R rounds an image
  !! receives one buf_ready post per send peer per round but consumes them
  !! only R-1 times, since the first nbsend skips its wait (nothing has been
  !! credited yet) -- so every instance leaves size(send_pe) credits behind.
  !! Left alone they accumulate by one round per instance until a later wait
  !! is satisfied by stale credits: the sender then runs further ahead than
  !! the double buffer tolerates and can overwrite a half the receiver has
  !! not unpacked. data_ready is balanced (every nbsend is paired with an
  !! nbwait) but is drained too, so the invariant is simply "both events are
  !! zero when an instance starts".
  !!
  !! event_query reports what has already landed, so neither wait can block.
  !! The caller brackets this with the sync all of gs_caf_init: the earlier
  !! one has delivered every post of the previous instance, the later one
  !! keeps any image from signalling before all of them have drained.
  subroutine gs_caf_event_drain()
    integer :: pending

    call event_query(gs_caf_buf_ready_ev, pending)
    if (pending .gt. 0) then
       event wait(gs_caf_buf_ready_ev, until_count=pending)
    end if

    call event_query(gs_caf_data_ready_ev, pending)
    if (pending .gt. 0) then
       event wait(gs_caf_data_ready_ev, until_count=pending)
    end if

  end subroutine gs_caf_event_drain

#endif
#endif

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
    integer :: me_rank, thrdid

    thrdid = 0
    !$ thrdid = omp_get_thread_num()

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

    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          half_off = gs_caf_peer_parity(dst) * gs_caf_buf_size
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
          if (thrdid .eq. 0) then
             gs_caf_recv_buf(half_off + doff + 1 : half_off + doff + ndst)[dimg] &
                  = this%send_buf(off + 1 : off + ndst)
          end if
       end do
       !$omp barrier
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       ! Event-mode signalling is a sequence of image-control statements and
       ! is executed by the master thread alone.
       if (thrdid .eq. 0) then
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
       end if
       ! The back-pressure wait above must complete before any put below
       ! overwrites the receivers' buffers.
       !$omp barrier

       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          half_off = gs_caf_peer_parity(dst) * gs_caf_buf_size
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
          if (thrdid .eq. 0) then
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
          end if
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
       if (thrdid .eq. 0) then
          me_rank = this_image() - 1
          do i = 1, size(this%send_pe)
             dst = this%send_pe(i)
             off = this%send_offset(i)
             ndst = this%send_len(i)
             dimg = this%send_img(i)
             doff = this%dest_offset(i)
             half_off = gs_caf_peer_parity(dst) * gs_caf_buf_size

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
       end if
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
    integer :: me_rank, thrdid

    thrdid = 0
    !$ thrdid = omp_get_thread_num()

    ! All incoming data is awaited by the master (see gs_nbsend_caf on why
    ! coarray operations are funnelled); the barrier then releases the team
    ! into the unpack once every slab has landed.
    if (thrdid .eq. 0) then
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
    end if
    !$omp barrier

    ! Reduce each received slab into u. The loop over peers stays serial: a
    ! dof shared by 3+ ranks appears in several recv_dof lists, so reducing
    ! two slabs concurrently would race on that dof. Parallelism is taken
    ! within each slab instead.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       off = this%recv_offset(i)
       nsrc = this%recv_len(i)
       half_off = gs_caf_peer_parity(src) * gs_caf_buf_size
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
    if (thrdid .eq. 0) then
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

       ! Flip the double-buffer parity of every peer this round involved, so
       ! our next exchange with each of them uses the other half.
       do i = 1, size(this%peer)
          gs_caf_peer_parity(this%peer(i)) = &
               1 - gs_caf_peer_parity(this%peer(i))
       end do
    end if
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
    integer :: me_rank, thrdid

    thrdid = 0
    !$ thrdid = omp_get_thread_num()

    ! Same threading split as the scalar path: the nc-component pack is
    ! work-shared over the dofs of one peer at a time, every coarray
    ! operation is issued by the master alone.

    if (gs_caf_mode .eq. GS_CAF_SIGNAL_SYNC) then
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          half_off = gs_caf_peer_parity(dst) * gs_caf_buf_size
          sp => this%send_dof(dst)%array()
          !$omp do
          do j = 1, ndst
             do c = 1, nc
                this%send_buf(nc*off + (c-1)*ndst + j) = u((c-1)*n + sp(j))
             end do
          end do
          !$omp end do
          if (thrdid .eq. 0) then
             gs_caf_recv_buf(half_off + nc*doff + 1 : half_off + nc*doff + nc*ndst) &
                  [dimg] = this%send_buf(nc*off + 1 : nc*off + nc*ndst)
          end if
       end do
       !$omp barrier
#ifdef HAVE_COARRAY_EVENTS
    else if (gs_caf_mode .eq. GS_CAF_SIGNAL_EVENT) then
       if (thrdid .eq. 0) then
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
       end if
       !$omp barrier

       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          dimg = this%send_img(i)
          doff = this%dest_offset(i)
          half_off = gs_caf_peer_parity(dst) * gs_caf_buf_size
          sp => this%send_dof(dst)%array()
          !$omp do
          do j = 1, ndst
             do c = 1, nc
                this%send_buf(nc*off + (c-1)*ndst + j) = u((c-1)*n + sp(j))
             end do
          end do
          !$omp end do
          if (thrdid .eq. 0) then
             gs_caf_recv_buf(half_off + nc*doff + 1 : half_off + nc*doff + nc*ndst) &
                  [dimg] = this%send_buf(nc*off + 1 : nc*off + nc*ndst)
             ! event post is meant to act as an image-control statement
             ! that establishes segment ordering with the matching event
             ! wait, but real-world coarray runtimes can let a small event
             ! message race past a still-in-flight RDMA put -- the
             ! receiver's wait then completes before the data has landed.
             ! sync memory forces the put to commit locally before the post.
             sync memory
             event post(gs_caf_data_ready_ev[dimg])
          end if
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

       if (thrdid .eq. 0) then
          me_rank = this_image() - 1
          do i = 1, size(this%send_pe)
             dst = this%send_pe(i)
             off = this%send_offset(i)
             ndst = this%send_len(i)
             dimg = this%send_img(i)
             doff = this%dest_offset(i)
             half_off = gs_caf_peer_parity(dst) * gs_caf_buf_size

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
       end if
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
    integer :: me_rank, thrdid

    thrdid = 0
    !$ thrdid = omp_get_thread_num()
    if (thrdid .eq. 0) then
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
    end if
    !$omp barrier

    ! Serial over peers (a dof shared by 3+ ranks appears in several recv
    ! lists); parallelism is taken within each slab.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       off = this%recv_offset(i)
       nsrc = this%recv_len(i)
       half_off = gs_caf_peer_parity(src) * gs_caf_buf_size
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

    if (thrdid .eq. 0) then
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

       do i = 1, size(this%peer)
          gs_caf_peer_parity(this%peer(i)) = 1 - gs_caf_peer_parity(this%peer(i))
       end do
    end if
    !$omp barrier
#else
    call neko_error("Coarray Fortran support not built")
#endif
  end subroutine gs_nbwait_vec_caf

end module gs_caf
