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
!> Defines OpenSHMEM gather-scatter communication
module gs_shmem
  use num_types, only : rp, c_rp, i8
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use comm, only : NEKO_COMM, pe_rank, pe_size
  use mpi_f08, only : MPI_Allreduce, MPI_Alltoall, MPI_INTEGER, MPI_MAX
  use utils, only : neko_error
#ifdef HAVE_OPENSHMEM
  use shmem, only : shmem_malloc, shmem_calloc, shmem_free, &
       shmem_putmem_signal_nbi, shmem_signal_wait_until, &
       shmem_uint64_atomic_set, shmem_uint64_wait_until, &
       shmem_barrier_all, shmem_query_thread, SHMEM_SIGNAL_SET, &
       SHMEM_CMP_GE, SHMEM_THREAD_MULTIPLE
#endif
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_loc, &
       c_f_pointer, c_associated, c_sizeof, c_size_t, c_int64_t, c_int
  implicit none
  private

  !> True when the OpenSHMEM library provides SHMEM_THREAD_MULTIPLE. The
  !! send path then lets each OpenMP thread run the back-pressure wait, the
  !! pack and the put for its own share of the send peers; otherwise every
  !! SHMEM call is funnelled through the master thread, mirroring how gs_mpi
  !! keys off NEKO_MPI_THREAD_PROVIDED. Queried on every gs_shmem_t init
  !! (the level is a property of the library, so all instances agree).
  logical :: gs_shmem_thread_multiple = .false.

  !> Whether a native OpenSHMEM library was built into this Neko
  !! (--with-openshmem). Lets callers (e.g. the gs comm. autotuner) skip
  !! the backend rather than aborting in init on builds without it.
#ifdef HAVE_OPENSHMEM
  logical, parameter, public :: GS_SHMEM_AVAIL = .true.
#else
  logical, parameter, public :: GS_SHMEM_AVAIL = .false.
#endif

  !> Symmetric buffer for one direction of OpenSHMEM communication
  type :: gs_shmem_buf_t
     !> Number of dofs per neighbor
     integer, allocatable :: ndofs(:)
     !> Local offset in this buffer for each neighbor
     integer, allocatable :: offset(:)
     !> For send_buf: offset in the remote PE's recv buffer where our
     !! data should land. Unused for recv_buf.
     integer, allocatable :: remote_offset(:)
     !> Total number of dofs in this buffer on this PE
     integer :: total = 0
     !> Maximum total across all PEs (size of the symmetric allocation)
     integer :: max_total = 0
     !> Raw pointer to the symmetric buffer returned by shmem_malloc.
     type(c_ptr) :: buf_ptr = C_NULL_PTR
   contains
     procedure, pass(this) :: init => gs_shmem_buf_init
     procedure, pass(this) :: free => gs_shmem_buf_free
  end type gs_shmem_buf_t

  !> Gather-scatter communication using OpenSHMEM one-sided puts with
  !! per-rank signaling for completion (OpenSHMEM 1.5).
  !!
  !! Each PE allocates symmetric send and recv buffers sized to the
  !! global maximum number of dofs (so puts can target a well-defined
  !! offset on the remote PE's recv buffer) plus two symmetric uint64
  !! arrays of size pe_size, indexed by the remote PE's rank:
  !! data_signals[r] holds an "iter" value written by PE r when it has
  !! put data into our recv buffer (via shmem_putmem_signal_nbi);
  !! ack_signals[r] holds an "iter" value written by PE r after it has
  !! consumed our most recent put (via shmem_uint64_atomic_set). The
  !! sender waits on its own ack_signals[dst] before each put so a fast
  !! sender cannot overwrite a buffer the receiver hasn't consumed.
  !!
  !! Per-rank slot indexing avoids the slot-exchange handshake that a
  !! neighbor-list-position scheme would require, which is fragile
  !! because send_pe and recv_pe are pushed independently in
  !! gs_schedule and may differ in size and ordering.
  !!
  !! Under OpenMP the pack and unpack loops are work-shared across the
  !! calling team. The send path additionally farms the whole per-peer
  !! protocol out to the threads when the library provides
  !! SHMEM_THREAD_MULTIPLE (see gs_shmem_thread_multiple); otherwise, and
  !! always on the receive side, the SHMEM calls are issued by the master
  !! thread alone.
  type, public, extends(gs_comm_t) :: gs_shmem_t
     type(gs_shmem_buf_t) :: send_buf
     type(gs_shmem_buf_t) :: recv_buf
     ! Symmetric data-arrival signals; data_signals[r] is set to `iter`
     ! by PE r via shmem_putmem_signal_nbi when it has put data into
     ! our recv buffer. Sized to pe_size on every PE.
     type(c_ptr) :: data_signals_ptr = C_NULL_PTR
     ! Symmetric ack signals; ack_signals[r] is set to `iter` by PE r
     ! via shmem_uint64_atomic_set after it has consumed our most
     ! recent put. Sized to pe_size on every PE.
     type(c_ptr) :: ack_signals_ptr = C_NULL_PTR
     ! Monotonically increasing iteration counter; sender writes this
     ! value into the receiver's data signal slot via SHMEM_SIGNAL_SET,
     ! and the receiver writes the same value into the sender's ack
     ! slot after consumption.
     integer(kind=i8) :: iter = 0
   contains
     procedure, pass(this) :: init => gs_shmem_init
     procedure, pass(this) :: free => gs_shmem_free
     procedure, pass(this) :: nbsend => gs_shmem_nbsend
     procedure, pass(this) :: nbrecv => gs_shmem_nbrecv
     procedure, pass(this) :: nbwait => gs_shmem_nbwait
     procedure, pass(this) :: nbsend_vec => gs_shmem_nbsend_vec
     procedure, pass(this) :: nbrecv_vec => gs_shmem_nbrecv_vec
     procedure, pass(this) :: nbwait_vec => gs_shmem_nbwait_vec
  end type gs_shmem_t

contains

  !> Allocate symmetric memory and per-neighbor bookkeeping for one
  !! direction of communication.
  !! @param pe_order  ranks in send_pe / recv_pe order (1-based)
  !! @param dof_stack indexed by rank, lower bound = 0 (per init_dofs)
  subroutine gs_shmem_buf_init(this, pe_order, dof_stack)
    class(gs_shmem_buf_t), intent(inout) :: this
    integer, intent(in) :: pe_order(:)
    type(stack_i4_t), intent(inout) :: dof_stack(0:)
    integer :: i, ierr, n
    integer(c_size_t) :: sz
    real(c_rp) :: rp_dummy
#ifndef HAVE_OPENSHMEM
    call neko_error('Neko was not built with OpenSHMEM support')
#else

    n = size(pe_order)

    allocate(this%ndofs(n))
    allocate(this%offset(n))
    allocate(this%remote_offset(n))

    do i = 1, n
       this%remote_offset(i) = -1
    end do

    this%total = 0
    do i = 1, n
       this%ndofs(i) = dof_stack(pe_order(i))%size()
       this%offset(i) = this%total
       this%total = this%total + this%ndofs(i)
    end do

    ! OpenSHMEM symmetric memory must be allocated with the same size on
    ! every PE.
    call MPI_Allreduce(this%total, this%max_total, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)

    ! Sized for up to GS_VEC_NC components so the fused vector path can
    ! reuse the same symmetric buffer; the scalar path uses the first
    ! max_total elements.
    sz = c_sizeof(rp_dummy) * int(max(GS_VEC_NC*this%max_total, 1), c_size_t)
    this%buf_ptr = shmem_malloc(sz)
    if (.not. c_associated(this%buf_ptr)) then
       call neko_error('shmem_malloc failed for gs_shmem buffer')
    end if
#endif
  end subroutine gs_shmem_buf_init

  !> Release symmetric memory and bookkeeping.
  subroutine gs_shmem_buf_free(this)
    class(gs_shmem_buf_t), intent(inout) :: this

    if (allocated(this%ndofs)) deallocate(this%ndofs)
    if (allocated(this%offset)) deallocate(this%offset)
    if (allocated(this%remote_offset)) deallocate(this%remote_offset)

#ifdef HAVE_OPENSHMEM
    if (c_associated(this%buf_ptr)) then
       ! shmem_free is collective; gs_shmem_free issues a barrier before
       ! tearing down buffers.
       call shmem_free(this%buf_ptr)
    end if
#endif
    this%buf_ptr = C_NULL_PTR
    this%total = 0
    this%max_total = 0

  end subroutine gs_shmem_buf_free

  !> Initialise OpenSHMEM based communication method
  subroutine gs_shmem_init(this, send_pe, recv_pe)
    class(gs_shmem_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer :: i, ierr
    integer(c_size_t) :: i64_size
    integer(c_int64_t) :: i64_dummy
    integer, allocatable :: local_offsets(:), remote_offsets(:)
#ifdef HAVE_OPENSHMEM
    integer(c_int) :: thread_level
#endif
#ifndef HAVE_OPENSHMEM
    call neko_error('Neko was not built with OpenSHMEM support')
#else

    ! No self-free here: gs_schedule has just filled send_dof/recv_dof and
    ! gs_shmem_free would deallocate them again, leaving the buffer setup
    ! below reading from a dangling descriptor. Like every other comm
    ! backend, init assumes a freshly allocated object.
    call this%init_order(send_pe, recv_pe)

    ! Decide once whether the OpenSHMEM library tolerates concurrent calls
    ! from several OpenMP threads; the send path branches on it.
    call shmem_query_thread(thread_level)
    gs_shmem_thread_multiple = (thread_level .eq. SHMEM_THREAD_MULTIPLE)

    call this%send_buf%init(this%send_pe, this%send_dof)
    call this%recv_buf%init(this%recv_pe, this%recv_dof)

    ! Allocate the per-rank symmetric signal arrays. Size pe_size on
    ! every PE -> no Allreduce needed and no slot handshake required;
    ! every PE writes/reads at offset = remote PE's rank.
    i64_size = c_sizeof(i64_dummy)
    this%data_signals_ptr = shmem_calloc(int(pe_size, c_size_t), i64_size)
    if (.not. c_associated(this%data_signals_ptr)) then
       call neko_error('shmem_calloc failed for gs_shmem data signals')
    end if
    this%ack_signals_ptr = shmem_calloc(int(pe_size, c_size_t), i64_size)
    if (.not. c_associated(this%ack_signals_ptr)) then
       call neko_error('shmem_calloc failed for gs_shmem ack signals')
    end if

    ! Offset exchange via Alltoall. send_pe and recv_pe are built in
    ! gs_schedule with a shifted-modulo neighbor pattern; pairwise
    ! Isend/Irecv keyed by neighbor rank turned out to deadlock at
    ! certain PE counts. Alltoall is a single collective call that
    ! cannot mismatch, and the payload is just pe_size ints per PE.
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

    this%iter = 0
    this%vec_supported = .true.

    ! Ensure all PEs have completed symmetric allocation before any
    ! one-sided communication is issued.
    call shmem_barrier_all()
#endif
  end subroutine gs_shmem_init

  !> Deallocate OpenSHMEM based communication method
  subroutine gs_shmem_free(this)
    class(gs_shmem_t), intent(inout) :: this

#ifdef HAVE_OPENSHMEM
    ! shmem_free is collective; synchronize first to ensure no in-flight
    ! puts target the buffers we are about to release.
    call shmem_barrier_all()

    if (c_associated(this%data_signals_ptr)) then
       call shmem_free(this%data_signals_ptr)
    end if
    if (c_associated(this%ack_signals_ptr)) then
       call shmem_free(this%ack_signals_ptr)
    end if
    this%data_signals_ptr = C_NULL_PTR
    this%ack_signals_ptr = C_NULL_PTR
    this%iter = 0
#endif

    call this%send_buf%free()
    call this%recv_buf%free()

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_shmem_free

  !> Pack the gathered shared dofs into the symmetric send buffer and
  !! issue non-blocking puts with signaling to each neighbor's recv
  !! buffer. Before each put, wait for the receiver's ack of our
  !! previous round so we never overwrite a buffer the receiver hasn't
  !! consumed yet.
  subroutine gs_shmem_nbsend(this, u, n, tag, deps, strm)
    class(gs_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, dst, base, ndst
    integer(c_size_t) :: nbytes
    integer , pointer :: sp(:)
    real(kind=rp), pointer :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:), ack_signals(:)
    real(c_rp) :: rp_dummy
#ifdef HAVE_OPENSHMEM

    ! Entered from inside the gs_op_vector OpenMP parallel region. Every
    ! thread binds its own pointers to the symmetric buffers (the c_f_pointer
    ! targets are private locals), but the iteration counter is shared state
    ! and must be bumped exactly once per round.
    !$omp master
    ! Each gs op gets a fresh signal value so receivers can distinguish
    ! puts from one call vs. the next.
    this%iter = this%iter + 1
    !$omp end master

    call c_f_pointer(this%send_buf%buf_ptr, send_data, &
         [max(this%send_buf%max_total, 1)])
    call c_f_pointer(this%recv_buf%buf_ptr, recv_data, &
         [max(this%recv_buf%max_total, 1)])
    call c_f_pointer(this%data_signals_ptr, data_signals, [pe_size])
    call c_f_pointer(this%ack_signals_ptr, ack_signals, [pe_size])

    ! Publish the new iter before anyone reads it below.
    !$omp barrier

    if (gs_shmem_thread_multiple) then
       ! Each thread drives the whole protocol for its share of the peers,
       ! keeping the per-peer wait/pack/put order of the serial code.
       !$omp do
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)

          ! Wait for dst to have consumed our previous round's data before
          ! we overwrite their recv buffer. ack_signals[dst] is calloc'd
          ! to 0 and dst writes the iter value (via atomic_set) after
          ! consumption, so for iter == 1 the wait succeeds immediately.
          call shmem_uint64_wait_until(c_loc(ack_signals(dst + 1)), &
               SHMEM_CMP_GE, this%iter - 1_8)

          sp => this%send_dof(dst)%array()
          base = this%send_buf%offset(i)
          ndst = this%send_buf%ndofs(i)
          do concurrent (j = 1:ndst)
             send_data(base + j) = u(sp(j))
          end do

          nbytes = int(ndst, c_size_t) * c_sizeof(rp_dummy)
          ! Put data + signal dst's data_signals[my_rank].
          call shmem_putmem_signal_nbi( &
               c_loc(recv_data(this%send_buf%remote_offset(i) + 1)), &
               c_loc(send_data(base + 1)), &
               nbytes, &
               c_loc(data_signals(pe_rank + 1)), &
               this%iter, SHMEM_SIGNAL_SET, dst)
       end do
       !$omp end do
    else
       ! Funnelled: every SHMEM call is left to the master and only the pack
       ! is work-shared, one peer at a time. The per-peer wait/pack/put order
       ! of the serial code is kept -- the ack is what makes a slab safe to
       ! repack as well, since it means the previous round's non-blocking put
       ! has been consumed and no longer reads our send buffer. The barrier
       ! after the ack wait releases the team into the pack, and the implicit
       ! barrier at end do completes the slab before the master puts it.
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          base = this%send_buf%offset(i)
          ndst = this%send_buf%ndofs(i)

          !$omp master
          call shmem_uint64_wait_until(c_loc(ack_signals(dst + 1)), &
               SHMEM_CMP_GE, this%iter - 1_8)
          !$omp end master
          !$omp barrier

          sp => this%send_dof(dst)%array()
          !OCL NORECURRENCE, NOVREC, NOALIAS
          !DIR$ CONCURRENT
          !DIR$ IVDEP
          !GCC$ ivdep
          !NEC$ IVDEP
          !$omp do
          do j = 1, ndst
             send_data(base + j) = u(sp(j))
          end do
          !$omp end do

          !$omp master
          nbytes = int(ndst, c_size_t) * c_sizeof(rp_dummy)
          call shmem_putmem_signal_nbi( &
               c_loc(recv_data(this%send_buf%remote_offset(i) + 1)), &
               c_loc(send_data(base + 1)), &
               nbytes, &
               c_loc(data_signals(pe_rank + 1)), &
               this%iter, SHMEM_SIGNAL_SET, dst)
          !$omp end master
       end do
       !$omp barrier
    end if
#endif
  end subroutine gs_shmem_nbsend

  !> No-op: receives are completed via remote put-with-signal.
  subroutine gs_shmem_nbrecv(this, tag)
    class(gs_shmem_t), intent(inout) :: this
    integer, intent(in) :: tag

  end subroutine gs_shmem_nbrecv

  !> Wait per-neighbor for the signal indicating that data has landed,
  !! apply the gather-scatter operation from the recv buffer into u, and
  !! ack the sender so they may overwrite the buffer in the next round.
  subroutine gs_shmem_nbwait(this, u, n, op, strm)
    class(gs_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer :: i, j, src, base, nsrc
    integer(c_int64_t) :: dummy
    integer , pointer :: sp(:)
    real(kind=rp), pointer :: recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:), ack_signals(:)
#ifdef HAVE_OPENSHMEM

    call c_f_pointer(this%recv_buf%buf_ptr, recv_data, &
         [max(this%recv_buf%max_total, 1)])
    call c_f_pointer(this%data_signals_ptr, data_signals, [pe_size])
    call c_f_pointer(this%ack_signals_ptr, ack_signals, [pe_size])

    ! The loop over peers stays serial even when the library is thread
    ! multiple: a dof shared by 3+ ranks appears in several recv_dof lists,
    ! so reducing two slabs concurrently would race on that dof. The master
    ! awaits and acks each slab in turn and the team reduces it, taking the
    ! parallelism within the slab instead.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       base = this%recv_buf%offset(i)
       nsrc = this%recv_buf%ndofs(i)

       ! Wait for data from src; data_signals[src] is set by src's put.
       !$omp master
       dummy = shmem_signal_wait_until(c_loc(data_signals(src + 1)), &
            SHMEM_CMP_GE, this%iter)
       !$omp end master
       !$omp barrier

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
             u(sp(j)) = u(sp(j)) + recv_data(base + j)
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
             u(sp(j)) = u(sp(j)) * recv_data(base + j)
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
             u(sp(j)) = min(u(sp(j)), recv_data(base + j))
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
             u(sp(j)) = max(u(sp(j)), recv_data(base + j))
          end do
          !$omp end do
       case default
          call neko_error("Unknown operation in gs_nbwait_shmem")
       end select

       ! Tell src we are done with this round's buffer so they may
       ! overwrite it on the next round. We set our own rank's slot
       ! in src's ack_signals; src waits there on the next nbsend. The
       ! implicit barrier at the end do above guarantees the whole team is
       ! finished with the slab before it is released.
       !$omp master
       call shmem_uint64_atomic_set( &
            c_loc(ack_signals(pe_rank + 1)), &
            this%iter, src)
       !$omp end master
    end do
#endif
  end subroutine gs_shmem_nbwait

  !> Fused nc-component send: pack nc contiguous component blocks per peer
  !! slab and put nc*ndofs reals with a single signal. Buffer indexing and
  !! put size scale by nc; the per-rank signalling is unchanged.
  !! @param u compact shared buffer, component-outer: u((c-1)*n + idx).
  subroutine gs_shmem_nbsend_vec(this, u, n, nc, tag, deps, strm)
    class(gs_shmem_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, c, dst, base, ndst
    integer(c_size_t) :: nbytes
    integer, pointer :: sp(:)
    real(kind=rp), pointer :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:), ack_signals(:)
    real(c_rp) :: rp_dummy
#ifdef HAVE_OPENSHMEM

    !$omp master
    this%iter = this%iter + 1
    !$omp end master

    call c_f_pointer(this%send_buf%buf_ptr, send_data, &
         [max(nc*this%send_buf%max_total, 1)])
    call c_f_pointer(this%recv_buf%buf_ptr, recv_data, &
         [max(nc*this%recv_buf%max_total, 1)])
    call c_f_pointer(this%data_signals_ptr, data_signals, [pe_size])
    call c_f_pointer(this%ack_signals_ptr, ack_signals, [pe_size])

    !$omp barrier

    ! Same threading split as the scalar path; see gs_shmem_nbsend.
    if (gs_shmem_thread_multiple) then
       !$omp do
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)

          call shmem_uint64_wait_until(c_loc(ack_signals(dst + 1)), &
               SHMEM_CMP_GE, this%iter - 1_8)

          sp => this%send_dof(dst)%array()
          base = this%send_buf%offset(i)
          ndst = this%send_buf%ndofs(i)
          do c = 1, nc
             do concurrent (j = 1:ndst)
                send_data(nc*base + (c-1)*ndst + j) = u((c-1)*n + sp(j))
             end do
          end do

          nbytes = int(nc*ndst, c_size_t) * c_sizeof(rp_dummy)
          call shmem_putmem_signal_nbi( &
               c_loc(recv_data(nc*this%send_buf%remote_offset(i) + 1)), &
               c_loc(send_data(nc*base + 1)), &
               nbytes, &
               c_loc(data_signals(pe_rank + 1)), &
               this%iter, SHMEM_SIGNAL_SET, dst)
       end do
       !$omp end do
    else
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          base = this%send_buf%offset(i)
          ndst = this%send_buf%ndofs(i)

          !$omp master
          call shmem_uint64_wait_until(c_loc(ack_signals(dst + 1)), &
               SHMEM_CMP_GE, this%iter - 1_8)
          !$omp end master
          !$omp barrier

          sp => this%send_dof(dst)%array()
          !$omp do
          do j = 1, ndst
             do c = 1, nc
                send_data(nc*base + (c-1)*ndst + j) = u((c-1)*n + sp(j))
             end do
          end do
          !$omp end do

          !$omp master
          nbytes = int(nc*ndst, c_size_t) * c_sizeof(rp_dummy)
          call shmem_putmem_signal_nbi( &
               c_loc(recv_data(nc*this%send_buf%remote_offset(i) + 1)), &
               c_loc(send_data(nc*base + 1)), &
               nbytes, &
               c_loc(data_signals(pe_rank + 1)), &
               this%iter, SHMEM_SIGNAL_SET, dst)
          !$omp end master
       end do
       !$omp barrier
    end if
#endif
  end subroutine gs_shmem_nbsend_vec

  !> No-op: receives are completed via remote put-with-signal.
  subroutine gs_shmem_nbrecv_vec(this, tag, nc)
    class(gs_shmem_t), intent(inout) :: this
    integer, intent(in) :: tag, nc
  end subroutine gs_shmem_nbrecv_vec

  !> Fused nc-component wait/reduce: per peer, wait on the data signal and
  !! reduce nc component blocks into u, then ack the sender.
  subroutine gs_shmem_nbwait_vec(this, u, n, nc, op, strm)
    class(gs_shmem_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer :: i, j, c, src, base, nsrc
    integer(c_int64_t) :: dummy
    integer, pointer :: sp(:)
    real(kind=rp), pointer :: recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:), ack_signals(:)
#ifdef HAVE_OPENSHMEM

    call c_f_pointer(this%recv_buf%buf_ptr, recv_data, &
         [max(nc*this%recv_buf%max_total, 1)])
    call c_f_pointer(this%data_signals_ptr, data_signals, [pe_size])
    call c_f_pointer(this%ack_signals_ptr, ack_signals, [pe_size])

    ! Serial over peers (a dof shared by 3+ ranks appears in several recv
    ! lists); parallelism is taken within each slab. See gs_shmem_nbwait.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       base = this%recv_buf%offset(i)
       nsrc = this%recv_buf%ndofs(i)

       !$omp master
       dummy = shmem_signal_wait_until(c_loc(data_signals(src + 1)), &
            SHMEM_CMP_GE, this%iter)
       !$omp end master
       !$omp barrier

       sp => this%recv_dof(src)%array()
       select case (op)
       case (GS_OP_ADD)
          !$omp do
          do j = 1, nsrc
             do c = 1, nc
                u((c-1)*n + sp(j)) = u((c-1)*n + sp(j)) + &
                     recv_data(nc*base + (c-1)*nsrc + j)
             end do
          end do
          !$omp end do
       case (GS_OP_MUL)
          !$omp do
          do j = 1, nsrc
             do c = 1, nc
                u((c-1)*n + sp(j)) = u((c-1)*n + sp(j)) * &
                     recv_data(nc*base + (c-1)*nsrc + j)
             end do
          end do
          !$omp end do
       case (GS_OP_MIN)
          !$omp do
          do j = 1, nsrc
             do c = 1, nc
                u((c-1)*n + sp(j)) = min(u((c-1)*n + sp(j)), &
                     recv_data(nc*base + (c-1)*nsrc + j))
             end do
          end do
          !$omp end do
       case (GS_OP_MAX)
          !$omp do
          do j = 1, nsrc
             do c = 1, nc
                u((c-1)*n + sp(j)) = max(u((c-1)*n + sp(j)), &
                     recv_data(nc*base + (c-1)*nsrc + j))
             end do
          end do
          !$omp end do
       case default
          call neko_error("Unknown operation in gs_shmem_nbwait_vec")
       end select

       !$omp master
       call shmem_uint64_atomic_set( &
            c_loc(ack_signals(pe_rank + 1)), &
            this%iter, src)
       !$omp end master
    end do
#endif
  end subroutine gs_shmem_nbwait_vec

end module gs_shmem
