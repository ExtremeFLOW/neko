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
  use num_types, only : rp, c_rp
  use gs_comm, only : gs_comm_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use comm, only : NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_Sendrecv, MPI_INTEGER, MPI_MAX, &
       MPI_STATUS_IGNORE
  use utils, only : neko_error
#ifdef HAVE_OPENSHMEM
  use shmem, only : shmem_malloc, shmem_calloc, shmem_free, &
       shmem_putmem_signal_nbi, shmem_signal_wait_until, &
       shmem_uint64_atomic_set, shmem_barrier_all, &
       SHMEM_SIGNAL_SET, SHMEM_CMP_GE
#endif
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_loc, &
       c_f_pointer, c_associated, c_sizeof, c_size_t, c_int64_t
  implicit none
  private

  !> Symmetric buffer for one direction of OpenSHMEM communication
  type :: gs_shmem_buf_t
     !> Number of dofs per neighbor
     integer, allocatable :: ndofs(:)
     !> Local offset in this buffer for each neighbor
     integer, allocatable :: offset(:)
     !> For send_buf: offset in the remote PE's recv buffer where our
     !! data should land (= remote PE's recv-side offset for our data).
     !! Unused for recv_buf.
     integer, allocatable :: remote_offset(:)
     !> Slot index on the remote PE for our signal puts.
     !! For send_buf(i), this is the index in send_pe(i)'s data_signals
     !! array that we should set when delivering data to send_pe(i)
     !! (= send_pe(i)'s recv index for us).
     !! For recv_buf(i), this is the index in recv_pe(i)'s ack_signals
     !! array that we should set when acking recv_pe(i) (= recv_pe(i)'s
     !! send index for us). Both populated by MPI_Sendrecv during init.
     integer, allocatable :: remote_signal_idx(:)
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
  !! per-neighbor signaling for completion (OpenSHMEM 1.5).
  !!
  !! Each PE allocates symmetric send and recv buffers sized to the
  !! global maximum number of dofs (so puts can target a well-defined
  !! offset on the remote PE's recv buffer) plus two symmetric uint64
  !! signal arrays sized to the global maximum number of neighbors:
  !! data_signals (incoming "data has arrived" notifications) and
  !! ack_signals (incoming "buffer consumed, you may overwrite"
  !! notifications). The sender waits on ack_signals before each put so
  !! a fast sender can never overwrite a buffer the receiver hasn't
  !! consumed yet.
  type, public, extends(gs_comm_t) :: gs_shmem_t
     type(gs_shmem_buf_t) :: send_buf
     type(gs_shmem_buf_t) :: recv_buf
     !> Symmetric data-arrival signals; data_signals[i] is set to `iter`
     !! by recv_pe(i) when it has put data into our recv buffer.
     !! Sized to n_signals (global max neighbors) on every PE.
     type(c_ptr) :: data_signals_ptr = C_NULL_PTR
     !> Symmetric ack signals; ack_signals[i] is set to `iter` by
     !! send_pe(i) after it has consumed our most recent put. The sender
     !! on this PE waits on its own ack_signals[i] before issuing the
     !! next put to send_pe(i), preventing a buffer-overwrite race.
     !! Sized to n_signals (global max neighbors) on every PE.
     type(c_ptr) :: ack_signals_ptr = C_NULL_PTR
     !> Global maximum of max(size(send_pe), size(recv_pe)) across all
     !! PEs; this is the size of the symmetric signal allocations.
     integer :: n_signals = 0
     !> Monotonically increasing iteration counter; sender writes this
     !! value into the receiver's data signal slot via SHMEM_SIGNAL_SET,
     !! and the receiver waits for it. Receiver writes the same value
     !! into the sender's ack slot after consumption.
     integer(kind=8) :: iter = 0
   contains
     procedure, pass(this) :: init => gs_shmem_init
     procedure, pass(this) :: free => gs_shmem_free
     procedure, pass(this) :: nbsend => gs_shmem_nbsend
     procedure, pass(this) :: nbrecv => gs_shmem_nbrecv
     procedure, pass(this) :: nbwait => gs_shmem_nbwait
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
    allocate(this%remote_signal_idx(n))

    do i = 1, n
       this%remote_offset(i) = -1
       this%remote_signal_idx(i) = -1
    end do

    this%total = 0
    do i = 1, n
       this%ndofs(i) = dof_stack(pe_order(i))%size()
       this%offset(i) = this%total
       this%total = this%total + this%ndofs(i)
    end do

    call MPI_Allreduce(this%total, this%max_total, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)

    sz = c_sizeof(rp_dummy) * int(max(this%max_total, 1), c_size_t)
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
    if (allocated(this%remote_signal_idx)) deallocate(this%remote_signal_idx)

#ifdef HAVE_OPENSHMEM
    if (c_associated(this%buf_ptr)) then
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
    integer :: i, ierr, n_local, slot, n_alloc
    integer(c_size_t) :: i64_size
    integer(c_int64_t) :: i64_dummy
#ifndef HAVE_OPENSHMEM
    call neko_error('Neko was not built with OpenSHMEM support')
#else

    call this%init_order(send_pe, recv_pe)

    call this%send_buf%init(this%send_pe, this%send_dof)
    call this%recv_buf%init(this%recv_pe, this%recv_dof)

    ! Slot-exchange variant: signal arrays are sized to the global
    ! maximum of max(n_send, n_recv) so they scale with neighbor count
    ! rather than total PE count. Each PE tells its peers which slot
    ! to target via two MPI_Sendrecv exchanges (data slots on tag 1,
    ! ack slots on tag 2). The (recv_pe(i), send_pe(i)) pair in those
    ! Sendrecvs is generally a pair of *different* PEs because
    ! gs_schedule pushes send_pe and recv_pe in shifted-modulo order;
    ! MPI matches by source/dest + tag, so this still works.
    n_local = max(size(this%send_pe), size(this%recv_pe))
    call MPI_Allreduce(n_local, this%n_signals, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)
    n_alloc = max(this%n_signals, 1)
    i64_size = c_sizeof(i64_dummy)
    this%data_signals_ptr = shmem_calloc(int(n_alloc, c_size_t), i64_size)
    if (.not. c_associated(this%data_signals_ptr)) then
       call neko_error('shmem_calloc failed for gs_shmem data signals')
    end if
    this%ack_signals_ptr = shmem_calloc(int(n_alloc, c_size_t), i64_size)
    if (.not. c_associated(this%ack_signals_ptr)) then
       call neko_error('shmem_calloc failed for gs_shmem ack signals')
    end if

    ! Exchange recv-buffer offsets (existing): tell recv_pe(i) where in
    ! our recv buffer their data should land; learn from send_pe(i)
    ! where in their recv buffer our data should land.
    do i = 1, max(size(this%send_pe), size(this%recv_pe))
       if (i .le. size(this%send_pe) .and. i .le. size(this%recv_pe)) then
          call MPI_Sendrecv( &
               this%recv_buf%offset(i), 1, MPI_INTEGER, this%recv_pe(i), 0, &
               this%send_buf%remote_offset(i), 1, MPI_INTEGER, &
               this%send_pe(i), 0, &
               NEKO_COMM, MPI_STATUS_IGNORE, ierr)
       end if
    end do

    ! Exchange data-signal slots (tag 1): tell recv_pe(i) "my recv slot
    ! for incoming from you = i" so they target my data_signals[i]
    ! when sending; learn from send_pe(i) where to target their
    ! data_signals when we send to them -> send_buf%remote_signal_idx(i).
    do i = 1, max(size(this%send_pe), size(this%recv_pe))
       if (i .le. size(this%send_pe) .and. i .le. size(this%recv_pe)) then
          slot = i
          call MPI_Sendrecv( &
               slot, 1, MPI_INTEGER, this%recv_pe(i), 1, &
               this%send_buf%remote_signal_idx(i), 1, MPI_INTEGER, &
               this%send_pe(i), 1, &
               NEKO_COMM, MPI_STATUS_IGNORE, ierr)
       end if
    end do

    ! Exchange ack-signal slots (tag 2): tell send_pe(i) "my send slot
    ! for you = i" so they target my ack_signals[i] when acking; learn
    ! from recv_pe(i) where to target their ack_signals when we ack
    ! them -> recv_buf%remote_signal_idx(i).
    do i = 1, max(size(this%send_pe), size(this%recv_pe))
       if (i .le. size(this%send_pe) .and. i .le. size(this%recv_pe)) then
          slot = i
          call MPI_Sendrecv( &
               slot, 1, MPI_INTEGER, this%send_pe(i), 2, &
               this%recv_buf%remote_signal_idx(i), 1, MPI_INTEGER, &
               this%recv_pe(i), 2, &
               NEKO_COMM, MPI_STATUS_IGNORE, ierr)
       end if
    end do

    this%iter = 0

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
    this%n_signals = 0
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
  subroutine gs_shmem_nbsend(this, u, n, deps, strm)
    class(gs_shmem_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, dst, base
    integer(c_size_t) :: nbytes
    integer , pointer :: sp(:)
    real(kind=rp), pointer :: send_data(:), recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:), ack_signals(:)
    integer(c_int64_t) :: dummy
    real(c_rp) :: rp_dummy
#ifdef HAVE_OPENSHMEM

    ! Each gs op gets a fresh signal value so receivers can distinguish
    ! puts from one call vs. the next.
    this%iter = this%iter + 1

    call c_f_pointer(this%send_buf%buf_ptr, send_data, &
         [max(this%send_buf%max_total, 1)])
    call c_f_pointer(this%recv_buf%buf_ptr, recv_data, &
         [max(this%recv_buf%max_total, 1)])
    call c_f_pointer(this%data_signals_ptr, data_signals, &
         [max(this%n_signals, 1)])
    call c_f_pointer(this%ack_signals_ptr, ack_signals, &
         [max(this%n_signals, 1)])

    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)

       ! Wait for the receiver to have consumed our previous round's
       ! data before we overwrite their recv buffer. ack_signals is
       ! calloc'd to 0 and the receiver writes the iter value after
       ! consumption, so for iter == 1 the wait succeeds immediately.
       ! Local index i corresponds to send_pe(i) by construction.
       dummy = shmem_signal_wait_until(c_loc(ack_signals(i)), &
            SHMEM_CMP_GE, this%iter - 1_8)

       sp => this%send_dof(dst)%array()
       base = this%send_buf%offset(i)
        do concurrent (j = 1:this%send_dof(dst)%size())
          send_data(base + j) = u(sp(j))
       end do

       nbytes = int(this%send_buf%ndofs(i), c_size_t) * c_sizeof(rp_dummy)
       call shmem_putmem_signal_nbi( &
            c_loc(recv_data(this%send_buf%remote_offset(i) + 1)), &
            c_loc(send_data(base + 1)), &
            nbytes, &
            c_loc(data_signals(this%send_buf%remote_signal_idx(i))), &
            this%iter, SHMEM_SIGNAL_SET, dst)
    end do
#endif
  end subroutine gs_shmem_nbsend

  !> No-op: receives are completed via remote put-with-signal.
  subroutine gs_shmem_nbrecv(this)
    class(gs_shmem_t), intent(inout) :: this

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
    integer :: i, j, src, base
    integer(c_int64_t) :: dummy
    integer , pointer :: sp(:)
    real(kind=rp), pointer :: recv_data(:)
    integer(c_int64_t), pointer :: data_signals(:), ack_signals(:)
#ifdef HAVE_OPENSHMEM

    call c_f_pointer(this%recv_buf%buf_ptr, recv_data, &
         [max(this%recv_buf%max_total, 1)])
    call c_f_pointer(this%data_signals_ptr, data_signals, &
         [max(this%n_signals, 1)])
    call c_f_pointer(this%ack_signals_ptr, ack_signals, &
         [max(this%n_signals, 1)])

    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)

       ! Local index i corresponds to recv_pe(i) by construction.
       dummy = shmem_signal_wait_until(c_loc(data_signals(i)), &
            SHMEM_CMP_GE, this%iter)

       sp => this%recv_dof(src)%array()
       base = this%recv_buf%offset(i)
       select case (op)
       case (GS_OP_ADD)
          !NEC$ IVDEP
          do concurrent (j = 1:this%recv_dof(src)%size())
             u(sp(j)) = u(sp(j)) + recv_data(base + j)
          end do
       case (GS_OP_MUL)
          !NEC$ IVDEP
          do concurrent (j = 1:this%recv_dof(src)%size())
             u(sp(j)) = u(sp(j)) * recv_data(base + j)
          end do
       case (GS_OP_MIN)
          !NEC$ IVDEP
          do concurrent (j = 1:this%recv_dof(src)%size())
             u(sp(j)) = min(u(sp(j)), recv_data(base + j))
          end do
       case (GS_OP_MAX)
          !NEC$ IVDEP
          do concurrent (j = 1:this%recv_dof(src)%size())
             u(sp(j)) = max(u(sp(j)), recv_data(base + j))
          end do
       case default
          call neko_error("Unknown operation in gs_nbwait_shmem")
       end select

       ! Tell the sender we are done with this round's buffer so they
       ! may overwrite it on the next round.
       call shmem_uint64_atomic_set( &
            c_loc(ack_signals(this%recv_buf%remote_signal_idx(i))), &
            this%iter, src)
    end do
#endif
  end subroutine gs_shmem_nbwait

end module gs_shmem
