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
!> Defines MPI one-sided (RMA) gather-scatter communication
module gs_mpi_rma
  use num_types, only : rp, i8
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use comm, only : NEKO_COMM, pe_rank, pe_size, MPI_REAL_PRECISION
  use utils, only : neko_error
  use mpi_f08, only : MPI_Win, MPI_Info, MPI_Win_allocate, MPI_Win_free, &
       MPI_Win_lock_all, MPI_Win_unlock_all, MPI_Win_flush, &
       MPI_Win_flush_all, &
       MPI_Win_sync, MPI_Win_get_attr, MPI_Put, MPI_Accumulate, &
       MPI_Alltoall, MPI_Barrier, MPI_Iprobe, MPI_Status, &
       MPI_Info_create, MPI_Info_set, MPI_Info_free, &
       MPI_INTEGER, MPI_INTEGER8, MPI_REPLACE, MPI_MODE_NOCHECK, &
       MPI_ADDRESS_KIND, MPI_WIN_MODEL, MPI_WIN_UNIFIED, &
       MPI_ANY_SOURCE, MPI_ANY_TAG
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_f_pointer, &
       c_associated, c_sizeof
  implicit none
  private

  !> MPI RMA needs nothing beyond MPI-3, so the backend is always built.
  !! Kept as a parameter for symmetry with the other one-sided backends,
  !! which the gs comm. autotuner queries before benchmarking them.
  logical, parameter, public :: GS_MPI_RMA_AVAIL = .true.

  !> Number of spin iterations between the progress pokes in gs_mpi_rma_wait_ge
  integer, parameter :: GS_MPI_RMA_POKE = 64

  !> Whether to complete the payload puts with a single MPI_Win_flush_all
  !! rather than one MPI_Win_flush per peer. Batched is the default,
  !! set NEKO_GS_RMA_FLUSH_ALL=0 for the per-peer form.
  logical, save :: gs_mpi_rma_flush_all = .true.
  !> Whether the flush strategy has been read from the environment. It is a
  !! program-wide binding, read once on the first init.
  logical, save :: gs_mpi_rma_flush_bound = .false.

  !> Gather-scatter communication using MPI one-sided puts into a passive
  !! target window, with per-rank signalling for completion.
  type, public, extends(gs_comm_t) :: gs_mpi_rma_t
     !> Origin buffer for the puts. Plain local memory, no window needed.
     real(kind=rp), allocatable :: send_buf(:)
     !> Number of dofs to send to each peer, in send_pe order
     integer, allocatable :: send_ndofs(:)
     !> Offset of each peer's slab in send_buf (scalar path units)
     integer, allocatable :: send_offset(:)
     !> Offset of our slab in each peer's receive window (scalar path
     !! units), exchanged at init
     integer, allocatable :: send_rdisp(:)

     !> Receive window, holding GS_VEC_NC*recv_total reals
     type(MPI_Win) :: win_data
     !> Base of the receive window
     type(c_ptr) :: recv_ptr = C_NULL_PTR
     !> Number of dofs received from each peer, in recv_pe order
     integer, allocatable :: recv_ndofs(:)
     !> Offset of each peer's slab in the receive window (scalar path units)
     integer, allocatable :: recv_offset(:)
     !> Total number of received dofs on this rank
     integer :: recv_total = 0

     !> Signal window, 2*pe_size counters (data signals then acks)
     type(MPI_Win) :: win_sig
     !> Base of the signal window
     type(c_ptr) :: sig_ptr = C_NULL_PTR

     !> Monotonically increasing round counter. The sender replaces the
     !! receiver's data signal with it; the receiver replaces the sender's
     !! ack with the same value once consumed.
     integer(kind=i8) :: iter = 0

     !> Whether the windows use MPI_WIN_UNIFIED. In the separate memory
     !! model a remote update is only guaranteed visible to a local load
     !! after MPI_Win_sync, so the spin waits call it every iteration.
     logical :: unified = .true.

     !> Whether the windows have been allocated, so free can tell a live
     !! backend from one that never got past init
     logical :: win_alloc = .false.
   contains
     procedure, pass(this) :: init => gs_mpi_rma_init
     procedure, pass(this) :: free => gs_mpi_rma_free
     procedure, pass(this) :: nbsend => gs_mpi_rma_nbsend
     procedure, pass(this) :: nbrecv => gs_mpi_rma_nbrecv
     procedure, pass(this) :: nbwait => gs_mpi_rma_nbwait
     procedure, pass(this) :: nbsend_vec => gs_mpi_rma_nbsend_vec
     procedure, pass(this) :: nbrecv_vec => gs_mpi_rma_nbrecv_vec
     procedure, pass(this) :: nbwait_vec => gs_mpi_rma_nbwait_vec
  end type gs_mpi_rma_t

contains

  !> Initialise MPI RMA based communication method
  subroutine gs_mpi_rma_init(this, send_pe, recv_pe)
    class(gs_mpi_rma_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    type(MPI_Info) :: info
    integer(kind=MPI_ADDRESS_KIND) :: win_size, attr_val
    integer, allocatable :: local_disp(:), remote_disp(:)
    integer(kind=i8), pointer :: sig(:)
    integer :: i, n, send_total, ierr, disp_unit, env_len
    integer(kind=i8) :: i8_dummy
    real(kind=rp) :: rp_dummy
    character(len=255) :: env_val
    logical :: flag

    ! No self-free here: gs_schedule has just filled send_dof/recv_dof and
    ! gs_mpi_rma_free would deallocate them again. Like every other comm
    ! backend, init assumes a freshly allocated object.
    call this%init_order(send_pe, recv_pe)

    ! Bind the flush strategy once; every instance shares it, so that an
    ! A/B run measures one protocol rather than a mixture.
    if (.not. gs_mpi_rma_flush_bound) then
       call get_environment_variable("NEKO_GS_RMA_FLUSH_ALL", env_val, env_len)
       if (env_len .gt. 0) then
          env_len = min(env_len, len(env_val))
          gs_mpi_rma_flush_all = (env_val(1:env_len) .ne. '0')
       end if
       gs_mpi_rma_flush_bound = .true.
    end if

    ! Send side: exact local sizing, plain memory
    n = size(this%send_pe)
    allocate(this%send_ndofs(n), this%send_offset(n), this%send_rdisp(n))
    send_total = 0
    do i = 1, n
       this%send_ndofs(i) = this%send_dof(this%send_pe(i))%size()
       this%send_offset(i) = send_total
       send_total = send_total + this%send_ndofs(i)
    end do
    this%send_rdisp = -1
    allocate(this%send_buf(max(GS_VEC_NC*send_total, 1)))

    ! Receive side: exact local sizing, exposed in a window
    n = size(this%recv_pe)
    allocate(this%recv_ndofs(n), this%recv_offset(n))
    this%recv_total = 0
    do i = 1, n
       this%recv_ndofs(i) = this%recv_dof(this%recv_pe(i))%size()
       this%recv_offset(i) = this%recv_total
       this%recv_total = this%recv_total + this%recv_ndofs(i)
    end do

    ! same_disp_unit lets the implementation skip the per-target lookup on
    ! every operation.
    !
    ! The two accumulate assertions matter more than they look. By default
    ! MPI must preserve rar/raw/war/waw ordering between accumulates, which
    ! can rule out a hardware-atomic path for the signals in favour of a
    ! general, ordered one -- and the signals are two of the four RMA calls
    ! this backend makes per peer per round. Neither assertion costs
    ! anything here: a round issues at most one accumulate to any given
    ! location and consecutive rounds are separated by a flush, so no two
    ! are ever outstanding against the same slot, and every accumulate in
    ! the protocol is MPI_REPLACE.
    call MPI_Info_create(info, ierr)
    call MPI_Info_set(info, 'same_disp_unit', 'true', ierr)
    call MPI_Info_set(info, 'accumulate_ordering', 'none', ierr)
    call MPI_Info_set(info, 'accumulate_ops', 'same_op', ierr)

    disp_unit = int(c_sizeof(rp_dummy))
    win_size = int(max(GS_VEC_NC*this%recv_total, 1), MPI_ADDRESS_KIND) * &
         int(disp_unit, MPI_ADDRESS_KIND)
    call MPI_Win_allocate(win_size, disp_unit, info, NEKO_COMM, &
         this%recv_ptr, this%win_data, ierr)
    if (.not. c_associated(this%recv_ptr)) then
       call neko_error('MPI_Win_allocate failed for gs_mpi_rma data window')
    end if

    disp_unit = int(c_sizeof(i8_dummy))
    win_size = int(2*pe_size, MPI_ADDRESS_KIND) * &
         int(disp_unit, MPI_ADDRESS_KIND)
    call MPI_Win_allocate(win_size, disp_unit, info, NEKO_COMM, &
         this%sig_ptr, this%win_sig, ierr)
    if (.not. c_associated(this%sig_ptr)) then
       call neko_error('MPI_Win_allocate failed for gs_mpi_rma signal window')
    end if

    call MPI_Info_free(info, ierr)
    this%win_alloc = .true.

    ! MPI_Win_allocate does not zero the memory, and the protocol reads the
    ! counters before anyone has written them (iter == 1 waits for an ack
    ! of >= 0).
    call c_f_pointer(this%sig_ptr, sig, [2*pe_size])
    sig = 0_i8

    ! Both windows come from the same call sequence, so one query answers
    ! for both.
    call MPI_Win_get_attr(this%win_sig, MPI_WIN_MODEL, attr_val, flag, ierr)
    this%unified = flag .and. (int(attr_val) .eq. MPI_WIN_UNIFIED)

    call MPI_Win_lock_all(MPI_MODE_NOCHECK, this%win_data, ierr)
    call MPI_Win_lock_all(MPI_MODE_NOCHECK, this%win_sig, ierr)

    ! Tell each peer where in our receive window its slab belongs. Alltoall
    ! rather than pairwise Isend/Irecv for the same reason as in gs_shmem:
    ! send_pe and recv_pe follow a shifted-modulo pattern that deadlocked
    ! pairwise exchanges at certain rank counts, and the payload is one int
    ! per rank.
    allocate(local_disp(0:pe_size - 1), remote_disp(0:pe_size - 1))
    local_disp = -1
    do i = 1, size(this%recv_pe)
       local_disp(this%recv_pe(i)) = this%recv_offset(i)
    end do
    call MPI_Alltoall(local_disp, 1, MPI_INTEGER, &
         remote_disp, 1, MPI_INTEGER, NEKO_COMM, ierr)
    do i = 1, size(this%send_pe)
       this%send_rdisp(i) = remote_disp(this%send_pe(i))
    end do
    deallocate(local_disp, remote_disp)

    this%iter = 0
    this%vec_supported = .true.

    ! The zeroing above is a local store into window memory; no rank may
    ! issue RMA until every rank has done it.
    call MPI_Barrier(NEKO_COMM, ierr)

  end subroutine gs_mpi_rma_init

  !> Deallocate MPI RMA based communication method
  subroutine gs_mpi_rma_free(this)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer :: ierr

    if (this%win_alloc) then
       ! MPI_Win_free is collective and completes the epoch; barrier first
       ! so no rank is still spinning on, or putting into, a window that
       ! another rank is about to tear down.
       call MPI_Barrier(NEKO_COMM, ierr)
       call MPI_Win_unlock_all(this%win_sig, ierr)
       call MPI_Win_unlock_all(this%win_data, ierr)
       call MPI_Win_free(this%win_sig, ierr)
       call MPI_Win_free(this%win_data, ierr)
       this%win_alloc = .false.
    end if
    this%recv_ptr = C_NULL_PTR
    this%sig_ptr = C_NULL_PTR
    this%recv_total = 0
    this%iter = 0

    if (allocated(this%send_buf)) deallocate(this%send_buf)
    if (allocated(this%send_ndofs)) deallocate(this%send_ndofs)
    if (allocated(this%send_offset)) deallocate(this%send_offset)
    if (allocated(this%send_rdisp)) deallocate(this%send_rdisp)
    if (allocated(this%recv_ndofs)) deallocate(this%recv_ndofs)
    if (allocated(this%recv_offset)) deallocate(this%recv_offset)

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_mpi_rma_free

  !> Reload a signal counter through a VOLATILE dummy, so the spin in
  !! gs_mpi_rma_wait_ge reads memory on every iteration rather than a value the
  !! compiler hoisted out of the loop.
  !! @param s the signal array to read from
  !! @param slot 1-based slot to read
  !! @return its current value
  function gs_mpi_rma_load(s, slot) result(v)
    ! The whole array is the dummy rather than the one element that is read:
    ! a VOLATILE dummy fed from a pointer array has to be assumed-shape or a
    ! pointer array itself. No INTENT on it either, since gfortran rejects
    ! INTENT(IN) together with VOLATILE. Read only.
    integer(kind=i8), volatile :: s(:)
    integer, intent(in) :: slot
    integer(kind=i8) :: v

    v = s(slot)

  end function gs_mpi_rma_load

  !> Spin until the local signal counter in slot @a slot has reached @a val.
  !! Called by the master thread only.
  !! @param this the backend whose signal window is polled
  !! @param slot 1-based slot in the signal window: rank+1 for a data
  !! signal, pe_size+rank+1 for an ack
  !! @param val the value to wait for, compared with >=
  subroutine gs_mpi_rma_wait_ge(this, slot, val)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer, intent(in) :: slot
    integer(kind=i8), intent(in) :: val
    integer(kind=i8), pointer :: sig(:)
    type(MPI_Status) :: status
    integer :: spin, ierr
    logical :: flag

    call c_f_pointer(this%sig_ptr, sig, [2*pe_size])

    spin = 0
    do
       ! In the separate memory model the public window copy is only made
       ! visible to local loads by a sync; in the unified model this is a
       ! cheap memory barrier.
       if (.not. this%unified) call MPI_Win_sync(this%win_sig, ierr)

       if (gs_mpi_rma_load(sig, slot) .ge. val) exit

       ! Poke MPI so implementations without asynchronous progress still
       ! advance incoming RMA. Not needed on hardware-driven one-sided
       ! components, hence only every GS_MPI_RMA_POKE iterations.
       spin = spin + 1
       if (mod(spin, GS_MPI_RMA_POKE) .eq. 0) then
          call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, NEKO_COMM, flag, &
               status, ierr)
       end if
    end do

  end subroutine gs_mpi_rma_wait_ge

  !> Pack the gathered shared dofs and put them into each neighbour's
  !! receive window, then announce the puts. See the type comment for why
  !! the puts, the flush and the signals are separate phases.
  subroutine gs_mpi_rma_nbsend(this, u, n, tag, deps, strm)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer(kind=MPI_ADDRESS_KIND) :: rdisp
    integer, pointer :: sp(:)
    integer :: i, j, dst, base, ndst, ierr

    ! Entered from inside the gs_op_vector OpenMP parallel region. Every MPI
    ! call is funnelled through the master: RMA under MPI_THREAD_MULTIPLE is
    ! the least tuned path in most implementations.
    !$omp master
    this%iter = this%iter + 1
    !$omp end master

    ! Pack and put one peer at a time.
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       base = this%send_offset(i)
       ndst = this%send_ndofs(i)
       sp => this%send_dof(dst)%array()
       !OCL NORECURRENCE, NOVREC, NOALIAS
       !DIR$ CONCURRENT
       !DIR$ IVDEP
       !GCC$ ivdep
       !NEC$ IVDEP
       !$omp do
       do j = 1, ndst
          this%send_buf(base + j) = u(sp(j))
       end do
       !$omp end do

       ! Wait for dst to have consumed our previous round before we
       ! overwrite its receive window. The ack slots start at 0 and dst
       ! writes iter after consuming, so round 1 passes immediately.
       !$omp master
       call gs_mpi_rma_wait_ge(this, pe_size + dst + 1, this%iter - 1_i8)

       rdisp = int(this%send_rdisp(i), MPI_ADDRESS_KIND)
       call MPI_Put(this%send_buf(base + 1), ndst, MPI_REAL_PRECISION, &
            dst, rdisp, ndst, MPI_REAL_PRECISION, this%win_data, ierr)
       !$omp end master
    end do

    ! Announce each peer as soon as its own put is remotely complete. MPI
    ! orders accumulates only against the same location, so the flush is
    ! what stands in for a put-with-notify.
    !$omp master
    if (gs_mpi_rma_flush_all) call MPI_Win_flush_all(this%win_data, ierr)
    rdisp = int(pe_rank, MPI_ADDRESS_KIND)
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       if (.not. gs_mpi_rma_flush_all) then
          call MPI_Win_flush(dst, this%win_data, ierr)
       end if
       call MPI_Accumulate(this%iter, 1, MPI_INTEGER8, dst, rdisp, 1, &
            MPI_INTEGER8, MPI_REPLACE, this%win_sig, ierr)
    end do
    call MPI_Win_flush_all(this%win_sig, ierr)
    !$omp end master
    !$omp barrier

  end subroutine gs_mpi_rma_nbsend

  !> No-op: receives are completed by the remote put and its signal.
  subroutine gs_mpi_rma_nbrecv(this, tag)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer, intent(in) :: tag

  end subroutine gs_mpi_rma_nbrecv

  !> Wait per neighbour for the signal that its data has landed, reduce the
  !! slab into u, and ack the sender so it may overwrite the slab next round.
  subroutine gs_mpi_rma_nbwait(this, u, n, op, strm)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer(kind=MPI_ADDRESS_KIND) :: rdisp
    integer, pointer :: sp(:)
    real(kind=rp), pointer :: recv_data(:)
    integer :: i, j, src, base, nsrc, ierr

    call c_f_pointer(this%recv_ptr, recv_data, [max(this%recv_total, 1)])

    ! Serial over peers: a dof shared by three or more ranks appears in
    ! several recv_dof lists, so reducing two slabs concurrently would race
    ! on it. The parallelism is taken within each slab instead.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       base = this%recv_offset(i)
       nsrc = this%recv_ndofs(i)

       !$omp master
       call gs_mpi_rma_wait_ge(this, src + 1, this%iter)
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
          call neko_error("Unknown operation in gs_mpi_rma_nbwait")
       end select

       ! The implicit barrier at end do guarantees the whole team is done
       ! with the slab before it is released back to src.
       !$omp master
       rdisp = int(pe_size + pe_rank, MPI_ADDRESS_KIND)
       call MPI_Accumulate(this%iter, 1, MPI_INTEGER8, src, rdisp, 1, &
            MPI_INTEGER8, MPI_REPLACE, this%win_sig, ierr)
       !$omp end master
    end do

    ! An accumulate is only guaranteed to be issued at a flush, so the acks
    ! have to be pushed out before we go back to computing.
    !$omp master
    call MPI_Win_flush_all(this%win_sig, ierr)
    !$omp end master
    !$omp barrier

  end subroutine gs_mpi_rma_nbwait

  !> Fused nc-component send: pack nc contiguous component blocks per peer
  !! slab and put nc*ndofs reals. Buffer indexing and put sizes scale by nc;
  !! the signalling is unchanged.
  !! @param u compact shared buffer, component-outer: u((c-1)*n + idx).
  subroutine gs_mpi_rma_nbsend_vec(this, u, n, nc, tag, deps, strm)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer(kind=MPI_ADDRESS_KIND) :: rdisp
    integer, pointer :: sp(:)
    integer :: i, j, c, dst, base, ndst, ierr

    !$omp master
    this%iter = this%iter + 1
    !$omp end master

    ! Pack and put per peer, as in the scalar path; see gs_mpi_rma_nbsend.
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       base = this%send_offset(i)
       ndst = this%send_ndofs(i)
       sp => this%send_dof(dst)%array()
       !$omp do
       do j = 1, ndst
          do c = 1, nc
             this%send_buf(nc*base + (c-1)*ndst + j) = u((c-1)*n + sp(j))
          end do
       end do
       !$omp end do

       !$omp master
       call gs_mpi_rma_wait_ge(this, pe_size + dst + 1, this%iter - 1_i8)

       rdisp = int(nc*this%send_rdisp(i), MPI_ADDRESS_KIND)
       call MPI_Put(this%send_buf(nc*base + 1), nc*ndst, &
            MPI_REAL_PRECISION, dst, rdisp, nc*ndst, MPI_REAL_PRECISION, &
            this%win_data, ierr)
       !$omp end master
    end do

    !$omp master
    if (gs_mpi_rma_flush_all) call MPI_Win_flush_all(this%win_data, ierr)
    rdisp = int(pe_rank, MPI_ADDRESS_KIND)
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       if (.not. gs_mpi_rma_flush_all) then
          call MPI_Win_flush(dst, this%win_data, ierr)
       end if
       call MPI_Accumulate(this%iter, 1, MPI_INTEGER8, dst, rdisp, 1, &
            MPI_INTEGER8, MPI_REPLACE, this%win_sig, ierr)
    end do
    call MPI_Win_flush_all(this%win_sig, ierr)
    !$omp end master
    !$omp barrier

  end subroutine gs_mpi_rma_nbsend_vec

  !> No-op: receives are completed by the remote put and its signal.
  subroutine gs_mpi_rma_nbrecv_vec(this, tag, nc)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer, intent(in) :: tag, nc

  end subroutine gs_mpi_rma_nbrecv_vec

  !> Fused nc-component wait/reduce: per peer, wait on the data signal and
  !! reduce nc component blocks into u, then ack the sender.
  subroutine gs_mpi_rma_nbwait_vec(this, u, n, nc, op, strm)
    class(gs_mpi_rma_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer(kind=MPI_ADDRESS_KIND) :: rdisp
    integer, pointer :: sp(:)
    real(kind=rp), pointer :: recv_data(:)
    integer :: i, j, c, src, base, nsrc, ierr

    call c_f_pointer(this%recv_ptr, recv_data, [max(nc*this%recv_total, 1)])

    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       base = this%recv_offset(i)
       nsrc = this%recv_ndofs(i)

       !$omp master
       call gs_mpi_rma_wait_ge(this, src + 1, this%iter)
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
          call neko_error("Unknown operation in gs_mpi_rma_nbwait_vec")
       end select

       !$omp master
       rdisp = int(pe_size + pe_rank, MPI_ADDRESS_KIND)
       call MPI_Accumulate(this%iter, 1, MPI_INTEGER8, src, rdisp, 1, &
            MPI_INTEGER8, MPI_REPLACE, this%win_sig, ierr)
       !$omp end master
    end do

    !$omp master
    call MPI_Win_flush_all(this%win_sig, ierr)
    !$omp end master
    !$omp barrier

  end subroutine gs_mpi_rma_nbwait_vec

end module gs_mpi_rma
