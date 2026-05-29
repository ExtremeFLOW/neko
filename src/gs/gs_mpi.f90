! Copyright (c) 2020-2026, The Neko Authors
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
!> Defines MPI gather-scatter communication
module gs_mpi
  use num_types, only : rp
  use gs_comm, only : gs_comm_t, GS_COMM_MPI, GS_COMM_MPIGPU
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use mpi_f08, only : MPI_STATUSES_IGNORE, MPI_Status, &
       MPI_Request, MPI_Isend, MPI_IRecv, MPI_Testsome, MPI_Testall, &
       MPI_THREAD_MULTIPLE
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, NEKO_MPI_THREAD_PROVIDED
  use, intrinsic :: iso_c_binding
  use utils, only : neko_error
  implicit none
  private

  !> Gather-scatter communication using MPI.
  type, public, extends(gs_comm_t) :: gs_mpi_t
     !> Concatenated send slabs, one per peer.
     real(kind=rp), allocatable :: send_buf(:)
     !> Concatenated recv slabs, one per peer.
     real(kind=rp), allocatable :: recv_buf(:)
     !> Number of dofs to send to / receive from each peer.
     integer, allocatable :: send_len(:), recv_len(:)
     !> 0-based offsets into send_buf / recv_buf for each peer.
     integer, allocatable :: send_offset(:), recv_offset(:)
     !> Per-peer non-blocking MPI requests.
     type(MPI_Request), allocatable :: send_request(:), recv_request(:)
     !> Scratch arrays for MPI_Testsome on the recv side: indices of the
     !! requests completed by the call and their statuses.
     integer, allocatable :: recv_indices(:)
     type(MPI_Status), allocatable :: recv_statuses(:)
     integer :: ncompleted
   contains
     procedure, pass(this) :: init => gs_mpi_init
     procedure, pass(this) :: free => gs_mpi_free
     !> See gs_comm.f90 for more details on these routines
     procedure, pass(this) :: nbsend => gs_nbsend_mpi
     procedure, pass(this) :: nbrecv => gs_nbrecv_mpi
     procedure, pass(this) :: nbwait => gs_nbwait_mpi
  end type gs_mpi_t

contains

  !> Initialise MPI based communication method
  !! See gs_comm.f90 for details
  subroutine gs_mpi_init(this, send_pe, recv_pe)
    class(gs_mpi_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer :: i, nsend, nrecv, send_total, recv_total

    call this%init_order(send_pe, recv_pe)

    nsend = size(this%send_pe)
    nrecv = size(this%recv_pe)

    allocate(this%send_len(nsend), this%send_offset(nsend))
    allocate(this%send_request(nsend))

    allocate(this%recv_len(nrecv), this%recv_offset(nrecv))
    allocate(this%recv_request(nrecv))
    allocate(this%recv_indices(nrecv), this%recv_statuses(nrecv))

    send_total = 0
    do i = 1, nsend
       this%send_len(i) = this%send_dof(this%send_pe(i))%size()
       this%send_offset(i) = send_total
       send_total = send_total + this%send_len(i)
    end do
    allocate(this%send_buf(max(1, send_total)))

    recv_total = 0
    do i = 1, nrecv
       this%recv_len(i) = this%recv_dof(this%recv_pe(i))%size()
       this%recv_offset(i) = recv_total
       recv_total = recv_total + this%recv_len(i)
    end do
    allocate(this%recv_buf(max(1, recv_total)))

  end subroutine gs_mpi_init

  !> Deallocate MPI based communication method
  subroutine gs_mpi_free(this)
    class(gs_mpi_t), intent(inout) :: this

    if (allocated(this%send_buf)) then
       deallocate(this%send_buf)
    end if

    if (allocated(this%recv_buf)) then
       deallocate(this%recv_buf)
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

    if (allocated(this%send_request)) then
       deallocate(this%send_request)
    end if

    if (allocated(this%recv_request)) then
       deallocate(this%recv_request)
    end if

    if (allocated(this%recv_indices)) then
       deallocate(this%recv_indices)
    end if

    if (allocated(this%recv_statuses)) then
       deallocate(this%recv_statuses)
    end if

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_mpi_free

  !> Post non-blocking send operations
  subroutine gs_nbsend_mpi(this, u, n, tag, deps, strm)
    class(gs_mpi_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, ierr, dst, off, ndst

    ! Gather data from u into the per-peer slab of send_buf according
    ! to indices in send_dof. Slabs are contiguous so each MPI_Isend
    ! sends a single contiguous block.

    ! If the MPI library doesn't support MULTIPLE, use all threads to
    ! pack the send buffer. Otherwise, let each thread pack and send
    ! to a different neighbour
    if (NEKO_MPI_THREAD_PROVIDED .lt. MPI_THREAD_MULTIPLE) then
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          select type (sp => this%send_dof(dst)%data)
          type is (integer)
             !$omp do
             do j = 1, ndst
                this%send_buf(off + j) = u(sp(j))
             end do
             !$omp end do
          end select
          !$omp master
          call MPI_Isend(this%send_buf(off + 1), ndst, &
               MPI_REAL_PRECISION, dst, tag, &
               NEKO_COMM, this%send_request(i), ierr)
          !$omp end master
       end do
       !$omp barrier
    else
       !$omp do
       do i = 1, size(this%send_pe)
          dst = this%send_pe(i)
          off = this%send_offset(i)
          ndst = this%send_len(i)
          select type (sp => this%send_dof(dst)%data)
          type is (integer)
             !$omp simd
             do j = 1, ndst
                this%send_buf(off + j) = u(sp(j))
             end do
          end select
          call MPI_Isend(this%send_buf(off + 1), ndst, &
               MPI_REAL_PRECISION, dst, tag, &
               NEKO_COMM, this%send_request(i), ierr)
       end do
       !$omp end do
    end if

  end subroutine gs_nbsend_mpi

  !> Post non-blocking receive operations
  subroutine gs_nbrecv_mpi(this, tag)
    class(gs_mpi_t), intent(inout) :: this
    integer, intent(in) :: tag
    integer :: i, ierr, off, nsrc

    ! Issue recv requests, we will later check that these have finished
    ! in nbwait


    ! If the MPI library doesn't support MULTIPLE, the master thread
    ! will issue all Irecv's. Otherwise, threads will issue Irecv's
    ! concurrently.
    if (NEKO_MPI_THREAD_PROVIDED .lt. MPI_THREAD_MULTIPLE) then
       !$omp master
       do i = 1, size(this%recv_pe)
          off = this%recv_offset(i)
          nsrc = this%recv_len(i)
          call MPI_IRecv(this%recv_buf(off + 1), nsrc, &
               MPI_REAL_PRECISION, this%recv_pe(i), tag, &
               NEKO_COMM, this%recv_request(i), ierr)
       end do
       !$omp end master
       !$omp barrier
    else
       !$omp do
       do i = 1, size(this%recv_pe)
          off = this%recv_offset(i)
          nsrc = this%recv_len(i)
          call MPI_IRecv(this%recv_buf(off + 1), nsrc, &
               MPI_REAL_PRECISION, this%recv_pe(i), tag, &
               NEKO_COMM, this%recv_request(i), ierr)
       end do
       !$omp end do
    end if
  end subroutine gs_nbrecv_mpi

  !> Wait for non-blocking operations
  subroutine gs_nbwait_mpi(this, u, n, op, strm)
    class(gs_mpi_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, k, src, off, nsrc, ierr
    integer :: op
    integer, pointer :: sp(:)
    integer :: nreqs
    logical :: sends_done

    ! Poll for any subset of the outstanding recv requests to complete,
    ! reduce each completed slab into u, and repeat until all are done.
    nreqs = size(this%recv_pe)
    do while (nreqs .gt. 0)
       !$omp master
       call MPI_Testsome(size(this%recv_request), this%recv_request, &
            this%ncompleted, this%recv_indices, this%recv_statuses, ierr)
       !$omp end master
       !$omp barrier
       do k = 1, this%ncompleted
          i = this%recv_indices(k)
          !> @todo Check size etc against status
          src = this%recv_pe(i)
          off = this%recv_offset(i)
          nsrc = this%recv_len(i)
          select type (sp => this%recv_dof(src)%data)
          type is (integer)
             ! Do operation with data in buffer on dof specified by recv_dof
             select case (op)
             case (GS_OP_ADD)
                !OCL NORECURRENCE, NOVREC, NOALIAS
                !DIR$ CONCURRENT
                !GCC$ ivdep
                !NEC$ IVDEP
                !$omp do
                do j = 1, nsrc
                   u(sp(j)) = u(sp(j)) + this%recv_buf(off + j)
                end do
                !$omp end do
             case (GS_OP_MUL)
                !OCL NORECURRENCE, NOVREC, NOALIAS
                !DIR$ CONCURRENT
                !GCC$ ivdep
                !NEC$ IVDEP
                !$omp do
                do j = 1, nsrc
                   u(sp(j)) = u(sp(j)) * this%recv_buf(off + j)
                end do
                !$omp end do
             case (GS_OP_MIN)
                !OCL NORECURRENCE, NOVREC, NOALIAS
                !DIR$ CONCURRENT
                !GCC$ ivdep
                !NEC$ IVDEP
                !$omp do
                do j = 1, nsrc
                   u(sp(j)) = min(u(sp(j)), this%recv_buf(off + j))
                end do
                !$omp end do
             case (GS_OP_MAX)
                !OCL NORECURRENCE, NOVREC, NOALIAS
                !DIR$ CONCURRENT
                !GCC$ ivdep
                !NEC$ IVDEP
                !$omp do
                do j = 1, nsrc
                   u(sp(j)) = max(u(sp(j)), this%recv_buf(off + j))
                end do
                !$omp end do
             case default
                call neko_error("Unknown operation in gs_nbwait_mpi")
             end select
          end select
       end do
       nreqs = nreqs - this%ncompleted
       ! Synchronise before master can re-enter MPI_Testsome and overwrite
       ! this%ncompleted / this%recv_indices, which non-master threads are
       ! still reading in the unpack above.
       !$omp barrier
    end do
    !$omp master
    ! Finally, poll until all outstanding non-blocking sends have drained.
    if (size(this%send_request) .gt. 0) then
       sends_done = .false.
       do while (.not. sends_done)
          call MPI_Testall(size(this%send_request), this%send_request, &
               sends_done, MPI_STATUSES_IGNORE, ierr)
       end do
    end if
    !$omp end master
    !$omp barrier

  end subroutine gs_nbwait_mpi

end module gs_mpi
