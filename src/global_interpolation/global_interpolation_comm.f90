! Copyright (c) 2020-2025, The Neko Authors
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
!> Defines global interpolation communication
!! Based on the MPI based gather-scatter kernel
module glb_intrp_comm
  use num_types, only : rp
  use comm, only : pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use stack, only : stack_i4_t
  use, intrinsic :: iso_c_binding
  use mpi_f08, only : MPI_Test, MPI_STATUS_IGNORE, MPI_Status, &
       MPI_Request, MPI_Isend, MPI_IRecv, MPI_Comm
  implicit none
  private


  !> MPI buffer for non-blocking operations
  type, private :: glb_intrp_comm_mpi_t
     !> status of this operation
     type(MPI_Status) :: status
     !> unique request for this buffer
     type(MPI_Request) :: request
     !> flag = false when operation is in process
     !! Set to true when operation has succeeded
     logical :: flag
     !> Buffer with data to send/recieve
     real(kind=rp), allocatable :: data(:)
  end type glb_intrp_comm_mpi_t

  !> Global interpolation communication method
  type, public :: glb_intrp_comm_t
     !> A list of stacks of dof indices local to this process to send to rank_i
     type(stack_i4_t), allocatable :: send_dof(:)
     !> recv_dof(rank_i) is a stack of dof indices local to this process to
     !! receive from rank_i. size(recv_dof) == pe_size
     type(stack_i4_t), allocatable :: recv_dof(:)
     !> Size of communicator
     integer :: pe_size
     !> Array of ranks that this process should send to
     !! @note: this will usually be fewer than the total number of ranks
     !! size(send_pe) <= pe_size
     integer, allocatable :: send_pe(:)
     !> array of ranks that this process will receive messages from
     integer, allocatable :: recv_pe(:)
     !> Comm. buffers for send operations
     type(glb_intrp_comm_mpi_t), allocatable :: send_buf(:)
     !> Comm. buffers for recv operations
     type(glb_intrp_comm_mpi_t), allocatable :: recv_buf(:)
     !> Communicator
     type(MPI_Comm) :: comm
   contains
     procedure, pass(this) :: init => glb_intrp_comm_init
     procedure, pass(this) :: free => glb_intrp_comm_free
     procedure, pass(this) :: init_dofs => glb_intrp_comm_init_dofs
     procedure, pass(this) :: free_dofs => glb_intrp_comm_free_dofs
     procedure, pass(this) :: init_order => glb_intrp_comm_init_order
     procedure, pass(this) :: free_order => glb_intrp_comm_free_order
     procedure, pass(this) :: nbwait_no_op => glb_intrp_comm_nbwait_no_op
     procedure, pass(this) :: sendrecv => glb_intrp_comm_sendrecv
  end type glb_intrp_comm_t

contains

  !> Initialise MPI based communication method
  subroutine glb_intrp_comm_init(this, send_pe, recv_pe, comm)
    class(glb_intrp_comm_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    type(MPI_Comm), intent(inout), optional :: comm
    integer, pointer :: sp(:), rp(:)
    integer :: i
    if (present(comm)) then
       this%comm = comm
    else
       this%comm = NEKO_COMM
    end if

    call this%init_order(send_pe, recv_pe)

    allocate(this%send_buf(send_pe%size()))

    sp => send_pe%array()
    do i = 1, send_pe%size()
       allocate(this%send_buf(i)%data(this%send_dof(sp(i))%size()))
    end do

    allocate(this%recv_buf(recv_pe%size()))

    rp => recv_pe%array()
    do i = 1, recv_pe%size()
       allocate(this%recv_buf(i)%data(this%recv_dof(rp(i))%size()))
    end do

  end subroutine glb_intrp_comm_init

  !> Deallocate MPI based communication method
  subroutine glb_intrp_comm_free(this)
    class(glb_intrp_comm_t), intent(inout) :: this
    integer :: i

    if (allocated(this%send_buf)) then
       do i = 1, size(this%send_buf)
          if (allocated(this%send_buf(i)%data)) then
             deallocate(this%send_buf(i)%data)
          end if
       end do
       deallocate(this%send_buf)
    end if

    if (allocated(this%recv_buf)) then
       do i = 1, size(this%recv_buf)
          if (allocated(this%recv_buf(i)%data)) then
             deallocate(this%recv_buf(i)%data)
          end if
       end do
       deallocate(this%recv_buf)
    end if

    call this%free_order()
    call this%free_dofs()

  end subroutine glb_intrp_comm_free

  !Initalize stacks for each rank of dof indices to send/recv
  subroutine glb_intrp_comm_init_dofs(this, comm_size)
    class(glb_intrp_comm_t), intent(inout) :: this
    integer, optional, intent(in) :: comm_size
    integer :: i

    if (present(comm_size)) then
       this%pe_size = comm_size
    else
       this%pe_size = pe_size
    end if

    call this%free_dofs()

    allocate(this%send_dof(0:this%pe_size-1))
    allocate(this%recv_dof(0:this%pe_size-1))

    do i = 0, this%pe_size -1
       call this%send_dof(i)%init()
       call this%recv_dof(i)%init()
    end do

  end subroutine glb_intrp_comm_init_dofs

  subroutine glb_intrp_comm_free_dofs(this)
    class(glb_intrp_comm_t), intent(inout) :: this
    integer :: i

    if (allocated(this%send_dof)) then
       do i = 0, this%pe_size - 1
          call this%send_dof(i)%free()
       end do
       deallocate(this%send_dof)
    end if

    if (allocated(this%recv_dof)) then
       do i = 0, this%pe_size - 1
          call this%recv_dof(i)%free()
       end do
       deallocate(this%recv_dof)
    end if

  end subroutine glb_intrp_comm_free_dofs

  !>Obtains which ranks to send and receive data from
  !! @param send_pe, only contains rank ids this process should send to
  !! @param recv_pe, only the ranks this process should receive from
  subroutine glb_intrp_comm_init_order(this, send_pe, recv_pe)
    class(glb_intrp_comm_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer, contiguous, pointer :: sp(:), rp(:)
    integer :: i

    allocate(this%send_pe(send_pe%size()))

    sp => send_pe%array()
    do i = 1, send_pe%size()
       this%send_pe(i) = sp(i)
    end do

    allocate(this%recv_pe(recv_pe%size()))

    rp => recv_pe%array()
    do i = 1, recv_pe%size()
       this%recv_pe(i) = rp(i)
    end do

  end subroutine glb_intrp_comm_init_order

  subroutine glb_intrp_comm_free_order(this)
    class(glb_intrp_comm_t), intent(inout) :: this

    if (allocated(this%send_pe)) then
       deallocate(this%send_pe)
    end if

    if (allocated(this%recv_pe)) then
       deallocate(this%recv_pe)
    end if

  end subroutine glb_intrp_comm_free_order

  !> Non-blocking sendrecv
  subroutine glb_intrp_comm_sendrecv(this, send, recv, n_send, n_recv)
    class(glb_intrp_comm_t), intent(inout) :: this
    integer, intent(in) :: n_send, n_recv
    real(kind=rp), dimension(n_send), intent(inout) :: send
    real(kind=rp), dimension(n_recv), intent(inout) :: recv
    type(c_ptr) :: null_ptr = c_null_ptr
    integer :: i, j, ierr, src, dst, thrdid
    integer, pointer :: sp(:)
    integer :: nreqs

    thrdid = 0
    !$ thrdid = omp_get_thread_num()

    !
    ! Issue non-blocking receives
    !
    do i = 1, size(this%recv_pe)
       ! We should not need this extra associate block, ant it works
       ! great without it for GNU, Intel, NEC and Cray, but throws an
       ! ICE with NAG.
       ! Issue recv requests, we will later check that these have finished
       ! in nbwait
       associate(recv_data => this%recv_buf(i)%data)
         call MPI_IRecv(recv_data, size(recv_data), &
              MPI_REAL_PRECISION, this%recv_pe(i), thrdid, &
              this%comm, this%recv_buf(i)%request, ierr)
       end associate
       this%recv_buf(i)%flag = .false.
    end do

    !
    ! Issue non-blocking sends
    !
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       ! Gather data from u into buffers according to indices in send_dof
       ! We want to send contigous data to each process in send_pe
       sp => this%send_dof(dst)%array()
       do concurrent (j = 1:this%send_dof(dst)%size())
          this%send_buf(i)%data(j) = send(sp(j))
       end do
       ! We should not need this extra associate block, ant it works
       ! great without it for GNU, Intel, NEC and Cray, but throws an
       ! ICE with NAG.
       associate(send_data => this%send_buf(i)%data)
         call MPI_Isend(send_data, this%send_dof(dst)%size(), &
              MPI_REAL_PRECISION, this%send_pe(i), thrdid, &
              this%comm, this%send_buf(i)%request, ierr)
       end associate
       this%send_buf(i)%flag = .false.
    end do

    !
    ! Wait for non-blocking operations
    !

    nreqs = size(this%recv_pe)

    do while (nreqs .gt. 0)
       do i = 1, size(this%recv_pe)
          if (.not. this%recv_buf(i)%flag) then
             ! Check if we have recieved the data we want
             call MPI_Test(this%recv_buf(i)%request, this%recv_buf(i)%flag, &
                  this%recv_buf(i)%status, ierr)
             ! If it has been received
             if (this%recv_buf(i)%flag) then
                ! One more request has been succesful
                nreqs = nreqs - 1
                !> @todo Check size etc against status
                src = this%recv_pe(i)
                sp => this%recv_dof(src)%array()
                !NEC$ IVDEP
                do concurrent (j = 1:this%recv_dof(src)%size())
                   recv(sp(j)) = this%recv_buf(i)%data(j)
                end do
             end if
          end if
       end do
    end do
    ! Finally, check that the non-blocking sends this rank have issued have also
    ! completed successfully

    nreqs = size(this%send_pe)
    do while (nreqs .gt. 0)
       do i = 1, size(this%send_pe)
          if (.not. this%send_buf(i)%flag) then
             call MPI_Test(this%send_buf(i)%request, this%send_buf(i)%flag, &
                  MPI_STATUS_IGNORE, ierr)
             if (this%send_buf(i)%flag) nreqs = nreqs - 1
          end if
       end do
    end do

  end subroutine glb_intrp_comm_sendrecv

  !> Wait for non-blocking operations
  subroutine glb_intrp_comm_nbwait_no_op(this)
    class(glb_intrp_comm_t), intent(inout) :: this
    integer :: i, j, src, ierr
    integer , pointer :: sp(:)
    integer :: nreqs

    nreqs = size(this%recv_pe)

    do while (nreqs .gt. 0)
       do i = 1, size(this%recv_pe)
          if (.not. this%recv_buf(i)%flag) then
             ! Check if we have recieved the data we want
             call MPI_Test(this%recv_buf(i)%request, this%recv_buf(i)%flag, &
                  this%recv_buf(i)%status, ierr)
             ! If it has been received
             if (this%recv_buf(i)%flag) then
                ! One more request has been succesful
                nreqs = nreqs - 1
             end if
          end if
       end do
    end do
    ! Finally, check that the non-blocking sends this rank have issued have also
    ! completed successfully

    nreqs = size(this%send_pe)
    do while (nreqs .gt. 0)
       do i = 1, size(this%send_pe)
          if (.not. this%send_buf(i)%flag) then
             call MPI_Test(this%send_buf(i)%request, this%send_buf(i)%flag, &
                  MPI_STATUS_IGNORE, ierr)
             if (this%send_buf(i)%flag) nreqs = nreqs - 1
          end if
       end do
    end do

  end subroutine glb_intrp_comm_nbwait_no_op


end module glb_intrp_comm
