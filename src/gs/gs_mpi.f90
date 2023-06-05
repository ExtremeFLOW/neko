! Copyright (c) 2020-2022, The Neko Authors
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
  use neko_config
  use num_types
  use gs_comm
  use gs_ops
  use stack
  use mpi_f08
  use comm
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none

  !> MPI buffer for non-blocking operations
  type, private :: gs_comm_mpi_t
     type(MPI_Status) :: status
     type(MPI_Request) :: request
     logical :: flag
     real(kind=rp), allocatable :: data(:)
  end type  gs_comm_mpi_t

  !> Gather-scatter communication using MPI
  type, extends(gs_comm_t) :: gs_mpi_t
     type(gs_comm_mpi_t), allocatable :: send_buf(:)     !< Comm. buffers
     type(gs_comm_mpi_t), allocatable :: recv_buf(:)     !< Comm. buffers
   contains
     procedure, pass(this) :: init => gs_mpi_init
     procedure, pass(this) :: free => gs_mpi_free
     procedure, pass(this) :: nbsend => gs_nbsend_mpi
     procedure, pass(this) :: nbrecv => gs_nbrecv_mpi
     procedure, pass(this) :: nbwait => gs_nbwait_mpi
  end type gs_mpi_t

contains

  !> Initialise MPI based communication method
  subroutine gs_mpi_init(this, send_pe, recv_pe)
    class(gs_mpi_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe       
    integer, pointer :: sp(:), rp(:)
    integer :: i
    
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

  end subroutine gs_mpi_init

  !> Deallocate MPI based communication method
  subroutine gs_mpi_free(this)
    class(gs_mpi_t), intent(inout) :: this
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
    
  end subroutine gs_mpi_free

  !> Post non-blocking send operations
  subroutine gs_nbsend_mpi(this, u, n, deps)
    class(gs_mpi_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: deps
    integer ::  i, j, ierr, dst
    integer , pointer :: sp(:)

    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       sp => this%send_dof(dst)%array()
       do j = 1, this%send_dof(dst)%size()
          this%send_buf(i)%data(j) = u(sp(j))
       end do
       ! We should not need this extra associate block, ant it works
       ! great without it for GNU, Intel, NEC and Cray, but throws an
       ! ICE with NAG.
       associate(send_data => this%send_buf(i)%data)
         call MPI_Isend(send_data, size(send_data), &
              MPI_REAL_PRECISION, this%send_pe(i), 0, &
              NEKO_COMM, this%send_buf(i)%request, ierr)
       end associate
       this%send_buf(i)%flag = .false.
    end do
  end subroutine gs_nbsend_mpi

  !> Post non-blocking receive operations
  subroutine gs_nbrecv_mpi(this)
    class(gs_mpi_t), intent(inout) :: this
    integer :: i, ierr

    do i = 1, size(this%recv_pe)
       ! We should not need this extra associate block, ant it works
       ! great without it for GNU, Intel, NEC and Cray, but throws an
       ! ICE with NAG.
       associate(recv_data => this%recv_buf(i)%data)
         call MPI_IRecv(recv_data, size(recv_data), &
              MPI_REAL_PRECISION, this%recv_pe(i), 0, &
              NEKO_COMM, this%recv_buf(i)%request, ierr)
       end associate
       this%recv_buf(i)%flag = .false.
    end do
    
  end subroutine gs_nbrecv_mpi

  !> Wait for non-blocking operations
  subroutine gs_nbwait_mpi(this, u, n, op)
    class(gs_mpi_t), intent(inout) :: this
    integer, intent(in) :: n    
    real(kind=rp), dimension(n), intent(inout) :: u
    integer :: i, j, src, ierr
    integer :: op
    integer , pointer :: sp(:)
    integer :: nreqs

    nreqs = size(this%recv_pe)

    do while (nreqs .gt. 0) 
       do i = 1, size(this%recv_pe)
          if (.not. this%recv_buf(i)%flag) then
             call MPI_Test(this%recv_buf(i)%request, this%recv_buf(i)%flag, &
                  this%recv_buf(i)%status, ierr)
             if (this%recv_buf(i)%flag) then
                nreqs = nreqs - 1
                !> @todo Check size etc against status
                src = this%recv_pe(i)
                sp => this%recv_dof(src)%array()
                select case(op)
                case (GS_OP_ADD)
                   !NEC$ IVDEP
                   do j = 1, this%send_dof(src)%size()
                      u(sp(j)) = u(sp(j)) + this%recv_buf(i)%data(j)
                   end do
                case (GS_OP_MUL)
                   !NEC$ IVDEP
                   do j = 1, this%send_dof(src)%size()
                      u(sp(j)) = u(sp(j)) * this%recv_buf(i)%data(j)
                   end do
                case (GS_OP_MIN)
                   !NEC$ IVDEP
                   do j = 1, this%send_dof(src)%size()
                      u(sp(j)) = min(u(sp(j)), this%recv_buf(i)%data(j))
                   end do
                case (GS_OP_MAX)
                   !NEC$ IVDEP
                   do j = 1, this%send_dof(src)%size()
                      u(sp(j)) = max(u(sp(j)), this%recv_buf(i)%data(j))
                   end do
                end select
             end if
          end if
       end do
    end do

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

  end subroutine gs_nbwait_mpi
  
end module gs_mpi
