! Copyright (c) 2022, The Neko Authors
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
!> Defines a gather-scatter communication method
module gs_comm
  use num_types
  use comm
  use stack
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none  

  integer, public, parameter :: GS_COMM_MPI = 1, GS_COMM_MPIGPU = 2

  !> Gather-scatter communication method
  type, public, abstract :: gs_comm_t
     type(stack_i4_t), allocatable :: send_dof(:)     !< Send dof to shared-gs
     type(stack_i4_t), allocatable :: recv_dof(:)     !< Recv dof to shared-gs
     integer, allocatable :: send_pe(:)               !< Send order
     integer, allocatable :: recv_pe(:)               !< Recv order
   contains
     procedure(gs_comm_init), pass(this), deferred :: init
     procedure(gs_comm_free), pass(this), deferred :: free
     procedure(gs_nbsend), pass(this), deferred :: nbsend
     procedure(gs_nbrecv), pass(this), deferred :: nbrecv
     procedure(gs_nbwait), pass(this), deferred :: nbwait
     procedure, pass(this) :: init_dofs
     procedure, pass(this) :: free_dofs
     procedure, pass(this) :: init_order
     procedure, pass(this) :: free_order
  end type gs_comm_t

  !> Abstract interface for initialising a Gather-scatter communication method
  abstract interface
     subroutine gs_comm_init(this, send_pe, recv_pe)
       import gs_comm_t
       import stack_i4_t       
       class(gs_comm_t), intent(inout) :: this
       type(stack_i4_t), intent(inout) :: send_pe
       type(stack_i4_t), intent(inout) :: recv_pe       
     end subroutine gs_comm_init
  end interface

  !> Abstract interface for deallocating a Gather-scatter communication method
  abstract interface
     subroutine gs_comm_free(this)
       import gs_comm_t
       class(gs_comm_t), intent(inout) :: this
     end subroutine gs_comm_free
  end interface

  !> Abstract interface for initiating non-blocking send operations
  abstract interface
     subroutine gs_nbsend(this, u, n, deps)
       import gs_comm_t
       import stack_i4_t
       import c_ptr
       import rp
       class(gs_comm_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: u
       type(c_ptr), intent(inout) :: deps
     end subroutine gs_nbsend
  end interface

  !> Abstract interface for initiating non-blocking receive operations
  abstract interface
     subroutine gs_nbrecv(this)
       import gs_comm_t
       class(gs_comm_t), intent(inout) :: this
     end subroutine gs_nbrecv
  end interface

  !> Abstract interface for watining on non-blocking operations
  abstract interface
     subroutine gs_nbwait(this, u, n, op)
       import gs_comm_t
       import stack_i4_t
       import c_ptr
       import rp
       class(gs_comm_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: u
       integer :: op
     end subroutine gs_nbwait
  end interface

contains

  subroutine init_dofs(this)
    class(gs_comm_t), intent(inout) :: this
    integer :: i

    call this%free_dofs()

    allocate(this%send_dof(0:pe_size-1))
    allocate(this%recv_dof(0:pe_size-1))

    do i = 0, pe_size -1
       call this%send_dof(i)%init()
       call this%recv_dof(i)%init()
    end do
    
  end subroutine init_dofs

  subroutine free_dofs(this)
    class(gs_comm_t), intent(inout) :: this
    integer :: i

    if (allocated(this%send_dof)) then
       do i = 0, pe_size - 1
          call this%send_dof(i)%free()
       end do
       deallocate(this%send_dof)
    end if

    if (allocated(this%recv_dof)) then
       do i = 0, pe_size - 1
          call this%recv_dof(i)%free()
       end do
       deallocate(this%recv_dof)
    end if
    
  end subroutine free_dofs

  subroutine init_order(this, send_pe, recv_pe)
    class(gs_comm_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer, pointer :: sp(:)
    integer :: i
    
    allocate(this%send_pe(send_pe%size()))

    sp => send_pe%array()
    do i = 1, send_pe%size()
       this%send_pe(i) = sp(i)
    end do

    allocate(this%recv_pe(recv_pe%size()))

    sp => recv_pe%array()
    do i = 1, recv_pe%size()
       this%recv_pe(i) = sp(i)
    end do
    
  end subroutine init_order

  subroutine free_order(this)
    class(gs_comm_t), intent(inout) :: this

    if (allocated(this%send_pe)) then
       deallocate(this%send_pe)
    end if

    if (allocated(this%recv_pe)) then
       deallocate(this%recv_pe)
    end if
    
  end subroutine free_order
     
end module gs_comm
