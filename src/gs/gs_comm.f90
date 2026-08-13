! Copyright (c) 2022-2026, The Neko Authors
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
  use num_types, only : rp
  use comm, only : pe_size
  use stack, only : stack_i4_t
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer, public, parameter :: GS_COMM_MPI = 1, GS_COMM_MPIGPU = 2, &
       GS_COMM_NCCL = 3, GS_COMM_NVSHMEM = 4, GS_COMM_OPENSHMEM = 5, &
       GS_COMM_CAF = 6, GS_COMM_NEIGHBOUR = 7, GS_COMM_UTOFU = 8

  !> Maximum number of components handled by the fused vector (multi-component)
  !! halo exchange used by gs_op_r3. Sizes the backend vector buffers.
  integer, public, parameter :: GS_VEC_NC = 3

  !> Gather-scatter communication method
  type, public, abstract :: gs_comm_t
     !> A list of stacks of dof indices local to this process to send to rank_i
     type(stack_i4_t), allocatable :: send_dof(:)
     !> recv_dof(rank_i) is a stack of dof indices local to this process to
     !! receive from rank_i. size(recv_dof) == pe_size
     type(stack_i4_t), allocatable :: recv_dof(:)
     !> Array of ranks that this process should send to
     !! @note: this will usually be fewer than the total number of ranks
     !! size(send_pe) <= pe_size
     integer, allocatable :: send_pe(:)
     !> array of ranks that this process will receive messages from
     integer, allocatable :: recv_pe(:)
     !> Whether this backend implements the fused vector (multi-component)
     !! halo exchange (nbsend_vec/nbrecv_vec/nbwait_vec). When .false., the
     !! gs_op_r3 caller falls back to nc independent scalar exchanges.
     logical :: vec_supported = .false.
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
     procedure, pass(this) :: take_schedule
     procedure, pass(this) :: init_schedule
     !> Fused vector halo exchange. Default implementations abort; backends
     !! that set vec_supported = .true. override them.
     procedure, pass(this) :: nbsend_vec => gs_nbsend_vec
     procedure, pass(this) :: nbrecv_vec => gs_nbrecv_vec
     procedure, pass(this) :: nbwait_vec => gs_nbwait_vec
  end type gs_comm_t

  !> Abstract interface for initializing a Gather-scatter communication method
  !! @param send_pe, stack of ranks this process will send messages to
  !! @param recv_pe, stack of ranks this process will receive messages from
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
  !! Sends the values in u(send_dof(send_pe(i))) to each rank send_pe(i) for all
  !! ranks in send_pe
  !! @param n, length of u (redundant)
  !! @param deps, gather_event (for device aware mpi)
  !! @param strm, device stream to execute operation on
  abstract interface
     subroutine gs_nbsend(this, u, n, tag, deps, strm)
       import gs_comm_t
       import stack_i4_t
       import c_ptr
       import rp
       class(gs_comm_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: u
       integer, intent(in) :: tag
       type(c_ptr), intent(inout) :: deps
       type(c_ptr), intent(inout) :: strm
     end subroutine gs_nbsend
  end interface


  !> Abstract interface for initiating non-blocking recieve operations
  !! Posts non-blocking recieve of values and puts the values into buffers
  abstract interface
     subroutine gs_nbrecv(this, tag)
       import gs_comm_t
       class(gs_comm_t), intent(inout) :: this
       integer, intent(in) :: tag
     end subroutine gs_nbrecv
  end interface

  !> Abstract interface for waiting on non-blocking operations
  !! Waits and checks that data is in buffers and unpacks buffers
  !! into correct location in u
  !! u(recv_dof(recv_pe(i))) = gs_op(recieve_buffers(recv_pe) for this dof)
  !! @param u, data to store operation into
  !! @param n, length of u (redundant)
  !! @param op, gather scatter operation to carry out
  !! @param strm, device stream to execute this operation on
  abstract interface
     subroutine gs_nbwait(this, u, n, op, strm)
       import gs_comm_t
       import stack_i4_t
       import c_ptr
       import rp
       class(gs_comm_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: u
       integer :: op
       type(c_ptr), intent(inout) :: strm
     end subroutine gs_nbwait
  end interface

  public :: gs_comm_init, gs_comm_free, gs_nbsend, gs_nbrecv, gs_nbwait
contains
  !Initalize stacks for each rank of dof indices to send/recv
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

  !>Obtains which ranks to send and receive data from
  !! @param send_pe, only contains rank ids this porcesss should send to
  !! @param recv_pe, only the ranks this process should receive from
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

  !> Take over the gather-scatter schedule (dof lists and peer order) of
  !! @a src, avoiding a second (expensive) pass over the connectivity. The
  !! data is moved rather than copied, so @a src is left without a schedule
  !! and must not be used for communication afterwards (it can still be
  !! freed). No communication resources are set up here; complete the
  !! handover with @a init_schedule once @a src has been freed, so that the
  !! two backends never hold their resources at the same time.
  !! @param src the backend to take the schedule from, left without one
  subroutine take_schedule(this, src)
    class(gs_comm_t), intent(inout) :: this
    class(gs_comm_t), intent(inout) :: src

    if (.not. allocated(src%send_dof) .or. .not. allocated(src%recv_dof) .or. &
         .not. allocated(src%send_pe) .or. .not. allocated(src%recv_pe)) then
       call neko_error('Gather-scatter comm. method has no schedule')
    end if

    call this%free_dofs()
    call this%free_order()

    call move_alloc(src%send_dof, this%send_dof)
    call move_alloc(src%recv_dof, this%recv_dof)
    call move_alloc(src%send_pe, this%send_pe)
    call move_alloc(src%recv_pe, this%recv_pe)

  end subroutine take_schedule

  !> Set up this communication method for the schedule taken over by
  !! @a take_schedule. Collective, as @a init is.
  subroutine init_schedule(this)
    class(gs_comm_t), intent(inout) :: this
    type(stack_i4_t) :: send_pe, recv_pe
    integer, allocatable :: sp(:), rp(:)
    integer :: i, pe

    if (.not. allocated(this%send_pe) .or. .not. allocated(this%recv_pe)) then
       call neko_error('Gather-scatter comm. method has no schedule')
    end if

    ! init_order allocates send_pe/recv_pe from the stacks, so hand the
    ! peer lists back as stacks and leave the arrays deallocated
    call move_alloc(this%send_pe, sp)
    call move_alloc(this%recv_pe, rp)

    call send_pe%init(max(size(sp), 1))
    do i = 1, size(sp)
       pe = sp(i)
       call send_pe%push(pe)
    end do

    call recv_pe%init(max(size(rp), 1))
    do i = 1, size(rp)
       pe = rp(i)
       call recv_pe%push(pe)
    end do

    deallocate(sp, rp)

    call this%init(send_pe, recv_pe)

    call send_pe%free()
    call recv_pe%free()

  end subroutine init_schedule

  !> Default fused vector send. Abort unless a backend overrides it.
  !! @param u compact shared buffer, component-outer: u((c-1)*n + idx)
  !! @param n number of shared dofs (per component)
  !! @param nc number of components
  subroutine gs_nbsend_vec(this, u, n, nc, tag, deps, strm)
    class(gs_comm_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    call neko_error('Vector gather-scatter not supported by this comm backend')
  end subroutine gs_nbsend_vec

  !> Default fused vector receive. Abort unless a backend overrides it.
  subroutine gs_nbrecv_vec(this, tag, nc)
    class(gs_comm_t), intent(inout) :: this
    integer, intent(in) :: tag, nc
    call neko_error('Vector gather-scatter not supported by this comm backend')
  end subroutine gs_nbrecv_vec

  !> Default fused vector wait/reduce. Abort unless a backend overrides it.
  subroutine gs_nbwait_vec(this, u, n, nc, op, strm)
    class(gs_comm_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer :: op
    type(c_ptr), intent(inout) :: strm
    call neko_error('Vector gather-scatter not supported by this comm backend')
  end subroutine gs_nbwait_vec

end module gs_comm
