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
!> Defines GPU aware MPI gather-scatter communication
module gs_device_mpi
  use neko_config
  use num_types
  use gs_comm
  use gs_ops
  use stack
  use mpi_f08
  use comm
  use device
  implicit none

  !> MPI buffer for non-blocking operations
  type, private :: gs_device_mpi_buf_t
     integer, allocatable :: ndofs(:)
     integer, allocatable :: offset(:)
     integer :: total
     type(c_ptr) :: reqs = C_NULL_PTR
     type(c_ptr) :: buf_d = C_NULL_PTR
     type(c_ptr) :: dof_d = C_NULL_PTR
   contains
     procedure, pass(this) :: init => gs_device_mpi_buf_init
     procedure, pass(this) :: free => gs_device_mpi_buf_free
  end type gs_device_mpi_buf_t

  !> Gather-scatter communication using device MPI.
  !! The arrays are indexed per PE like @a send_pe and @ recv_pe.
  type, extends(gs_comm_t) :: gs_device_mpi_t
     type(gs_device_mpi_buf_t) :: send_buf
     type(gs_device_mpi_buf_t) :: recv_buf
   contains
     procedure, pass(this) :: init => gs_device_mpi_init
     procedure, pass(this) :: free => gs_device_mpi_free
     procedure, pass(this) :: nbsend => gs_device_mpi_nbsend
     procedure, pass(this) :: nbrecv => gs_device_mpi_nbrecv
     procedure, pass(this) :: nbwait => gs_device_mpi_nbwait
  end type gs_device_mpi_t

  interface
    subroutine hip_gs_pack(u_d, buf_d, dof_d, n) &
          bind(c, name='hip_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n
       type(c_ptr), value :: u_d, buf_d, dof_d
     end subroutine hip_gs_pack
  end interface

  interface
    subroutine cuda_gs_pack(u_d, buf_d, dof_d, n) &
          bind(c, name='cuda_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n
       type(c_ptr), value :: u_d, buf_d, dof_d
     end subroutine cuda_gs_pack
  end interface

  interface
    subroutine hip_gs_unpack(u_d, op, buf_d, dof_d, n) &
          bind(c, name='hip_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, n
       type(c_ptr), value :: u_d, buf_d, dof_d
     end subroutine hip_gs_unpack
  end interface

  interface
    subroutine cuda_gs_unpack(u_d, op, buf_d, dof_d, n) &
          bind(c, name='cuda_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, n
       type(c_ptr), value :: u_d, buf_d, dof_d
     end subroutine cuda_gs_unpack
  end interface

  interface
    subroutine device_mpi_init_reqs(n, reqs) &
          bind(c, name='device_mpi_init_reqs')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n
       type(c_ptr) :: reqs
     end subroutine device_mpi_init_reqs
  end interface

  interface
    subroutine device_mpi_free_reqs(n, reqs) &
          bind(c, name='device_mpi_free_reqs')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n
       type(c_ptr) :: reqs
     end subroutine device_mpi_free_reqs
  end interface

  interface
    subroutine device_mpi_isend(buf_d, offset, nbytes, rank, reqs, i) &
          bind(c, name='device_mpi_isend')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, nbytes, rank, i
       type(c_ptr), value :: buf_d, reqs
     end subroutine device_mpi_isend
  end interface

  interface
    subroutine device_mpi_irecv(buf_d, offset, nbytes, rank, reqs, i) &
          bind(c, name='device_mpi_irecv')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, nbytes, rank, i
       type(c_ptr), value :: buf_d, reqs
     end subroutine device_mpi_irecv
  end interface

  interface
    integer(c_int) function device_mpi_test(req) &
          bind(c, name='device_mpi_test')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: req
     end function device_mpi_test
  end interface

  interface
    subroutine device_mpi_waitall(n, reqs) &
          bind(c, name='device_mpi_waitall')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n
       type(c_ptr), value :: reqs
    end subroutine device_mpi_waitall
  end interface

contains

  subroutine gs_device_mpi_buf_init(this, pe_order, dof_stack, mark_dupes)
    class(gs_device_mpi_buf_t), intent(inout) :: this
    integer, allocatable, intent(inout) :: pe_order(:)
    type(stack_i4_t), allocatable, intent(inout) :: dof_stack(:)
    logical, intent(in) :: mark_dupes
    integer, allocatable :: dofs(:)
    integer :: i, j, total
    integer(c_size_t) :: sz
    type(htable_i4_t) :: doftable
    integer :: dupe, marked, k

    call device_mpi_init_reqs(size(pe_order), this%reqs)

    allocate(this%ndofs(size(pe_order)))
    allocate(this%offset(size(pe_order)))

    total = 0
    do i = 1, size(pe_order)
       this%ndofs(i) = dof_stack(pe_order(i))%size()
       this%offset(i) = total
       total = total + this%ndofs(i)
    end do

    this%total = total

    sz = rp * total
    call device_alloc(this%buf_d, sz)

    sz = 4 * total
    call device_alloc(this%dof_d, sz)

    if (mark_dupes) call doftable%init(2*total)
    allocate(dofs(total))

    marked = 0
    do i = 1, size(pe_order)
       ! %array() breaks on cray
       select type (arr => dof_stack(pe_order(i))%data)
       type is (integer)
         do j = 1, this%ndofs(i)
            k = this%offset(i) + j
            if (mark_dupes) then
               if (doftable%get(arr(j), dupe) .eq. 0) then
                  if (dofs(dupe) .gt. 0) then
                     dofs(dupe) = -dofs(dupe)
                     marked = marked + 1
                  end if
                  dofs(k) = -arr(j)
                  marked = marked + 1
               else
                  call doftable%set(arr(j), k)
                  dofs(k) = arr(j)
               end if
            else
               dofs(k) = arr(j)
            end if
         end do
       end select
    end do

    call device_memcpy(dofs, this%dof_d, total, HOST_TO_DEVICE)

    deallocate(dofs)
    call doftable%free()

  end subroutine gs_device_mpi_buf_init

  subroutine gs_device_mpi_buf_free(this)
    class(gs_device_mpi_buf_t), intent(inout) :: this
  end subroutine

  !> Initialise MPI based communication method
  subroutine gs_device_mpi_init(this, send_pe, recv_pe)
    class(gs_device_mpi_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe

    call this%init_order(send_pe, recv_pe)

    call this%send_buf%init(this%send_pe, this%send_dof, .false.)
    call this%recv_buf%init(this%recv_pe, this%recv_dof, .true.)

  end subroutine gs_device_mpi_init

  !> Deallocate MPI based communication method
  subroutine gs_device_mpi_free(this)
    class(gs_device_mpi_t), intent(inout) :: this
    integer :: i

    call this%send_buf%free()
    call this%recv_buf%free()

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_device_mpi_free

  !> Post non-blocking send operations
  subroutine gs_device_mpi_nbsend(this, u, n)
    class(gs_device_mpi_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer ::  i
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u, n)

#ifdef HAVE_HIP
    call hip_gs_pack(u_d, &
                     this%send_buf%buf_d, &
                     this%send_buf%dof_d, &
                     this%send_buf%total)
#elif HAVE_CUDA
    call cuda_gs_pack(u_d, &
                      this%send_buf%buf_d, &
                      this%send_buf%dof_d, &
                      this%send_buf%total)
#else
    call neko_error('gs_device_mpi: no backend')
#endif

    call device_sync()

    do i = 1, size(this%send_pe)
       call device_mpi_isend(this%send_buf%buf_d, rp*this%send_buf%offset(i), &
                             rp*this%send_buf%ndofs(i), this%send_pe(i), &
                             this%send_buf%reqs, i)
    end do

  end subroutine gs_device_mpi_nbsend

  !> Post non-blocking receive operations
  subroutine gs_device_mpi_nbrecv(this)
    class(gs_device_mpi_t), intent(inout) :: this
    integer :: i

    do i = 1, size(this%recv_pe)
       call device_mpi_irecv(this%recv_buf%buf_d, rp*this%recv_buf%offset(i), &
                             rp*this%recv_buf%ndofs(i), this%recv_pe(i), &
                             this%recv_buf%reqs, i)
    end do

  end subroutine gs_device_mpi_nbrecv

  !> Wait for non-blocking operations
  subroutine gs_device_mpi_nbwait(this, u, n, op)
    class(gs_device_mpi_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer :: i, flag
    integer :: op
    integer :: nreqs
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u, n)

    call device_mpi_waitall(size(this%recv_pe), this%recv_buf%reqs)

#ifdef HAVE_HIP
    call hip_gs_unpack(u_d, op, &
                       this%recv_buf%buf_d, &
                       this%recv_buf%dof_d, &
                       this%recv_buf%total)
#elif HAVE_CUDA
    call cuda_gs_unpack(u_d, op, &
                        this%recv_buf%buf_d, &
                        this%recv_buf%dof_d, &
                        this%recv_buf%total)
#else
                call neko_error('gs_device_mpi: no backend')
#endif

    call device_mpi_waitall(size(this%send_pe), this%send_buf%reqs)

    ! ...
    call device_sync()

  end subroutine gs_device_mpi_nbwait

end module gs_device_mpi
