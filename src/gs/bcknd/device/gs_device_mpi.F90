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
  type, private :: gs_comm_device_mpi_t
     type(c_ptr) :: request = C_NULL_PTR        !< MPI_Request in C
     type(c_ptr) :: buf_d = C_NULL_PTR          !< MPI send/recv buffer
     type(c_ptr) :: dof_d = C_NULL_PTR          !< buf->dof mapping (1-indexed)
     integer :: ndofs
     logical :: flag                            !< Request done flag
  end type gs_comm_device_mpi_t

  !> Gather-scatter communication using device MPI.
  !! The arrays are indexed per PE like @a send_pe and @ recv_pe.
  type, extends(gs_comm_t) :: gs_device_mpi_t
     type(gs_comm_device_mpi_t), allocatable :: send_buf(:)
     type(gs_comm_device_mpi_t), allocatable :: recv_buf(:)

     type(c_ptr) :: send_buf_ptrs_d = C_NULL_PTR    !< Array of buf pointers
     type(c_ptr) :: send_dof_ptrs_d = C_NULL_PTR    !< Array of dof pointers
     type(c_ptr) :: send_ndofs_d = C_NULL_PTR       !< Array of ndofs
     integer :: send_ndofs_max                      !< Max ndofs
   contains
     procedure, pass(this) :: init => gs_device_mpi_init
     procedure, pass(this) :: free => gs_device_mpi_free
     procedure, pass(this) :: nbsend => gs_device_mpi_nbsend
     procedure, pass(this) :: nbrecv => gs_device_mpi_nbrecv
     procedure, pass(this) :: nbwait => gs_device_mpi_nbwait
  end type gs_device_mpi_t

  interface
     subroutine hip_gs_pack(dof_ptrs_d, buf_ptrs_d, ndofs_d, npe, u_d, n, &
                            ndofs_max) &
          bind(c, name='hip_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: npe, n, ndofs_max
       type(c_ptr), value :: dof_ptrs_d, buf_ptrs_d, ndofs_d, u_d
     end subroutine hip_gs_pack
  end interface

  interface
     subroutine cuda_gs_pack(dof_ptrs_d, buf_ptrs_d, ndofs_d, npe, u_d, n, &
                             ndofs_max) &
          bind(c, name='cuda_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: npe, n, ndofs_max
       type(c_ptr), value :: dof_ptrs_d, buf_ptrs_d, ndofs_d, u_d
     end subroutine cuda_gs_pack
  end interface

  interface
    subroutine hip_gs_unpack(buf_d, dof_d, ndofs, u_d, op) &
          bind(c, name='hip_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: ndofs, op
       type(c_ptr), value :: buf_d, dof_d, u_d
     end subroutine hip_gs_unpack
  end interface

  interface
    subroutine cuda_gs_unpack(buf_d, dof_d, ndofs, u_d, op) &
          bind(c, name='cuda_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: ndofs, op
       type(c_ptr), value :: buf_d, dof_d, u_d
     end subroutine cuda_gs_unpack
  end interface

  interface
    subroutine device_mpi_init_request(req) &
          bind(c, name='device_mpi_init_request')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: req
     end subroutine device_mpi_init_request
  end interface

  interface
    subroutine device_mpi_free_request(req) &
          bind(c, name='device_mpi_free_request')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: req
     end subroutine device_mpi_free_request
  end interface

  interface
    subroutine device_mpi_isend(buf_d, nbytes, rank, req) &
          bind(c, name='device_mpi_isend')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: nbytes, rank
       type(c_ptr), value :: buf_d, req
     end subroutine device_mpi_isend
  end interface

  interface
    subroutine device_mpi_irecv(buf_d, nbytes, rank, req) &
          bind(c, name='device_mpi_irecv')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: nbytes, rank
       type(c_ptr), value :: buf_d, req
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

contains

  !> Initialise MPI based communication method
  subroutine gs_device_mpi_init(this, send_pe, recv_pe)
    class(gs_device_mpi_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    type(c_ptr), allocatable, target :: buf_ptrs(:), dof_ptrs(:)
    type(c_ptr) :: host_ptr
    integer, allocatable :: ndofs(:)
    integer, pointer :: pe(:)
    integer :: i, ndof
    integer(c_size_t) :: sz

    call this%init_order(send_pe, recv_pe)

    allocate(this%send_buf(send_pe%size()))
    allocate(buf_ptrs(send_pe%size()))
    allocate(dof_ptrs(send_pe%size()))
    allocate(ndofs(send_pe%size()))

    this%send_ndofs_max = 0

    pe => send_pe%array()
    do i = 1, send_pe%size()
       ndof = this%send_dof(pe(i))%size()
       ndofs(i) = ndof

       if (ndof .gt. this%send_ndofs_max) then
          this%send_ndofs_max = ndof
       end if

       sz = rp * ndof
       call device_alloc(buf_ptrs(i), sz)

       sz = 4 * ndof
       call device_alloc(dof_ptrs(i), sz)

       ! %array() breaks on cray
       select type (arr=>this%send_dof(pe(i))%data)
       type is (integer)
       call device_memcpy(arr, dof_ptrs(i), ndof, HOST_TO_DEVICE)
       end select

       this%send_buf(i)%buf_d = buf_ptrs(i)
       this%send_buf(i)%dof_d = dof_ptrs(i)
       this%send_buf(i)%ndofs = ndof

       call device_mpi_init_request(this%send_buf(i)%request)
    end do

    sz = send_pe%size() * c_sizeof(C_NULL_PTR)
    call device_alloc(this%send_buf_ptrs_d, sz)
    host_ptr = c_loc(buf_ptrs)
    call device_memcpy(host_ptr, this%send_buf_ptrs_d, &
                       sz, HOST_TO_DEVICE)

    sz = send_pe%size() * c_sizeof(C_NULL_PTR)
    call device_alloc(this%send_dof_ptrs_d, sz)
    host_ptr = c_loc(dof_ptrs)
    call device_memcpy(host_ptr, this%send_dof_ptrs_d, &
                       sz, HOST_TO_DEVICE)

    sz = 4 * send_pe%size()
    call device_alloc(this%send_ndofs_d, sz)
    call device_memcpy(ndofs, this%send_ndofs_d, send_pe%size(), HOST_TO_DEVICE)

    deallocate(buf_ptrs)
    deallocate(dof_ptrs)
    deallocate(ndofs)

    allocate(this%recv_buf(recv_pe%size()))

    pe => recv_pe%array()
    do i = 1, recv_pe%size()
       ndof = this%recv_dof(pe(i))%size()

       sz = rp * ndof
       call device_alloc(this%recv_buf(i)%buf_d, sz)

       sz = 4 * ndof
       call device_alloc(this%recv_buf(i)%dof_d, sz)

       ! %array() breaks on cray
       select type (arr=>this%recv_dof(pe(i))%data)
       type is (integer)
         call device_memcpy(arr, this%recv_buf(i)%dof_d, ndof, HOST_TO_DEVICE)
       end select

       this%recv_buf(i)%ndofs = ndof

       call device_mpi_init_request(this%recv_buf(i)%request)
    end do

  end subroutine gs_device_mpi_init

  !> Deallocate MPI based communication method
  subroutine gs_device_mpi_free(this)
    class(gs_device_mpi_t), intent(inout) :: this
    integer :: i

    if (allocated(this%send_buf)) then
       do i = 1, size(this%send_buf)
          if (c_associated(this%send_buf(i)%request)) then
             call device_mpi_free_request(this%send_buf(i)%request)
          end if
          if (c_associated(this%send_buf(i)%buf_d)) then
             call device_free(this%send_buf(i)%buf_d)
          end if
          if (c_associated(this%send_buf(i)%dof_d)) then
             call device_free(this%send_buf(i)%dof_d)
          end if
       end do
       deallocate(this%send_buf)
    end if

    if (c_associated(this%send_buf_ptrs_d)) then
       call device_free(this%send_buf_ptrs_d)
    end if
    if (c_associated(this%send_dof_ptrs_d)) then
       call device_free(this%send_dof_ptrs_d)
    end if
    if (c_associated(this%send_ndofs_d)) then
       call device_free(this%send_ndofs_d)
    end if

    if (allocated(this%recv_buf)) then
       do i = 1, size(this%recv_buf)
          if (c_associated(this%recv_buf(i)%request)) then
             call device_mpi_free_request(this%recv_buf(i)%request)
          end if
          if (c_associated(this%recv_buf(i)%buf_d)) then
             call device_free(this%recv_buf(i)%buf_d)
          end if
          if (c_associated(this%recv_buf(i)%dof_d)) then
             call device_free(this%recv_buf(i)%dof_d)
          end if
       end do
       deallocate(this%recv_buf)
    end if

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
    call hip_gs_pack(this%send_dof_ptrs_d, this%send_buf_ptrs_d, &
                     this%send_ndofs_d, size(this%send_pe), u_d, n, &
                     this%send_ndofs_max)
#elif HAVE_CUDA
    call cuda_gs_pack(this%send_dof_ptrs_d, this%send_buf_ptrs_d, &
                      this%send_ndofs_d, size(this%send_pe), u_d, n, &
                      this%send_ndofs_max)
#else
    call neko_error('gs_device_mpi: no backend')
#endif

    call device_sync()

    do i = 1, size(this%send_pe)
       call device_mpi_isend(this%send_buf(i)%buf_d, rp*this%send_buf(i)%ndofs, &
            this%send_pe(i), this%send_buf(i)%request)
       this%send_buf(i)%flag = .false.
    end do

  end subroutine gs_device_mpi_nbsend

  !> Post non-blocking receive operations
  subroutine gs_device_mpi_nbrecv(this)
    class(gs_device_mpi_t), intent(inout) :: this
    integer :: i

    do i = 1, size(this%recv_pe)
       call device_mpi_irecv(this%recv_buf(i)%buf_d, rp*this%recv_buf(i)%ndofs, &
            this%recv_pe(i), this%recv_buf(i)%request)
       this%recv_buf(i)%flag = .false.
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

    nreqs = size(this%recv_pe)
    do while (nreqs .gt. 0)
       do i = 1, size(this%recv_pe)
          if (.not. this%recv_buf(i)%flag) then
             flag = device_mpi_test(this%recv_buf(i)%request)
             this%recv_buf(i)%flag = flag .ne. 0

             if (this%recv_buf(i)%flag) then
                nreqs = nreqs - 1
#ifdef HAVE_HIP
                call hip_gs_unpack(this%recv_buf(i)%buf_d, &
                                   this%recv_buf(i)%dof_d, &
                                   this%recv_buf(i)%ndofs, &
                                   u_d, op)
#elif HAVE_CUDA
                call cuda_gs_unpack(this%recv_buf(i)%buf_d, &
                                    this%recv_buf(i)%dof_d, &
                                    this%recv_buf(i)%ndofs, &
                                    u_d, op)
#else
                call neko_error('gs_device_mpi: no backend')
#endif
             end if
          end if
       end do
    end do

    nreqs = size(this%send_pe)
    do while (nreqs .gt. 0)
       do i = 1, size(this%send_pe)
          if (.not. this%send_buf(i)%flag) then
             flag = device_mpi_test(this%send_buf(i)%request)
             this%send_buf(i)%flag = flag .ne. 0
             if (this%send_buf(i)%flag) nreqs = nreqs - 1
          end if
       end do
    end do

  end subroutine gs_device_mpi_nbwait

end module gs_device_mpi
