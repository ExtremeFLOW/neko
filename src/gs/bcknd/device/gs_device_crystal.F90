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
!> Defines GPU aware crystal router gather-scatter communication
!! @details
!! The device counterpart of @a gs_crystal: the same routing plan, with the
!! halo never leaving the device. Every word movement the routing asks for is
!! an indexed gather, which is what the existing gs pack kernel already does
!! (`buf[j] = u[dof[j]-1]`), so forwarding a stage costs one kernel launch and
!! no new device code -- the working buffer simply takes the place of the
!! shared vector as the kernel's source.
!!
!! The two working buffer columns are separate allocations, so each stage
!! addresses its source and destination from offset zero and the received
!! words land straight behind the ones that stay. As on the host, only the
!! first stage overlaps the local gather-scatter.
module gs_device_crystal
  use num_types, only : rp, c_rp
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use gs_crystal_plan, only : gs_crystal_plan_t
  use stack, only : stack_i4_t
  use htable, only : htable_i4_t
  use device, only : device_alloc, device_free, device_memcpy, &
       device_get_ptr, device_sync, HOST_TO_DEVICE
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated, &
       c_sizeof, c_size_t, c_int32_t
  implicit none
  private

  !> Gather-scatter communication using a crystal router on the device.
  type, public, extends(gs_comm_t) :: gs_device_crystal_t
     !> The routing plan, worked out once from the schedule
     type(gs_crystal_plan_t) :: plan
     !> Working buffer columns, plan%nwrk words each
     type(c_ptr) :: buf_d(2) = C_NULL_PTR
     !> Send buffer, plan%nsmax words
     type(c_ptr) :: sbuf_d = C_NULL_PTR
     !> Fused vector working and send buffers, GS_VEC_NC times as long.
     !! Positions are interleaved, matching the gs vector pack kernels
     type(c_ptr) :: buf_v_d(2) = C_NULL_PTR
     type(c_ptr) :: sbuf_v_d = C_NULL_PTR
     !> Per-stage index lists of the words that stay and the words that go,
     !! as positions in the stage's source column
     type(c_ptr), allocatable :: keep_idx_d(:), send_idx_d(:)
     !> The same lists expanded over the components of the fused vector
     !! exchange, built on first use for the number of components seen
     type(c_ptr), allocatable :: keep_idx_v_d(:), send_idx_v_d(:)
     !> Shared vector indices for the first stage and for the reduction
     type(c_ptr) :: pack_keep_d = C_NULL_PTR
     type(c_ptr) :: pack_send_d = C_NULL_PTR
     type(c_ptr) :: unpack_d = C_NULL_PTR
     !> Requests of the stage in flight, at most one send and two receives
     type(c_ptr) :: sreqs = C_NULL_PTR
     type(c_ptr) :: rreqs = C_NULL_PTR
     integer :: nsreq = 0, nrreq = 0
     !> Tag of the operation in flight, taken from nbrecv since the later
     !! stages are driven from nbwait, which is not given one
     integer :: tag = 0
     !> Number of components the expanded index lists were built for, 0
     !! until the first fused vector exchange
     integer :: vec_nc = 0
   contains
     procedure, pass(this) :: init => gs_device_crystal_init
     procedure, pass(this) :: free => gs_device_crystal_free
     procedure, pass(this) :: nbsend => gs_device_crystal_nbsend
     procedure, pass(this) :: nbrecv => gs_device_crystal_nbrecv
     procedure, pass(this) :: nbwait => gs_device_crystal_nbwait
     procedure, pass(this) :: nbsend_vec => gs_device_crystal_nbsend_vec
     procedure, pass(this) :: nbrecv_vec => gs_device_crystal_nbrecv_vec
     procedure, pass(this) :: nbwait_vec => gs_device_crystal_nbwait_vec
  end type gs_device_crystal_t

#ifdef HAVE_HIP
  interface
     subroutine hip_gs_pack(u_d, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'hip_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_pack
  end interface

  interface
     subroutine hip_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'hip_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_unpack
  end interface

  interface
     subroutine hip_gs_pack_vec(u_d, buf_d, dof_d, offset, n, nc, ns, stream) &
          bind(c, name = 'hip_gs_pack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_pack_vec
  end interface

  interface
     subroutine hip_gs_unpack_vec(u_d, op, buf_d, dof_d, offset, n, nc, ns, &
          stream) bind(c, name = 'hip_gs_unpack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_unpack_vec
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_gs_pack(u_d, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'cuda_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_pack
  end interface

  interface
     subroutine cuda_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name = 'cuda_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack
  end interface

  interface
     subroutine cuda_gs_pack_vec(u_d, buf_d, dof_d, offset, n, nc, ns, stream) &
          bind(c, name = 'cuda_gs_pack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_pack_vec
  end interface

  interface
     subroutine cuda_gs_unpack_vec(u_d, op, buf_d, dof_d, offset, n, nc, ns, &
          stream) bind(c, name = 'cuda_gs_unpack_vec')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n, nc, ns
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack_vec
  end interface
#endif

  interface
     subroutine device_mpi_init_reqs(n, reqs) &
          bind(c, name = 'device_mpi_init_reqs')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n
       type(c_ptr) :: reqs
     end subroutine device_mpi_init_reqs
  end interface

  interface
     subroutine device_mpi_free_reqs(reqs) &
          bind(c, name = 'device_mpi_free_reqs')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: reqs
     end subroutine device_mpi_free_reqs
  end interface

  interface
     subroutine device_mpi_isend(buf_d, offset, nbytes, rank, tag, reqs, i) &
          bind(c, name = 'device_mpi_isend')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, nbytes, rank, tag, i
       type(c_ptr), value :: buf_d, reqs
     end subroutine device_mpi_isend
  end interface

  interface
     subroutine device_mpi_irecv(buf_d, offset, nbytes, rank, tag, reqs, i) &
          bind(c, name = 'device_mpi_irecv')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: offset, nbytes, rank, tag, i
       type(c_ptr), value :: buf_d, reqs
     end subroutine device_mpi_irecv
  end interface

  interface
     subroutine device_mpi_waitall(n, reqs) &
          bind(c, name = 'device_mpi_waitall')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n
       type(c_ptr), value :: reqs
     end subroutine device_mpi_waitall
  end interface

contains

  !> Initialise crystal router based device communication
  !! See gs_comm.f90 for details
  subroutine gs_device_crystal_init(this, send_pe, recv_pe)
    class(gs_device_crystal_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer, allocatable :: dofs(:)
    integer :: i, nst
    integer(c_size_t) :: sz
    real(c_rp) :: rp_dummy

    call this%init_order(send_pe, recv_pe)

    call this%plan%init(this%send_pe, this%recv_pe, this%send_dof, &
         this%recv_dof)

    call device_mpi_init_reqs(1, this%sreqs)
    call device_mpi_init_reqs(2, this%rreqs)

    sz = c_sizeof(rp_dummy) * this%plan%nwrk
    call device_alloc(this%buf_d(1), sz)
    call device_alloc(this%buf_d(2), sz)
    sz = c_sizeof(rp_dummy) * this%plan%nsmax
    call device_alloc(this%sbuf_d, sz)

    sz = c_sizeof(rp_dummy) * GS_VEC_NC * this%plan%nwrk
    call device_alloc(this%buf_v_d(1), sz)
    call device_alloc(this%buf_v_d(2), sz)
    sz = c_sizeof(rp_dummy) * GS_VEC_NC * this%plan%nsmax
    call device_alloc(this%sbuf_v_d, sz)

    ! The first stage gathers straight out of the shared vector
    if (this%plan%nstage .gt. 0) then
       call cr_upload(this%pack_keep_d, this%plan%pack_keep_dof, &
            this%plan%stage(1)%nkw)
       call cr_upload(this%pack_send_d, this%plan%pack_send_dof, &
            this%plan%stage(1)%nsw)
    end if

    ! A shared dof may be delivered by several peers, so the reduction
    ! marks the repeats for the unpack kernel to handle atomically
    allocate(dofs(max(this%plan%nfinal, 1)))
    call cr_mark_dupes(this%plan%unpack_dof, dofs, this%plan%nfinal)
    call cr_upload(this%unpack_d, dofs, this%plan%nfinal)
    deallocate(dofs)

    nst = max(this%plan%nstage, 1)
    allocate(this%keep_idx_d(nst), this%send_idx_d(nst))
    allocate(this%keep_idx_v_d(nst), this%send_idx_v_d(nst))
    this%keep_idx_d = C_NULL_PTR
    this%send_idx_d = C_NULL_PTR
    this%keep_idx_v_d = C_NULL_PTR
    this%send_idx_v_d = C_NULL_PTR

    ! Stage one has no index lists: its words come from the shared vector
    do i = 2, this%plan%nstage
       associate (st => this%plan%stage(i))
         if (.not. st%inplace) then
            call cr_upload(this%keep_idx_d(i), st%keep_idx, st%nkw)
         end if
         if (st%dst .ge. 0) then
            call cr_upload(this%send_idx_d(i), st%send_idx, st%nsw)
         end if
       end associate
    end do

    this%vec_supported = .true.

  end subroutine gs_device_crystal_init

  !> Deallocate crystal router based device communication
  subroutine gs_device_crystal_free(this)
    class(gs_device_crystal_t), intent(inout) :: this
    integer :: i

    if (c_associated(this%sreqs)) call device_mpi_free_reqs(this%sreqs)
    if (c_associated(this%rreqs)) call device_mpi_free_reqs(this%rreqs)

    do i = 1, 2
       if (c_associated(this%buf_d(i))) call device_free(this%buf_d(i))
       if (c_associated(this%buf_v_d(i))) call device_free(this%buf_v_d(i))
    end do
    if (c_associated(this%sbuf_d)) call device_free(this%sbuf_d)
    if (c_associated(this%sbuf_v_d)) call device_free(this%sbuf_v_d)

    if (c_associated(this%pack_keep_d)) call device_free(this%pack_keep_d)
    if (c_associated(this%pack_send_d)) call device_free(this%pack_send_d)
    if (c_associated(this%unpack_d)) call device_free(this%unpack_d)

    call cr_free_ptrs(this%keep_idx_d)
    call cr_free_ptrs(this%send_idx_d)
    call cr_free_ptrs(this%keep_idx_v_d)
    call cr_free_ptrs(this%send_idx_v_d)

    if (allocated(this%keep_idx_d)) deallocate(this%keep_idx_d)
    if (allocated(this%send_idx_d)) deallocate(this%send_idx_d)
    if (allocated(this%keep_idx_v_d)) deallocate(this%keep_idx_v_d)
    if (allocated(this%send_idx_v_d)) deallocate(this%send_idx_v_d)

    this%vec_nc = 0

    call this%plan%free()

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_device_crystal_free

  !> Post the receives of the first routing stage
  subroutine gs_device_crystal_nbrecv(this, tag)
    class(gs_device_crystal_t), intent(inout) :: this
    integer, intent(in) :: tag

    this%tag = tag
    this%nrreq = 0
    this%nsreq = 0
    if (this%plan%nstage .eq. 0) return

    associate (st => this%plan%stage(1))
      if (st%src .ge. 0) then
         this%nrreq = this%nrreq + 1
         call device_mpi_irecv(this%buf_d(st%dst_sel), rp*st%nkw, &
              rp*st%nrw, st%src, tag, this%rreqs, this%nrreq)
      end if
      if (st%src2 .ge. 0) then
         this%nrreq = this%nrreq + 1
         call device_mpi_irecv(this%buf_d(st%dst_sel), &
              rp*(st%nkw + st%nrw), rp*st%nr2w, st%src2, tag, &
              this%rreqs, this%nrreq)
      end if
    end associate

  end subroutine gs_device_crystal_nbrecv

  !> Pack the shared vector and post the send of the first routing stage
  subroutine gs_device_crystal_nbsend(this, u, n, tag, deps, strm)
    class(gs_device_crystal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    type(c_ptr) :: u_d

    if (this%plan%nstage .eq. 0) return

    u_d = device_get_ptr(u)

    associate (st => this%plan%stage(1))
      if (st%dst .ge. 0) then
         call cr_gather(u_d, this%sbuf_d, this%pack_send_d, st%nsw, strm)
         call device_sync(strm)
         this%nsreq = 1
         call device_mpi_isend(this%sbuf_d, 0, rp*st%nsw, st%dst, tag, &
              this%sreqs, 1)
      end if

      ! The words that stay put, into the column the arriving ones are
      ! being received behind
      call cr_gather(u_d, this%buf_d(st%dst_sel), this%pack_keep_d, &
           st%nkw, strm)
    end associate

  end subroutine gs_device_crystal_nbsend

  !> Drive the remaining routing stages and reduce what is delivered into
  !! the shared vector
  subroutine gs_device_crystal_nbwait(this, u, n, op, strm)
    class(gs_device_crystal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer :: s
    type(c_ptr) :: u_d

    if (this%plan%nstage .eq. 0) return

    u_d = device_get_ptr(u)

    call device_mpi_waitall(this%nrreq, this%rreqs)
    call device_mpi_waitall(this%nsreq, this%sreqs)

    do s = 2, this%plan%nstage
       associate (st => this%plan%stage(s))
         this%nrreq = 0
         this%nsreq = 0
         if (st%src .ge. 0) then
            this%nrreq = this%nrreq + 1
            call device_mpi_irecv(this%buf_d(st%dst_sel), rp*st%nkw, &
                 rp*st%nrw, st%src, this%tag, this%rreqs, this%nrreq)
         end if
         if (st%src2 .ge. 0) then
            this%nrreq = this%nrreq + 1
            call device_mpi_irecv(this%buf_d(st%dst_sel), &
                 rp*(st%nkw + st%nrw), rp*st%nr2w, st%src2, this%tag, &
                 this%rreqs, this%nrreq)
         end if

         if (st%dst .ge. 0) then
            call cr_gather(this%buf_d(st%src_sel), this%sbuf_d, &
                 this%send_idx_d(s), st%nsw, strm)
            call device_sync(strm)
            this%nsreq = 1
            call device_mpi_isend(this%sbuf_d, 0, rp*st%nsw, st%dst, &
                 this%tag, this%sreqs, 1)
         end if

         ! A stage that sends nothing keeps every word where it is, and is
         ! receiving into the column it already occupies
         if (.not. st%inplace) then
            call cr_gather(this%buf_d(st%src_sel), this%buf_d(st%dst_sel), &
                 this%keep_idx_d(s), st%nkw, strm)
         end if

         call device_mpi_waitall(this%nrreq, this%rreqs)
         call device_mpi_waitall(this%nsreq, this%sreqs)
         call device_sync(strm)
       end associate
    end do

    call cr_scatter(u_d, op, this%buf_d(this%plan%final_sel), &
         this%unpack_d, this%plan%nfinal, strm)

    call device_sync(strm)

  end subroutine gs_device_crystal_nbwait

  !> Post the receives of the first routing stage, fused nc-component
  subroutine gs_device_crystal_nbrecv_vec(this, tag, nc)
    class(gs_device_crystal_t), intent(inout) :: this
    integer, intent(in) :: tag, nc

    if (nc .gt. GS_VEC_NC) then
       call neko_error('gs_device_crystal: too many components in ' // &
            'vector exchange')
    end if

    call cr_vec_index(this, nc)

    this%tag = tag
    this%nrreq = 0
    this%nsreq = 0
    if (this%plan%nstage .eq. 0) return

    associate (st => this%plan%stage(1))
      if (st%src .ge. 0) then
         this%nrreq = this%nrreq + 1
         call device_mpi_irecv(this%buf_v_d(st%dst_sel), rp*nc*st%nkw, &
              rp*nc*st%nrw, st%src, tag, this%rreqs, this%nrreq)
      end if
      if (st%src2 .ge. 0) then
         this%nrreq = this%nrreq + 1
         call device_mpi_irecv(this%buf_v_d(st%dst_sel), &
              rp*nc*(st%nkw + st%nrw), rp*nc*st%nr2w, st%src2, tag, &
              this%rreqs, this%nrreq)
      end if
    end associate

  end subroutine gs_device_crystal_nbrecv_vec

  !> Pack the shared vector and post the send of the first routing stage,
  !! fused nc-component
  !! @param u compact shared device buffer, component-outer with stride n
  subroutine gs_device_crystal_nbsend_vec(this, u, n, nc, tag, deps, strm)
    class(gs_device_crystal_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    type(c_ptr) :: u_d

    if (this%plan%nstage .eq. 0) return

    u_d = device_get_ptr(u)

    associate (st => this%plan%stage(1))
      if (st%dst .ge. 0) then
         call cr_gather_vec(u_d, this%sbuf_v_d, this%pack_send_d, st%nsw, &
              nc, n, strm)
         call device_sync(strm)
         this%nsreq = 1
         call device_mpi_isend(this%sbuf_v_d, 0, rp*nc*st%nsw, st%dst, &
              tag, this%sreqs, 1)
      end if

      call cr_gather_vec(u_d, this%buf_v_d(st%dst_sel), this%pack_keep_d, &
           st%nkw, nc, n, strm)
    end associate

  end subroutine gs_device_crystal_nbsend_vec

  !> Drive the remaining routing stages and reduce what is delivered into
  !! the shared vector, fused nc-component
  subroutine gs_device_crystal_nbwait_vec(this, u, n, nc, op, strm)
    class(gs_device_crystal_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer :: s
    type(c_ptr) :: u_d

    if (this%plan%nstage .eq. 0) return

    u_d = device_get_ptr(u)

    call device_mpi_waitall(this%nrreq, this%rreqs)
    call device_mpi_waitall(this%nsreq, this%sreqs)

    do s = 2, this%plan%nstage
       associate (st => this%plan%stage(s))
         this%nrreq = 0
         this%nsreq = 0
         if (st%src .ge. 0) then
            this%nrreq = this%nrreq + 1
            call device_mpi_irecv(this%buf_v_d(st%dst_sel), rp*nc*st%nkw, &
                 rp*nc*st%nrw, st%src, this%tag, this%rreqs, this%nrreq)
         end if
         if (st%src2 .ge. 0) then
            this%nrreq = this%nrreq + 1
            call device_mpi_irecv(this%buf_v_d(st%dst_sel), &
                 rp*nc*(st%nkw + st%nrw), rp*nc*st%nr2w, st%src2, &
                 this%tag, this%rreqs, this%nrreq)
         end if

         ! The interleaved layout makes a stage's word movement the same
         ! indexed gather as in the scalar case, over nc times as many
         ! words, so the plain pack kernel drives it here too
         if (st%dst .ge. 0) then
            call cr_gather(this%buf_v_d(st%src_sel), this%sbuf_v_d, &
                 this%send_idx_v_d(s), nc*st%nsw, strm)
            call device_sync(strm)
            this%nsreq = 1
            call device_mpi_isend(this%sbuf_v_d, 0, rp*nc*st%nsw, st%dst, &
                 this%tag, this%sreqs, 1)
         end if

         if (.not. st%inplace) then
            call cr_gather(this%buf_v_d(st%src_sel), &
                 this%buf_v_d(st%dst_sel), this%keep_idx_v_d(s), &
                 nc*st%nkw, strm)
         end if

         call device_mpi_waitall(this%nrreq, this%rreqs)
         call device_mpi_waitall(this%nsreq, this%sreqs)
         call device_sync(strm)
       end associate
    end do

    call cr_scatter_vec(u_d, op, this%buf_v_d(this%plan%final_sel), &
         this%unpack_d, this%plan%nfinal, nc, n, strm)

    call device_sync(strm)

  end subroutine gs_device_crystal_nbwait_vec

  !> Gather @a n words of @a src into @a dst through the 1-based index list
  !! @a idx_d, i.e. `dst(j) = src(idx(j))`
  subroutine cr_gather(src_d, dst_d, idx_d, n, strm)
    type(c_ptr), intent(in) :: src_d, dst_d, idx_d
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: strm

    if (n .le. 0) return

#ifdef HAVE_HIP
    call hip_gs_pack(src_d, dst_d, idx_d, 0, n, strm)
#elif HAVE_CUDA
    call cuda_gs_pack(src_d, dst_d, idx_d, 0, n, strm)
#else
    call neko_error('gs_device_crystal: no backend')
#endif

  end subroutine cr_gather

  !> Gather @a n positions of the component-outer @a src (stride @a ns) into
  !! the interleaved @a dst
  subroutine cr_gather_vec(src_d, dst_d, idx_d, n, nc, ns, strm)
    type(c_ptr), intent(in) :: src_d, dst_d, idx_d
    integer, intent(in) :: n, nc, ns
    type(c_ptr), intent(inout) :: strm

    if (n .le. 0) return

#ifdef HAVE_HIP
    call hip_gs_pack_vec(src_d, dst_d, idx_d, 0, n, nc, ns, strm)
#elif HAVE_CUDA
    call cuda_gs_pack_vec(src_d, dst_d, idx_d, 0, n, nc, ns, strm)
#else
    call neko_error('gs_device_crystal: no backend')
#endif

  end subroutine cr_gather_vec

  !> Reduce @a n delivered words into the shared vector under @a op
  subroutine cr_scatter(u_d, op, buf_d, idx_d, n, strm)
    type(c_ptr), intent(in) :: u_d, buf_d, idx_d
    integer, intent(in) :: op, n
    type(c_ptr), intent(inout) :: strm

    if (n .le. 0) return

#ifdef HAVE_HIP
    call hip_gs_unpack(u_d, op, buf_d, idx_d, 0, n, strm)
#elif HAVE_CUDA
    call cuda_gs_unpack(u_d, op, buf_d, idx_d, 0, n, strm)
#else
    call neko_error('gs_device_crystal: no backend')
#endif

  end subroutine cr_scatter

  !> Reduce @a n delivered positions into the component-outer shared vector
  subroutine cr_scatter_vec(u_d, op, buf_d, idx_d, n, nc, ns, strm)
    type(c_ptr), intent(in) :: u_d, buf_d, idx_d
    integer, intent(in) :: op, n, nc, ns
    type(c_ptr), intent(inout) :: strm

    if (n .le. 0) return

#ifdef HAVE_HIP
    call hip_gs_unpack_vec(u_d, op, buf_d, idx_d, 0, n, nc, ns, strm)
#elif HAVE_CUDA
    call cuda_gs_unpack_vec(u_d, op, buf_d, idx_d, 0, n, nc, ns, strm)
#else
    call neko_error('gs_device_crystal: no backend')
#endif

  end subroutine cr_scatter_vec

  !> Build the per-stage index lists of the fused vector exchange, which
  !! address the same words as the scalar ones with the components spelled
  !! out: position p of a stage becomes nc*(p-1) + 1 .. nc*p.
  !! @param nc number of components to build the lists for. Rebuilt if a
  !! later exchange asks for a different number, which nothing does today
  subroutine cr_vec_index(this, nc)
    class(gs_device_crystal_t), intent(inout) :: this
    integer, intent(in) :: nc
    integer, allocatable :: idx(:)
    integer :: i

    if (this%vec_nc .eq. nc) return

    call cr_free_ptrs(this%keep_idx_v_d)
    call cr_free_ptrs(this%send_idx_v_d)

    do i = 2, this%plan%nstage
       associate (st => this%plan%stage(i))
         if (.not. st%inplace) then
            call cr_expand(st%keep_idx, st%nkw, nc, idx)
            call cr_upload(this%keep_idx_v_d(i), idx, nc*st%nkw)
            deallocate(idx)
         end if
         if (st%dst .ge. 0) then
            call cr_expand(st%send_idx, st%nsw, nc, idx)
            call cr_upload(this%send_idx_v_d(i), idx, nc*st%nsw)
            deallocate(idx)
         end if
       end associate
    end do

    this%vec_nc = nc

  end subroutine cr_vec_index

  !> Spell out the components of the 1-based index list @a idx over @a nc
  subroutine cr_expand(idx, n, nc, out)
    integer, intent(in) :: idx(:)
    integer, intent(in) :: n, nc
    integer, allocatable, intent(out) :: out(:)
    integer :: j, c

    allocate(out(max(nc*n, 1)))
    do j = 1, n
       do c = 1, nc
          out(nc*(j-1) + c) = nc*(idx(j) - 1) + c
       end do
    end do

  end subroutine cr_expand

  !> Copy the first @a n entries of an index list to the device, leaving
  !! @a ptr null when there are none
  subroutine cr_upload(ptr, idx, n)
    type(c_ptr), intent(inout) :: ptr
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: n
    integer(c_size_t) :: sz
    integer(c_int32_t) :: i4_dummy

    if (c_associated(ptr)) call device_free(ptr)
    ptr = C_NULL_PTR

    if (n .le. 0) return

    sz = c_sizeof(i4_dummy) * n
    call device_alloc(ptr, sz)
    call device_memcpy(idx, ptr, n, HOST_TO_DEVICE, sync = .true.)

  end subroutine cr_upload

  !> Release a list of device pointers, leaving them null
  subroutine cr_free_ptrs(ptrs)
    type(c_ptr), allocatable, intent(inout) :: ptrs(:)
    integer :: i

    if (.not. allocated(ptrs)) return

    do i = 1, size(ptrs)
       if (c_associated(ptrs(i))) call device_free(ptrs(i))
       ptrs(i) = C_NULL_PTR
    end do

  end subroutine cr_free_ptrs

  !> Copy @a dof into @a out, negating every index that appears more than
  !! once so that the unpack kernel reduces those atomically
  subroutine cr_mark_dupes(dof, out, n)
    integer, intent(in) :: dof(:)
    integer, intent(out) :: out(:)
    integer, intent(in) :: n
    type(htable_i4_t) :: doftable
    integer :: j, dupe, key, val

    if (n .le. 0) return

    ! The table takes its key and data as intent(inout), so neither the
    ! index list nor the loop counter can be handed to it directly
    call doftable%init(2*n)
    do j = 1, n
       key = dof(j)
       if (doftable%get(key, dupe) .eq. 0) then
          if (out(dupe) .gt. 0) out(dupe) = -out(dupe)
          out(j) = -dof(j)
       else
          key = dof(j)
          val = j
          call doftable%set(key, val)
          out(j) = dof(j)
       end if
    end do
    call doftable%free()

  end subroutine cr_mark_dupes

end module gs_device_crystal
