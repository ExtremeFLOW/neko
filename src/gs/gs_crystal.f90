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
!> Defines crystal router gather-scatter communication
!! @details
!! Runs the halo exchange as a crystal router rather than as a message to
!! every peer: the routing plan worked out once by @a gs_crystal_plan says,
!! for each stage, which words leave for a single partner, which stay, and
!! how many arrive. One gs operation therefore costs at most one message per
!! active stage instead of one per peer, at the price of forwarding the words
!! that are not yet home -- see gs_crystal_plan.f90 for when that trade pays.
!!
!! Only the first stage overlaps the local gather-scatter: it is posted from
!! @a nbsend and completed in @a nbwait, which then drives the remaining
!! stages back to back. Stages are dependent, so there is nothing to hide the
!! later ones behind.
!!
!! All MPI is issued from the master thread, so the backend needs nothing
!! beyond MPI_THREAD_FUNNELED; the packing, forwarding and reduction loops
!! are shared out over the team.
module gs_crystal
  use num_types, only : rp
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use gs_crystal_plan, only : gs_crystal_plan_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Request, MPI_Isend, MPI_Irecv, MPI_Waitall, &
       MPI_STATUSES_IGNORE
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Gather-scatter communication using a crystal router.
  type, public, extends(gs_comm_t) :: gs_crystal_t
     !> The routing plan, worked out once from the schedule
     type(gs_crystal_plan_t) :: plan
     !> Working buffer, two columns of plan%nwrk laid end to end. Column q
     !! starts at (q-1)*plan%nwrk
     real(kind=rp), allocatable :: buf(:)
     !> Send buffer, plan%nsmax words
     real(kind=rp), allocatable :: sbuf(:)
     !> Fused vector working and send buffers, the same layout scaled by
     !! GS_VEC_NC components. Positions are interleaved: position p holds
     !! its nc components at nc*(p-1) + 1 .. nc*p
     real(kind=rp), allocatable :: buf_v(:)
     real(kind=rp), allocatable :: sbuf_v(:)
     !> Requests of the stage in flight, at most one send and two receives
     type(MPI_Request) :: sreq(1), rreq(2)
     integer :: nsreq = 0, nrreq = 0
     !> Tag of the operation in flight, taken from nbrecv since the later
     !! stages are driven from nbwait, which is not given one
     integer :: tag = 0
   contains
     procedure, pass(this) :: init => gs_crystal_init
     procedure, pass(this) :: free => gs_crystal_free
     !> See gs_comm.f90 for more details on these routines
     procedure, pass(this) :: nbsend => gs_crystal_nbsend
     procedure, pass(this) :: nbrecv => gs_crystal_nbrecv
     procedure, pass(this) :: nbwait => gs_crystal_nbwait
     procedure, pass(this) :: nbsend_vec => gs_crystal_nbsend_vec
     procedure, pass(this) :: nbrecv_vec => gs_crystal_nbrecv_vec
     procedure, pass(this) :: nbwait_vec => gs_crystal_nbwait_vec
  end type gs_crystal_t

contains

  !> Initialise crystal router based communication method
  !! See gs_comm.f90 for details
  subroutine gs_crystal_init(this, send_pe, recv_pe)
    class(gs_crystal_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe

    call this%init_order(send_pe, recv_pe)

    call this%plan%init(this%send_pe, this%recv_pe, this%send_dof, &
         this%recv_dof)

    allocate(this%buf(2*this%plan%nwrk))
    allocate(this%sbuf(this%plan%nsmax))
    allocate(this%buf_v(2*GS_VEC_NC*this%plan%nwrk))
    allocate(this%sbuf_v(GS_VEC_NC*this%plan%nsmax))

    this%vec_supported = .true.

  end subroutine gs_crystal_init

  !> Deallocate crystal router based communication method
  subroutine gs_crystal_free(this)
    class(gs_crystal_t), intent(inout) :: this

    if (allocated(this%buf)) deallocate(this%buf)
    if (allocated(this%sbuf)) deallocate(this%sbuf)
    if (allocated(this%buf_v)) deallocate(this%buf_v)
    if (allocated(this%sbuf_v)) deallocate(this%sbuf_v)

    call this%plan%free()

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_crystal_free

  !> Post the receives of the first routing stage
  subroutine gs_crystal_nbrecv(this, tag)
    class(gs_crystal_t), intent(inout) :: this
    integer, intent(in) :: tag
    integer :: co, ierr

    !$omp master
    this%tag = tag
    this%nrreq = 0
    this%nsreq = 0
    if (this%plan%nstage .gt. 0) then
       associate (st => this%plan%stage(1))
         ! The words that stay are packed into [1, nkw] of the same column
         ! by nbsend, so the arriving ones can land straight behind them
         co = (st%dst_sel - 1) * this%plan%nwrk
         if (st%src .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf(co + st%nkw + 1), st%nrw, &
                 MPI_REAL_PRECISION, st%src, tag, NEKO_COMM, &
                 this%rreq(this%nrreq), ierr)
         end if
         if (st%src2 .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf(co + st%nkw + st%nrw + 1), st%nr2w, &
                 MPI_REAL_PRECISION, st%src2, tag, NEKO_COMM, &
                 this%rreq(this%nrreq), ierr)
         end if
       end associate
    end if
    !$omp end master
    !$omp barrier

  end subroutine gs_crystal_nbrecv

  !> Pack the shared vector and post the send of the first routing stage
  subroutine gs_crystal_nbsend(this, u, n, tag, deps, strm)
    class(gs_crystal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: j, co, ierr

    if (this%plan%nstage .eq. 0) return

    associate (st => this%plan%stage(1))
      co = (st%dst_sel - 1) * this%plan%nwrk

      ! Gather the words leaving on the first stage straight out of the
      ! shared vector, so the routing never sees a separate packed copy
      !$omp do
      do j = 1, st%nsw
         this%sbuf(j) = u(this%plan%pack_send_dof(j))
      end do
      !$omp end do

      !$omp master
      if (st%dst .ge. 0) then
         this%nsreq = 1
         call MPI_Isend(this%sbuf(1), st%nsw, MPI_REAL_PRECISION, &
              st%dst, tag, NEKO_COMM, this%sreq(1), ierr)
      end if
      !$omp end master

      ! The words that stay put on this stage, into the same column the
      ! arriving ones are being received behind
      !$omp do
      do j = 1, st%nkw
         this%buf(co + j) = u(this%plan%pack_keep_dof(j))
      end do
      !$omp end do
    end associate

  end subroutine gs_crystal_nbsend

  !> Drive the remaining routing stages and reduce what is delivered into
  !! the shared vector
  subroutine gs_crystal_nbwait(this, u, n, op, strm)
    class(gs_crystal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer :: s, i, j, co, so, fo, ierr

    if (this%plan%nstage .eq. 0) return

    !$omp master
    call MPI_Waitall(this%nrreq, this%rreq, MPI_STATUSES_IGNORE, ierr)
    call MPI_Waitall(this%nsreq, this%sreq, MPI_STATUSES_IGNORE, ierr)
    !$omp end master
    !$omp barrier

    do s = 2, this%plan%nstage
       associate (st => this%plan%stage(s))
         so = (st%src_sel - 1) * this%plan%nwrk
         co = (st%dst_sel - 1) * this%plan%nwrk

         !$omp master
         this%nrreq = 0
         this%nsreq = 0
         if (st%src .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf(co + st%nkw + 1), st%nrw, &
                 MPI_REAL_PRECISION, st%src, this%tag, NEKO_COMM, &
                 this%rreq(this%nrreq), ierr)
         end if
         if (st%src2 .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf(co + st%nkw + st%nrw + 1), st%nr2w, &
                 MPI_REAL_PRECISION, st%src2, this%tag, NEKO_COMM, &
                 this%rreq(this%nrreq), ierr)
         end if
         !$omp end master

         if (st%dst .ge. 0) then
            !$omp do
            do j = 1, st%nsw
               this%sbuf(j) = this%buf(so + st%send_idx(j))
            end do
            !$omp end do
            !$omp master
            this%nsreq = 1
            call MPI_Isend(this%sbuf(1), st%nsw, MPI_REAL_PRECISION, &
                 st%dst, this%tag, NEKO_COMM, this%sreq(1), ierr)
            !$omp end master
         end if

         ! A stage that sends nothing keeps every word where it is, and is
         ! receiving into the column it already occupies
         if (.not. st%inplace) then
            !$omp do
            do j = 1, st%nkw
               this%buf(co + j) = this%buf(so + st%keep_idx(j))
            end do
            !$omp end do
         end if

         !$omp master
         call MPI_Waitall(this%nrreq, this%rreq, MPI_STATUSES_IGNORE, ierr)
         call MPI_Waitall(this%nsreq, this%sreq, MPI_STATUSES_IGNORE, ierr)
         !$omp end master
         !$omp barrier
       end associate
    end do

    ! Reduce record by record: a shared dof may be delivered by several
    ! peers, so its contributions are separated by the record barriers
    fo = (this%plan%final_sel - 1) * this%plan%nwrk
    do i = 1, this%plan%nfinal_rec
       associate (ro => this%plan%final_off(i), rl => this%plan%final_len(i))
         select case (op)
         case (GS_OP_ADD)
            !OCL NORECURRENCE, NOVREC, NOALIAS
            !DIR$ CONCURRENT
            !DIR$ IVDEP
            !GCC$ ivdep
            !NEC$ IVDEP
            !$omp do
            do j = 1, rl
               u(this%plan%unpack_dof(ro + j)) = &
                    u(this%plan%unpack_dof(ro + j)) + this%buf(fo + ro + j)
            end do
            !$omp end do
         case (GS_OP_MUL)
            !OCL NORECURRENCE, NOVREC, NOALIAS
            !DIR$ CONCURRENT
            !DIR$ IVDEP
            !GCC$ ivdep
            !NEC$ IVDEP
            !$omp do
            do j = 1, rl
               u(this%plan%unpack_dof(ro + j)) = &
                    u(this%plan%unpack_dof(ro + j)) * this%buf(fo + ro + j)
            end do
            !$omp end do
         case (GS_OP_MIN)
            !OCL NORECURRENCE, NOVREC, NOALIAS
            !DIR$ CONCURRENT
            !DIR$ IVDEP
            !GCC$ ivdep
            !NEC$ IVDEP
            !$omp do
            do j = 1, rl
               u(this%plan%unpack_dof(ro + j)) = &
                    min(u(this%plan%unpack_dof(ro + j)), &
                    this%buf(fo + ro + j))
            end do
            !$omp end do
         case (GS_OP_MAX)
            !OCL NORECURRENCE, NOVREC, NOALIAS
            !DIR$ CONCURRENT
            !DIR$ IVDEP
            !GCC$ ivdep
            !NEC$ IVDEP
            !$omp do
            do j = 1, rl
               u(this%plan%unpack_dof(ro + j)) = &
                    max(u(this%plan%unpack_dof(ro + j)), &
                    this%buf(fo + ro + j))
            end do
            !$omp end do
         case default
            call neko_error("Unknown operation in gs_crystal_nbwait")
         end select
       end associate
    end do

  end subroutine gs_crystal_nbwait

  !> Post the receives of the first routing stage, fused nc-component
  subroutine gs_crystal_nbrecv_vec(this, tag, nc)
    class(gs_crystal_t), intent(inout) :: this
    integer, intent(in) :: tag, nc
    integer :: co, ierr

    if (nc .gt. GS_VEC_NC) then
       call neko_error('gs_crystal: too many components in vector exchange')
    end if

    !$omp master
    this%tag = tag
    this%nrreq = 0
    this%nsreq = 0
    if (this%plan%nstage .gt. 0) then
       associate (st => this%plan%stage(1))
         co = (st%dst_sel - 1) * nc * this%plan%nwrk
         if (st%src .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf_v(co + nc*st%nkw + 1), nc*st%nrw, &
                 MPI_REAL_PRECISION, st%src, tag, NEKO_COMM, &
                 this%rreq(this%nrreq), ierr)
         end if
         if (st%src2 .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf_v(co + nc*(st%nkw + st%nrw) + 1), &
                 nc*st%nr2w, MPI_REAL_PRECISION, st%src2, tag, NEKO_COMM, &
                 this%rreq(this%nrreq), ierr)
         end if
       end associate
    end if
    !$omp end master
    !$omp barrier

  end subroutine gs_crystal_nbrecv_vec

  !> Pack the shared vector and post the send of the first routing stage,
  !! fused nc-component
  !! @param u compact shared buffer, component-outer: u((c-1)*n + idx)
  subroutine gs_crystal_nbsend_vec(this, u, n, nc, tag, deps, strm)
    class(gs_crystal_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: j, c, co, ierr

    if (this%plan%nstage .eq. 0) return

    associate (st => this%plan%stage(1))
      co = (st%dst_sel - 1) * nc * this%plan%nwrk

      !$omp do
      do j = 1, st%nsw
         do c = 1, nc
            this%sbuf_v(nc*(j-1) + c) = &
                 u((c-1)*n + this%plan%pack_send_dof(j))
         end do
      end do
      !$omp end do

      !$omp master
      if (st%dst .ge. 0) then
         this%nsreq = 1
         call MPI_Isend(this%sbuf_v(1), nc*st%nsw, MPI_REAL_PRECISION, &
              st%dst, tag, NEKO_COMM, this%sreq(1), ierr)
      end if
      !$omp end master

      !$omp do
      do j = 1, st%nkw
         do c = 1, nc
            this%buf_v(co + nc*(j-1) + c) = &
                 u((c-1)*n + this%plan%pack_keep_dof(j))
         end do
      end do
      !$omp end do
    end associate

  end subroutine gs_crystal_nbsend_vec

  !> Drive the remaining routing stages and reduce what is delivered into
  !! the shared vector, fused nc-component
  subroutine gs_crystal_nbwait_vec(this, u, n, nc, op, strm)
    class(gs_crystal_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op
    integer :: s, i, j, c, co, so, fo, ierr

    if (this%plan%nstage .eq. 0) return

    !$omp master
    call MPI_Waitall(this%nrreq, this%rreq, MPI_STATUSES_IGNORE, ierr)
    call MPI_Waitall(this%nsreq, this%sreq, MPI_STATUSES_IGNORE, ierr)
    !$omp end master
    !$omp barrier

    do s = 2, this%plan%nstage
       associate (st => this%plan%stage(s))
         so = (st%src_sel - 1) * nc * this%plan%nwrk
         co = (st%dst_sel - 1) * nc * this%plan%nwrk

         !$omp master
         this%nrreq = 0
         this%nsreq = 0
         if (st%src .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf_v(co + nc*st%nkw + 1), nc*st%nrw, &
                 MPI_REAL_PRECISION, st%src, this%tag, NEKO_COMM, &
                 this%rreq(this%nrreq), ierr)
         end if
         if (st%src2 .ge. 0) then
            this%nrreq = this%nrreq + 1
            call MPI_Irecv(this%buf_v(co + nc*(st%nkw + st%nrw) + 1), &
                 nc*st%nr2w, MPI_REAL_PRECISION, st%src2, this%tag, &
                 NEKO_COMM, this%rreq(this%nrreq), ierr)
         end if
         !$omp end master

         if (st%dst .ge. 0) then
            !$omp do
            do j = 1, st%nsw
               do c = 1, nc
                  this%sbuf_v(nc*(j-1) + c) = &
                       this%buf_v(so + nc*(st%send_idx(j) - 1) + c)
               end do
            end do
            !$omp end do
            !$omp master
            this%nsreq = 1
            call MPI_Isend(this%sbuf_v(1), nc*st%nsw, MPI_REAL_PRECISION, &
                 st%dst, this%tag, NEKO_COMM, this%sreq(1), ierr)
            !$omp end master
         end if

         if (.not. st%inplace) then
            !$omp do
            do j = 1, st%nkw
               do c = 1, nc
                  this%buf_v(co + nc*(j-1) + c) = &
                       this%buf_v(so + nc*(st%keep_idx(j) - 1) + c)
               end do
            end do
            !$omp end do
         end if

         !$omp master
         call MPI_Waitall(this%nrreq, this%rreq, MPI_STATUSES_IGNORE, ierr)
         call MPI_Waitall(this%nsreq, this%sreq, MPI_STATUSES_IGNORE, ierr)
         !$omp end master
         !$omp barrier
       end associate
    end do

    fo = (this%plan%final_sel - 1) * nc * this%plan%nwrk
    do i = 1, this%plan%nfinal_rec
       associate (ro => this%plan%final_off(i), rl => this%plan%final_len(i))
         select case (op)
         case (GS_OP_ADD)
            !$omp do
            do j = 1, rl
               do c = 1, nc
                  u((c-1)*n + this%plan%unpack_dof(ro + j)) = &
                       u((c-1)*n + this%plan%unpack_dof(ro + j)) + &
                       this%buf_v(fo + nc*(ro + j - 1) + c)
               end do
            end do
            !$omp end do
         case (GS_OP_MUL)
            !$omp do
            do j = 1, rl
               do c = 1, nc
                  u((c-1)*n + this%plan%unpack_dof(ro + j)) = &
                       u((c-1)*n + this%plan%unpack_dof(ro + j)) * &
                       this%buf_v(fo + nc*(ro + j - 1) + c)
               end do
            end do
            !$omp end do
         case (GS_OP_MIN)
            !$omp do
            do j = 1, rl
               do c = 1, nc
                  u((c-1)*n + this%plan%unpack_dof(ro + j)) = &
                       min(u((c-1)*n + this%plan%unpack_dof(ro + j)), &
                       this%buf_v(fo + nc*(ro + j - 1) + c))
               end do
            end do
            !$omp end do
         case (GS_OP_MAX)
            !$omp do
            do j = 1, rl
               do c = 1, nc
                  u((c-1)*n + this%plan%unpack_dof(ro + j)) = &
                       max(u((c-1)*n + this%plan%unpack_dof(ro + j)), &
                       this%buf_v(fo + nc*(ro + j - 1) + c))
               end do
            end do
            !$omp end do
         case default
            call neko_error("Unknown operation in gs_crystal_nbwait_vec")
         end select
       end associate
    end do

  end subroutine gs_crystal_nbwait_vec

end module gs_crystal
