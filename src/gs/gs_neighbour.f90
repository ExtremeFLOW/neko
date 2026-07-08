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
!> Defines gather-scatter communication using MPI neighbourhood collectives
module gs_neighbour
  use num_types, only : rp
  use gs_comm, only : gs_comm_t, GS_VEC_NC
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use stack, only : stack_i4_t
  use mpi_f08, only : MPI_Comm, MPI_Request, MPI_STATUS_IGNORE, &
       MPI_INFO_NULL, MPI_Dist_graph_create_adjacent, &
       MPI_Ineighbor_alltoallv, MPI_Wait, MPI_Comm_free
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION
  use, intrinsic :: iso_c_binding
  use utils, only : neko_error
  implicit none
  private

  !> Gather-scatter communication using an MPI neighbourhood collective.
  !! The whole halo exchange is carried out by a single
  !! `MPI_Ineighbor_alltoallv` over a distributed-graph communicator built
  !! once from the send/recv peer lists. Compared to the point-to-point
  !! `gs_mpi_t` this trades N `Isend`/`Irecv` pairs for one collective call:
  !! a single entry into the (serialised, on Fugaku) MPI runtime, leaving
  !! the library free to stripe the transfers across all network interfaces
  !! internally.
  type, public, extends(gs_comm_t) :: gs_neighbour_t
     !> Concatenated send slabs, one per peer, in send_pe order.
     real(kind=rp), allocatable :: send_buf(:)
     !> Concatenated recv slabs, one per peer, in recv_pe order.
     real(kind=rp), allocatable :: recv_buf(:)
     !> Per-peer counts and 0-based element displacements for the
     !! neighbourhood collective. send_* follow send_pe (the graph
     !! destinations), recv_* follow recv_pe (the graph sources).
     integer, allocatable :: sendcounts(:), sdispls(:)
     integer, allocatable :: recvcounts(:), rdispls(:)
     !> Distributed-graph communicator over the halo neighbourhood.
     type(MPI_Comm) :: neigh_comm
     !> In-flight non-blocking neighbourhood collective.
     type(MPI_Request) :: request
     !> Fused vector (multi-component) exchange buffers, sized for up to
     !! GS_VEC_NC components. Each peer slab holds nc consecutive component
     !! blocks. sendcounts_v/sdispls_v etc. are the nc-scaled collective
     !! descriptors, filled per call.
     real(kind=rp), allocatable :: send_buf_v(:)
     real(kind=rp), allocatable :: recv_buf_v(:)
     integer, allocatable :: sendcounts_v(:), sdispls_v(:)
     integer, allocatable :: recvcounts_v(:), rdispls_v(:)
   contains
     procedure, pass(this) :: init => gs_neighbour_init
     procedure, pass(this) :: free => gs_neighbour_free
     !> See gs_comm.f90 for more details on these routines
     procedure, pass(this) :: nbsend => gs_nbsend_neighbour
     procedure, pass(this) :: nbrecv => gs_nbrecv_neighbour
     procedure, pass(this) :: nbwait => gs_nbwait_neighbour
     procedure, pass(this) :: nbsend_vec => gs_nbsend_vec_neighbour
     procedure, pass(this) :: nbrecv_vec => gs_nbrecv_vec_neighbour
     procedure, pass(this) :: nbwait_vec => gs_nbwait_vec_neighbour
  end type gs_neighbour_t

contains

  !> Initialise the neighbourhood-collective communication method
  !! See gs_comm.f90 for details
  subroutine gs_neighbour_init(this, send_pe, recv_pe)
    class(gs_neighbour_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer :: i, nsend, nrecv, send_total, recv_total, ierr
    !> Unit weights for the graph. MPI_UNWEIGHTED is the natural choice but
    !! trips frt's mpi_f08 generic resolution; with reorder = .false. the
    !! weights never affect the neighbour ordering, so all-ones is equivalent.
    integer, allocatable :: src_weights(:), dst_weights(:)

    call this%init_order(send_pe, recv_pe)

    nsend = size(this%send_pe)
    nrecv = size(this%recv_pe)

    ! Allocate the count/displacement arrays with at least one element so
    ! that a peerless rank (nsend or nrecv == 0) still passes a valid array
    ! to the collective.
    allocate(this%sendcounts(max(1, nsend)), this%sdispls(max(1, nsend)))
    allocate(this%recvcounts(max(1, nrecv)), this%rdispls(max(1, nrecv)))
    this%sendcounts = 0
    this%sdispls = 0
    this%recvcounts = 0
    this%rdispls = 0

    send_total = 0
    do i = 1, nsend
       this%sendcounts(i) = this%send_dof(this%send_pe(i))%size()
       this%sdispls(i) = send_total
       send_total = send_total + this%sendcounts(i)
    end do
    allocate(this%send_buf(max(1, send_total)))

    recv_total = 0
    do i = 1, nrecv
       this%recvcounts(i) = this%recv_dof(this%recv_pe(i))%size()
       this%rdispls(i) = recv_total
       recv_total = recv_total + this%recvcounts(i)
    end do
    allocate(this%recv_buf(max(1, recv_total)))

    ! Fused vector exchange buffers/descriptors, sized for GS_VEC_NC comps.
    allocate(this%send_buf_v(max(1, GS_VEC_NC*send_total)))
    allocate(this%recv_buf_v(max(1, GS_VEC_NC*recv_total)))
    allocate(this%sendcounts_v(max(1, nsend)), this%sdispls_v(max(1, nsend)))
    allocate(this%recvcounts_v(max(1, nrecv)), this%rdispls_v(max(1, nrecv)))
    this%sendcounts_v = 0
    this%sdispls_v = 0
    this%recvcounts_v = 0
    this%rdispls_v = 0
    this%vec_supported = .true.

    ! Build a distributed-graph communicator over the halo neighbourhood.
    ! With MPI_Dist_graph_create_adjacent and reorder = .false. the order of
    ! sources/destinations is preserved, so the j-th block of the collective
    ! corresponds to recv_pe(j)/send_pe(j) and the count/displacement arrays
    ! above can be used directly.
    allocate(src_weights(max(1, nrecv)), dst_weights(max(1, nsend)))
    src_weights = 1
    dst_weights = 1
    call MPI_Dist_graph_create_adjacent(NEKO_COMM, &
         nrecv, this%recv_pe, src_weights, &
         nsend, this%send_pe, dst_weights, &
         MPI_INFO_NULL, .false., this%neigh_comm, ierr)
    deallocate(src_weights, dst_weights)

  end subroutine gs_neighbour_init

  !> Deallocate the neighbourhood-collective communication method
  subroutine gs_neighbour_free(this)
    class(gs_neighbour_t), intent(inout) :: this
    integer :: ierr

    if (allocated(this%send_buf)) deallocate(this%send_buf)
    if (allocated(this%recv_buf)) deallocate(this%recv_buf)
    if (allocated(this%sendcounts)) deallocate(this%sendcounts)
    if (allocated(this%sdispls)) deallocate(this%sdispls)
    if (allocated(this%recvcounts)) deallocate(this%recvcounts)
    if (allocated(this%rdispls)) deallocate(this%rdispls)

    if (allocated(this%send_buf_v)) deallocate(this%send_buf_v)
    if (allocated(this%recv_buf_v)) deallocate(this%recv_buf_v)
    if (allocated(this%sendcounts_v)) deallocate(this%sendcounts_v)
    if (allocated(this%sdispls_v)) deallocate(this%sdispls_v)
    if (allocated(this%recvcounts_v)) deallocate(this%recvcounts_v)
    if (allocated(this%rdispls_v)) deallocate(this%rdispls_v)

    call MPI_Comm_free(this%neigh_comm, ierr)

    call this%free_order()
    call this%free_dofs()

  end subroutine gs_neighbour_free

  !> Pack the send buffer and initiate the neighbourhood collective.
  !! The collective covers both directions, so all of the communication is
  !! launched here; gs_nbrecv_neighbour is a no-op.
  subroutine gs_nbsend_neighbour(this, u, n, tag, deps, strm)
    class(gs_neighbour_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, dst, off, ndst, ierr

    ! Gather data from u into the per-peer slab of send_buf according to the
    ! indices in send_dof. Each thread packs a different peer's slab; the
    ! implicit barrier at end-do guarantees the buffer is complete before
    ! the master fires the collective.
    !$omp do
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       off = this%sdispls(i)
       ndst = this%sendcounts(i)
       select type (dofs => this%send_dof(dst)%data)
       type is (integer)
          !$omp simd
          do j = 1, ndst
             this%send_buf(off + j) = u(dofs(j))
          end do
       end select
    end do
    !$omp end do

    ! A neighbourhood collective is one call over neigh_comm and must be
    ! issued by a single thread.
    !$omp master
    call MPI_Ineighbor_alltoallv(this%send_buf, this%sendcounts, &
         this%sdispls, MPI_REAL_PRECISION, this%recv_buf, this%recvcounts, &
         this%rdispls, MPI_REAL_PRECISION, this%neigh_comm, this%request, &
         ierr)
    !$omp end master

  end subroutine gs_nbsend_neighbour

  !> No-op: the neighbourhood collective issued in nbsend handles both the
  !! send and receive directions, so there is nothing to pre-post here.
  subroutine gs_nbrecv_neighbour(this, tag)
    class(gs_neighbour_t), intent(inout) :: this
    integer, intent(in) :: tag
    ! Nothing to do: the collective is launched in nbsend.
  end subroutine gs_nbrecv_neighbour

  !> Wait for the neighbourhood collective and reduce the received slabs
  subroutine gs_nbwait_neighbour(this, u, n, op, strm)
    class(gs_neighbour_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, src, off, nsrc
    integer :: op
    integer :: ierr

    ! The collective is master-issued, so wait on it from the master and
    ! barrier before anyone touches recv_buf.
    !$omp master
    call MPI_Wait(this%request, MPI_STATUS_IGNORE, ierr)
    !$omp end master
    !$omp barrier

    ! Reduce each received slab into u. The outer loop over peers stays
    ! serial: a dof shared by 3+ ranks appears in several recv_dof lists, so
    ! reducing two slabs concurrently would race on that dof. Parallelism is
    ! taken within each slab instead.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       off = this%rdispls(i)
       nsrc = this%recvcounts(i)
       select type (dofs => this%recv_dof(src)%data)
       type is (integer)
          select case (op)
          case (GS_OP_ADD)
             !OCL NORECURRENCE, NOVREC, NOALIAS
             !DIR$ CONCURRENT
             !DIR$ IVDEP
             !GCC$ ivdep
             !NEC$ IVDEP
             !$omp do
             do j = 1, nsrc
                u(dofs(j)) = u(dofs(j)) + this%recv_buf(off + j)
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
                u(dofs(j)) = u(dofs(j)) * this%recv_buf(off + j)
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
                u(dofs(j)) = min(u(dofs(j)), this%recv_buf(off + j))
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
                u(dofs(j)) = max(u(dofs(j)), this%recv_buf(off + j))
             end do
             !$omp end do
          case default
             call neko_error("Unknown operation in gs_nbwait_neighbour")
          end select
       end select
    end do

  end subroutine gs_nbwait_neighbour

  !> Pack the nc components and launch a single neighbourhood collective.
  !! @param u compact shared buffer, component-outer: u((c-1)*n + idx).
  !! Peer i's slab in send_buf_v holds nc consecutive blocks of sendcounts(i)
  !! values, so the collective descriptors are simply nc times the scalar ones.
  subroutine gs_nbsend_vec_neighbour(this, u, n, nc, tag, deps, strm)
    class(gs_neighbour_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    integer, intent(in) :: tag
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, c, dst, off, ndst, ierr

    !$omp do
    do i = 1, size(this%send_pe)
       dst = this%send_pe(i)
       off = this%sdispls(i)
       ndst = this%sendcounts(i)
       select type (dofs => this%send_dof(dst)%data)
       type is (integer)
          do c = 1, nc
             !$omp simd
             do j = 1, ndst
                this%send_buf_v(nc*off + (c-1)*ndst + j) = u((c-1)*n + dofs(j))
             end do
          end do
       end select
    end do
    !$omp end do

    ! A neighbourhood collective is one call over neigh_comm and must be
    ! issued by a single thread. The end-do barrier above guarantees the
    ! send buffer is fully packed first.
    !$omp master
    do i = 1, size(this%send_pe)
       this%sendcounts_v(i) = nc * this%sendcounts(i)
       this%sdispls_v(i) = nc * this%sdispls(i)
    end do
    do i = 1, size(this%recv_pe)
       this%recvcounts_v(i) = nc * this%recvcounts(i)
       this%rdispls_v(i) = nc * this%rdispls(i)
    end do
    call MPI_Ineighbor_alltoallv(this%send_buf_v, this%sendcounts_v, &
         this%sdispls_v, MPI_REAL_PRECISION, this%recv_buf_v, &
         this%recvcounts_v, this%rdispls_v, MPI_REAL_PRECISION, &
         this%neigh_comm, this%request, ierr)
    !$omp end master

  end subroutine gs_nbsend_vec_neighbour

  !> No-op: the collective is launched in nbsend_vec.
  subroutine gs_nbrecv_vec_neighbour(this, tag, nc)
    class(gs_neighbour_t), intent(inout) :: this
    integer, intent(in) :: tag, nc
    ! Nothing to do: the collective is launched in nbsend_vec.
  end subroutine gs_nbrecv_vec_neighbour

  !> Wait for the vector collective and reduce each received slab into u.
  !! Peer i's received slab holds nc consecutive blocks of recvcounts(i)
  !! values, matching the sender's packing convention.
  subroutine gs_nbwait_vec_neighbour(this, u, n, nc, op, strm)
    class(gs_neighbour_t), intent(inout) :: this
    integer, intent(in) :: n, nc
    real(kind=rp), dimension(nc*n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: i, j, c, src, off, nsrc, ierr
    integer :: op

    !$omp master
    call MPI_Wait(this%request, MPI_STATUS_IGNORE, ierr)
    !$omp end master
    !$omp barrier

    ! Serial over peers (a dof shared by 3+ ranks appears in several recv
    ! lists); parallelism is taken within each slab.
    do i = 1, size(this%recv_pe)
       src = this%recv_pe(i)
       off = this%rdispls(i)
       nsrc = this%recvcounts(i)
       select type (dofs => this%recv_dof(src)%data)
       type is (integer)
          select case (op)
          case (GS_OP_ADD)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + dofs(j)) = u((c-1)*n + dofs(j)) + &
                        this%recv_buf_v(nc*off + (c-1)*nsrc + j)
                end do
             end do
             !$omp end do
          case (GS_OP_MUL)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + dofs(j)) = u((c-1)*n + dofs(j)) * &
                        this%recv_buf_v(nc*off + (c-1)*nsrc + j)
                end do
             end do
             !$omp end do
          case (GS_OP_MIN)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + dofs(j)) = min(u((c-1)*n + dofs(j)), &
                        this%recv_buf_v(nc*off + (c-1)*nsrc + j))
                end do
             end do
             !$omp end do
          case (GS_OP_MAX)
             !$omp do
             do j = 1, nsrc
                do c = 1, nc
                   u((c-1)*n + dofs(j)) = max(u((c-1)*n + dofs(j)), &
                        this%recv_buf_v(nc*off + (c-1)*nsrc + j))
                end do
             end do
             !$omp end do
          case default
             call neko_error("Unknown operation in gs_nbwait_vec_neighbour")
          end select
       end select
    end do

  end subroutine gs_nbwait_vec_neighbour

end module gs_neighbour
