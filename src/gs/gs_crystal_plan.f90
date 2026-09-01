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
!> Routing plan for the crystal router gather-scatter comm. backends
!! @details
!! The crystal router (see @a crystal_router) delivers records to arbitrary
!! destinations in @f$ \lceil \log_2 P \rceil @f$ recursive-bisection stages,
!! talking to at most one partner per stage instead of to every peer at once.
!! For the halo exchange the peer set is fixed by the gather-scatter schedule,
!! so the routing is the same for every gs operation and can be worked out
!! once, here, rather than re-derived from the payload on every call.
!!
!! The algorithm is the crystal router of Fox et al., "Solving Problems on
!! Concurrent Processors, Volume 1", Prentice-Hall, 1988.
!!
!! What this module computes is that plan: for each stage, the partner ranks,
!! the exact word counts to send and receive, and the index lists that move
!! the words which stay put. The runtime is then a fixed sequence of
!! `Isend`/`Irecv` of known size plus indexed gathers -- no destination scan,
!! no size negotiation, and no allocation, none of which the setup-phase
!! router avoids because it cannot.
!!
!! ### What the aggregation buys, and what it costs
!! One gs operation goes from @f$ k @f$ messages to @f$ k @f$ peers to at most
!! one message per active stage. Because a partition's peers are close in rank
!! index, most stages have nothing to send or receive and are dropped from the
!! plan entirely, so the active stage count is set by the largest rank distance
!! among the peers rather than by @f$ \log_2 P @f$. The price is store and
!! forward: a word bound for a peer several stages away crosses the network
!! once per stage it survives, and is copied locally each time. This trades
!! bandwidth for message count and is worth it only where per-message overhead
!! dominates -- few elements per rank, coarse multigrid levels, or a runtime
!! that handles many concurrent peers badly. Which is why the backends built on
!! it are candidates for the runtime autotuning rather than a default.
!!
!! ### Record identity
!! A record is the whole slab of shared dofs one rank sends to one peer, and it
!! is routed intact. Its payload length and its originating rank are all a
!! relay needs; the destination's local dof indices are never shipped, since
!! the destination already knows them as `recv_dof(origin)`. Contributions are
!! therefore reduced at the destination only, exactly as in the pairwise
!! backends. Folding contributions to the same dof at intermediate hops would
!! cut the @f$ O(k^2) @f$ traffic of a multiplicity @f$ k @f$ vertex to
!! @f$ O(k \log k) @f$, but needs the destination's dof identity at the relay
!! and is deliberately left out here.
module gs_crystal_plan
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use stack, only : stack_i4_t
  use utils, only : neko_error
  use mpi_f08, only : MPI_Sendrecv, MPI_Status, MPI_INTEGER, MPI_PROC_NULL
  implicit none
  private

  !> Metadata words per in-flight record in the symbolic pass: destination
  !! rank, originating rank, payload length.
  integer, parameter :: CR_REC = 3

  !> Message tag for the symbolic pass. Setup only, and every exchange in it
  !! is a fully synchronised `MPI_Sendrecv`, so one tag suffices.
  integer, parameter :: CR_PLAN_TAG = 0

  !> One active communication stage of the routing plan.
  !! @details The working buffer is held in one of two columns; a stage reads
  !! column @a src_sel and leaves the result in column @a dst_sel, laid out as
  !! @verbatim
  !!   [ words that stay (nkw) | from src (nrw) | from src2 (nr2w) ]
  !! @endverbatim
  !! A stage that sends nothing keeps every word it holds, so it needs no copy
  !! and receives straight into the column it already occupies (@a inplace).
  type, public :: gs_crystal_stage_t
     !> Rank the packed send buffer goes to, -1 if this stage sends nothing
     integer :: dst = -1
     !> Paired partner to receive from, -1 if this stage receives nothing
     integer :: src = -1
     !> Unpaired top rank of an odd sized range offloading to us, -1 if none
     integer :: src2 = -1
     integer :: nsw = 0 !< Words sent to @a dst
     integer :: nrw = 0 !< Words received from @a src
     integer :: nr2w = 0 !< Words received from @a src2
     integer :: nkw = 0 !< Words that stay local over this stage
     integer :: src_sel = 1 !< Buffer column read by this stage
     integer :: dst_sel = 1 !< Buffer column written by this stage
     !> Whether the stage sends nothing, and so needs no copy and can
     !! receive into the column it already occupies
     logical :: inplace = .false.
     !> 1-based positions in column @a src_sel of the words that stay, in
     !! the order they take in column @a dst_sel. Unallocated on the first
     !! stage, which materialises them straight from the shared vector
     integer, allocatable :: keep_idx(:)
     !> 1-based positions in column @a src_sel of the words that are sent,
     !! in the order they take in the send buffer. Unallocated on the first
     !! stage, see @a keep_idx
     integer, allocatable :: send_idx(:)
  end type gs_crystal_stage_t

  !> The full routing plan for one gather-scatter schedule.
  type, public :: gs_crystal_plan_t
     !> The active stages, in order. Stages that neither send nor receive
     !! leave the buffer untouched and are not represented
     type(gs_crystal_stage_t), allocatable :: stage(:)
     integer :: nstage = 0 !< Number of active stages
     integer :: nwrk = 0 !< Words each working buffer column must hold
     integer :: nsmax = 0 !< Words the send buffer must hold
     integer :: ntotal = 0 !< Words packed out of the shared vector
     integer :: nfinal = 0 !< Words delivered to this rank
     integer :: final_sel = 1 !< Buffer column holding the delivered words
     !> 1-based indices into the shared vector of the words that stay local
     !! over the first stage, in the order they take in the working buffer
     integer, allocatable :: pack_keep_dof(:)
     !> 1-based indices into the shared vector of the words sent by the
     !! first stage, in the order they take in the send buffer
     integer, allocatable :: pack_send_dof(:)
     !> 1-based indices into the shared vector that the delivered words
     !! reduce into, in the order the words arrive
     integer, allocatable :: unpack_dof(:)
     !> Offsets into @a unpack_dof of each delivered record, 0-based
     integer, allocatable :: final_off(:)
     !> Lengths of each delivered record
     integer, allocatable :: final_len(:)
     integer :: nfinal_rec = 0 !< Number of delivered records
   contains
     procedure, pass(this) :: init => gs_crystal_plan_init
     procedure, pass(this) :: free => gs_crystal_plan_free
  end type gs_crystal_plan_t

contains

  !> Work out the routing plan for a gather-scatter schedule.
  !! @details Runs the recursive bisection of the crystal router over record
  !! metadata alone -- three integers per slab rather than the payload -- so
  !! the partner ranks, word counts and index lists of every stage are known
  !! before a single gs operation runs.
  !! @param send_pe ranks this process sends to
  !! @param recv_pe ranks this process receives from
  !! @param send_dof per-rank lists of shared dofs to send
  !! @param recv_dof per-rank lists of shared dofs to reduce into
  !! @note Collective over NEKO_COMM, as the bisection is.
  subroutine gs_crystal_plan_init(this, send_pe, recv_pe, send_dof, recv_dof)
    class(gs_crystal_plan_t), intent(inout) :: this
    integer, intent(in) :: send_pe(:), recv_pe(:)
    type(stack_i4_t), intent(inout) :: send_dof(0:), recv_dof(0:)
    type(gs_crystal_stage_t), allocatable :: st(:)
    logical, allocatable :: active(:)
    integer, allocatable :: rdest(:), rorig(:), rlen(:)
    integer, allocatable :: ndest(:), norig(:), nlen(:)
    integer, allocatable :: smeta(:), rmeta(:), rmeta2(:)
    integer, allocatable :: kidx(:), sidx(:), pack_dof(:)
    integer :: lo, hi, m, half, mid, r, partner
    integer :: s, nst, nstmax, sel, nrec, nnew, off, t
    integer :: i, j, k, nkw, nsw, nrw, nr2w, nsrec, nrrec, nr2rec
    integer :: dst, src, src2, curlen
    logical :: lower, keepit

    call this%free()

    ! The records this rank starts with: one per peer, holding that peer's
    ! whole slab of shared dofs
    nrec = size(send_pe)
    allocate(rdest(max(nrec, 1)), rorig(max(nrec, 1)), rlen(max(nrec, 1)))
    this%ntotal = 0
    do i = 1, nrec
       rdest(i) = send_pe(i)
       rorig(i) = pe_rank
       rlen(i) = send_dof(send_pe(i))%size()
       this%ntotal = this%ntotal + rlen(i)
    end do
    this%nwrk = this%ntotal

    ! The bisection visits ceil(log2(pe_size)) ranges; size the scratch plan
    ! for that many stages before any of them are dropped
    nstmax = 0
    m = pe_size
    do while (m .gt. 1)
       nstmax = nstmax + 1
       m = m - m/2
    end do
    allocate(st(max(nstmax, 1)), active(max(nstmax, 1)))
    active = .false.

    lo = 0
    hi = pe_size
    s = 0

    do while (hi - lo .gt. 1)
       s = s + 1
       m = hi - lo
       half = m / 2
       mid = lo + half
       lower = (pe_rank .lt. mid)

       ! Partner in the opposite half. Lower ranks always have one; in an
       ! odd sized range the last upper rank is unpaired and offloads its
       ! cross-half records to lo, which absorbs them with an extra receive
       if (lower) then
          partner = mid + (pe_rank - lo)
       else
          r = pe_rank - mid
          if (r .lt. half) then
             partner = lo + r
          else
             partner = -1
          end if
       end if

       dst = partner
       src = partner
       src2 = -1
       if (iand(m, 1) .eq. 1) then
          if ((.not. lower) .and. ((pe_rank - mid) .eq. half)) then
             dst = lo
          else if (pe_rank .eq. lo) then
             src2 = hi - 1
          end if
       end if

       ! Split the records into those already in our half and those bound
       ! for the other one, and note where their words sit in the buffer
       nkw = 0
       nsw = 0
       nsrec = 0
       do i = 1, nrec
          keepit = ((rdest(i) .lt. mid) .eqv. lower)
          if (keepit) then
             nkw = nkw + rlen(i)
          else
             nsw = nsw + rlen(i)
             nsrec = nsrec + 1
          end if
       end do

       allocate(kidx(max(nkw, 1)), sidx(max(nsw, 1)))
       allocate(smeta(max(CR_REC*nsrec, 1)))
       off = 0
       j = 0
       k = 0
       nsrec = 0
       do i = 1, nrec
          keepit = ((rdest(i) .lt. mid) .eqv. lower)
          if (keepit) then
             do t = 1, rlen(i)
                j = j + 1
                kidx(j) = off + t
             end do
          else
             do t = 1, rlen(i)
                k = k + 1
                sidx(k) = off + t
             end do
             smeta(CR_REC*nsrec + 1) = rdest(i)
             smeta(CR_REC*nsrec + 2) = rorig(i)
             smeta(CR_REC*nsrec + 3) = rlen(i)
             nsrec = nsrec + 1
          end if
          off = off + rlen(i)
       end do

       ! Hand the outgoing records' metadata to the partner and take in what
       ! it routes to us, then, in an odd sized range, the unpaired rank's
       call cr_meta_exchange(smeta, nsrec, dst, rmeta, nrrec, src)
       nr2rec = 0
       if (src2 .ge. 0) then
          call cr_meta_exchange(smeta, 0, -1, rmeta2, nr2rec, src2)
       else
          allocate(rmeta2(1))
       end if

       nrw = 0
       do i = 1, nrrec
          nrw = nrw + rmeta(CR_REC*(i-1) + 3)
       end do
       nr2w = 0
       do i = 1, nr2rec
          nr2w = nr2w + rmeta2(CR_REC*(i-1) + 3)
       end do

       ! Rebuild the record list as (kept) ++ (from src) ++ (from src2), the
       ! order the words take in the buffer after this stage
       nnew = (nrec - nsrec) + nrrec + nr2rec
       allocate(ndest(max(nnew, 1)), norig(max(nnew, 1)), nlen(max(nnew, 1)))
       j = 0
       do i = 1, nrec
          if ((rdest(i) .lt. mid) .eqv. lower) then
             j = j + 1
             ndest(j) = rdest(i)
             norig(j) = rorig(i)
             nlen(j) = rlen(i)
          end if
       end do
       do i = 1, nrrec
          j = j + 1
          ndest(j) = rmeta(CR_REC*(i-1) + 1)
          norig(j) = rmeta(CR_REC*(i-1) + 2)
          nlen(j) = rmeta(CR_REC*(i-1) + 3)
       end do
       do i = 1, nr2rec
          j = j + 1
          ndest(j) = rmeta2(CR_REC*(i-1) + 1)
          norig(j) = rmeta2(CR_REC*(i-1) + 2)
          nlen(j) = rmeta2(CR_REC*(i-1) + 3)
       end do

       ! A stage with nothing on the wire in either direction keeps every
       ! word where it is, and is dropped from the plan below
       if (nsw .eq. 0) dst = -1
       if (nrw .eq. 0) src = -1
       if (nr2w .eq. 0) src2 = -1

       st(s)%dst = dst
       st(s)%src = src
       st(s)%src2 = src2
       st(s)%nsw = nsw
       st(s)%nrw = nrw
       st(s)%nr2w = nr2w
       st(s)%nkw = nkw
       st(s)%inplace = (nsw .eq. 0)
       call move_alloc(kidx, st(s)%keep_idx)
       call move_alloc(sidx, st(s)%send_idx)
       active(s) = (dst .ge. 0) .or. (src .ge. 0) .or. (src2 .ge. 0)

       this%nsmax = max(this%nsmax, nsw)
       curlen = nkw + nrw + nr2w
       this%nwrk = max(this%nwrk, curlen)

       deallocate(smeta, rmeta, rmeta2)
       call move_alloc(ndest, rdest)
       call move_alloc(norig, rorig)
       call move_alloc(nlen, rlen)
       nrec = nnew

       if (lower) then
          hi = mid
       else
          lo = mid
       end if
    end do

    ! Everything left is addressed to this rank; its origin says which of
    ! the recv_dof lists it reduces into
    call cr_final_dofs(this, rdest, rorig, rlen, nrec, recv_pe, recv_dof)

    ! Keep only the stages that move something, and pin down which buffer
    ! column each of them reads and writes
    nst = 0
    do i = 1, s
       if (active(i)) nst = nst + 1
    end do
    allocate(this%stage(max(nst, 1)))
    this%nstage = nst
    j = 0
    do i = 1, s
       if (.not. active(i)) cycle
       j = j + 1
       this%stage(j)%dst = st(i)%dst
       this%stage(j)%src = st(i)%src
       this%stage(j)%src2 = st(i)%src2
       this%stage(j)%nsw = st(i)%nsw
       this%stage(j)%nrw = st(i)%nrw
       this%stage(j)%nr2w = st(i)%nr2w
       this%stage(j)%nkw = st(i)%nkw
       this%stage(j)%inplace = st(i)%inplace
       call move_alloc(st(i)%keep_idx, this%stage(j)%keep_idx)
       call move_alloc(st(i)%send_idx, this%stage(j)%send_idx)
    end do
    deallocate(st, active)

    ! The first stage materialises its words out of the shared vector, so it
    ! writes column 1 whether or not it would otherwise have swapped
    sel = 1
    if (nst .gt. 0) then
       this%stage(1)%src_sel = 1
       this%stage(1)%dst_sel = 1
       do j = 2, nst
          this%stage(j)%src_sel = sel
          if (this%stage(j)%inplace) then
             this%stage(j)%dst_sel = sel
          else
             this%stage(j)%dst_sel = 3 - sel
             sel = 3 - sel
          end if
       end do
    end if
    this%final_sel = sel

    ! Turn the first stage's buffer positions into shared vector indices, so
    ! its words are gathered straight out of the shared vector
    if (nst .gt. 0) then
       allocate(pack_dof(max(this%ntotal, 1)))
       j = 0
       do i = 1, size(send_pe)
          select type (sp => send_dof(send_pe(i))%data)
          type is (integer)
             do t = 1, send_dof(send_pe(i))%size()
                j = j + 1
                pack_dof(j) = sp(t)
             end do
          end select
       end do

       allocate(this%pack_keep_dof(max(this%stage(1)%nkw, 1)))
       do j = 1, this%stage(1)%nkw
          this%pack_keep_dof(j) = pack_dof(this%stage(1)%keep_idx(j))
       end do
       allocate(this%pack_send_dof(max(this%stage(1)%nsw, 1)))
       do j = 1, this%stage(1)%nsw
          this%pack_send_dof(j) = pack_dof(this%stage(1)%send_idx(j))
       end do
       deallocate(pack_dof)
       deallocate(this%stage(1)%keep_idx, this%stage(1)%send_idx)
    else
       allocate(this%pack_keep_dof(1), this%pack_send_dof(1))
    end if

    this%nwrk = max(this%nwrk, 1)
    this%nsmax = max(this%nsmax, 1)

    deallocate(rdest, rorig, rlen)

  end subroutine gs_crystal_plan_init

  !> Turn the records left at the end of the routing into the reduction
  !! order of the delivered words.
  !! @param rdest, rorig, rlen the delivered records
  !! @param nrec number of delivered records
  !! @param recv_pe ranks this process receives from
  !! @param recv_dof per-rank lists of shared dofs to reduce into
  subroutine cr_final_dofs(this, rdest, rorig, rlen, nrec, recv_pe, recv_dof)
    class(gs_crystal_plan_t), intent(inout) :: this
    integer, intent(in) :: rdest(:), rorig(:), rlen(:), nrec
    integer, intent(in) :: recv_pe(:)
    type(stack_i4_t), intent(inout) :: recv_dof(0:)
    integer :: i, j, t, nfin

    if (nrec .ne. size(recv_pe)) then
       call neko_error('gs_crystal_plan: routed record count does not ' // &
            'match the gather-scatter schedule')
    end if

    nfin = 0
    do i = 1, nrec
       if (rdest(i) .ne. pe_rank) then
          call neko_error('gs_crystal_plan: record left undelivered')
       end if
       if (rorig(i) .lt. 0 .or. rorig(i) .ge. pe_size) then
          call neko_error('gs_crystal_plan: record origin out of range')
       end if
       if (rlen(i) .ne. recv_dof(rorig(i))%size()) then
          call neko_error('gs_crystal_plan: routed record length does ' // &
               'not match the gather-scatter schedule')
       end if
       nfin = nfin + rlen(i)
    end do

    this%nfinal = nfin
    this%nfinal_rec = nrec
    this%nwrk = max(this%nwrk, nfin)
    allocate(this%unpack_dof(max(nfin, 1)))
    allocate(this%final_off(max(nrec, 1)), this%final_len(max(nrec, 1)))

    j = 0
    do i = 1, nrec
       this%final_off(i) = j
       this%final_len(i) = rlen(i)
       select type (rp_dof => recv_dof(rorig(i))%data)
       type is (integer)
          do t = 1, rlen(i)
             j = j + 1
             this%unpack_dof(j) = rp_dof(t)
          end do
       end select
    end do

  end subroutine cr_final_dofs

  !> Size-negotiated bidirectional exchange of record metadata with a
  !! partner. Either @a dst or @a src may be -1 to make that direction a
  !! no-op, which is how the unpaired ranks of an odd sized range are kept
  !! matched and deadlock free.
  !! @param sbuf packed metadata of the records to hand over
  !! @param nsrec number of records in @a sbuf
  !! @param dst rank to send to, or -1
  !! @param rbuf allocated here, holds the metadata taken in
  !! @param nrrec number of records in @a rbuf
  !! @param src rank to receive from, or -1
  subroutine cr_meta_exchange(sbuf, nsrec, dst, rbuf, nrrec, src)
    integer, intent(in) :: sbuf(:)
    integer, intent(in) :: nsrec, dst, src
    integer, allocatable, intent(out) :: rbuf(:)
    integer, intent(out) :: nrrec
    type(MPI_Status) :: status
    integer :: ierr, d, s

    d = MPI_PROC_NULL
    if (dst .ge. 0) d = dst
    s = MPI_PROC_NULL
    if (src .ge. 0) s = src

    nrrec = 0
    call MPI_Sendrecv(nsrec, 1, MPI_INTEGER, d, CR_PLAN_TAG, &
         nrrec, 1, MPI_INTEGER, s, CR_PLAN_TAG, NEKO_COMM, status, ierr)

    allocate(rbuf(max(CR_REC*nrrec, 1)))
    call MPI_Sendrecv(sbuf, CR_REC*nsrec, MPI_INTEGER, d, CR_PLAN_TAG, &
         rbuf, CR_REC*nrrec, MPI_INTEGER, s, CR_PLAN_TAG, NEKO_COMM, &
         status, ierr)

  end subroutine cr_meta_exchange

  !> Release the routing plan
  subroutine gs_crystal_plan_free(this)
    class(gs_crystal_plan_t), intent(inout) :: this
    integer :: i

    if (allocated(this%stage)) then
       do i = 1, size(this%stage)
          if (allocated(this%stage(i)%keep_idx)) &
               deallocate(this%stage(i)%keep_idx)
          if (allocated(this%stage(i)%send_idx)) &
               deallocate(this%stage(i)%send_idx)
       end do
       deallocate(this%stage)
    end if

    if (allocated(this%pack_keep_dof)) deallocate(this%pack_keep_dof)
    if (allocated(this%pack_send_dof)) deallocate(this%pack_send_dof)
    if (allocated(this%unpack_dof)) deallocate(this%unpack_dof)
    if (allocated(this%final_off)) deallocate(this%final_off)
    if (allocated(this%final_len)) deallocate(this%final_len)

    this%nstage = 0
    this%nwrk = 0
    this%nsmax = 0
    this%ntotal = 0
    this%nfinal = 0
    this%nfinal_rec = 0
    this%final_sel = 1

  end subroutine gs_crystal_plan_free

end module gs_crystal_plan
