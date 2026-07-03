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
!> Crystal router: scalable all-to-some personalized exchange
!! @details
!! Routes a set of variable-length records, each tagged with a destination
!! rank, from any source to any destination in @f$ \lceil \log_2 P \rceil @f$
!! communication stages. This replaces dense @f$ O(P) @f$-round neighbour
!! discovery / personalized exchange (which is @f$ O(P^2) @f$ globally and,
!! under flat MPI at high rank counts, exhausts the MPI runtime's internal
!! unexpected-message buffers) with a recursive-bisection routing that talks
!! to at most one or two partners per stage.
!!
!! Records are passed in a single packed `integer(i8)` buffer as a sequence of
!! self-describing entries
!! @verbatim
!!     [ dest, len, payload(1), ..., payload(len) ]
!! @endverbatim
!! where `dest` is the destination rank in @f$ [0, P) @f$ and `len` is the
!! payload length. On return the buffer holds exactly the records whose `dest`
!! equals `pe_rank`, in the same packed format (the `dest`/`len` headers are
!! preserved so the caller can unpack uniformly).
!!
!! ### Algorithm and correctness
!! The active rank range `[lo, hi)` always contains `pe_rank`; it starts as
!! `[0, P)` and is halved each stage at `mid = lo + (hi-lo)/2`. A record is
!! kept if its destination lies in the local half and otherwise shipped to a
!! partner in the opposite half. Since a record is always moved into the half
!! containing its destination, the range invariant "`dest in [lo, hi)`" is
!! preserved, and when the range collapses to a single rank that rank is the
!! destination -- so every record is delivered exactly once and none is
!! dropped.
!!
!! For odd-sized ranges the upper half has exactly one more rank than the
!! lower half, leaving one unpaired upper rank (`hi-1`). That rank still holds
!! records destined for the lower half; it offloads them to rank `lo`, which
!! performs one extra receive. All other ranks exchange with a single partner
!! via `MPI_Sendrecv`, and the unpaired transfers use `MPI_PROC_NULL` for the
!! idle direction, so every stage is fully matched and deadlock-free for any
!! `P` (not only powers of two). Message sizes are negotiated per stage (a
!! one-integer count exchange precedes each data exchange), so no global
!! maximum buffer is ever allocated.
module crystal_router
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use num_types, only : i8
  use stack, only : stack_i8_t
  use utils, only : neko_error
  use mpi_f08, only : MPI_Sendrecv, MPI_Status, MPI_INTEGER, MPI_INTEGER8, &
       MPI_PROC_NULL
  implicit none
  private

  !> Message tag used for all crystal-router exchanges. Stages are fully
  !! synchronised (blocking `MPI_Sendrecv`), so a single tag suffices.
  integer, parameter :: CR_TAG = 0

  public :: crystal_router_transfer, crystal_router_pack

contains

  !> Append one record to a packed crystal-router buffer.
  !! @details Writes the self-describing entry `[dest, size(body), body...]`
  !! that crystal_router_transfer consumes. Owning the record format here means
  !! call sites only supply a destination rank and a payload, and never encode
  !! the `[dest, len, ...]` envelope themselves.
  !! @param out   Stack accumulating the packed records.
  !! @param dest  Destination rank for this record, in @f$ [0, P) @f$.
  !! @param body  The record payload (its length becomes the `len` field).
  subroutine crystal_router_pack(out, dest, body)
    type(stack_i8_t), intent(inout) :: out
    integer, intent(in) :: dest
    integer(i8), intent(in) :: body(:)
    integer(i8) :: w
    integer :: i

    w = int(dest, i8)
    call out%push(w)
    w = int(size(body), i8)
    call out%push(w)
    do i = 1, size(body)
       w = body(i)
       call out%push(w)
    end do
  end subroutine crystal_router_pack

  !> Route packed records to their destination ranks.
  !! @param buf  Packed records `[dest, len, payload...]` repeated. On entry,
  !!             arbitrary destinations in @f$ [0, P) @f$; on exit, only the
  !!             records addressed to `pe_rank` (may be reallocated).
  !! @param n    Used length of @a buf (number of `integer(i8)` words). Updated
  !!             on exit to the length of the received set.
  subroutine crystal_router_transfer(buf, n)
    integer(i8), allocatable, intent(inout) :: buf(:)
    integer, intent(inout) :: n
    integer(i8), allocatable :: keep(:), snd(:), rcv(:), extra(:)
    integer :: nkeep, nsnd, nrcv, nextra
    integer :: lo, hi, m, half, mid, r, partner
    logical :: lower

    ! Validate destinations up front so a packing bug surfaces here rather
    ! than as a silently misrouted (lost) record at scale.
    call cr_check_dest(buf, n)

    lo = 0
    hi = pe_size

    do while (hi - lo > 1)
       m = hi - lo
       half = m / 2
       mid = lo + half
       lower = (pe_rank < mid)

       ! Split the local buffer: records staying in our half vs. those bound
       ! for the opposite half.
       call cr_partition(buf, n, mid, lower, keep, nkeep, snd, nsnd)

       ! Pick the partner in the opposite half. Lower ranks always have one;
       ! in an odd-sized range the last upper rank (r == half) is unpaired.
       if (lower) then
          partner = mid + (pe_rank - lo)
       else
          r = pe_rank - mid
          if (r .lt. half) then
             partner = lo + r
          else
             partner = MPI_PROC_NULL
          end if
       end if

       ! Bidirectional exchange with the partner, then rebuild the local
       ! buffer as (kept records) ++ (received records).
       call cr_exchange(snd, nsnd, partner, rcv, nrcv, partner)
       call cr_concat(buf, n, keep, nkeep, rcv, nrcv)

       ! Odd-sized range: the unpaired upper rank (hi-1) hands its lower-half
       ! records to rank lo, which absorbs them with one extra receive.
       if (iand(m, 1) .eq. 1) then
          if ((.not. lower) .and. ((pe_rank - mid) .eq. half)) then
             call cr_exchange(snd, nsnd, lo, extra, nextra, MPI_PROC_NULL)
          else if (pe_rank .eq. lo) then
             call cr_exchange(snd, 0, MPI_PROC_NULL, extra, nextra, hi - 1)
             call cr_append(buf, n, extra, nextra)
          end if
       end if

       ! Narrow the active range to the half containing pe_rank.
       if (lower) then
          hi = mid
       else
          lo = mid
       end if
    end do

    ! Release the per-stage scratch buffers (buf is returned to the caller).
    if (allocated(keep)) deallocate(keep)
    if (allocated(snd)) deallocate(snd)
    if (allocated(rcv)) deallocate(rcv)
    if (allocated(extra)) deallocate(extra)

  end subroutine crystal_router_transfer

  !> Abort if any record destination falls outside @f$ [0, P) @f$.
  subroutine cr_check_dest(buf, n)
    integer(i8), allocatable, intent(in) :: buf(:)
    integer, intent(in) :: n
    integer :: p, dest, rlen

    p = 1
    do while (p .le. n)
       dest = int(buf(p))
       rlen = int(buf(p + 1))
       if (dest .lt. 0 .or. dest .ge. pe_size) then
          call neko_error('crystal_router: record destination out of range')
       end if
       p = p + 2 + rlen
    end do
  end subroutine cr_check_dest

  !> Partition @a buf into records kept locally and records to be sent.
  !! A record is kept when the side of its destination (relative to @a mid)
  !! matches @a lower, i.e. it already belongs to our half.
  subroutine cr_partition(buf, n, mid, lower, keep, nkeep, snd, nsnd)
    integer(i8), allocatable, intent(in) :: buf(:)
    integer, intent(in) :: n, mid
    logical, intent(in) :: lower
    integer(i8), allocatable, intent(out) :: keep(:), snd(:)
    integer, intent(out) :: nkeep, nsnd
    integer :: p, dest, rlen, reclen
    logical :: in_lower

    allocate(keep(max(n, 1)))
    allocate(snd(max(n, 1)))
    nkeep = 0
    nsnd = 0

    p = 1
    do while (p .le. n)
       dest = int(buf(p))
       rlen = int(buf(p + 1))
       reclen = 2 + rlen
       in_lower = dest .lt. mid
       if (in_lower .eqv. lower) then
          keep(nkeep + 1:nkeep + reclen) = buf(p:p + reclen - 1)
          nkeep = nkeep + reclen
       else
          snd(nsnd + 1:nsnd + reclen) = buf(p:p + reclen - 1)
          nsnd = nsnd + reclen
       end if
       p = p + reclen
    end do
  end subroutine cr_partition

  !> Size-negotiated bidirectional exchange with a partner. Either @a dst or
  !! @a src may be `MPI_PROC_NULL` to make that direction a no-op.
  subroutine cr_exchange(sbuf, sn, dst, rbuf, rn, src)
    integer(i8), intent(in) :: sbuf(:)
    integer, intent(in) :: sn, dst, src
    integer(i8), allocatable, intent(out) :: rbuf(:)
    integer, intent(out) :: rn
    type(MPI_Status) :: status
    integer :: ierr

    rn = 0
    call MPI_Sendrecv(sn, 1, MPI_INTEGER, dst, CR_TAG, &
         rn, 1, MPI_INTEGER, src, CR_TAG, NEKO_COMM, status, ierr)

    allocate(rbuf(max(rn, 1)))
    call MPI_Sendrecv(sbuf, sn, MPI_INTEGER8, dst, CR_TAG, &
         rbuf, rn, MPI_INTEGER8, src, CR_TAG, NEKO_COMM, status, ierr)
  end subroutine cr_exchange

  !> Rebuild @a buf as the concatenation @a a(1:na) ++ @a b(1:nb).
  subroutine cr_concat(buf, n, a, na, b, nb)
    integer(i8), allocatable, intent(inout) :: buf(:)
    integer, intent(out) :: n
    integer(i8), allocatable, intent(in) :: a(:), b(:)
    integer, intent(in) :: na, nb

    if (allocated(buf)) deallocate(buf)
    allocate(buf(max(na + nb, 1)))
    if (na .gt. 0) buf(1:na) = a(1:na)
    if (nb .gt. 0) buf(na + 1:na + nb) = b(1:nb)
    n = na + nb
  end subroutine cr_concat

  !> Append @a e(1:ne) to @a buf(1:n) in place.
  subroutine cr_append(buf, n, e, ne)
    integer(i8), allocatable, intent(inout) :: buf(:)
    integer, intent(inout) :: n
    integer(i8), allocatable, intent(in) :: e(:)
    integer, intent(in) :: ne
    integer(i8), allocatable :: tmp(:)

    if (ne .le. 0) return
    allocate(tmp(n + ne))
    if (n .gt. 0) tmp(1:n) = buf(1:n)
    tmp(n + 1:n + ne) = e(1:ne)
    call move_alloc(tmp, buf)
    n = n + ne
  end subroutine cr_append

end module crystal_router
