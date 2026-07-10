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
!     _  __  ____  __ __  ____
!    / |/ / / __/ / //_/ / __ \
!   /    / / _/  / ,<   / /_/ /
!  /_/|_/ /___/ /_/|_|  \____/
!
!> genmap_light -- a lightweight, CPU-only spectral-bisection mesh partitioner
!! for Neko `.nmsh` files.
!!
!! This is a standalone reimplementation of the partitioning algorithm in
!! Nek5000's `genmap`: recursive spectral bisection (RSB) of the element dual
!! graph using the Fiedler vector of the graph Laplacian, computed with a
!! restarted Lanczos iteration. Like `contrib/prepart` it then reorders the mesh
!! so that a later *linear* read reproduces the partition -- it reads an `.nmsh`,
!! partitions into `nparts` parts, and writes a new `.nmsh` whose elements are
!! grouped into contiguous per-partition blocks (curves and zones carried
!! through). It does **no Neko initialisation** and needs only a Fortran compiler
!! (no MPI, no Neko library, no ParMETIS), so -- like the other *_light tools --
!! it runs anywhere and needs no separate CPU build of a GPU-configured Neko.
!!
!! genmap vs ParMETIS: genmap's spectral bisection is lighter weight and has no
!! external dependency; ParMETIS (used by Neko's `prepart` / load balancer) is
!! multilevel k-way and often gives a slightly lower edge cut, but pulls in an
!! MPI library. genmap_light is the drop-in, dependency-free alternative.
!!
!! Usage (mirrors prepart):
!!   genmap_light mesh.nmsh nparts [out.nmsh]
!! Default output is <base>_<nparts>.nmsh (same convention as prepart).
!!
!! Tested against the meshes shipped with Neko and cross-checked against an
!! independent (SciPy) Fiedler eigensolver, but it remains your responsibility to
!! confirm the result is correct for your own mesh -- `make check` on a case you
!! trust first. See genmap_light.md for the algorithm, validation, and scope.
!!
!! 3D hex meshes only (genmap's common case); 2D exits with a clear error.

!> Compact open-addressing int32 -> int32 hash, used for the periodic point
!! merge (raw id -> merged id) and for compressing point ids to 1..nrnk.
module gm_hash
  use, intrinsic :: iso_fortran_env, only: i4 => int32, i8 => int64
  implicit none
  private
  public :: imap_t

  type :: imap_t
     integer(i4), allocatable :: key(:), val(:)
     integer(i8) :: cap = 0, mask = 0, n = 0
   contains
     procedure :: init => im_init
     procedure :: set => im_set !< insert/overwrite key -> val
     procedure :: getor => im_getor !< val if present, else the supplied default
     procedure :: intern => im_intern !< map key to a dense 1..n id, assigning on first sight
     procedure :: free => im_free
  end type imap_t

contains

  pure function mix(k) result(h)
    integer(i4), intent(in) :: k
    integer(i8) :: h
    h = int(k, i8)
    h = ieor(h, ishft(h, -33)) * (-49064778989728563_i8)
    h = ieor(h, ishft(h, -29)) * (-4417276706812531889_i8)
    h = ieor(h, ishft(h, -32))
  end function mix

  subroutine im_init(this, cap0)
    class(imap_t), intent(inout) :: this
    integer(i8), intent(in) :: cap0
    integer(i8) :: c
    c = 1024
    do while (c < cap0 * 2_i8)
       c = ishft(c, 1)
    end do
    this%cap = c; this%mask = c - 1; this%n = 0
    allocate(this%key(0:c-1), this%val(0:c-1)); this%key = 0
  end subroutine im_init

  subroutine im_set(this, k, v)
    class(imap_t), intent(inout) :: this
    integer(i4), intent(in) :: k, v
    integer(i8) :: h
    if (this%cap == 0) call im_init(this, 1024_i8)
    if (this%n * 3_i8 > this%cap * 2_i8) call im_grow(this)
    h = iand(mix(k), this%mask)
    do
       if (this%key(h) == 0) then
          this%key(h) = k; this%val(h) = v; this%n = this%n + 1; return
       else if (this%key(h) == k) then
          this%val(h) = v; return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end subroutine im_set

  function im_getor(this, k, default) result(v)
    class(imap_t), intent(in) :: this
    integer(i4), intent(in) :: k, default
    integer(i4) :: v
    integer(i8) :: h
    if (this%cap == 0) then; v = default; return; end if
    h = iand(mix(k), this%mask)
    do
       if (this%key(h) == 0) then
          v = default; return
       else if (this%key(h) == k) then
          v = this%val(h); return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end function im_getor

  !> Return a dense id in 1..(#distinct keys) for k, assigning the next free id
  !! the first time k is seen. `next` is the running counter (caller-owned).
  function im_intern(this, k, next) result(id)
    class(imap_t), intent(inout) :: this
    integer(i4), intent(in) :: k
    integer(i4), intent(inout) :: next
    integer(i4) :: id
    integer(i8) :: h
    if (this%cap == 0) call im_init(this, 1024_i8)
    if (this%n * 3_i8 > this%cap * 2_i8) call im_grow(this)
    h = iand(mix(k), this%mask)
    do
       if (this%key(h) == 0) then
          next = next + 1
          this%key(h) = k; this%val(h) = next; this%n = this%n + 1
          id = next; return
       else if (this%key(h) == k) then
          id = this%val(h); return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end function im_intern

  subroutine im_grow(this)
    class(imap_t), intent(inout) :: this
    integer(i4), allocatable :: ok(:), ov(:)
    integer(i8) :: i, h, oldcap
    oldcap = this%cap
    call move_alloc(this%key, ok); call move_alloc(this%val, ov)
    this%cap = ishft(this%cap, 1); this%mask = this%cap - 1
    allocate(this%key(0:this%cap-1), this%val(0:this%cap-1)); this%key = 0
    do i = 0, oldcap - 1
       if (ok(i) /= 0) then
          h = iand(mix(ok(i)), this%mask)
          do while (this%key(h) /= 0)
             h = iand(h + 1_i8, this%mask)
          end do
          this%key(h) = ok(i); this%val(h) = ov(i)
       end if
    end do
    deallocate(ok, ov)
  end subroutine im_grow

  subroutine im_free(this)
    class(imap_t), intent(inout) :: this
    if (allocated(this%key)) deallocate(this%key, this%val)
    this%cap = 0; this%mask = 0; this%n = 0
  end subroutine im_free

end module gm_hash


!> The genmap numerics: element dual-graph Laplacian (weight = # shared
!! vertices, exactly as Nek5000 genmap's c2c), the Fiedler vector via a restarted
!! Lanczos iteration on that Laplacian (glanczos/calcvec/lanczos2), a
!! connectivity (DFS) test, and the recursive spectral bisection with
!! connectivity-preserving fallbacks. Faithful to genmap.f.
module gm_rsb
  use, intrinsic :: iso_fortran_env, only: dp => real64, i4 => int32, i8 => int64, &
       error_unit
  implicit none
  private
  public :: rsb_partition

  integer, parameter :: MLANC = 100 ! Krylov dimension for the Fiedler Lanczos
  logical :: dbg = .false. ! GENMAP_DEBUG diagnostic prints

contains

  !> Partition elements 1..nelv into `nparts` parts via recursive spectral
  !! bisection. cell(nv,nelv) holds the (periodic-merged, compressed) vertex ids
  !! of each element; nrnk = number of distinct vertex ids. part(1:nelv) receives
  !! 0-based partition ids.
  subroutine rsb_partition(cell, nv, nelv, nrnk, nparts, part)
    integer(i4), intent(in) :: cell(nv, nelv)
    integer, intent(in) :: nv, nelv, nrnk, nparts
    integer(i4), intent(out) :: part(nelv)
    integer(i4), allocatable :: loc(:)
    integer :: i, elen
    character(len=8) :: dbgs
    call get_environment_variable('GENMAP_DEBUG', dbgs, elen)
    dbg = (elen > 0)
    allocate(loc(nelv))
    do i = 1, nelv
       loc(i) = i
    end do
    part = 0
    call rsb_rec(loc, nelv, nparts, 0, cell, nv, nelv, nrnk, part)
    deallocate(loc)
  end subroutine rsb_partition

  !> Recursively bisect the element subset loc(1:n) into q parts, numbering them
  !! base .. base+q-1.
  recursive subroutine rsb_rec(loc, n, q, base, cell, nv, nelv, nrnk, part)
    integer(i4), intent(inout) :: loc(n)
    integer, intent(in) :: n, q, base, nv, nelv, nrnk
    integer(i4), intent(in) :: cell(nv, nelv)
    integer(i4), intent(inout) :: part(nelv)
    integer :: q1, q2, n1, i
    if (q <= 1 .or. n == 0) then
       do i = 1, n
          part(loc(i)) = base
       end do
       return
    end if
    q1 = q / 2
    q2 = q - q1
    ! proportional target so all final parts differ by <= 1 element
    n1 = nint(real(n, dp) * real(q1, dp) / real(q, dp))
    if (n1 < 1) n1 = 1
    if (n1 > n - 1) n1 = n - 1
    ! bisect: reorder loc(1:n) so the first n1 form one side, the rest the other
    call bisect(loc, n, n1, cell, nv, nelv, nrnk)
    call rsb_rec(loc(1), n1, q1, base, cell, nv, nelv, nrnk, part)
    call rsb_rec(loc(n1+1), n - n1, q2, base + q1, cell, nv, nelv, nrnk, part)
  end subroutine rsb_rec

  !> Reorder loc(1:n) so loc(1:n1) is one half and loc(n1+1:n) the other, using
  !! spectral bisection. It computes the KWANT lowest spectral modes of the
  !! subset dual-graph Laplacian and, for each, forms the balanced split at rank
  !! n1 and measures the resulting edge cut; it keeps the lowest-cut split, giving
  !! preference to splits whose two halves are both connected. If none is
  !! connected it region-grows the best split to restore connectivity. Trying the
  !! few lowest modes (not just the single Fiedler vector) makes the bisection
  !! robust on meshes with a degenerate/near-degenerate Fiedler space (cylinders,
  !! symmetric boxes) as well as high-aspect-ratio/periodic meshes.
  subroutine bisect(loc, n, n1, cell, nv, nelv, nrnk)
    integer(i4), intent(inout) :: loc(n)
    integer, intent(in) :: n, n1, nv, nelv, nrnk
    integer(i4), intent(in) :: cell(nv, nelv)
    integer(i4), allocatable :: ia(:), ja(:)
    real(dp), allocatable :: va(:), fk(:,:)
    integer(i4), allocatable :: perm(:), side(:), bside(:), tmp(:)
    integer :: i, c, kwant, kgot
    real(dp) :: cut, score, bscore, bcut
    logical :: okconn, bconn

    if (n <= 2) return ! loc order is fine; first n1 vs rest

    call build_dual(loc, n, cell, nv, nelv, nrnk, ia, ja, va)

    kwant = min(4, n - 1)
    allocate(fk(n, kwant), perm(n), side(n), bside(n), tmp(n))
    call fiedler_k(ia, ja, va, n, kwant, fk, kgot)

    bscore = huge(1.0_dp); bcut = 0.0_dp; bconn = .false.; bside = 2
    do c = 1, max(kgot, 1)
       call rsort(fk(1, c), perm, n)
       side = 2
       do i = 1, n1
          side(perm(i)) = 1
       end do
       cut = split_cut(side, ia, ja, va, n)
       okconn = conn_side(side, 1, n, ia, ja) .and. conn_side(side, 2, n, ia, ja)
       ! prefer connected splits, then lowest cut
       score = cut
       if (.not. okconn) score = cut + 1.0e12_dp
       if (score < bscore) then
          bscore = score; bcut = cut; bside = side; bconn = okconn
       end if
    end do

    if (.not. bconn) then
       call region_grow(n, n1, bside, ia, ja) ! guarantee connectivity
    end if
    ! (print bcut, never bscore: the +1e12 disconnected penalty would overflow
    ! a default integer)
    if (dbg) write(error_unit,'(a,i0,a,i0,a,l1,a,i0)') &
         '   [dbg] bisect n=', n, ' n1=', n1, ' connected=', bconn, &
         ' cut=', int(bcut, i8)

    call apply_side(loc, n, bside, tmp)

    deallocate(ia, ja, va, fk, perm, side, bside, tmp)
  end subroutine bisect

  !> Weighted edge cut of a bipartition `side` (1/2) over the CSR graph (ia,ja,va);
  !! sums -va (the shared-vertex weight) over off-diagonal edges that cross, /2.
  function split_cut(side, ia, ja, va, n) result(cut)
    integer, intent(in) :: n
    integer(i4), intent(in) :: side(n), ia(0:n), ja(:)
    real(dp), intent(in) :: va(:)
    real(dp) :: cut
    integer :: i, j, jc
    cut = 0.0_dp
    do i = 1, n
       do j = int(ia(i-1)) + 1, int(ia(i))
          jc = ja(j)
          if (jc /= i .and. side(jc) /= side(i)) cut = cut - va(j)
       end do
    end do
    cut = cut * 0.5_dp
  end function split_cut

  !> Rearrange loc so all side==1 come first (in original relative order), then
  !! side==2. Two-pass stable partition via scratch tmp(n).
  subroutine apply_side(loc, n, side, tmp)
    integer(i4), intent(inout) :: loc(n)
    integer, intent(in) :: n
    integer(i4), intent(in) :: side(n)
    integer(i4), intent(inout) :: tmp(n)
    integer :: i, k
    k = 0
    do i = 1, n
       if (side(i) == 1) then
          k = k + 1; tmp(k) = loc(i)
       end if
    end do
    do i = 1, n
       if (side(i) == 2) then
          k = k + 1; tmp(k) = loc(i)
       end if
    end do
    loc = tmp
  end subroutine apply_side

  !> Build the dual-graph Laplacian (CSR: ia(0:n), ja, va) for the subset
  !! loc(1:n). Off-diagonal va = -(number of shared vertices); diagonal va =
  !! +sum of shared-vertex counts to neighbours -- the weighted graph Laplacian,
  !! exactly Nek5000 genmap's c2c weighting.
  subroutine build_dual(loc, n, cell, nv, nelv, nrnk, ia, ja, va)
    integer(i4), intent(in) :: loc(n)
    integer, intent(in) :: n, nv, nelv, nrnk
    integer(i4), intent(in) :: cell(nv, nelv)
    integer(i4), allocatable, intent(out) :: ia(:), ja(:)
    real(dp), allocatable, intent(out) :: va(:)
    integer(i4), allocatable :: v2c_cnt(:), v2c_ptr(:), v2c(:)
    integer(i4), allocatable :: touched(:), wcnt(:), nbr(:)
    integer :: i, k, e, vtx, jj, j0, j1, ne, nnz, deg, off, nb
    integer(i4) :: rowsum

    ! --- vertex -> element (local) CSR over the subset ---
    allocate(v2c_cnt(nrnk)); v2c_cnt = 0
    do i = 1, n
       e = loc(i)
       do k = 1, nv
          vtx = cell(k, e)
          v2c_cnt(vtx) = v2c_cnt(vtx) + 1
       end do
    end do
    allocate(v2c_ptr(nrnk + 1))
    v2c_ptr(1) = 1
    do i = 1, nrnk
       v2c_ptr(i + 1) = v2c_ptr(i) + v2c_cnt(i)
    end do
    allocate(v2c(v2c_ptr(nrnk + 1) - 1))
    v2c_cnt = 0
    do i = 1, n
       e = loc(i)
       do k = 1, nv
          vtx = cell(k, e)
          v2c(v2c_ptr(vtx) + v2c_cnt(vtx)) = int(i, i4) ! store LOCAL index
          v2c_cnt(vtx) = v2c_cnt(vtx) + 1
       end do
    end do

    ! --- accumulate adjacency; two passes (count, then fill) ---
    allocate(touched(n)); touched = 0
    allocate(wcnt(n)); wcnt = 0
    allocate(nbr(n))
    allocate(ia(0:n))
    ia(0) = 0
    nnz = 0
    ! pass 1: count degree (incl. self) per row
    do i = 1, n
       e = loc(i)
       nb = 0
       do k = 1, nv
          vtx = cell(k, e)
          j0 = v2c_ptr(vtx); j1 = v2c_ptr(vtx + 1) - 1
          do jj = j0, j1
             if (touched(v2c(jj)) /= i) then
                touched(v2c(jj)) = i
                nb = nb + 1
             end if
          end do
       end do
       nnz = nnz + nb
       ia(i) = nnz
    end do
    allocate(ja(nnz), va(nnz))
    ! pass 2: fill neighbours + shared-vertex weights
    touched = 0
    do i = 1, n
       e = loc(i)
       off = int(ia(i-1))
       deg = 0
       do k = 1, nv
          vtx = cell(k, e)
          j0 = v2c_ptr(vtx); j1 = v2c_ptr(vtx + 1) - 1
          do jj = j0, j1
             ne = v2c(jj)
             if (touched(ne) /= i) then
                touched(ne) = i
                deg = deg + 1
                nbr(deg) = int(ne, i4)
                wcnt(ne) = 1
             else
                wcnt(ne) = wcnt(ne) + 1
             end if
          end do
       end do
       rowsum = 0
       do k = 1, deg
          ne = nbr(k)
          ja(off + k) = int(ne, i4)
          if (ne == i) then
             va(off + k) = 0.0_dp ! placeholder; set to rowsum below
          else
             rowsum = rowsum + wcnt(ne)
             va(off + k) = -real(wcnt(ne), dp)
          end if
       end do
       ! set the diagonal entry to +rowsum
       do k = 1, deg
          if (ja(off + k) == i) va(off + k) = real(rowsum, dp)
       end do
    end do

    deallocate(v2c_cnt, v2c_ptr, v2c, touched, wcnt, nbr)
  end subroutine build_dual

  !> The KWANT lowest non-trivial eigenvectors (Ritz vectors) of the graph
  !! Laplacian (ia,ja,va) into columns of F(n,KWANT); kgot returns how many were
  !! produced. Uses a Lanczos iteration with FULL reorthogonalisation that
  !! deflates the constant null vector each step (ortho1) -- genmap's Lanczos ->
  !! spectral idea, but numerically robust (plain Lanczos loses orthogonality and
  !! under-resolves the low modes on hard spectra). The smallest mode is the
  !! Fiedler vector; the caller picks the lowest-cut split among the KWANT modes.
  subroutine fiedler_k(ia, ja, va, n, kwant, F, kgot)
    integer(i4), intent(in) :: ia(0:n), ja(:)
    real(dp), intent(in) :: va(:)
    integer, intent(in) :: n, kwant
    real(dp), intent(out) :: F(n, kwant)
    integer, intent(out) :: kgot
    real(dp), allocatable :: qv(:,:), alpha(:), bet(:), w(:), f0(:)
    real(dp), allocatable :: dd(:), ee(:), zz(:,:)
    integer(i4), allocatable :: ord(:)
    real(dp) :: nrm, s, bk
    integer :: m, i, j, k, mm, c, idum
    integer(i8) :: hh

    m = min(MLANC, n - 1)
    allocate(qv(n, m + 1), alpha(m), bet(m), w(n), f0(n))

    ! Deterministic generic start with components along every eigenmode (a fixed
    ! integer hash -> [-0.5,0.5]), so Lanczos resolves the low modes regardless
    ! of how the elements are ordered in the file.
    do i = 1, n
       hh = iand(int(i, i8) * 2654435761_i8, 2147483647_i8)
       hh = iand(ieor(hh, ishft(hh, -15)), 1048575_i8)
       f0(i) = real(hh, dp) / 1048575.0_dp - 0.5_dp
    end do
    call ortho1(f0, n)
    nrm = sqrt(dot(f0, f0, n))
    if (nrm <= 0.0_dp) then
       do i = 1, n
          f0(i) = real(i, dp)
       end do
       call ortho1(f0, n); nrm = sqrt(dot(f0, f0, n))
    end if
    do i = 1, n
       qv(i, 1) = f0(i) / nrm
    end do

    mm = m
    do k = 1, m
       call ax(w, qv(1, k), ia, ja, va, n)
       call ortho1(w, n) ! stay orthogonal to 1 (null space)
       alpha(k) = dot(w, qv(1, k), n)
       if (k == 1) then
          do i = 1, n
             w(i) = w(i) - alpha(k) * qv(i, k)
          end do
       else
          do i = 1, n
             w(i) = w(i) - alpha(k) * qv(i, k) - bet(k-1) * qv(i, k-1)
          end do
       end if
       do j = 1, k ! full reorthogonalisation
          s = dot(w, qv(1, j), n)
          do i = 1, n
             w(i) = w(i) - s * qv(i, j)
          end do
       end do
       bk = sqrt(dot(w, w, n))
       if (bk <= 1.0e-12_dp) then
          mm = k; exit
       end if
       bet(k) = bk
       do i = 1, n
          qv(i, k+1) = w(i) / bk
       end do
    end do

    allocate(dd(mm), ee(mm), zz(mm, mm), ord(mm))
    call tql2_eig(alpha, bet, dd, ee, zz, mm, idum) ! all eigenpairs of T
    call rsort(dd, ord, mm) ! ord(1) = smallest eigenvalue
    kgot = min(kwant, mm)
    do c = 1, kgot
       do i = 1, n
          F(i, c) = 0.0_dp
       end do
       do k = 1, mm
          s = zz(k, ord(c))
          do i = 1, n
             F(i, c) = F(i, c) + qv(i, k) * s
          end do
       end do
    end do
    deallocate(qv, alpha, bet, w, f0, dd, ee, zz, ord)
  end subroutine fiedler_k

  !> y = L x, L the graph Laplacian in CSR (diagonal stored explicitly in va).
  subroutine ax(y, x, ia, ja, va, n)
    integer, intent(in) :: n
    real(dp), intent(out) :: y(n)
    real(dp), intent(in) :: x(n), va(:)
    integer(i4), intent(in) :: ia(0:n), ja(:)
    integer :: i, j
    do i = 1, n
       y(i) = 0.0_dp
       do j = int(ia(i-1)) + 1, int(ia(i))
          y(i) = y(i) + va(j) * x(ja(j))
       end do
    end do
  end subroutine ax

  !> Symmetric tridiagonal QL with implicit shifts (Numerical Recipes / genmap
  !! calcvec). diag(1:m), upper(1:m-1) input; eigenvalues in d, eigenvectors in
  !! columns of z; imin = index of the smallest eigenvalue. Eigenvectors are
  !! orthonormalised as in genmap.
  subroutine tql2_eig(diag, upper, d, e, z, n, imin)
    integer, intent(in) :: n
    real(dp), intent(in) :: diag(*), upper(*)
    real(dp), intent(inout) :: d(n), e(n), z(n, n)
    integer, intent(out) :: imin
    integer :: i, k, l, m, iter
    real(dp) :: dd, g, r, s, c, p, f, b, dmin, sc, enorm
    do i = 1, n
       d(i) = diag(i)
    end do
    do i = 1, n - 1
       e(i) = upper(i)
    end do
    e(n) = 0.0_dp
    z = 0.0_dp
    do i = 1, n
       z(i, i) = 1.0_dp
    end do
    do l = 1, n
       iter = 0
1      continue
       do m = l, n - 1
          dd = abs(d(m)) + abs(d(m + 1))
          if (abs(e(m)) + dd == dd) go to 2
       end do
       m = n
2      continue
       if (m /= l) then
          if (iter == 30) then
             write(error_unit,*) 'genmap_light: tql2 too many iterations'
             exit
          end if
          iter = iter + 1
          g = (d(l + 1) - d(l)) / (2.0_dp * e(l))
          r = sqrt(g**2 + 1.0_dp)
          g = d(m) - d(l) + e(l) / (g + sign(r, g))
          s = 1.0_dp; c = 1.0_dp; p = 0.0_dp
          do i = m - 1, l, -1
             f = s * e(i)
             b = c * e(i)
             if (abs(f) >= abs(g)) then
                c = g / f
                r = sqrt(c**2 + 1.0_dp)
                e(i + 1) = f * r
                s = 1.0_dp / r
                c = c * s
             else
                s = f / g
                r = sqrt(s**2 + 1.0_dp)
                e(i + 1) = g * r
                c = 1.0_dp / r
                s = s * c
             end if
             g = d(i + 1) - p
             r = (d(i) - g) * s + 2.0_dp * c * b
             p = s * r
             d(i + 1) = g + p
             g = c * r - b
             do k = 1, n
                f = z(k, i + 1)
                z(k, i + 1) = s * z(k, i) + c * f
                z(k, i) = c * z(k, i) - s * f
             end do
          end do
          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0_dp
          go to 1
       end if
    end do
    ! smallest eigenvalue
    imin = 1
    dmin = d(1)
    do i = 1, n
       if (d(i) < dmin) then
          dmin = d(i); imin = i
       end if
    end do
    ! orthonormalise eigenvectors (genmap: scale each column by 1/||.||)
    do k = 1, n
       enorm = sqrt(abs(dot(z(1, k), z(1, k), n)))
       if (enorm /= 0.0_dp) then
          sc = 1.0_dp / enorm
          do i = 1, n
             z(i, k) = z(i, k) * sc
          end do
       end if
    end do
  end subroutine tql2_eig

  subroutine ortho1(p, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: p(n)
    real(dp) :: s
    integer :: i
    s = 0.0_dp
    do i = 1, n
       s = s + p(i)
    end do
    s = s / real(n, dp)
    do i = 1, n
       p(i) = p(i) - s
    end do
  end subroutine ortho1

  pure function dot(x, y, n) result(s)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n)
    real(dp) :: s
    integer :: i
    s = 0.0_dp
    do i = 1, n
       s = s + x(i) * y(i)
    end do
  end function dot

  !> Heap sort of real key a(1:n) ascending; returns permutation perm so that
  !! a(perm(1)) <= a(perm(2)) <= ...  (stable enough; ties keep heap order).
  subroutine rsort(a, perm, n)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n)
    integer(i4), intent(out) :: perm(n)
    integer :: i, l, ir, indx, j
    real(dp) :: q
    do i = 1, n
       perm(i) = int(i, i4)
    end do
    if (n <= 1) return
    l = n / 2 + 1
    ir = n
    do
       if (l > 1) then
          l = l - 1
          indx = perm(l)
          q = a(indx)
       else
          indx = perm(ir)
          q = a(indx)
          perm(ir) = perm(1)
          ir = ir - 1
          if (ir == 1) then
             perm(1) = indx
             return
          end if
       end if
       i = l
       j = l + l
       do while (j <= ir)
          if (j < ir) then
             if (a(perm(j)) < a(perm(j + 1))) j = j + 1
          end if
          if (q < a(perm(j))) then
             perm(i) = perm(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if
       end do
       perm(i) = indx
    end do
  end subroutine rsort

  !> DFS connectivity of the subgraph induced by {i : side(i)==s}, over local
  !! indices with adjacency (ia,ja). Empty set counts as connected.
  function conn_side(side, s, n, ia, ja) result(ok)
    integer, intent(in) :: s, n
    integer(i4), intent(in) :: side(n), ia(0:n), ja(:)
    logical :: ok
    integer(i4), allocatable :: stack(:), seen(:)
    integer :: i, top, cur, j, cnt, total, start
    allocate(stack(n), seen(n)); seen = 0
    total = 0; start = 0
    do i = 1, n
       if (side(i) == s) then
          total = total + 1
          if (start == 0) start = i
       end if
    end do
    if (total == 0) then
       ok = .true.; deallocate(stack, seen); return
    end if
    top = 1; stack(1) = int(start, i4); seen(start) = 1; cnt = 1
    do while (top > 0)
       cur = stack(top); top = top - 1
       do j = int(ia(cur-1)) + 1, int(ia(cur))
          i = ja(j)
          if (i >= 1 .and. i <= n) then
             if (side(i) == s .and. seen(i) == 0) then
                seen(i) = 1; cnt = cnt + 1
                top = top + 1; stack(top) = int(i, i4)
             end if
          end if
       end do
    end do
    ok = (cnt == total)
    deallocate(stack, seen)
  end function conn_side

  !> Region-growing fallback bisection: grow side 1 by BFS over the dual graph
  !! from a fixed seed until it reaches n1 elements; everything not reached is
  !! side 2. This guarantees a connected side 1 at the exact target size (side 2
  !! is the complement and may itself be disconnected) -- genmap's growth
  !! strategy in simplified form. If the input subgraph is disconnected and the
  !! BFS runs dry early, side 1 is topped up deterministically.
  subroutine region_grow(n, n1, side, ia, ja)
    integer, intent(in) :: n, n1
    integer(i4), intent(in) :: ia(0:n), ja(:)
    integer(i4), intent(inout) :: side(n)
    integer(i4), allocatable :: q(:)
    integer :: head, tail, cur, j, i, c1
    ! seed 1 = local index 1; seed 2 = the farthest-labelled index (n)
    allocate(q(n))
    side = 0
    side(1) = 1
    c1 = 1
    head = 1; tail = 1; q(1) = 1_i4
    do while (head <= tail .and. c1 < n1)
       cur = q(head); head = head + 1
       do j = int(ia(cur-1)) + 1, int(ia(cur))
          i = ja(j)
          if (i >= 1 .and. i <= n) then
             if (side(i) == 0) then
                side(i) = 1; c1 = c1 + 1
                tail = tail + 1; q(tail) = int(i, i4)
                if (c1 >= n1) exit
             end if
          end if
       end do
    end do
    ! any element not reached (still 0) goes to side 2
    do i = 1, n
       if (side(i) == 0) side(i) = 2
    end do
    ! if side 1 fell short (disconnected input), top it up with lowest-index
    ! side-2 elements to hit the exact target n1
    if (c1 < n1) then
       do i = 1, n
          if (c1 >= n1) exit
          if (side(i) == 2) then
             side(i) = 1; c1 = c1 + 1
          end if
       end do
    end if
    deallocate(q)
  end subroutine region_grow

end module gm_rsb


program genmap_light
  use, intrinsic :: iso_fortran_env, only: dp => real64, i4 => int32, i8 => int64, &
       error_unit
  use gm_hash, only: imap_t
  use gm_rsb, only: rsb_partition
  implicit none

  ! Neko facet (1..6) corners in nmsh vertex indices (= face_nodes composed with
  ! the nmsh->mesh read swap [1,2,4,3,5,6,8,7]); same table validated byte-exact
  ! for periodic glb_pt_ids in rea2nbin_light / mesh_checker_light.
  integer(i4), parameter :: FACE_RE2(4,6) = reshape((/ &
       1,5,8,4, 2,6,7,3, 1,2,6,5, 4,3,7,8, 1,2,3,4, 5,6,7,8 /), (/4,6/))

  character(len=1024) :: fin, fout, argbuf
  integer :: uin, uout, ios, i, k, argc, nparts
  integer(i4) :: nelv, gdim, elid, nz, ncur, ztype
  integer(i4) :: ze, zf, zpe, zpf, g1, g2, g3, g4
  integer(i4), allocatable :: elid_of_pos(:), pos_of_elid(:)
  integer(i4), allocatable :: ev(:,:) ! (8,nelv) original nmsh point ids
  real(dp), allocatable :: cx(:,:), cy(:,:), cz(:,:) ! (8,nelv) coords
  integer(i4), allocatable :: cell(:,:) ! (8,nelv) merged+compressed vtx ids
  integer(i4), allocatable :: part(:), newid(:), pos_of_newid(:)
  integer(i4), allocatable :: blk(:) ! partition prefix offsets
  ! zone storage
  integer(i4), allocatable :: pz_e(:), pz_f(:), pz_pe(:), pz_pf(:), pz_g(:,:)
  integer(i4), allocatable :: lz_e(:), lz_f(:), lz_lab(:)
  integer(i4), allocatable :: xz(:,:) ! legacy zone types, raw 9-int records
  integer :: npz, nlz, nxz
  ! curve storage
  integer(i4), allocatable :: cu_e(:), cu_type(:,:)
  real(dp), allocatable :: cu_data(:,:,:)
  integer :: ncu
  type(imap_t) :: rm, comp
  integer(i4) :: nrnk, nextid, vid, mid
  real(dp) :: xyz(3)
  integer(i8) :: cut, edges

  argc = command_argument_count()
  if (argc < 2) then
     write(*,*) 'Usage: genmap_light mesh.nmsh nparts [out.nmsh]'
     write(*,*) '  Partition a Neko .nmsh into <nparts> parts by recursive'
     write(*,*) '  spectral bisection (Nek5000 genmap algorithm) and write a'
     write(*,*) '  reordered .nmsh whose elements are grouped per partition'
     write(*,*) '  (a linear read then reproduces the partition, like prepart).'
     stop 1
  end if
  call get_command_argument(1, fin)
  call get_command_argument(2, argbuf)
  read(argbuf, *, iostat=ios) nparts
  if (ios /= 0 .or. nparts < 1) then
     write(error_unit,*) 'Error: nparts must be a positive integer'; stop 1
  end if
  if (argc >= 3) then
     call get_command_argument(3, fout)
  else
     block
       integer :: islash, idot
       islash = scan(trim(fin), '/', back=.true.) ! 0 if no directory part
       idot = scan(trim(fin), '.', back=.true.)
       if (idot > islash) then ! a real extension to strip
          fout = trim(fin(1:idot-1)) // '_' // trim(adjustl(argbuf)) // '.nmsh'
       else ! no extension (or '.' only in dir)
          fout = trim(fin) // '_' // trim(adjustl(argbuf)) // '.nmsh'
       end if
     end block
  end if

  write(*,'(a)') ''
  write(*,'(a)') '     _  __  ____  __ __  ____'
  write(*,'(a)') '    / |/ / / __/ / //_/ / __ \'
  write(*,'(a)') '   /    / / _/  / ,<   / /_/ /'
  write(*,'(a)') '  /_/|_/ /___/ /_/|_|  \____/   genmap_light  (spectral partitioner)'
  write(*,'(a)') ''
  write(*,'(a,a)') '  input     : ', trim(fin)
  write(*,'(a,a)') '  output    : ', trim(fout)
  write(*,'(a,i0)') '  nparts    : ', nparts

  ! ---- read header + elements ----
  open(newunit=uin, file=trim(fin), access='stream', form='unformatted', &
       status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(error_unit,*) 'Error: cannot open ', trim(fin); stop 1
  end if
  read(uin, iostat=ios) nelv, gdim; call chk_read(ios)
  if (gdim /= 3) then
     write(error_unit,*) 'Error: genmap_light supports 3D (hex) meshes only (gdim=', &
          gdim, ').'; stop 1
  end if
  if (nparts > nelv) then
     write(error_unit,*) 'Error: nparts (', nparts, ') > number of elements (', &
          nelv, ').'; stop 1
  end if
  write(*,'(a,i0,a)') '  mesh      : ', nelv, ' hex elements'

  allocate(ev(8, nelv), cx(8, nelv), cy(8, nelv), cz(8, nelv))
  allocate(elid_of_pos(nelv), pos_of_elid(nelv))
  write(*,'(a)') '  [1/5] reading elements ...'
  do i = 1, nelv
     read(uin, iostat=ios) elid; call chk_read(ios)
     elid_of_pos(i) = elid
     do k = 1, 8
        read(uin, iostat=ios) ev(k, i), xyz(1), xyz(2), xyz(3); call chk_read(ios)
        cx(k, i) = xyz(1); cy(k, i) = xyz(2); cz(k, i) = xyz(3)
     end do
  end do
  pos_of_elid = 0
  do i = 1, nelv
     if (elid_of_pos(i) < 1 .or. elid_of_pos(i) > nelv) then
        write(error_unit,*) 'Error: element id out of [1,nelv]: ', elid_of_pos(i); stop 1
     end if
     pos_of_elid(elid_of_pos(i)) = int(i, i4)
  end do
  do i = 1, nelv ! element ids must be a permutation of 1..nelv
     if (pos_of_elid(i) == 0) then
        write(error_unit,*) 'Error: element ids are not a permutation of 1..nelv', &
             ' (duplicate or missing id ', i, ').'; stop 1
     end if
  end do

  ! ---- read zones ----
  write(*,'(a)') '  [2/5] reading zones and curves ...'
  read(uin, iostat=ios) nz; call chk_read(ios)
  allocate(pz_e(nz), pz_f(nz), pz_pe(nz), pz_pf(nz), pz_g(4, nz))
  allocate(lz_e(nz), lz_f(nz), lz_lab(nz), xz(9, max(nz, 1)))
  npz = 0; nlz = 0; nxz = 0
  call rm%init(max(1024_i8, int(nelv, i8)))
  do i = 1, nz
     read(uin, iostat=ios) ze, zf, zpe, zpf, g1, g2, g3, g4, ztype; call chk_read(ios)
     if (ze < 1 .or. ze > nelv .or. zf < 1 .or. zf > 6) cycle
     if (ztype == 5) then
        npz = npz + 1
        pz_e(npz) = ze; pz_f(npz) = zf; pz_pe(npz) = zpe; pz_pf(npz) = zpf
        pz_g(1, npz) = g1; pz_g(2, npz) = g2; pz_g(3, npz) = g3; pz_g(4, npz) = g4
        ! periodic point merge, exactly as mesh_checker_light / Neko on read
        block
          integer :: pos
          pos = pos_of_elid(ze)
          call rm%set(ev(FACE_RE2(1,zf), pos), g1)
          call rm%set(ev(FACE_RE2(2,zf), pos), g2)
          call rm%set(ev(FACE_RE2(3,zf), pos), g3)
          call rm%set(ev(FACE_RE2(4,zf), pos), g4)
        end block
     else if (ztype == 7) then
        nlz = nlz + 1
        lz_e(nlz) = ze; lz_f(nlz) = zf; lz_lab(nlz) = zpf
     else
        ! legacy zone types (e.g. type 1/2 in older meshes): carry through
        ! verbatim, renumbering only the element id at write time
        nxz = nxz + 1
        xz(:, nxz) = (/ ze, zf, zpe, zpf, g1, g2, g3, g4, ztype /)
     end if
  end do
  if (nxz > 0) then
     write(*,'(a,i0,a)') '        note: carrying ', nxz, &
          ' zone records of legacy types through (element ids renumbered)'
  end if

  ! ---- read curves ----
  read(uin, iostat=ios) ncur; call chk_read(ios)
  ncu = ncur
  if (ncu > 0) then
     allocate(cu_e(ncu), cu_type(12, ncu), cu_data(5, 12, ncu))
     do i = 1, ncu
        read(uin, iostat=ios) cu_e(i), cu_data(:, :, i), cu_type(:, i); call chk_read(ios)
        if (cu_e(i) < 1 .or. cu_e(i) > nelv) then
           write(error_unit,*) 'Error: curve element id out of [1,nelv]: ', cu_e(i); stop 1
        end if
     end do
  end if
  close(uin)

  ! ---- build cell(): apply periodic merge, then compress to 1..nrnk ----
  write(*,'(a)') '  [3/5] building element dual graph ...'
  allocate(cell(8, nelv))
  call comp%init(max(1024_i8, int(nelv, i8)))
  nextid = 0
  do i = 1, nelv
     do k = 1, 8
        vid = rm%getor(ev(k, i), ev(k, i)) ! merged id (identity if not periodic)
        mid = comp%intern(vid, nextid) ! dense 1..nrnk
        cell(k, i) = mid
     end do
  end do
  nrnk = nextid

  ! ---- recursive spectral bisection ----
  write(*,'(a,i0,a)') '  [4/5] recursive spectral bisection into ', nparts, ' parts ...'
  allocate(part(nelv))
  call rsb_partition(cell, 8, nelv, nrnk, nparts, part)

  ! edge-cut diagnostic (on the merged dual graph)
  call count_cut(cell, 8, nelv, nrnk, part, cut, edges)
  write(*,'(a,i0,a,i0,a,f7.3,a)') '        edge cut: ', cut, ' / ', edges, &
       ' shared-vertex links cut (', 100.0_dp*real(cut,dp)/real(max(edges,1_i8),dp), '%)'

  ! ---- reorder (prepart style): contiguous per-partition blocks ----
  write(*,'(a)') '  [5/5] renumbering and writing reordered mesh ...'
  allocate(blk(0:nparts))
  blk = 0
  do i = 1, nelv
     blk(part(i) + 1) = blk(part(i) + 1) + 1
  end do
  do i = 1, nparts
     blk(i) = blk(i) + blk(i - 1) ! exclusive prefix sum -> block starts
  end do
  allocate(newid(nelv), pos_of_newid(nelv))
  block
    integer(i4) :: cnt(0:nparts-1)
    integer :: r
    cnt = 0
    do i = 1, nelv
       r = part(i)
       cnt(r) = cnt(r) + 1
       newid(i) = blk(r) + cnt(r) ! 1-based new global id
       pos_of_newid(newid(i)) = int(i, i4)
    end do
  end block

  ! ---- write output nmsh ----
  open(newunit=uout, file=trim(fout), access='stream', form='unformatted', &
       status='replace', action='write')
  write(uout) nelv, gdim
  do i = 1, nelv
     block
       integer :: pos
       pos = pos_of_newid(i)
       write(uout) int(i, i4) ! el_idx = new contiguous global id
       do k = 1, 8
          write(uout) ev(k, pos), cx(k, pos), cy(k, pos), cz(k, pos)
       end do
     end block
  end do

  ! zones: type-5 (periodic) first, then type-7 (labeled), matching Neko's
  ! writer, then any legacy zone types carried through
  write(uout) int(npz + nlz + nxz, i4)
  do i = 1, npz
     block
       integer(i4) :: enew, penew
       enew = newid(pos_of_elid(pz_e(i)))
       if (pz_pe(i) >= 1 .and. pz_pe(i) <= nelv) then
          penew = newid(pos_of_elid(pz_pe(i)))
       else
          penew = pz_pe(i)
       end if
       write(uout) enew, pz_f(i), penew, pz_pf(i), &
            pz_g(1,i), pz_g(2,i), pz_g(3,i), pz_g(4,i), 5_i4
     end block
  end do
  do i = 1, nlz
     write(uout) newid(pos_of_elid(lz_e(i))), lz_f(i), 0_i4, lz_lab(i), &
          0_i4, 0_i4, 0_i4, 0_i4, 7_i4
  end do
  do i = 1, nxz
     write(uout) newid(pos_of_elid(xz(1, i))), xz(2:9, i)
  end do

  ! curves: first-class -- remap element id, carry curve_data/type verbatim
  write(uout) int(ncu, i4)
  do i = 1, ncu
     write(uout) newid(pos_of_elid(cu_e(i))), cu_data(:, :, i), cu_type(:, i)
  end do
  close(uout)

  write(*,'(a,a)') '  done -> ', trim(fout)

  call rm%free(); call comp%free()

contains

  !> Abort cleanly if the .nmsh is truncated (reads happen before any output
  !! file is opened, so nothing partial is left behind).
  subroutine chk_read(ios_)
    integer, intent(in) :: ios_
    if (ios_ == 0) return
    write(error_unit,*) 'Error: truncated or corrupt .nmsh file (unexpected ', &
         'end of data)'
    stop 1
  end subroutine chk_read


  !> Weighted edge cut of the partition on the merged dual graph. Edge weight =
  !! number of shared vertices between two elements; `edges` is the total weight,
  !! `cut` the weight crossing partition boundaries. Computed by, for each vertex,
  !! counting element pairs that share it (a face-pair shares 4 vertices -> weight
  !! 4, matching genmap's connectivity weighting).
  subroutine count_cut(cell, nv, nelv, nrnk, part, cut, edges)
    integer(i4), intent(in) :: cell(nv, nelv), part(nelv)
    integer, intent(in) :: nv, nelv, nrnk
    integer(i8), intent(out) :: cut, edges
    integer(i4), allocatable :: cptr(:), ccnt(:), clist(:)
    integer :: i, k, vtx, a, b, j0, j1, ea, eb
    allocate(ccnt(nrnk)); ccnt = 0
    do i = 1, nelv
       do k = 1, nv
          ccnt(cell(k, i)) = ccnt(cell(k, i)) + 1
       end do
    end do
    allocate(cptr(nrnk + 1)); cptr(1) = 1
    do i = 1, nrnk
       cptr(i + 1) = cptr(i) + ccnt(i)
    end do
    allocate(clist(cptr(nrnk + 1) - 1)); ccnt = 0
    do i = 1, nelv
       do k = 1, nv
          vtx = cell(k, i)
          clist(cptr(vtx) + ccnt(vtx)) = int(i, i4)
          ccnt(vtx) = ccnt(vtx) + 1
       end do
    end do
    cut = 0_i8; edges = 0_i8
    do vtx = 1, nrnk
       j0 = cptr(vtx); j1 = cptr(vtx + 1) - 1
       do a = j0, j1
          ea = clist(a)
          do b = a + 1, j1
             eb = clist(b)
             if (ea == eb) cycle ! skip an element sharing a merged vertex
             ! with itself (degenerate periodic/thin mesh)
             edges = edges + 1_i8
             if (part(ea) /= part(eb)) cut = cut + 1_i8
          end do
       end do
    end do
    deallocate(ccnt, cptr, clist)
  end subroutine count_cut

end program genmap_light
