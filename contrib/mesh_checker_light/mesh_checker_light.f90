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
!> mesh_checker_light -- a lightweight, CPU-only diagnostics tool for Neko `.nmsh`
!! meshes.
!!
!! Neko's `contrib/mesh_checker` builds a full `mesh_t` using Neko's
!! general-purpose, production hash tables. Those are robust and reused across the
!! code, but they are sized generously and store polymorphic keys, so they need a
!! lot of memory for a large mesh. This tool is a lightweight, CPU-only
!! reimplementation of just the checks that are needed, using plain typed integer
!! arrays (a face hash and an edge hash keyed on the nmsh point ids) and no Neko
!! initialisation. Periodic facet-pairs are merged via the record's glb_pt_ids
!! exactly as Neko does on read, so the face/edge counts match glb_mfcs/glb_meds.
!! It gives the same result in a small fraction of the memory.
!!
!! Because it does not link Neko, it needs only a Fortran compiler (no MPI, no
!! Neko library) and runs anywhere -- in particular you do not need a separate
!! CPU build of Neko. The `contrib` tools are CPU-only, so on a Neko configured
!! for GPUs you would otherwise have to rebuild it for CPU just to run them; this
!! sidesteps that entirely.
!!
!! Reports, exactly like `mesh_checker`: element / point / face / edge counts,
!! bounding box, number of periodic faces, labeled-zone face counts, and -- the
!! key diagnostic -- the number of UNLABELED EXTERNAL faces (an external face
!! that is neither a labeled zone nor periodic; must be 0 for a valid mesh).
!! Exits non-zero if any unlabeled external face is found (as `mesh_checker` does).
!! Two optional checks are available: `--jacobian` (a streaming negative/zero
!! Jacobian scan, exact for straight-sided elements) and `--write-zone-indices`
!! (writes a field file marking boundary faces by zone index). 3D (hex) meshes
!! only.
!!
!! Tested against every 3D mesh shipped with Neko (cross-checked field for field
!! against an independent recomputation), but it remains your responsibility to
!! confirm the result is correct for your own mesh -- run `make check`, or
!! cross-check a real `mesh_checker` run, on a case you trust first. See
!! mesh_checker_light.md for the memory model, validation, and scope.
module mc_hash
  use, intrinsic :: iso_fortran_env, only: i4 => int32, i8 => int64, i1 => int8
  implicit none
  private
  public :: fhash_t, ehash_t, rmap_t

  !> Small remap of a raw nmsh point id -> its merged (periodic) point id.
  !! An id that is absent maps to itself (identity). This reproduces Neko's
  !! `apply_periodic_facet`, which overwrites a periodic facet's corner ids with
  !! the record's `glb_pt_ids` on read, before the face/edge tables are built.
  type :: rmap_t
     integer(i4), allocatable :: key(:), val(:)
     integer(i8) :: cap = 0, mask = 0, n = 0
   contains
     procedure :: init => r_init
     procedure :: set  => r_set
     procedure :: get  => r_get
     procedure :: free => r_free
  end type rmap_t

  !> Hash of unique faces: key = 4 sorted point ids; tracks an occurrence count
  !! (capped at 2 -- external = seen once, internal = seen twice).
  type :: fhash_t
     integer(i4), allocatable :: k1(:), k2(:), k3(:), k4(:)
     integer(i1), allocatable :: cnt(:)
     integer(i8) :: cap = 0, mask = 0, n = 0
   contains
     procedure :: init => f_init
     procedure :: add  => f_add
     procedure :: cnt_of => f_cnt_of
     procedure :: free => f_free
  end type fhash_t

  !> Hash of unique edges: key = 2 sorted point ids (membership only).
  type :: ehash_t
     integer(i4), allocatable :: k1(:), k2(:)
     integer(i8) :: cap = 0, mask = 0, n = 0
   contains
     procedure :: init => e_init
     procedure :: add  => e_add
     procedure :: free => e_free
  end type ehash_t

contains

  pure function mix4(a, b, c, d) result(h)
    integer(i4), intent(in) :: a, b, c, d
    integer(i8) :: h
    h = int(a, i8)
    h = ieor(h, ishft(h, -33)) * (-49064778989728563_i8)
    h = ieor(h, int(b, i8))
    h = ieor(h, ishft(h, -29)) * (-4417276706812531889_i8)
    h = ieor(h, int(c, i8))
    h = ieor(h, ishft(h, -32)) * (-49064778989728563_i8)
    h = ieor(h, int(d, i8))
    h = ieor(h, ishft(h, -32))
  end function mix4

  ! -------- face hash --------
  subroutine f_init(this, cap0)
    class(fhash_t), intent(inout) :: this
    integer(i8), intent(in) :: cap0
    integer(i8) :: c
    c = 1024
    do while (c < cap0)
       c = ishft(c, 1)
    end do
    this%cap = c; this%mask = c - 1; this%n = 0
    allocate(this%k1(0:c-1), this%k2(0:c-1), this%k3(0:c-1), this%k4(0:c-1))
    allocate(this%cnt(0:c-1))
    this%k1 = 0; this%cnt = 0_i1
  end subroutine f_init

  !> Insert a face (its 4 point ids, in any order) and bump its count.
  subroutine f_add(this, a, b, c, d)
    class(fhash_t), intent(inout) :: this
    integer(i4), intent(in) :: a, b, c, d
    integer(i4) :: s(4), t
    integer(i8) :: h
    s = (/ a, b, c, d /)
    ! sort4 ascending (5 compare-swaps)
    if (s(1) > s(2)) then; t = s(1); s(1) = s(2); s(2) = t; end if
    if (s(3) > s(4)) then; t = s(3); s(3) = s(4); s(4) = t; end if
    if (s(1) > s(3)) then; t = s(1); s(1) = s(3); s(3) = t; end if
    if (s(2) > s(4)) then; t = s(2); s(2) = s(4); s(4) = t; end if
    if (s(2) > s(3)) then; t = s(2); s(2) = s(3); s(3) = t; end if
    if (int(this%n, i8) * 3_i8 > this%cap * 2_i8) call f_grow(this)
    h = iand(mix4(s(1), s(2), s(3), s(4)), this%mask)
    do
       if (this%k1(h) == 0) then
          this%k1(h) = s(1); this%k2(h) = s(2)
          this%k3(h) = s(3); this%k4(h) = s(4)
          this%cnt(h) = 1_i1; this%n = this%n + 1; return
       else if (this%k1(h) == s(1) .and. this%k2(h) == s(2) .and. &
                this%k3(h) == s(3) .and. this%k4(h) == s(4)) then
          if (this%cnt(h) < 2_i1) this%cnt(h) = this%cnt(h) + 1_i1
          return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end subroutine f_add

  !> Occurrence count of a face (0 if absent).
  function f_cnt_of(this, a, b, c, d) result(cnt)
    class(fhash_t), intent(in) :: this
    integer(i4), intent(in) :: a, b, c, d
    integer(i4) :: s(4), t, cnt
    integer(i8) :: h
    s = (/ a, b, c, d /)
    if (s(1) > s(2)) then; t = s(1); s(1) = s(2); s(2) = t; end if
    if (s(3) > s(4)) then; t = s(3); s(3) = s(4); s(4) = t; end if
    if (s(1) > s(3)) then; t = s(1); s(1) = s(3); s(3) = t; end if
    if (s(2) > s(4)) then; t = s(2); s(2) = s(4); s(4) = t; end if
    if (s(2) > s(3)) then; t = s(2); s(2) = s(3); s(3) = t; end if
    h = iand(mix4(s(1), s(2), s(3), s(4)), this%mask)
    do
       if (this%k1(h) == 0) then
          cnt = 0; return
       else if (this%k1(h) == s(1) .and. this%k2(h) == s(2) .and. &
                this%k3(h) == s(3) .and. this%k4(h) == s(4)) then
          cnt = int(this%cnt(h), i4); return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end function f_cnt_of

  subroutine f_grow(this)
    class(fhash_t), intent(inout) :: this
    integer(i4), allocatable :: o1(:), o2(:), o3(:), o4(:)
    integer(i1), allocatable :: oc(:)
    integer(i8) :: i, h, oldcap
    oldcap = this%cap
    call move_alloc(this%k1, o1); call move_alloc(this%k2, o2)
    call move_alloc(this%k3, o3); call move_alloc(this%k4, o4)
    call move_alloc(this%cnt, oc)
    this%cap = ishft(this%cap, 1); this%mask = this%cap - 1
    allocate(this%k1(0:this%cap-1), this%k2(0:this%cap-1), &
             this%k3(0:this%cap-1), this%k4(0:this%cap-1))
    allocate(this%cnt(0:this%cap-1)); this%k1 = 0; this%cnt = 0_i1
    do i = 0, oldcap - 1
       if (o1(i) /= 0) then
          h = iand(mix4(o1(i), o2(i), o3(i), o4(i)), this%mask)
          do while (this%k1(h) /= 0)
             h = iand(h + 1_i8, this%mask)
          end do
          this%k1(h) = o1(i); this%k2(h) = o2(i)
          this%k3(h) = o3(i); this%k4(h) = o4(i); this%cnt(h) = oc(i)
       end if
    end do
    deallocate(o1, o2, o3, o4, oc)
  end subroutine f_grow

  subroutine f_free(this)
    class(fhash_t), intent(inout) :: this
    if (allocated(this%k1)) deallocate(this%k1, this%k2, this%k3, this%k4, this%cnt)
    this%cap = 0; this%mask = 0; this%n = 0
  end subroutine f_free

  ! -------- edge hash --------
  subroutine e_init(this, cap0)
    class(ehash_t), intent(inout) :: this
    integer(i8), intent(in) :: cap0
    integer(i8) :: c
    c = 1024
    do while (c < cap0)
       c = ishft(c, 1)
    end do
    this%cap = c; this%mask = c - 1; this%n = 0
    allocate(this%k1(0:c-1), this%k2(0:c-1)); this%k1 = 0
  end subroutine e_init

  subroutine e_add(this, a, b)
    class(ehash_t), intent(inout) :: this
    integer(i4), intent(in) :: a, b
    integer(i4) :: lo, hi
    integer(i8) :: h
    lo = min(a, b); hi = max(a, b)
    if (int(this%n, i8) * 3_i8 > this%cap * 2_i8) call e_grow(this)
    h = iand(mix4(lo, hi, 0, 0), this%mask)
    do
       if (this%k1(h) == 0) then
          this%k1(h) = lo; this%k2(h) = hi; this%n = this%n + 1; return
       else if (this%k1(h) == lo .and. this%k2(h) == hi) then
          return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end subroutine e_add

  subroutine e_grow(this)
    class(ehash_t), intent(inout) :: this
    integer(i4), allocatable :: o1(:), o2(:)
    integer(i8) :: i, h, oldcap
    oldcap = this%cap
    call move_alloc(this%k1, o1); call move_alloc(this%k2, o2)
    this%cap = ishft(this%cap, 1); this%mask = this%cap - 1
    allocate(this%k1(0:this%cap-1), this%k2(0:this%cap-1)); this%k1 = 0
    do i = 0, oldcap - 1
       if (o1(i) /= 0) then
          h = iand(mix4(o1(i), o2(i), 0, 0), this%mask)
          do while (this%k1(h) /= 0)
             h = iand(h + 1_i8, this%mask)
          end do
          this%k1(h) = o1(i); this%k2(h) = o2(i)
       end if
    end do
    deallocate(o1, o2)
  end subroutine e_grow

  subroutine e_free(this)
    class(ehash_t), intent(inout) :: this
    if (allocated(this%k1)) deallocate(this%k1, this%k2)
    this%cap = 0; this%mask = 0; this%n = 0
  end subroutine e_free

  ! -------- periodic point-id remap --------
  subroutine r_init(this, cap0)
    class(rmap_t), intent(inout) :: this
    integer(i8), intent(in) :: cap0
    integer(i8) :: c
    c = 1024
    do while (c < cap0)
       c = ishft(c, 1)
    end do
    this%cap = c; this%mask = c - 1; this%n = 0
    allocate(this%key(0:c-1), this%val(0:c-1)); this%key = 0
  end subroutine r_init

  !> Record raw id `k` -> merged id `v`. If `k` is already mapped, keep the
  !! smaller target (Neko merges to the minimum id); consistent records agree.
  subroutine r_set(this, k, v)
    class(rmap_t), intent(inout) :: this
    integer(i4), intent(in) :: k, v
    integer(i8) :: h
    if (this%cap == 0) call r_init(this, 1024_i8)
    if (this%n * 3_i8 > this%cap * 2_i8) call r_grow(this)
    h = iand(mix4(k, 0, 0, 0), this%mask)
    do
       if (this%key(h) == 0) then
          this%key(h) = k; this%val(h) = v; this%n = this%n + 1; return
       else if (this%key(h) == k) then
          if (v < this%val(h)) this%val(h) = v
          return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end subroutine r_set

  !> Merged id of `k` (returns `k` itself if it is not a periodic point).
  function r_get(this, k) result(v)
    class(rmap_t), intent(in) :: this
    integer(i4), intent(in) :: k
    integer(i4) :: v
    integer(i8) :: h
    if (this%cap == 0) then; v = k; return; end if
    h = iand(mix4(k, 0, 0, 0), this%mask)
    do
       if (this%key(h) == 0) then
          v = k; return
       else if (this%key(h) == k) then
          v = this%val(h); return
       end if
       h = iand(h + 1_i8, this%mask)
    end do
  end function r_get

  subroutine r_grow(this)
    class(rmap_t), intent(inout) :: this
    integer(i4), allocatable :: ok(:), ov(:)
    integer(i8) :: i, h, oldcap
    oldcap = this%cap
    call move_alloc(this%key, ok); call move_alloc(this%val, ov)
    this%cap = ishft(this%cap, 1); this%mask = this%cap - 1
    allocate(this%key(0:this%cap-1), this%val(0:this%cap-1)); this%key = 0
    do i = 0, oldcap - 1
       if (ok(i) /= 0) then
          h = iand(mix4(ok(i), 0, 0, 0), this%mask)
          do while (this%key(h) /= 0)
             h = iand(h + 1_i8, this%mask)
          end do
          this%key(h) = ok(i); this%val(h) = ov(i)
       end if
    end do
    deallocate(ok, ov)
  end subroutine r_grow

  subroutine r_free(this)
    class(rmap_t), intent(inout) :: this
    if (allocated(this%key)) deallocate(this%key, this%val)
    this%cap = 0; this%mask = 0; this%n = 0
  end subroutine r_free

end module mc_hash


!> Straight-sided hex geometry helpers (trilinear map at the 3-point GLL nodes
!! r,s,t in {-1,0,1}). Exact for non-curved elements -- for those it reproduces
!! mesh_checker's dofmap geometry and Jacobian, because a trilinear field lies
!! inside the lx=3 (quadratic) tensor space, so the spectral derivative is exact.
!! Curve records ('C'/'m') are NOT deformed here (see the tool's scope note).
module mc_geom
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: hex_min_jac, hex_gll_xyz, GLL3

  real(dp), parameter :: GLL3(3) = (/ -1.0_dp, 0.0_dp, 1.0_dp /)
  ! reference-cube corner coords for nmsh/re2 vertex order v1..v8
  real(dp), parameter :: RC(8) = (/ -1, 1, 1,-1,-1, 1, 1,-1 /)
  real(dp), parameter :: SC(8) = (/ -1,-1, 1, 1,-1,-1, 1, 1 /)
  real(dp), parameter :: TC(8) = (/ -1,-1,-1,-1, 1, 1, 1, 1 /)

contains

  !> Minimum Jacobian determinant of the trilinear map over the 27 GLL points.
  !! <= 0 flags an inverted/degenerate (tangled) element.
  function hex_min_jac(cx, cy, cz) result(jmin)
    real(dp), intent(in) :: cx(8), cy(8), cz(8)
    real(dp) :: jmin, r, s, t, dr, ds, dt
    real(dp) :: xr, xs, xt, yr, ys, yt, zr, zs, zt, jac
    integer :: ir, is, it, k
    jmin = huge(1.0_dp)
    do it = 1, 3
       t = GLL3(it)
       do is = 1, 3
          s = GLL3(is)
          do ir = 1, 3
             r = GLL3(ir)
             xr = 0; xs = 0; xt = 0; yr = 0; ys = 0; yt = 0; zr = 0; zs = 0; zt = 0
             do k = 1, 8
                dr = 0.125_dp * RC(k) * (1 + s*SC(k)) * (1 + t*TC(k))
                ds = 0.125_dp * SC(k) * (1 + r*RC(k)) * (1 + t*TC(k))
                dt = 0.125_dp * TC(k) * (1 + r*RC(k)) * (1 + s*SC(k))
                xr = xr + dr*cx(k); xs = xs + ds*cx(k); xt = xt + dt*cx(k)
                yr = yr + dr*cy(k); ys = ys + ds*cy(k); yt = yt + dt*cy(k)
                zr = zr + dr*cz(k); zs = zs + ds*cz(k); zt = zt + dt*cz(k)
             end do
             jac = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr)
             if (jac < jmin) jmin = jac
          end do
       end do
    end do
  end function hex_min_jac

  !> The 27 GLL-point coordinates of the trilinear map, in nek fld order
  !! (k outer, j, i inner -> index = i + 3*(j-1) + 9*(k-1)).
  subroutine hex_gll_xyz(cx, cy, cz, gx, gy, gz)
    real(dp), intent(in) :: cx(8), cy(8), cz(8)
    real(dp), intent(out) :: gx(27), gy(27), gz(27)
    real(dp) :: r, s, t, nshp
    integer :: ir, is, it, k, idx
    idx = 0
    do it = 1, 3
       t = GLL3(it)
       do is = 1, 3
          s = GLL3(is)
          do ir = 1, 3
             r = GLL3(ir); idx = idx + 1
             gx(idx) = 0; gy(idx) = 0; gz(idx) = 0
             do k = 1, 8
                nshp = 0.125_dp * (1 + r*RC(k)) * (1 + s*SC(k)) * (1 + t*TC(k))
                gx(idx) = gx(idx) + nshp*cx(k)
                gy(idx) = gy(idx) + nshp*cy(k)
                gz(idx) = gz(idx) + nshp*cz(k)
             end do
          end do
       end do
    end do
  end subroutine hex_gll_xyz

end module mc_geom


program mesh_checker_light
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, &
       i4 => int32, i8 => int64, i1 => int8, error_unit
  use mc_hash, only: fhash_t, ehash_t, rmap_t
  use mc_geom, only: hex_min_jac, hex_gll_xyz
  implicit none

  ! Neko facet (1..6) corners in nmsh vertex indices (= face_nodes composed with
  ! the nmsh->mesh read swap [1,2,4,3,5,6,8,7]); same table validated byte-exact
  ! for periodic glb_pt_ids in rea2nbin_light.
  integer(i4), parameter :: FACE_RE2(4,6) = reshape((/ &
       1,5,8,4,   2,6,7,3,   1,2,6,5,   4,3,7,8,   1,2,3,4,   5,6,7,8 /), (/4,6/))
  ! The 12 hex edges in nmsh vertex indices (edge_nodes composed with the swap).
  integer(i4), parameter :: EDGE_RE2(2,12) = reshape((/ &
       1,2,  3,4,  5,6,  7,8,  1,4,  2,3,  5,8,  6,7,  1,5,  2,6,  4,8,  3,7 /), (/2,12/))

  character(len=1024) :: fin
  character(len=64) :: argbuf
  integer :: uin, ios, i, f, argc
  integer :: iphase, nphase                    ! user-facing progress counters
  integer(i4) :: nelv, gdim, elid, vidx(8), nz, ztype, zlabel
  integer(i4) :: ze, zf
  integer(i4) :: max_vidx, periodic_size, n_unlabeled
  integer(i4), allocatable :: elid_of_pos(:), pos_of_elid(:)
  integer(i4), allocatable :: ev(:,:)          ! (8,nelv): the 8 nmsh point ids per element
  integer(i1), allocatable :: ftype(:,:)       ! (6,nelv): 0 none, 1 labeled, 2 periodic
  integer(i1), allocatable :: flabel(:,:)      ! (6,nelv): labeled-zone index 1..20 or 0
  integer(i4) :: labeled_cnt(20)
  real(dp) :: xyz(3), lo(3), hi(3)
  real(dp) :: cx8(8), cy8(8), cz8(8), jmin
  integer(i8) :: cnt
  integer(i4) :: n_neg_jac, first_bad
  real(dp) :: min_jac
  logical :: do_jac, do_zone
  type(fhash_t) :: fh
  type(ehash_t) :: eh
  type(rmap_t) :: rm
  logical :: failed

  argc = command_argument_count()
  if (argc < 1) then
     write(*,*) 'Usage: mesh_checker_light mesh.nmsh [--jacobian] [--write-zone-indices]'
     write(*,*) '  --jacobian            : also check for negative/zero Jacobians'
     write(*,*) '                          (streaming; exact for straight-sided elements).'
     write(*,*) '  --write-zone-indices  : also write zone_indices.fld (boundary faces'
     write(*,*) '                          marked by zone index) -- streamed, low memory.'
     write(*,*) '  Reports mesh size, bounding box, periodic/labeled zones, and'
     write(*,*) '  the number of unlabeled external faces (0 for a valid mesh).'
     stop 1
  end if
  call get_command_argument(1, fin)
  do_jac = .false.; do_zone = .false.
  do i = 2, argc
     call get_command_argument(i, argbuf)
     select case (trim(argbuf))
     case ('--jacobian');            do_jac = .true.
     case ('--write-zone-indices');  do_zone = .true.; do_jac = .true.
     case default
        write(error_unit,*) 'Unknown option: ', trim(argbuf); stop 1
     end select
  end do

  nphase = 3
  if (do_zone) nphase = 4
  iphase = 0

  write(*,'(a)') ''
  write(*,'(a)') '     _  __  ____  __ __  ____'
  write(*,'(a)') '    / |/ / / __/ / //_/ / __ \'
  write(*,'(a)') '   /    / / _/  / ,<   / /_/ /'
  write(*,'(a)') '  /_/|_/ /___/ /_/|_|  \____/   mesh_checker_light'
  write(*,'(a)') ''

  open(newunit=uin, file=trim(fin), access='stream', form='unformatted', &
       status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(error_unit,*) 'Error: cannot open ', trim(fin); stop 1
  end if
  read(uin, iostat=ios) nelv, gdim
  if (ios /= 0) then
     write(error_unit,*) 'Error: ', trim(fin), ' is not a Neko .nmsh file ', &
          '(short header)'; stop 1
  end if
  if (gdim /= 3) then
     write(error_unit,*) 'Error: mesh_checker_light supports 3D meshes only (gdim=', &
          gdim, '). Use the full mesh_checker for 2D.'; stop 1
  end if

  allocate(ftype(6, nelv)); ftype = 0_i1
  if (do_zone) then
     allocate(flabel(6, nelv)); flabel = 0_i1
  end if
  allocate(elid_of_pos(nelv), pos_of_elid(nelv))
  labeled_cnt = 0; max_vidx = 0; periodic_size = 0; n_unlabeled = 0; failed = .false.
  n_neg_jac = 0; first_bad = 0; min_jac = huge(1.0_dp)
  lo = huge(1.0_dp); hi = -huge(1.0_dp)
  call fh%init(max(1024_i8, int(nelv, i8)))     ! merged unique faces, grown on demand
  call eh%init(max(1024_i8, int(nelv, i8)))
  allocate(ev(8, nelv))

  ! ---- pass 1 (streaming): elements -> connectivity, bbox, Jacobians ----
  ! The face/edge tables are NOT built here: for a periodic mesh, Neko's
  ! mesh_checker reports the counts AFTER merging periodic point ids
  ! (glb_mfcs/glb_meds), so we must first read the zones (which carry the merged
  ! glb_pt_ids) before we can hash faces/edges. We stash the connectivity in `ev`
  ! and hash it once the periodic remap is known.
  iphase = iphase + 1
  if (do_jac) then
     write(*,'(a,i0,a,i0,a)') '  [', iphase, '/', nphase, &
          '] reading elements (bounding box, connectivity, Jacobians) ...'
  else
     write(*,'(a,i0,a,i0,a)') '  [', iphase, '/', nphase, &
          '] reading elements (bounding box, connectivity) ...'
  end if
  do i = 1, nelv
     read(uin, iostat=ios) elid
     if (ios /= 0) then
        write(error_unit,*) 'Error: truncated element section (element ', i, &
             ' of ', nelv, ')'; stop 1
     end if
     do f = 1, 8
        read(uin, iostat=ios) vidx(f), xyz(1), xyz(2), xyz(3)
        if (ios /= 0) then
           write(error_unit,*) 'Error: truncated element section (element ', i, &
                ' of ', nelv, ')'; stop 1
        end if
        lo = min(lo, xyz); hi = max(hi, xyz)
        if (vidx(f) > max_vidx) max_vidx = vidx(f)
        cx8(f) = xyz(1); cy8(f) = xyz(2); cz8(f) = xyz(3)
     end do
     elid_of_pos(i) = elid
     ev(:, i) = vidx
     if (do_jac) then
        jmin = hex_min_jac(cx8, cy8, cz8)
        if (jmin < min_jac) min_jac = jmin
        if (jmin <= 0.0_dp) then
           n_neg_jac = n_neg_jac + 1
           if (first_bad == 0) first_bad = i
        end if
     end if
  end do

  ! inverse element-id map (nmsh zone %e is a global element id)
  pos_of_elid = 0
  do i = 1, nelv
     if (elid_of_pos(i) < 1 .or. elid_of_pos(i) > nelv) then
        write(error_unit,*) 'Error: element id out of [1,nelv]: ', elid_of_pos(i); stop 1
     end if
     pos_of_elid(elid_of_pos(i)) = i
  end do
  do i = 1, nelv                 ! element ids must be a permutation of 1..nelv
     if (pos_of_elid(i) == 0) then
        write(error_unit,*) 'Error: element ids are not a permutation of ', &
             '1..nelv (duplicate or missing id ', i, ')'; stop 1
     end if
  end do

  ! ---- zones: mark labeled (type 7) and periodic (type 5) facets ----
  iphase = iphase + 1
  write(*,'(a,i0,a,i0,a)') '  [', iphase, '/', nphase, &
       '] reading and marking boundary zones ...'
  block
    integer(i4) :: zpe, zpf, g1, g2, g3, g4, pos
    read(uin, iostat=ios) nz
    if (ios /= 0) then
       write(error_unit,*) 'Error: truncated file (missing zone count)'; stop 1
    end if
    do i = 1, nz
       ! nmsh_zone_t: e, f, p_e, p_f, glb_pt_ids(4), type   (9 x i4)
       read(uin, iostat=ios) ze, zf, zpe, zpf, g1, g2, g3, g4, ztype
       if (ios /= 0) then
          write(error_unit,*) 'Error: truncated zone section (record ', i, &
               ' of ', nz, ')'; stop 1
       end if
       if (ze < 1 .or. ze > nelv .or. zf < 1 .or. zf > 6) cycle
       pos = pos_of_elid(ze)
       if (ztype == 5) then
          periodic_size = periodic_size + 1
          ftype(zf, pos) = 2_i1
          ! Record the periodic point merge exactly as Neko's apply_periodic_facet
          ! does: corner k of facet zf (nmsh slot FACE_RE2(k,zf)) takes the merged
          ! id glb_pt_ids(k). The file stores glb_pt_ids in that same face_nodes
          ! order (see mesh_reset_periodic_ids/hex_facet_order).
          call rm%set(ev(FACE_RE2(1,zf), pos), g1)
          call rm%set(ev(FACE_RE2(2,zf), pos), g2)
          call rm%set(ev(FACE_RE2(3,zf), pos), g3)
          call rm%set(ev(FACE_RE2(4,zf), pos), g4)
       else if (ztype == 7) then
          zlabel = zpf                       ! label lives in the p_f field
          if (zlabel >= 1 .and. zlabel <= 20) labeled_cnt(zlabel) = labeled_cnt(zlabel) + 1
          ftype(zf, pos) = 1_i1
          if (do_zone .and. zlabel >= 1 .and. zlabel <= 20) flabel(zf, pos) = int(zlabel, i1)
       end if
    end do
  end block
  close(uin)

  ! ---- apply the periodic merge to the stored connectivity (once) ----
  ! After this, every point carries its merged id, so a periodic facet-pair
  ! collapses to a single shared (internal) face/edge -- reproducing the
  ! glb_mfcs/glb_meds that Neko's mesh_checker reports.
  if (rm%n > 0) then
     do i = 1, nelv
        do f = 1, 8
           ev(f, i) = rm%get(ev(f, i))
        end do
     end do
  end if

  ! ---- build the (merged) face/edge tables and scan unlabeled external faces ----
  ! An external face appears exactly once across the mesh (merged face count == 1);
  ! a periodic facet-pair now shares one internal face (count == 2), so it is
  ! excluded automatically. Unlabeled external = external AND facet not labeled.
  iphase = iphase + 1
  write(*,'(a,i0,a,i0,a)') '  [', iphase, '/', nphase, &
       '] building face/edge tables and scanning unlabeled external faces ...'
  do i = 1, nelv
     do f = 1, 6
        call fh%add(ev(FACE_RE2(1,f),i), ev(FACE_RE2(2,f),i), &
                    ev(FACE_RE2(3,f),i), ev(FACE_RE2(4,f),i))
     end do
     do f = 1, 12
        call eh%add(ev(EDGE_RE2(1,f),i), ev(EDGE_RE2(2,f),i))
     end do
  end do
  do i = 1, nelv
     do f = 1, 6
        cnt = int(fh%cnt_of(ev(FACE_RE2(1,f),i), ev(FACE_RE2(2,f),i), &
                            ev(FACE_RE2(3,f),i), ev(FACE_RE2(4,f),i)), i8)
        if (cnt == 1 .and. ftype(f, i) == 0_i1) n_unlabeled = n_unlabeled + 1
     end do
  end do

  ! ---- report (mirrors mesh_checker) ----
  write(*,*) ''
  write(*,*) '--------------Size-------------'
  write(*,'(A,I0)') ' Number of elements: ', nelv
  write(*,'(A,I0)') ' Number of points:   ', max_vidx
  write(*,'(A,I0)') ' Number of faces:    ', fh%n
  write(*,'(A,I0)') ' Number of edges:    ', eh%n
  write(*,*) 'Bounding box:'
  write(*,'(A,2(g14.6,1X))') '    x', lo(1), hi(1)
  write(*,'(A,2(g14.6,1X))') '    y', lo(2), hi(2)
  write(*,'(A,2(g14.6,1X))') '    z', lo(3), hi(3)
  write(*,*) ''
  write(*,*) '--------------Zones------------'
  write(*,'(A,I0)') ' Number of periodic faces: ', periodic_size
  write(*,*) ''
  write(*,*) 'Labeled zones: '
  do i = 1, 20
     if (labeled_cnt(i) > 0) &
          write(*,'(A,I2,A,I0,A)') '    Zone ', i, ': ', labeled_cnt(i), ' faces'
  end do

  if (do_jac) then
     write(*,*) ''
     write(*,*) '------------Jacobian----------'
     write(*,'(A,g14.6)') ' Min Jacobian (straight-sided): ', min_jac
     if (n_neg_jac > 0) then
        failed = .true.
        write(*,'(A,I0,A,I0,A)') ' Error: Found ', n_neg_jac, &
             ' element(s) with a negative/zero Jacobian (first at element ', first_bad, ').'
     else
        write(*,*) ' No negative/zero Jacobians.'
     end if
  end if

  if (n_unlabeled > 0) then
     failed = .true.
     write(*,'(A,I0,A)') ' Error: Found ', n_unlabeled, ' unlabeled external faces.'
  end if

  if (do_zone) then
     iphase = iphase + 1
     write(*,'(a,i0,a,i0,a)') '  [', iphase, '/', nphase, &
          '] writing zone_indices.fld ...'
     call write_zone_indices_fld(trim(fin), nelv, flabel)
  end if

  write(*,*) 'Done'
  call fh%free(); call eh%free(); call rm%free()
  deallocate(ftype, elid_of_pos, pos_of_elid)
  if (allocated(ev)) deallocate(ev)
  if (allocated(flabel)) deallocate(flabel)

  if (failed) then
     write(error_unit,*) 'Mesh check failed with one or several errors.'
     stop 1
  end if

contains

  !> Write zone_indices.fld -- a single-precision Neko/NEKTON field file with the
  !! mesh (straight-sided trilinear geometry at the 3x3x3 GLL nodes) and one
  !! scalar field: the labeled-zone index at each GLL point (a point on a labeled
  !! facet takes that zone's index; the highest index wins where facets meet; 0
  !! in the interior). Streamed element-by-element -- O(1) memory, no dofmap.
  !! Geometry is straight-sided (curve records are not deformed); the zone labels
  !! -- the point of the file -- are exact. Also writes the .nek5000 companion.
  subroutine write_zone_indices_fld(meshfile, nel, flab)
    character(len=*), intent(in) :: meshfile
    integer(i4), intent(in) :: nel
    integer(i1), intent(in) :: flab(:,:)
    integer :: u, uf, kk, e, ff
    integer(i4) :: eidx
    real(sp) :: gx(27), gy(27), gz(27), sc(27), bb(6), test_pattern
    real(sp), allocatable :: scmn(:), scmx(:)
    real(dp) :: c8x(8), c8y(8), c8z(8), gdx(27), gdy(27), gdz(27)
    character(len=132) :: hdr

    write(hdr, 1) 4, 3, 3, 3, nel, nel, 0.0_dp, 1, 1, 1, 'XS01      '
1   format('#std', 1x, i1, 1x, i2, 1x, i2, 1x, i2, 1x, i10, 1x, i10, &
         1x, e20.13, 1x, i9, 1x, i6, 1x, i6, 1x, a10)

    open(newunit=u, file='zone_indices.fld', access='stream', form='unformatted', &
         status='replace', action='write')
    write(u) hdr
    test_pattern = 6.54321_sp
    write(u) test_pattern
    do e = 1, nel                                    ! element global-id list
       eidx = e; write(u) eidx
    end do

    ! mesh block: for each element, gx(1:27), gy(1:27), gz(1:27) (single prec)
    open(newunit=uf, file=meshfile, access='stream', form='unformatted', &
         status='old', action='read')
    read(uf) eidx, eidx                              ! skip nelv, gdim
    do e = 1, nel
       call read_corners(uf, c8x, c8y, c8z)
       call hex_gll_xyz(c8x, c8y, c8z, gdx, gdy, gdz)
       gx = real(gdx, sp); gy = real(gdy, sp); gz = real(gdz, sp)
       write(u) gx, gy, gz
    end do
    close(uf)

    ! scalar block: labeled-zone index at each of the 27 GLL points
    allocate(scmn(nel), scmx(nel))
    do e = 1, nel
       do kk = 1, 27
          sc(kk) = 0.0_sp
       end do
       do ff = 1, 6
          if (flab(ff, e) /= 0_i1) call mark_facet_points(sc, ff, real(flab(ff, e), sp))
       end do
       write(u) sc
       scmn(e) = minval(sc); scmx(e) = maxval(sc)
    end do

    ! per-element bounding-box metadata (min/max of x, y, z), single precision
    open(newunit=uf, file=meshfile, access='stream', form='unformatted', &
         status='old', action='read')
    read(uf) eidx, eidx
    do e = 1, nel
       call read_corners(uf, c8x, c8y, c8z)
       call hex_gll_xyz(c8x, c8y, c8z, gdx, gdy, gdz)
       bb(1) = real(minval(gdx), sp); bb(2) = real(maxval(gdx), sp)
       bb(3) = real(minval(gdy), sp); bb(4) = real(maxval(gdy), sp)
       bb(5) = real(minval(gdz), sp); bb(6) = real(maxval(gdz), sp)
       write(u) bb
    end do
    close(uf)

    ! per-element (min, max) metadata for the scalar field: Neko's fld writer
    ! appends this block after the geometry bounding boxes (see
    ! fld_file_write_metadata_scalar in src/io/fld_file.f90), so match it
    do e = 1, nel
       write(u) scmn(e), scmx(e)
    end do
    deallocate(scmn, scmx)
    close(u)

    ! ParaView companion
    open(newunit=uf, file='zone_indices.nek5000', status='replace', action='write')
    write(uf, '(A)') 'filetemplate: zone_indices.fld'
    write(uf, '(A)') 'firsttimestep: 1'
    write(uf, '(A)') 'numtimesteps: 1'
    close(uf)
    write(*,*) 'Wrote zone_indices.fld (+ .nek5000 companion)'
  end subroutine write_zone_indices_fld

  !> Read one nmsh hex element (el_idx + 8*(v_idx,x,y,z)); return the 8 corners.
  subroutine read_corners(u, c8x, c8y, c8z)
    integer, intent(in) :: u
    real(dp), intent(out) :: c8x(8), c8y(8), c8z(8)
    integer(i4) :: eidx, vi
    integer :: j
    read(u) eidx
    do j = 1, 8
       read(u) vi, c8x(j), c8y(j), c8z(j)
    end do
  end subroutine read_corners

  !> Set the 9 GLL points lying on facet @a ff to value @a v (max-wins).
  subroutine mark_facet_points(sc, ff, v)
    real(sp), intent(inout) :: sc(27)
    integer, intent(in) :: ff
    real(sp), intent(in) :: v
    integer :: ir, is, it, idx
    logical :: on
    do it = 1, 3
       do is = 1, 3
          do ir = 1, 3
             idx = ir + 3*(is-1) + 9*(it-1)
             on = (ff == 1 .and. ir == 1) .or. (ff == 2 .and. ir == 3) .or. &
                  (ff == 3 .and. is == 1) .or. (ff == 4 .and. is == 3) .or. &
                  (ff == 5 .and. it == 1) .or. (ff == 6 .and. it == 3)
             if (on .and. v > sc(idx)) sc(idx) = v
          end do
       end do
    end do
  end subroutine mark_facet_points

end program mesh_checker_light
