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
!> rea2nbin_light -- a lightweight, CPU-only re2 -> nmsh mesh converter.
!!
!! Neko's `contrib/rea2nbin` builds a full `mesh_t` using Neko's general-purpose,
!! production hash tables. Those are robust and reused everywhere in the code, but
!! they are sized generously and store polymorphic keys, so they need a lot of
!! memory for a large mesh. This tool is a lightweight reimplementation of just the
!! part needed for a serial re2 -> nmsh conversion: the point de-duplication and
!! the record streaming, using plain typed arrays and no Neko initialisation. The
!! result is the same `.nmsh`, produced in a small fraction of the memory.
!!
!! Because it does not link Neko, it needs only a Fortran compiler (no MPI, no
!! Neko library) and runs anywhere -- in particular you do not need a separate
!! CPU build of Neko. The `contrib` tools are CPU-only and serial, so on a Neko
!! configured for GPUs you would otherwise have to rebuild it for CPU just to
!! convert a mesh; this sidesteps that entirely.
!!
!! Supports labeled/named BCs, periodic BCs (any number of periodic directions:
!! channel = 2, TGV/HIT = 3), and curved elements ('C' circle / 'm' midside
!! records, passed through to nmsh curve records exactly as Neko does).
!!
!! Tested against the meshes shipped with Neko and a real wall-resolved pipe, but
!! it remains your responsibility to confirm the result is correct for your own
!! mesh -- run one of the `make check` suites on a case you trust first. See
!! rea2nbin_light.md for the memory model, validation, and scope.
module coord_hash
  use, intrinsic :: iso_fortran_env, only: dp => real64, i4 => int32, i8 => int64
  implicit none
  private
  type, public :: coord_hash_t
     real(dp), allocatable :: kx(:), ky(:), kz(:)
     integer(i4), allocatable :: id(:)
     integer(i8) :: cap = 0
     integer(i8) :: mask = 0
     integer(i4) :: n = 0
   contains
     procedure :: init  => ch_init
     procedure :: get_or_add => ch_get_or_add
     procedure :: free  => ch_free
  end type coord_hash_t
contains

  subroutine ch_init(this, init_cap)
    class(coord_hash_t), intent(inout) :: this
    integer(i8), intent(in) :: init_cap
    integer(i8) :: c
    c = 1024
    do while (c < init_cap)
       c = ishft(c, 1)
    end do
    this%cap = c; this%mask = c - 1; this%n = 0
    allocate(this%kx(0:c-1), this%ky(0:c-1), this%kz(0:c-1))
    allocate(this%id(0:c-1)); this%id = 0
  end subroutine ch_init
  pure function ch_hash(this, x, y, z) result(h)
    class(coord_hash_t), intent(in) :: this
    real(dp), intent(in) :: x, y, z
    integer(i8) :: h, a, b, c
    a = transfer(x, a); b = transfer(y, b); c = transfer(z, c)
    h = ieor(a, ishft(a, -33)) * (-49064778989728563_i8)
    h = ieor(h, b) + (-7046029254386353131_i8)
    h = ieor(h, ishft(h, -29)) * (-4417276706812531889_i8)
    h = ieor(h, c)
    h = ieor(h, ishft(h, -32))
    h = iand(h, this%mask)
  end function ch_hash
  subroutine ch_get_or_add(this, x, y, z, idx, is_new)
    class(coord_hash_t), intent(inout) :: this
    real(dp), intent(in) :: x, y, z
    integer(i4), intent(out) :: idx
    logical, intent(out) :: is_new
    integer(i8) :: s
    if (int(this%n, i8) * 3_i8 > this%cap * 2_i8) call ch_grow(this)
    s = ch_hash(this, x, y, z)
    do
       if (this%id(s) == 0) then
          this%n = this%n + 1
          this%kx(s) = x; this%ky(s) = y; this%kz(s) = z
          this%id(s) = this%n; idx = this%n; is_new = .true.; return
       ! exact REAL(dp) equality is intentional (-Wcompare-reals): coordinates are
       ! copied from the file, never computed, and Neko's own point dedup
       ! (htable_pt_t) also compares them bit-exactly
       else if (this%kx(s) == x .and. this%ky(s) == y .and. this%kz(s) == z) then
          idx = this%id(s); is_new = .false.; return
       end if
       s = iand(s + 1_i8, this%mask)
    end do
  end subroutine ch_get_or_add
  subroutine ch_grow(this)
    class(coord_hash_t), intent(inout) :: this
    real(dp), allocatable :: ox(:), oy(:), oz(:)
    integer(i4), allocatable :: oid(:)
    integer(i8) :: i, s, oldcap
    oldcap = this%cap
    call move_alloc(this%kx, ox); call move_alloc(this%ky, oy)
    call move_alloc(this%kz, oz); call move_alloc(this%id, oid)
    this%cap = ishft(this%cap, 1); this%mask = this%cap - 1
    allocate(this%kx(0:this%cap-1), this%ky(0:this%cap-1), this%kz(0:this%cap-1))
    allocate(this%id(0:this%cap-1)); this%id = 0
    do i = 0, oldcap - 1
       if (oid(i) /= 0) then
          s = ch_hash(this, ox(i), oy(i), oz(i))
          do while (this%id(s) /= 0)
             s = iand(s + 1_i8, this%mask)
          end do
          this%kx(s) = ox(i); this%ky(s) = oy(i); this%kz(s) = oz(i)
          this%id(s) = oid(i)
       end if
    end do
    deallocate(ox, oy, oz, oid)
  end subroutine ch_grow
  subroutine ch_free(this)
    class(coord_hash_t), intent(inout) :: this
    if (allocated(this%kx)) deallocate(this%kx, this%ky, this%kz, this%id)
    this%cap = 0; this%mask = 0; this%n = 0
  end subroutine ch_free
end module coord_hash


program rea2nbin_light
  use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, &
       i4 => int32, i8 => int64, error_unit
  use coord_hash, only: coord_hash_t
  implicit none

  ! re2 face -> Neko symmetric facet
  integer(i4), parameter :: FACET_MAP(6) = (/3, 2, 4, 1, 5, 6/)
  ! the 4 corners of each Neko facet (facet_order), in re2-vertex indices
  integer(i4), parameter :: FACE_RE2(4,6) = reshape((/ &
       1,5,8,4,   2,6,7,3,   1,2,6,5,   4,3,7,8,   1,2,3,4,   5,6,7,8 /), (/4,6/))
  integer(i4), parameter :: MAXZ = 20
  real(dp) :: PTOL = 1.0e-7_dp   ! periodic match tol; overridable via NEKO_PERIODIC_TOL

  character(len=1024) :: fin, fout
  character(len=80)   :: hdr
  character(len=5)    :: ver
  integer(i4) :: nel, ndim, nelv, nbc4
  logical :: v2
  real(sp) :: endian
  real(sp), parameter :: RE2_ENDIAN_TEST = 6.54321_sp
  integer :: uin, uout, ios, i, j
  type(coord_hash_t) :: ch
  real(dp) :: rg, x8(8), y8(8), z8(8)
  real(sp) :: rg4, x4(8), y4(8), z4(8)
  integer(i4) :: vidx(8)
  logical :: isnew
  integer(i4) :: gdim_out
  ! per-element vertex ids (re2 order) + id -> coordinate map (for periodic merge)
  integer(i4), allocatable :: elem_vid(:,:)
  real(dp), allocatable :: cx(:), cy(:), cz(:)
  integer(i4) :: cap_id

  write(*,'(a)') ''
  write(*,'(a)') '     _  __  ____  __ __  ____'
  write(*,'(a)') '    / |/ / / __/ / //_/ / __ \'
  write(*,'(a)') '   /    / / _/  / ,<   / /_/ /'
  write(*,'(a)') '  /_/|_/ /___/ /_/|_|  \____/   rea2nbin_light  (re2 -> nmsh)'
  write(*,'(a)') ''

  if (command_argument_count() < 1 .or. command_argument_count() > 2) then
     write(*,*) 'Usage: rea2nbin_light <mesh.re2> [<mesh.nmsh>]'; stop 1
  end if
  call get_command_argument(1, fin)
  if (command_argument_count() == 2) then
     call get_command_argument(2, fout)
  else
     fout = fin(1:len_trim(fin)-3) // 'nmsh'
  end if

  ! Periodic match tolerance: honour NEKO_PERIODIC_TOL like Neko (default 1e-7).
  block
    character(len=255) :: tol_str
    integer :: tol_len
    integer :: tios
    call get_environment_variable('NEKO_PERIODIC_TOL', tol_str, tol_len)
    if (tol_len > len(tol_str)) then
       write(error_unit,*) 'Error: NEKO_PERIODIC_TOL value too long'; stop 1
    end if
    if (tol_len > 0) then
       read(tol_str(1:tol_len), *, iostat=tios) PTOL
       if (tios /= 0 .or. PTOL <= 0.0_dp) then
          write(error_unit,*) 'Error: invalid NEKO_PERIODIC_TOL value: ', &
               tol_str(1:tol_len)
          stop 1
       end if
    end if
  end block

  open(newunit=uin, file=trim(fin), access='stream', form='unformatted', &
       status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(error_unit,*) 'Error: cannot open ', trim(fin); stop 1
  end if
  read(uin, iostat=ios) hdr; call chk_read(ios)
  ver = hdr(1:5); v2 = .false.
  if (ver == '#v004') then
     read(hdr, '(a5,i16,i3,i16,i4)') ver, nel, ndim, nelv, nbc4; v2 = .true.
  else if (ver == '#v002' .or. ver == '#v003') then
     read(hdr, '(a5,i9,i3,i9)') ver, nel, ndim, nelv; v2 = .true.
  else if (ver == '#v001') then
     read(hdr, '(a5,i9,i3,i9)') ver, nel, ndim, nelv; v2 = .false.
  else
     write(error_unit,*) 'Error: unknown re2 version ', ver; stop 1
  end if
  read(uin, iostat=ios) endian; call chk_read(ios)
  if (abs(endian - RE2_ENDIAN_TEST) > 1.0e-4_sp) then
     write(error_unit,*) 'Error: byte-swapped re2 not supported'; stop 1
  end if
  if (ndim /= 3) then
     write(error_unit,*) 'Error: this tool supports 3D (hex) meshes only'; stop 1
  end if

  write(*,'(a,a)')            '  input     : ', trim(fin)
  write(*,'(a,a)')            '  output    : ', trim(fout)
  write(*,'(a,i0,a,a,a)')     '  mesh      : ', nelv, ' hex elements  (format ', ver, ')'
  write(*,'(a)')              '  [1/3] reading elements and de-duplicating points ...'

  call ch%init(max(1024_i8, int(nelv, i8) * 2_i8))
  allocate(elem_vid(8, nelv))
  cap_id = max(1024, 2*nelv)
  allocate(cx(cap_id), cy(cap_id), cz(cap_id))

  open(newunit=uout, file=trim(fout), access='stream', form='unformatted', &
       status='replace', action='write')
  gdim_out = ndim
  write(uout) nelv, gdim_out

  ! Stream elements: read re2 record, dedup 8 vertices, write nmsh record,
  ! and record per-element vertex ids + id->coordinate for periodic merging.
  do i = 1, nelv
     if (v2) then
        read(uin, iostat=ios) rg, x8, y8, z8; call chk_read(ios)
     else
        read(uin, iostat=ios) rg4, x4, y4, z4; call chk_read(ios)
        x8 = real(x4, dp); y8 = real(y4, dp); z8 = real(z4, dp)
     end if
     do j = 1, 8
        call ch%get_or_add(x8(j), y8(j), z8(j), vidx(j), isnew)
        elem_vid(j, i) = vidx(j)
        if (isnew) then
           if (vidx(j) > cap_id) call grow_coords(cx, cy, cz, cap_id)
           cx(vidx(j)) = x8(j); cy(vidx(j)) = y8(j); cz(vidx(j)) = z8(j)
        end if
     end do
     ! nmsh hex record: el_idx, then 8 x (v_idx, x, y, z). Vertex order is
     ! identity w.r.t. re2 (the add_element swap and vcyc_to_sym cancel).
     write(uout) i
     do j = 1, 8
        write(uout) vidx(j), x8(j), y8(j), z8(j)
     end do
  end do

  write(*,'(a,i0,a)')        '        ', ch%n, ' unique points'
  write(*,'(a)')             '  [2/3] reading curves and boundary conditions ...'

  ! ---- tail: re2 curves + re2 BCs -> nmsh periodic & labeled zones ----
  block
    integer(i4) :: ncurv, nbcs, k, lbl, sym_facet, el, iface, nz
    real(dp) :: rc_elem, rc_zone, rc_pt(5), t2, bc_elem, bc_face, bc_data(5)
    real(sp) :: sc_pt(5), sbc_data(5)
    integer(i4) :: ci4, ei4, fi4
    character(len=8) :: btype8
    character(len=4) :: btype4
    character(len=8) :: bt
    ! curved-element records (stored compactly, re-emitted per element)
    integer(i4), allocatable :: rc_el(:), rc_ed(:), rc_ty(:)
    real(dp),    allocatable :: rc_d(:,:)
    integer(i4), allocatable :: ccnt(:), order(:)
    integer(i4) :: ncurves, acc, cc, cur_e, idx
    real(dp) :: cbuf(5,12)
    integer(i4) :: tbuf(12)
    logical :: curve_skip
    ! labeled zones
    integer(i4), allocatable :: z_e(:), z_f(:), z_lbl(:)
    integer(i4) :: nzl, cap_z, cur_internal, user_off, named_map(6)
    ! periodic facets
    integer(i4), allocatable :: p_e(:), p_f(:), p_pe(:), p_pf(:)
    integer(i4) :: npf, cap_p
    ! merged point ids (in-place min, exactly like Neko's create_periodic_ids)
    integer(i4), allocatable :: pid(:)
    integer(i4) :: gpid(4), pass

    if (v2) then
       read(uin, iostat=ios) t2; call chk_read(ios); ncurv = int(t2)
    else
       read(uin, iostat=ios) ncurv; call chk_read(ios)
    end if
    ! Read curve records into compact arrays. Type mapping mirrors Neko's
    ! re2_file_read_curve: 'C' (circle) -> 3, 'm' (midside) -> 4; 's'/'e' or any
    ! other type make Neko treat the whole mesh as non-curved, so a single
    ! unsupported record sets curve_skip and drops all curves (ncurves = 0).
    curve_skip = .false.
    if (ncurv > 0) then
       allocate(rc_el(ncurv), rc_ed(ncurv), rc_ty(ncurv), rc_d(5, ncurv))
       do k = 1, ncurv
          if (v2) then
             read(uin, iostat=ios) rc_elem, rc_zone, rc_pt, btype8; call chk_read(ios)
             rc_el(k) = int(rc_elem); rc_ed(k) = int(rc_zone)
             rc_d(1:5, k) = rc_pt; bt = btype8
          else
             read(uin, iostat=ios) ei4, ci4, sc_pt, btype4; call chk_read(ios)
             rc_el(k) = ei4; rc_ed(k) = ci4
             rc_d(1:5, k) = real(sc_pt, dp); bt = btype4
          end if
          if (rc_el(k) < 1 .or. rc_el(k) > nelv) then
             write(error_unit,'(a,i0)') &
                  'Error: curve record references element out of range: ', rc_el(k); stop 1
          end if
          if (rc_ed(k) < 1 .or. rc_ed(k) > 12) then
             write(error_unit,'(a,i0)') &
                  'Error: curve record edge index out of [1,12]: ', rc_ed(k); stop 1
          end if
          ! Neko decodes only the FIRST character (its chtemp is character(len=1)),
          ! so match on bt(1:1) -- robust to NUL/junk-padded type fields too.
          select case (bt(1:1))
          case ('C')
             rc_ty(k) = 3
          case ('m')
             rc_ty(k) = 4
          case default
             rc_ty(k) = 0; curve_skip = .true.
          end select
       end do
    end if

    if (v2) then
       read(uin, iostat=ios) t2; call chk_read(ios); nbcs = int(t2)
    else
       read(uin, iostat=ios) nbcs; call chk_read(ios)
    end if

    cap_z = max(64, nbcs); cap_p = max(64, nbcs)
    allocate(z_e(cap_z), z_f(cap_z), z_lbl(cap_z))
    allocate(p_e(cap_p), p_f(cap_p), p_pe(cap_p), p_pf(cap_p))
    nzl = 0; npf = 0; named_map = 0; cur_internal = 1; user_off = 0
    allocate(pid(ch%n))
    do i = 1, ch%n
       pid(i) = i
    end do

    block
      integer(i4), allocatable :: b_el(:), b_fc(:)
      character(len=8), allocatable :: b_ty(:)
      real(dp), allocatable :: b_d1(:), b_d2(:), b_d5(:)
      integer(i4) :: ib
      allocate(b_el(nbcs), b_fc(nbcs), b_ty(nbcs), b_d1(nbcs), b_d2(nbcs), b_d5(nbcs))
      do ib = 1, nbcs
         if (v2) then
            read(uin, iostat=ios) bc_elem, bc_face, bc_data, btype8; call chk_read(ios)
            b_el(ib) = int(bc_elem); b_fc(ib) = int(bc_face); b_ty(ib) = btype8
            b_d1(ib) = bc_data(1); b_d2(ib) = bc_data(2); b_d5(ib) = bc_data(5)
         else
            read(uin, iostat=ios) ei4, fi4, sbc_data, btype4; call chk_read(ios)
            b_el(ib) = ei4; b_fc(ib) = fi4; b_ty(ib) = btype4
            b_d1(ib) = real(sbc_data(1),dp); b_d2(ib) = real(sbc_data(2),dp)
            b_d5(ib) = real(sbc_data(5),dp)
         end if
         ! bounds (Neko neko_error's here; we stop cleanly rather than index OOB)
         if (b_el(ib) < 1 .or. b_el(ib) > nelv) then
            write(error_unit,'(a,i0)') &
                 'Error: BC record references element out of range: ', b_el(ib); stop 1
         end if
         if (b_fc(ib) < 1 .or. b_fc(ib) > 6) then
            write(error_unit,'(a,i0)') &
                 'Error: BC record face out of [1,6]: ', b_fc(ib); stop 1
         end if
      end do

      ! pass 1: MSH/EXO user labels (track the max user label)
      do ib = 1, nbcs
         bt = adjustl(b_ty(ib))
         if (is_msh(bt)) user_off = max(user_off, int(b_d5(ib)))
      end do
      do ib = 1, nbcs
         bt = adjustl(b_ty(ib))
         if (is_msh(bt)) call add_zone(z_e, z_f, z_lbl, nzl, b_el(ib), &
              FACET_MAP(b_fc(ib)), int(b_d5(ib)))
      end do

      ! pass 2: named BCs -> internal labels, and collect periodic facets
      do ib = 1, nbcs
         bt = adjustl(b_ty(ib)); el = b_el(ib); sym_facet = FACET_MAP(b_fc(ib))
         select case (trim(bt))
         case ('MSH','msh','EXO','exo')
         case ('W')
            call named_label(named_map, 1, cur_internal, user_off, lbl)
            call add_zone(z_e, z_f, z_lbl, nzl, el, sym_facet, lbl)
         case ('v','V')
            call named_label(named_map, 2, cur_internal, user_off, lbl)
            call add_zone(z_e, z_f, z_lbl, nzl, el, sym_facet, lbl)
         case ('O','o')
            call named_label(named_map, 3, cur_internal, user_off, lbl)
            call add_zone(z_e, z_f, z_lbl, nzl, el, sym_facet, lbl)
         case ('SYM','sym')
            call named_label(named_map, 4, cur_internal, user_off, lbl)
            call add_zone(z_e, z_f, z_lbl, nzl, el, sym_facet, lbl)
         case ('ON','on')
            call named_label(named_map, 5, cur_internal, user_off, lbl)
            call add_zone(z_e, z_f, z_lbl, nzl, el, sym_facet, lbl)
         case ('s','sl','sh','shl','S','SL','SH','SHL')
            call named_label(named_map, 6, cur_internal, user_off, lbl)
            call add_zone(z_e, z_f, z_lbl, nzl, el, sym_facet, lbl)
         case ('P')
            if (int(b_d1(ib)) < 1 .or. int(b_d1(ib)) > nelv .or. &
                 int(b_d2(ib)) < 1 .or. int(b_d2(ib)) > 6) then
               write(error_unit,'(a)') &
                    'Error: periodic BC references out-of-range partner element/face'; stop 1
            end if
            npf = npf + 1
            p_e(npf)  = el
            p_f(npf)  = sym_facet
            p_pe(npf) = int(b_d1(ib))
            p_pf(npf) = FACET_MAP(int(b_d2(ib)))
         case ('E','e')
            ! internal element connectivity: skip
         case default
            ! blank / unknown: skip
         end select
      end do
      deallocate(b_el, b_fc, b_ty, b_d1, b_d2, b_d5)
    end block

    ! ---- periodic merge: EXACT replica of Neko re2_file_read_bcs -- 3 sweeps
    !      (do j=1,3) over the P facets, each doing id = min(id_i, id_j) in place.
    !      Matching Neko's fixed 3 sweeps (not a full union-find) makes glb_pt_ids
    !      byte-identical for any mesh, incl. pathological >3-hop merge chains. ----
    do pass = 1, 3
       do i = 1, npf
          call merge_periodic_facet(pid, elem_vid, cx, cy, cz, &
               p_e(i), p_f(i), p_pe(i), p_pf(i))
       end do
    end do

    ! ---- write zones: periodic first (type 5), then labeled (type 7) ----
    write(*,'(a)')             '  [3/3] writing zones and curved-element records ...'
    nz = npf + nzl
    write(uout) nz
    do i = 1, npf
       do j = 1, 4
          gpid(j) = pid(elem_vid(FACE_RE2(j, p_f(i)), p_e(i)))
       end do
       ! nmsh_zone_t: e, f, p_e, p_f, glb_pt_ids(4), type
       write(uout) p_e(i), p_f(i), p_pe(i), p_pf(i), &
            gpid(1), gpid(2), gpid(3), gpid(4), 5_i4
    end do
    do lbl = 1, MAXZ
       do k = 1, nzl
          if (z_lbl(k) == lbl) then
             el = z_e(k); iface = z_f(k)
             write(uout) el, iface, 0_i4, lbl, 0_i4, 0_i4, 0_i4, 0_i4, 7_i4
          end if
       end do
    end do
    ! ---- curves: aggregate re2 curve records per element (edge -> slot in the
    !      nmsh curve_data(5,12)/type(12) buffers) and write one nmsh curve
    !      record per curved element, in ascending element order (as Neko's
    !      writer does). A counting sort groups the records by element in O(n).
    if (ncurv > 0 .and. .not. curve_skip) then
       allocate(ccnt(nelv)); ccnt = 0
       do k = 1, ncurv
          ccnt(rc_el(k)) = ccnt(rc_el(k)) + 1
       end do
       ncurves = 0
       do i = 1, nelv
          if (ccnt(i) > 0) ncurves = ncurves + 1
       end do
       acc = 1                               ! prefix-sum counts -> start offsets
       do i = 1, nelv
          cc = ccnt(i); ccnt(i) = acc; acc = acc + cc
       end do
       allocate(order(ncurv))                ! stable bucket by element
       do k = 1, ncurv
          order(ccnt(rc_el(k))) = k
          ccnt(rc_el(k)) = ccnt(rc_el(k)) + 1
       end do
       write(uout) ncurves
       cur_e = -1
       do idx = 1, ncurv                     ! walk records grouped by element
          k = order(idx)
          if (rc_el(k) /= cur_e) then
             if (cur_e > 0) write(uout) cur_e, cbuf, tbuf
             cbuf = 0.0_dp; tbuf = 0; cur_e = rc_el(k)
          end if
          cbuf(1:5, rc_ed(k)) = rc_d(1:5, k)
          tbuf(rc_ed(k)) = rc_ty(k)
       end do
       if (cur_e > 0) write(uout) cur_e, cbuf, tbuf
       deallocate(ccnt, order)
    else
       ncurves = 0
       write(uout) 0_i4                       ! ncurves
       if (curve_skip) write(*,'(a)') &
            '        note: unsupported curve type (s/e/other); mesh treated as '// &
            'non-curved (ncurves=0), as Neko also does'
    end if
    if (allocated(rc_el)) deallocate(rc_el, rc_ed, rc_ty, rc_d)

    write(*,'(a,i0,a,i0,a,i0,a)') '        ', npf, ' periodic + ', nzl, &
         ' labeled boundary facets, ', ncurves, ' curved elements'
    deallocate(z_e, z_f, z_lbl, p_e, p_f, p_pe, p_pf, pid)
  end block

  close(uin); close(uout)
  call ch%free()
  deallocate(elem_vid, cx, cy, cz)
  write(*,'(a)') ''
  write(*,'(a,a)') '  done -> ', trim(fout)
  write(*,'(a)') ''

contains

  !> Abort cleanly on a truncated/corrupt .re2, removing any partial output.
  subroutine chk_read(ios_)
    integer, intent(in) :: ios_
    logical :: op
    if (ios_ == 0) return
    write(error_unit,*) 'Error: truncated or corrupt .re2 file (unexpected ', &
         'end of data)'
    inquire(unit=uout, opened=op)
    if (op) close(uout, status='delete')
    stop 1
  end subroutine chk_read


  logical function is_msh(bt)
    character(len=*), intent(in) :: bt
    is_msh = (trim(bt)=='MSH' .or. trim(bt)=='msh' .or. &
              trim(bt)=='EXO' .or. trim(bt)=='exo')
  end function is_msh

  subroutine grow_coords(cx, cy, cz, cap)
    real(dp), allocatable, intent(inout) :: cx(:), cy(:), cz(:)
    integer(i4), intent(inout) :: cap
    real(dp), allocatable :: t(:)
    integer(i4) :: newcap
    newcap = cap * 2
    allocate(t(newcap)); t(1:cap) = cx; call move_alloc(t, cx)
    allocate(t(newcap)); t(1:cap) = cy; call move_alloc(t, cy)
    allocate(t(newcap)); t(1:cap) = cz; call move_alloc(t, cz)
    cap = newcap
  end subroutine grow_coords

  subroutine add_zone(ze, zf, zl, n, e, f, lbl)
    integer(i4), intent(inout) :: ze(:), zf(:), zl(:), n
    integer(i4), intent(in) :: e, f, lbl
    ! Neko aborts on labels outside [1, NEKO_MSH_MAX_ZLBLS=20]; match that rather
    ! than silently dropping them (the writer loop only emits labels 1..MAXZ).
    if (lbl < 1 .or. lbl > MAXZ) then
       write(error_unit,'(a,i0)') 'Error: boundary label out of [1,20]: ', lbl; stop 1
    end if
    n = n + 1; ze(n) = e; zf(n) = f; zl(n) = lbl
  end subroutine add_zone

  subroutine named_label(nm, which, cur, off, lbl)
    integer(i4), intent(inout) :: nm(6), cur
    integer(i4), intent(in) :: which, off
    integer(i4), intent(out) :: lbl
    if (nm(which) == 0) then
       nm(which) = cur; cur = cur + 1
    end if
    lbl = off + nm(which)
  end subroutine named_label

  !> Merge the ids of periodically-matching corners of facet (f,e) and (pf,pe),
  !! a line-for-line replica of Neko mesh_create_periodic_ids: translation
  !! L = mean facet offset; for each corner, the single partner matching within
  !! PTOL gets id = min(id_i, id_j) set in place on both. Neko aborts if a corner
  !! has 0 or >1 matches (a malformed periodic pairing); we do the same.
  subroutine merge_periodic_facet(pid, evid, cx, cy, cz, e, f, pe, pf)
    integer(i4), intent(inout) :: pid(:)
    integer(i4), intent(in) :: evid(:,:), e, f, pe, pf
    real(dp), intent(in) :: cx(:), cy(:), cz(:)
    real(dp) :: L(3), a(3), b(3)
    integer(i4) :: ia, ja, id_i, id_j, m, match
    L = 0.0_dp
    do ia = 1, 4
       id_i = evid(FACE_RE2(ia, f),  e)
       id_j = evid(FACE_RE2(ia, pf), pe)
       L(1) = L(1) + (cx(id_i) - cx(id_j))
       L(2) = L(2) + (cy(id_i) - cy(id_j))
       L(3) = L(3) + (cz(id_i) - cz(id_j))
    end do
    L = L / 4.0_dp
    do ia = 1, 4
       id_i = evid(FACE_RE2(ia, f), e)
       a = (/ cx(id_i), cy(id_i), cz(id_i) /)
       match = 0
       do ja = 1, 4
          id_j = evid(FACE_RE2(ja, pf), pe)
          b = (/ cx(id_j), cy(id_j), cz(id_j) /)
          if (norm2(a - b - L) < PTOL) then
             m = min(pid(id_i), pid(id_j))
             pid(id_i) = m; pid(id_j) = m
             match = match + 1
          end if
       end do
       if (match /= 1) then
          write(error_unit,'(a,i0,a)') 'Error: periodic facet corner has ', &
               match, ' matches (expected 1); malformed periodic pairing'; stop 1
       end if
    end do
  end subroutine merge_periodic_facet

end program rea2nbin_light
