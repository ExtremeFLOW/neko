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
!> genmeshbox_light -- a lightweight, CPU-only box-mesh generator for Neko.
!!
!! Neko's `contrib/genmeshbox` builds a full `mesh_t` using Neko's
!! general-purpose, production hash tables. Those are robust and reused across the
!! code, but they are sized generously and store polymorphic keys, so they need a
!! lot of memory for a large mesh. This tool is a lightweight, CPU-only
!! reimplementation of just the piece needed: it writes the `.nmsh` for a box mesh
!! by streaming the records directly, using plain typed arrays and no Neko
!! initialisation. Every point id and vertex is a closed-form function of the grid
!! indices, so it needs only O(1) memory beyond the three 1-D coordinate arrays,
!! giving the same mesh in a small fraction of the memory.
!!
!! Because it does not link Neko, it needs only a Fortran compiler (no MPI, no
!! Neko library) and runs anywhere -- in particular you do not need a separate
!! CPU build of Neko. The `contrib` tools are CPU-only, so on a Neko configured
!! for GPUs you would otherwise have to rebuild it for CPU just to generate a box;
!! this sidesteps that entirely.
!!
!! Output is byte-identical to genmeshbox except in the labeled-zone `p_e` and
!! `glb_pt_ids` fields, which genmeshbox leaves unset; this tool writes
!! deterministic zeros there instead.
!!
!! Usage matches genmeshbox:
!!   genmeshbox_light x0 x1 y0 y1 z0 z1 nelx nely nelz \
!!       [periodic_x periodic_y periodic_z] [dist_x dist_y dist_z] [out.nmsh]
!! periodic_* are .true./.false.; dist_* are 'uniform' or a CSV file of the nel+1
!! grid coordinates, read exactly as Neko does: a single line of comma/space
!! separated values, OR a header line followed by the values (when the file has
!! more than one line, the first line is treated as a header). Default out =
!! box.nmsh.
!!
!! Tested byte-for-byte against every genmeshbox box shipped with Neko and a
!! freshly-built current genmeshbox, but it remains your responsibility to confirm
!! the result is correct for your own box -- run `make check` on a case you trust
!! first. See genmeshbox_light.md for the memory model, validation, and scope.
program genmeshbox_light
  use, intrinsic :: iso_fortran_env, only: dp => real64, i4 => int32, i8 => int64, &
       error_unit
  implicit none

  ! Neko facet (1..6) corners as nmsh vertex slots (1..8); the table validated
  ! byte-exact for periodic ids in rea2nbin_light.
  integer, parameter :: FCORN(4,6) = reshape((/ &
       1,5,8,4, 2,6,7,3, 1,2,6,5, 4,3,7,8, 1,2,3,4, 5,6,7,8 /), (/4,6/))
  ! nmsh vertex slot (1..8) -> (ix,iy,iz) corner of the reference cube
  integer, parameter :: SIJK(3,8) = reshape((/ &
       0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 /), (/3,8/))

  character(len=1024) :: arg, fout
  real(dp) :: x0, x1, y0, y1, z0, z1
  integer :: nx, ny, nz, argc, uout, ex, ey, ez, s, lbl, na, i
  integer(i8) :: nel
  logical :: px, py, pz, direct
  character(len=1024) :: dxf, dyf, dzf
  character(len=1024), allocatable :: av(:) ! positional args (flags removed)
  real(dp), allocatable :: gx(:), gy(:), gz(:) ! grid-line coords (0..nx etc.)
  integer(i4) :: vid, elid
  real(dp) :: xx, yy, zz
  integer(i8) :: nzones, nper, nlab

  ! collect args, pulling out the --direct flag (may appear anywhere)
  argc = command_argument_count()
  allocate(av(max(1, argc))); na = 0; direct = .false.
  do i = 1, argc
     call get_command_argument(i, arg)
     if (trim(arg) == '--direct') then
        direct = .true.
     else
        na = na + 1; av(na) = arg
     end if
  end do
  if (na < 9) then
     write(*,*) 'Usage: genmeshbox_light x0 x1 y0 y1 z0 z1 nelx nely nelz \'
     write(*,*) '   [px py pz (.true./.false.)] [dist_x dist_y dist_z] [out.nmsh] [--direct]'
     write(*,*) 'px/py/pz: periodic in that direction. dist_*: "uniform" or a text file'
     write(*,*) 'of nel+1 grid coordinates. --direct: use x0+k*(x1-x0)/n coordinates'
     write(*,*) '(pre-2024 genmeshbox / exact grid values) instead of the default'
     write(*,*) 'cumulative-sum coordinates the current genmeshbox produces.'
     stop 1
  end if
  ! valid positional counts: 9 (basic), 12 (+periodic flags), 15 (+dist files),
  ! each optionally followed by a trailing output name (10, 13, 16)
  if (.not. (na == 9 .or. na == 10 .or. na == 12 .or. na == 13 .or. &
       na == 15 .or. na == 16)) then
     write(error_unit,*) 'Error: invalid number of positional arguments (', na, ')'
     write(error_unit,*) 'Valid: 9, 12 or 15 positional arguments, optionally'
     write(error_unit,*) 'followed by an output filename (10, 13, 16).'
     stop 1
  end if
  call parse_dp(av(1), x0); call parse_dp(av(2), x1); call parse_dp(av(3), y0)
  call parse_dp(av(4), y1); call parse_dp(av(5), z0); call parse_dp(av(6), z1)
  call parse_i4(av(7), nx); call parse_i4(av(8), ny); call parse_i4(av(9), nz)
  px = .false.; py = .false.; pz = .false.
  dxf = 'uniform'; dyf = 'uniform'; dzf = 'uniform'; fout = 'box.nmsh'
  if (na >= 12) then
     call parse_lg(av(10), px); call parse_lg(av(11), py); call parse_lg(av(12), pz)
  end if
  if (na >= 15) then
     dxf = av(13); dyf = av(14); dzf = av(15)
  end if
  if (na == 10 .or. na == 13 .or. na == 16) fout = av(na) ! trailing output name

  if (nx < 1 .or. ny < 1 .or. nz < 1) then
     write(error_unit,*) 'Error: nelx/nely/nelz must all be >= 1 (got ', &
          nx, ny, nz, ')'
     stop 1
  end if
  if (int(nx, i8) * int(ny, i8) * int(nz, i8) > int(huge(1_i4), i8) .or. &
      (int(nx, i8) + 1_i8) * (int(ny, i8) + 1_i8) * (int(nz, i8) + 1_i8) > &
      int(huge(1_i4), i8)) then
     write(error_unit,*) 'Error: mesh too large -- element/point counts must ', &
          'fit 32-bit .nmsh ids'
     stop 1
  end if

  write(*,'(a)') ''
  write(*,'(a)') '     _  __  ____  __ __  ____'
  write(*,'(a)') '    / |/ / / __/ / //_/ / __ \'
  write(*,'(a)') '   /    / / _/  / ,<   / /_/ /'
  write(*,'(a)') '  /_/|_/ /___/ /_/|_|  \____/   genmeshbox_light'
  write(*,'(a)') ''

  nel = int(nx, i8) * int(ny, i8) * int(nz, i8)

  ! grid-line coordinates in each direction (indices 0..n -> array 1..n+1)
  call build_axis(nx, x0, x1, dxf, direct, gx)
  call build_axis(ny, y0, y1, dyf, direct, gy)
  call build_axis(nz, z0, z1, dzf, direct, gz)

  open(newunit=uout, file=trim(fout), access='stream', form='unformatted', &
       status='replace', action='write', iostat=i)
  if (i /= 0) then
     write(error_unit,*) 'Error: cannot open output file ', trim(fout); stop 1
  end if
  write(uout) int(nel, i4), 3_i4 ! header: nelv, gdim

  ! ---- elements: genmeshbox loop order (ez outer, ey, ex inner) ----
  write(*,'(a,i0,a)') '  [1/2] writing ', nel, ' elements ...'
  elid = 0
  do ez = 0, nz-1
     do ey = 0, ny-1
        do ex = 0, nx-1
           elid = elid + 1
           write(uout) elid
           do s = 1, 8
              vid = pt_id(ex + SIJK(1,s), ey + SIJK(2,s), ez + SIJK(3,s))
              xx = gx(ex + SIJK(1,s) + 1)
              yy = gy(ey + SIJK(2,s) + 1)
              zz = gz(ez + SIJK(3,s) + 1)
              write(uout) vid, xx, yy, zz
           end do
        end do
     end do
  end do

  ! ---- zones ----
  write(*,'(a)') '  [2/2] writing zones ...'
  nper = 0; nlab = 0
  if (px) nper = nper + 2_i8 * ny * nz
  if (py) nper = nper + 2_i8 * nx * nz
  if (pz) nper = nper + 2_i8 * nx * ny
  if (.not. px) nlab = nlab + 2_i8 * ny * nz
  if (.not. py) nlab = nlab + 2_i8 * nx * nz
  if (.not. pz) nlab = nlab + 2_i8 * nx * ny
  nzones = nper + nlab
  write(uout) int(nzones, i4)

  ! periodic zones first, in genmeshbox marking order (x faces 1,2; y 3,4; z 5,6)
  if (px) then
     call emit_dir_periodic(uout, 1) ! x-low  (facet 1)
     call emit_dir_periodic(uout, 2) ! x-high (facet 2)
  end if
  if (py) then
     call emit_dir_periodic(uout, 3)
     call emit_dir_periodic(uout, 4)
  end if
  if (pz) then
     call emit_dir_periodic(uout, 5)
     call emit_dir_periodic(uout, 6)
  end if

  ! labeled zones, by label ascending (1..6); a face is labeled iff its
  ! direction is not periodic.
  do lbl = 1, 6
     if (lbl <= 2 .and. px) cycle
     if (lbl >= 3 .and. lbl <= 4 .and. py) cycle
     if (lbl >= 5 .and. pz) cycle
     call emit_dir_labeled(uout, lbl)
  end do

  write(uout) 0_i4 ! ncurves
  close(uout)

  write(*,'(a,i0,a,i0,a,i0,a,i0)') 'genmeshbox_light: ', nx, 'x', ny, 'x', nz, &
       ' = ', int(nel, i4)
  write(*,'(a,i0,a,i0,a)') 'Wrote ', nper, ' periodic + ', nlab, ' labeled zone facets'
  write(*,'(a,a)') 'Done. Wrote ', trim(fout)

contains

  !> Parse helpers: clean one-line diagnostics instead of runtime backtraces.
  subroutine parse_dp(s, v)
    character(len=*), intent(in) :: s
    real(dp), intent(out) :: v
    integer :: ios_
    read(s, *, iostat=ios_) v
    if (ios_ /= 0) then
       write(error_unit,*) 'Error: cannot parse "', trim(s), '" as a real number'
       stop 1
    end if
  end subroutine parse_dp

  subroutine parse_i4(s, v)
    character(len=*), intent(in) :: s
    integer, intent(out) :: v
    integer :: ios_
    read(s, *, iostat=ios_) v
    if (ios_ /= 0) then
       write(error_unit,*) 'Error: cannot parse "', trim(s), '" as an integer'
       stop 1
    end if
  end subroutine parse_i4

  subroutine parse_lg(s, v)
    character(len=*), intent(in) :: s
    logical, intent(out) :: v
    integer :: ios_
    read(s, *, iostat=ios_) v
    if (ios_ /= 0) then
       write(error_unit,*) 'Error: cannot parse "', trim(s), &
            '" as a logical (.true./.false.)'
       stop 1
    end if
  end subroutine parse_lg


  !> Raw lexicographic point id of grid node (px,py,pz), px in [0,nx] etc.
  pure integer(i4) function pt_id(a, b, c)
    integer, intent(in) :: a, b, c
    pt_id = int(1_i8 + a + int(b, i8)*(nx+1) + int(c, i8)*(nx+1)*(ny+1), i4)
  end function pt_id

  !> Periodic-merged id: wrap the high boundary of each periodic direction to 0
  !! (min-id merge), then take the lexicographic id -- matches create_periodic_ids.
  pure integer(i4) function merged_id(a, b, c)
    integer, intent(in) :: a, b, c
    integer :: ra, rb, rc
    ra = a; rb = b; rc = c
    if (px .and. a == nx) ra = 0
    if (py .and. b == ny) rb = 0
    if (pz .and. c == nz) rc = 0
    merged_id = pt_id(ra, rb, rc)
  end function merged_id

  !> Global element id from grid element indices (ex,ey,ez), 0-based.
  pure integer(i4) function el_id(a, b, c)
    integer, intent(in) :: a, b, c
    el_id = int(1_i8 + a + int(b, i8)*nx + int(c, i8)*nx*ny, i4)
  end function el_id

  !> Emit periodic zone records for face @a face (1..6), in genmeshbox marking order.
  subroutine emit_dir_periodic(u, face)
    integer, intent(in) :: u, face
    integer :: a, b, e_ex, e_ey, e_ez, pf, p_ex, p_ey, p_ez
    pf = opp_facet(face)
    ! marking order per direction (see genmeshbox.f90 zone loops)
    select case (face)
    case (1, 2) ! x: do e_y; do e_z
       do a = 0, ny-1
          do b = 0, nz-1
             e_ey = a; e_ez = b
             e_ex = merge(0, nx-1, face == 1)
             p_ex = merge(nx-1, 0, face == 1)
             p_ey = e_ey; p_ez = e_ez
             call wrec(u, face, e_ex, e_ey, e_ez, pf, p_ex, p_ey, p_ez)
          end do
       end do
    case (3, 4) ! y: do e_x; do e_z
       do a = 0, nx-1
          do b = 0, nz-1
             e_ex = a; e_ez = b
             e_ey = merge(0, ny-1, face == 3)
             p_ey = merge(ny-1, 0, face == 3)
             p_ex = e_ex; p_ez = e_ez
             call wrec(u, face, e_ex, e_ey, e_ez, pf, p_ex, p_ey, p_ez)
          end do
       end do
    case (5, 6) ! z: do e_x; do e_y
       do a = 0, nx-1
          do b = 0, ny-1
             e_ex = a; e_ey = b
             e_ez = merge(0, nz-1, face == 5)
             p_ez = merge(nz-1, 0, face == 5)
             p_ex = e_ex; p_ey = e_ey
             call wrec(u, face, e_ex, e_ey, e_ez, pf, p_ex, p_ey, p_ez)
          end do
       end do
    end select
  end subroutine emit_dir_periodic

  !> Write one periodic zone record (merged glb_pt_ids of the facet corners).
  subroutine wrec(u, face, e_ex, e_ey, e_ez, pf, p_ex, p_ey, p_ez)
    integer, intent(in) :: u, face, e_ex, e_ey, e_ez, pf, p_ex, p_ey, p_ez
    integer(i4) :: gp(4)
    integer :: k, sl
    do k = 1, 4
       sl = FCORN(k, face)
       gp(k) = merged_id(e_ex + SIJK(1,sl), e_ey + SIJK(2,sl), e_ez + SIJK(3,sl))
    end do
    ! nmsh_zone_t: e, f, p_e, p_f, glb_pt_ids(4), type=5
    write(u) el_id(e_ex, e_ey, e_ez), face, el_id(p_ex, p_ey, p_ez), pf, &
         gp(1), gp(2), gp(3), gp(4), 5_i4
  end subroutine wrec

  !> Emit labeled zone records for face @a face (label = face), marking order.
  !! p_e/glb_pt_ids are written 0 (Neko leaves them uninitialised).
  subroutine emit_dir_labeled(u, face)
    integer, intent(in) :: u, face
    integer :: a, b, e_ex, e_ey, e_ez
    select case (face)
    case (1, 2)
       do a = 0, ny-1
          do b = 0, nz-1
             e_ex = merge(0, nx-1, face == 1)
             write(u) el_id(e_ex, a, b), face, 0_i4, face, 0_i4, 0_i4, 0_i4, 0_i4, 7_i4
          end do
       end do
    case (3, 4)
       do a = 0, nx-1
          do b = 0, nz-1
             e_ey = merge(0, ny-1, face == 3)
             write(u) el_id(a, e_ey, b), face, 0_i4, face, 0_i4, 0_i4, 0_i4, 0_i4, 7_i4
          end do
       end do
    case (5, 6)
       do a = 0, nx-1
          do b = 0, ny-1
             e_ez = merge(0, nz-1, face == 5)
             write(u) el_id(a, b, e_ez), face, 0_i4, face, 0_i4, 0_i4, 0_i4, 0_i4, 7_i4
          end do
       end do
    end select
  end subroutine emit_dir_labeled

  pure integer function opp_facet(f)
    integer, intent(in) :: f
    select case (f)
    case (1); opp_facet = 2
    case (2); opp_facet = 1
    case (3); opp_facet = 4
    case (4); opp_facet = 3
    case (5); opp_facet = 6
    case (6); opp_facet = 5
    end select
  end function opp_facet

  !> Grid-line coordinates in one direction (n+1 values, indices 0..n).
  !!
  !! Default (accumulate) reproduces the CURRENT genmeshbox EXACTLY: element
  !! lengths are formed, then coordinates built by cumulative summation from a0
  !! (`cumm(1)=a0; cumm(k+1)=cumm(k)+el_len(k)`); for a distribution file the
  !! element lengths are the successive differences of its n+1 values.
  !!
  !! `--direct` uses the closed form `a0 + k*(a1-a0)/n` (what genmeshbox produced
  !! before the 2024 cumulative-sum change, and what every shipped example box
  !! stores), and takes distribution-file values verbatim as the grid lines.
  subroutine build_axis(n, a0, a1, fname, direct, g)
    integer, intent(in) :: n
    real(dp), intent(in) :: a0, a1
    character(len=*), intent(in) :: fname
    logical, intent(in) :: direct
    real(dp), allocatable, intent(out) :: g(:)
    real(dp), allocatable :: ellen(:), fv(:)
    integer :: k, u, ios, nlines, ios2
    character(len=256) :: hdr
    allocate(g(n+1))
    if (trim(fname) /= 'uniform') then
       allocate(fv(n+1))
       open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
       if (ios /= 0) then
          write(error_unit,*) 'Error: cannot open distribution file ', trim(fname); stop 1
       end if
       ! Match Neko's csv reader (csv_file_read_vector): count the records, and if
       ! the file has more than one line assume the first is a header and skip it,
       ! then list-directed read the n+1 grid values (commas/spaces/newlines all
       ! act as separators). A single-line CSV carries no header. This is exactly
       ! what genmeshbox does, so a distribution file that works with genmeshbox
       ! works here and produces the same grid.
       nlines = 0
       do
          read(u, *, iostat=ios2)
          if (ios2 /= 0) exit
          nlines = nlines + 1
       end do
       rewind(u)
       if (nlines > 1) read(u, '(A)', iostat=ios2) hdr
       read(u, *, iostat=ios) fv(1:n+1)
       if (ios /= 0) then
          write(error_unit,*) 'Error: expected ', n+1, ' grid values in ', trim(fname), &
               ' (Neko treats the first line as a header when the file has >1 line)'
          stop 1
       end if
       close(u)
    end if
    if (direct) then
       if (trim(fname) == 'uniform') then
          ! el_len formed first, then multiplied by k (the rounding order that
          ! reproduces the shipped example boxes exactly).
          do k = 0, n
             g(k+1) = a0 + ((a1 - a0) / real(n, dp)) * real(k, dp)
          end do
       else
          g(:) = fv(:) ! exact grid values, verbatim
       end if
    else
       allocate(ellen(n))
       if (trim(fname) == 'uniform') then
          ellen(:) = (a1 - a0) / real(n, dp)
       else
          do k = 1, n
             ellen(k) = fv(k+1) - fv(k)
          end do
       end if
       g(1) = a0
       do k = 1, n
          g(k+1) = g(k) + ellen(k)
       end do
       deallocate(ellen)
    end if
    if (allocated(fv)) deallocate(fv)
  end subroutine build_axis

end program genmeshbox_light
