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
!> Build periodic facet metadata between labeled zone pairs in an existing
!! `.nmsh` file.
!!
!! The element topology and any curved-element data are preserved. Only the
!! zone block is rebuilt so that requested labeled zone pairs become periodic
!! while untouched labeled zones remain unchanged.
program create_periodic_zones
  use neko
  use hex, only : hex_t
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS
  use utils, only : neko_error
  implicit none

  type :: facet_info_t
     integer :: f = 0
     integer :: e_local = 0
     integer :: e_global = 0
     real(kind=rp) :: center(3) = 0.0_rp
     ! Keep the facet corner coordinates so that a translated candidate can be
     ! validated geometrically and not only by its centroid.
     real(kind=rp) :: xyz(3, 4) = 0.0_rp
     integer :: ids(4) = 0
  end type facet_info_t

  character(len=NEKO_FNAME_LEN) :: input_fname, output_fname, arg
  character(len=4096) :: pair_spec
  type(file_t) :: input_file, output_file
  type(mesh_t) :: msh
  integer :: argc, i, npairs
  integer, allocatable :: zone_a(:), zone_b(:)
  integer :: n_existing_periodic
  integer, allocatable :: old_periodic_f(:), old_periodic_e(:)
  integer, allocatable :: old_partner_f(:), old_partner_e(:)
  integer, allocatable :: old_pids(:,:), old_org_ids(:,:)
  integer, allocatable :: old_label_counts(:)
  integer, allocatable :: old_labels_f(:,:), old_labels_e(:,:)
  integer :: label, j, next_arg
  integer :: ios
  real(kind=rp) :: tol
  logical :: tol_set
  character(len=512) :: err_msg

  argc = command_argument_count()
  if (argc < 3) then
     call print_usage()
     stop
  end if

  call get_command_argument(1, input_fname)
  call get_command_argument(2, output_fname)

  pair_spec = ''
  tol_set = .false.
  i = 3
  do while (i <= argc)
     call get_command_argument(i, arg)
     if (index(trim(arg), '--tol=') == 1) then
        read(arg(7:), *, iostat = ios) tol
        if (ios /= 0) then
           call neko_error('Invalid value passed to --tol=')
        end if
        tol_set = .true.
     else if (trim(arg) == '--tol') then
        next_arg = i + 1
        if (next_arg > argc) then
           call neko_error('Missing value after --tol')
        end if
        call get_command_argument(next_arg, arg)
        read(arg, *, iostat = ios) tol
        if (ios /= 0) then
           call neko_error('Invalid value passed after --tol')
        end if
        tol_set = .true.
        i = next_arg
     else if (index(trim(arg), '--') == 1) then
        call neko_error('Unknown option: ' // trim(arg))
     else
        ! Allow a pair specification split across several shell tokens.
        if (len_trim(pair_spec) > 0) pair_spec = trim(pair_spec) // ' '
        pair_spec = trim(pair_spec) // trim(arg)
     end if
     i = i + 1
  end do

  if (len_trim(pair_spec) == 0) then
     call print_usage()
     call neko_error('Missing periodic zone pair specification')
  end if

  call parse_zone_pairs(trim(pair_spec), zone_a, zone_b, npairs)

  call neko_init

  if (pe_size /= 1) then
     call neko_error('create_periodic_zones currently supports serial ' // &
          'execution only')
  end if

  call input_file%init(trim(input_fname))
  call input_file%read(msh)

  if (.not. tol_set) then
     tol = default_tolerance(msh)
  end if

  ! Snapshot the current zone tables before rebuilding them.
  call snapshot_existing_zones(msh, old_periodic_f, old_periodic_e, &
       old_partner_f, old_partner_e, old_pids, old_org_ids, old_label_counts, &
       old_labels_f, old_labels_e)

  ! Sanity checks
  do i = 1, npairs
     if (zone_a(i) < 1 .or. zone_a(i) > NEKO_MSH_MAX_ZLBLS .or. &
          zone_b(i) < 1 .or. zone_b(i) > NEKO_MSH_MAX_ZLBLS) then
        call neko_error('Zone indices must be in [1,20]')
     end if
     if (zone_a(i) == zone_b(i)) then
        call neko_error('Periodic zone pair labels must be distinct')
     end if
  end do

  ! Recreate the zone containers from scratch. This keeps the logic explicit:
  ! first restore preserved data, then append new periodic pairs.
  call msh%periodic%free()
  call msh%periodic%init(msh%nelv)
  do i = 1, NEKO_MSH_MAX_ZLBLS
     call msh%labeled_zones(i)%free()
     call msh%labeled_zones(i)%init(msh%nelv)
  end do
  msh%facet_type = 0

  ! Restore any periodic information already present in the input mesh.
  n_existing_periodic = size(old_periodic_f)
  do i = 1, n_existing_periodic
     call msh%periodic%add_periodic_facet(old_periodic_f(i), &
          old_periodic_e(i), old_partner_f(i), old_partner_e(i), &
          old_pids(:, i), old_org_ids(:, i))
  end do

  ! Copy back all labeled zones that are not being converted.
  do label = 1, NEKO_MSH_MAX_ZLBLS
     if (zone_is_converted(label, zone_a, zone_b)) cycle
     do i = 1, old_label_counts(label)
        call msh%mark_labeled_facet(old_labels_f(i, label), &
             old_labels_e(i, label), label)
     end do
  end do

  ! Convert each requested pair after the preserved data has been restored.
  do i = 1, npairs
     call build_periodic_pair(msh, zone_a(i), zone_b(i), old_label_counts, &
          old_labels_f, old_labels_e, tol)
  end do

  call msh%periodic%finalize()
  do i = 1, NEKO_MSH_MAX_ZLBLS
     call msh%labeled_zones(i)%finalize()
  end do

  if (pe_rank == 0) then
     write(*,'(A,1X,ES12.4)') 'Matching tolerance:', tol
  end if

  call output_file%init(trim(output_fname))
  call output_file%write(msh)

  deallocate(zone_a, zone_b)
  deallocate(old_periodic_f, old_periodic_e)
  deallocate(old_partner_f, old_partner_e)
  deallocate(old_pids, old_org_ids)
  deallocate(old_label_counts)
  deallocate(old_labels_f, old_labels_e)

  call msh%free()
  call neko_finalize

contains

  !> Print the command-line usage summary.
  subroutine print_usage()
    write(*,*) 'Usage: ./create_periodic_zones input.nmsh output.nmsh ' // &
         '"(1,2),(3,4)" [--tol=value]'
    write(*,*) 'The pair specification accepts any delimiter mix, e.g. ' // &
         '"(1,2),(3,4)" or "1:2 3:4".'
  end subroutine print_usage

  !> Parse the user-supplied periodic zone pair specification.
  !! @param spec Raw pair specification from the command line
  !! @param zone_a First zone index in each parsed pair
  !! @param zone_b Second zone index in each parsed pair
  !! @param npairs Number of parsed pairs
  subroutine parse_zone_pairs(spec, zone_a, zone_b, npairs)
    character(len=*), intent(in) :: spec
    integer, allocatable, intent(out) :: zone_a(:), zone_b(:)
    integer, intent(out) :: npairs
    character(len=len(spec)) :: cleaned
    integer, allocatable :: values(:)
    integer :: i, nvalues, ios
    logical :: in_token

    cleaned = spec
    ! Replace all separators by spaces so that list-directed parsing can handle
    ! formats such as "(1,2),(3,4)" or "1:2 3:4".
    do i = 1, len_trim(spec)
       if (.not. is_integer_char(spec(i:i))) cleaned(i:i) = ' '
    end do

    nvalues = 0
    in_token = .false.
    do i = 1, len_trim(cleaned)
       if (cleaned(i:i) /= ' ') then
          if (.not. in_token) then
             nvalues = nvalues + 1
             in_token = .true.
          end if
       else
          in_token = .false.
       end if
    end do

    if (nvalues == 0 .or. mod(nvalues, 2) /= 0) then
       call neko_error('Could not parse periodic zone pairs')
    end if

    allocate(values(nvalues))
    read(cleaned, *, iostat = ios) values
    if (ios /= 0) then
       call neko_error('Could not parse periodic zone pairs')
    end if

    npairs = nvalues / 2
    allocate(zone_a(npairs), zone_b(npairs))
    do i = 1, npairs
       zone_a(i) = values(2 * i - 1)
       zone_b(i) = values(2 * i)
    end do

    deallocate(values)
  end subroutine parse_zone_pairs

  !> Check whether a character can belong to a zone label token.
  !! @param ch Character to test
  logical function is_integer_char(ch)
    character(len=1), intent(in) :: ch
    is_integer_char = (ch >= '0' .and. ch <= '9')
  end function is_integer_char

  !> Test if a labeled zone is requested for periodic conversion.
  !! @param label Zone label to test
  !! @param zone_a First zone index in each requested pair
  !! @param zone_b Second zone index in each requested pair
  logical function zone_is_converted(label, zone_a, zone_b)
    integer, intent(in) :: label
    integer, intent(in) :: zone_a(:), zone_b(:)
    zone_is_converted = any(zone_a == label) .or. any(zone_b == label)
  end function zone_is_converted

  !> Copy the current periodic and labeled zone tables out of the mesh.
  !! @param msh Mesh whose zone metadata is being snapshotted
  !! @param periodic_f Facet id of each existing periodic entry
  !! @param periodic_e Local element id of each existing periodic entry
  !! @param partner_f Partner facet id of each existing periodic entry
  !! @param partner_e Partner element id of each existing periodic entry
  !! @param pids Stored periodic point ids for each periodic entry
  !! @param org_ids Stored original point ids for each periodic entry
  !! @param label_counts Number of facets in each labeled zone
  !! @param labels_f Facet ids grouped by labeled zone
  !! @param labels_e Local element ids grouped by labeled zone
  subroutine snapshot_existing_zones(msh, periodic_f, periodic_e, partner_f, &
       partner_e, pids, org_ids, label_counts, labels_f, labels_e)
    type(mesh_t), intent(in) :: msh
    integer, allocatable, intent(out) :: periodic_f(:), periodic_e(:)
    integer, allocatable, intent(out) :: partner_f(:), partner_e(:)
    integer, allocatable, intent(out) :: pids(:,:), org_ids(:,:)
    integer, allocatable, intent(out) :: label_counts(:)
    integer, allocatable, intent(out) :: labels_f(:,:), labels_e(:,:)
    integer :: i, label, max_label_size

    allocate(periodic_f(msh%periodic%size), periodic_e(msh%periodic%size))
    allocate(partner_f(msh%periodic%size), partner_e(msh%periodic%size))
    allocate(pids(4, msh%periodic%size), org_ids(4, msh%periodic%size))

    do i = 1, msh%periodic%size
       periodic_f(i) = msh%periodic%facet_el(i)%x(1)
       periodic_e(i) = msh%periodic%facet_el(i)%x(2)
       partner_f(i) = msh%periodic%p_facet_el(i)%x(1)
       partner_e(i) = msh%periodic%p_facet_el(i)%x(2)
       pids(:, i) = msh%periodic%p_ids(i)%x
       org_ids(:, i) = msh%periodic%org_ids(i)%x
    end do

    allocate(label_counts(NEKO_MSH_MAX_ZLBLS))
    label_counts = 0
    max_label_size = 0
    do label = 1, NEKO_MSH_MAX_ZLBLS
       label_counts(label) = msh%labeled_zones(label)%size
       max_label_size = max(max_label_size, label_counts(label))
    end do
    ! Keep the first extent non-zero to avoid zero-sized corner cases in the
    ! simple storage used below.
    allocate(labels_f(max(1, max_label_size), NEKO_MSH_MAX_ZLBLS))
    allocate(labels_e(max(1, max_label_size), NEKO_MSH_MAX_ZLBLS))
    labels_f = 0
    labels_e = 0

    do label = 1, NEKO_MSH_MAX_ZLBLS
       do i = 1, label_counts(label)
          labels_f(i, label) = msh%labeled_zones(label)%facet_el(i)%x(1)
          labels_e(i, label) = msh%labeled_zones(label)%facet_el(i)%x(2)
       end do
    end do
  end subroutine snapshot_existing_zones

  !> Convert one labeled zone pair into periodic facets.
  !! @param msh Mesh whose zone tables are being rebuilt
  !! @param label_a First labeled zone in the requested periodic pair
  !! @param label_b Second labeled zone in the requested periodic pair
  !! @param label_counts Number of facets in each labeled zone
  !! @param labels_f Facet ids grouped by labeled zone
  !! @param labels_e Local element ids grouped by labeled zone
  !! @param tol Matching tolerance for translated facets
  !! @details The repeated calls to `msh%create_periodic_ids` follow the same
  !! fixed-pass pattern used by the `re2` importer in
  !! `src/io/re2_file.f90:re2_file_read_bcs`, where periodic point ids are
  !! propagated over the matched facets several times to make shared vertices
  !! converge to consistent ids.
  subroutine build_periodic_pair(msh, label_a, label_b, label_counts, &
       labels_f, labels_e, tol)
    type(mesh_t), intent(inout) :: msh
    integer, intent(in) :: label_a, label_b
    integer, intent(in) :: label_counts(:)
    integer, intent(in) :: labels_f(:, :), labels_e(:, :)
    real(kind=rp), intent(in) :: tol
    type(facet_info_t), allocatable :: facets_a(:), facets_b(:)
    integer, allocatable :: match_ab(:), match_ba(:)
    real(kind=rp) :: offset(3)
    integer :: i, j

    if (label_counts(label_a) /= label_counts(label_b)) then
       call neko_error('Periodic zone pair has a different number of facets')
    end if
    if (label_counts(label_a) == 0) then
       call neko_error('Periodic zone pair refers to an empty labeled zone')
    end if

    allocate(facets_a(label_counts(label_a)))
    allocate(facets_b(label_counts(label_b)))

    do i = 1, size(facets_a)
       call fill_facet_info(msh, labels_f(i, label_a), labels_e(i, label_a), &
            facets_a(i))
       call fill_facet_info(msh, labels_f(i, label_b), labels_e(i, label_b), &
            facets_b(i))
    end do

    ! Infer one translation vector for the whole pair from the centroid shift.
    ! The subsequent geometric match verifies that the same translation works
    ! for every facet.
    offset = 0.0_rp
    do i = 1, size(facets_a)
       offset = offset + facets_b(i)%center - facets_a(i)%center
    end do
    offset = offset / real(size(facets_a), kind = rp)

    allocate(match_ab(size(facets_a)), match_ba(size(facets_b)))
    match_ab = 0
    match_ba = 0

    ! Enforce a strict one-to-one correspondence between the two labeled zones.
    do i = 1, size(facets_a)
       call find_matching_facet(facets_a(i), offset, facets_b, tol, j)
       if (match_ba(j) /= 0) then
          write(err_msg, '(A,I0,A,I0,A,3(1X,ES12.4),A,I0,A,I0,A)') &
               'Periodic zone mapping is not one-to-one for facet ', &
               facets_a(i)%f, ' on element ', facets_a(i)%e_global, &
               ' with centroid', facets_a(i)%center, '. It conflicts with ', &
               'a previous match to facet ', facets_b(j)%f, &
               ' on element ', facets_b(j)%e_global, '.'
          call neko_error(trim(err_msg))
       end if
       match_ab(i) = j
       match_ba(j) = i
    end do

    ! Neko stores one zone record per periodic facet, so add both directions.
    do i = 1, size(facets_a)
       j = match_ab(i)
       call msh%periodic%add_periodic_facet(facets_a(i)%f, &
            facets_a(i)%e_global, facets_b(j)%f, facets_b(j)%e_global, &
            facets_a(i)%ids, facets_a(i)%ids)
       call msh%periodic%add_periodic_facet(facets_b(j)%f, &
            facets_b(j)%e_global, facets_a(i)%f, facets_a(i)%e_global, &
            facets_b(j)%ids, facets_b(j)%ids)
    end do

    ! Repeat the periodic id propagation a few times, mirroring the existing
    ! importer workflow used for periodic faces read from other mesh formats.
    do j = 1, 3
       do i = 1, size(facets_a)
          call msh%create_periodic_ids(facets_a(i)%f, facets_a(i)%e_local, &
               facets_b(match_ab(i))%f, facets_b(match_ab(i))%e_local)
          call msh%create_periodic_ids(facets_b(match_ab(i))%f, &
               facets_b(match_ab(i))%e_local, facets_a(i)%f, &
               facets_a(i)%e_local)
       end do
    end do

    if (pe_rank == 0) then
       write(*,'(A,I0,A,I0,A,3(1X,ES12.4))') 'Periodic zones ', &
            label_a, ' <-> ', label_b, ', offset:', offset
    end if

    deallocate(facets_a, facets_b, match_ab, match_ba)
  end subroutine build_periodic_pair

  !> Collect the geometry and ids associated with one facet.
  !! @param msh Mesh containing the facet
  !! @param facet Facet index on the element
  !! @param el_local Local element index containing the facet
  !! @param info Gathered facet metadata
  subroutine fill_facet_info(msh, facet, el_local, info)
    type(mesh_t), intent(inout) :: msh
    integer, intent(in) :: facet, el_local
    type(facet_info_t), intent(out) :: info

    info%f = facet
    info%e_local = el_local
    info%e_global = msh%elements(el_local)%e%id()
    call get_hex_facet_geometry(msh, facet, el_local, info%xyz, info%center)
    call msh%get_facet_ids(facet, el_local, info%ids)
  end subroutine fill_facet_info

  !> Extract the corner coordinates and centroid of a hex facet.
  !! @param msh Mesh containing the facet
  !! @param facet Facet index on the element
  !! @param el_local Local element index containing the facet
  !! @param xyz Corner coordinates of the facet
  !! @param center Facet centroid
  subroutine get_hex_facet_geometry(msh, facet, el_local, xyz, center)
    type(mesh_t), intent(in) :: msh
    integer, intent(in) :: facet, el_local
    real(kind=rp), intent(out) :: xyz(3, 4), center(3)
    integer :: i
    integer, parameter :: face_nodes(4, 6) = reshape([ &
         1, 5, 7, 3, &
         2, 6, 8, 4, &
         1, 2, 6, 5, &
         3, 4, 8, 7, &
         1, 2, 4, 3, &
         5, 6, 8, 7], [4, 6])

    center = 0.0_rp
    select type (el => msh%elements(el_local)%e)
    type is (hex_t)
       do i = 1, 4
          xyz(:, i) = el%pts(face_nodes(i, facet))%p%x(1:3)
          center = center + xyz(:, i)
       end do
       center = center / 4.0_rp
    class default
       call neko_error('create_periodic_zones only supports hexahedral meshes')
    end select
  end subroutine get_hex_facet_geometry

  !> Find the unique candidate facet matching a translated source facet.
  !! @param facet Source facet metadata
  !! @param offset Translation vector from the source zone to the target zone
  !! @param candidates Candidate facets from the target zone
  !! @param tol Matching tolerance
  !! @param match_idx Index of the matching candidate facet
  subroutine find_matching_facet(facet, offset, candidates, tol, match_idx)
    type(facet_info_t), intent(in) :: facet
    real(kind=rp), intent(in) :: offset(3)
    type(facet_info_t), intent(in) :: candidates(:)
    real(kind=rp), intent(in) :: tol
    integer, intent(out) :: match_idx
    integer :: i, nmatch

    match_idx = 0
    nmatch = 0
    do i = 1, size(candidates)
       ! Use the centroid difference as a cheap pre-filter, then confirm the
       ! match by checking the translated facet corners.
       if (norm2(candidates(i)%center - facet%center - offset) .le. tol .and. &
            facets_match_under_translation(facet, candidates(i), offset, &
            tol)) then
          nmatch = nmatch + 1
          match_idx = i
       end if
    end do

    if (nmatch /= 1) then
       call neko_error('Could not determine a unique periodic facet match')
    end if
  end subroutine find_matching_facet

  !> Check whether two facets coincide after applying a translation.
  !! @param lhs Source facet metadata
  !! @param rhs Candidate target facet metadata
  !! @param offset Translation vector from `lhs` to `rhs`
  !! @param tol Matching tolerance
  logical function facets_match_under_translation(lhs, rhs, offset, tol)
    type(facet_info_t), intent(in) :: lhs, rhs
    real(kind=rp), intent(in) :: offset(3)
    real(kind=rp), intent(in) :: tol
    logical :: used(4)
    integer :: i, j, nmatch

    used = .false.
    facets_match_under_translation = .true.

    do i = 1, 4
       nmatch = 0
       do j = 1, 4
          if (used(j)) cycle
          if (norm2(rhs%xyz(:, j) - lhs%xyz(:, i) - offset) .le. tol) then
             used(j) = .true.
             nmatch = nmatch + 1
          end if
       end do
       if (nmatch /= 1) then
          facets_match_under_translation = .false.
          return
       end if
    end do
  end function facets_match_under_translation

  !> Compute a default geometric tolerance from the mesh extent.
  !! @details Ad hoc in terms of the numerical values: 1e-8 * lengthscale.
  !! @param msh Mesh used to estimate a characteristic length scale
  real(kind=rp) function default_tolerance(msh)
    type(mesh_t), intent(in) :: msh
    real(kind=rp) :: xyz(3), xmin(3), xmax(3), diag
    integer :: i

    xmin = huge(1.0_rp)
    xmax = -huge(1.0_rp)
    do i = 1, msh%mpts
       xyz = msh%points(i)%x(1:3)
       xmin = min(xmin, xyz)
       xmax = max(xmax, xyz)
    end do
    ! Scale the tolerance with the full bounding-box diagonal rather than with
    ! only the largest axis extent. This is more robust for anisotropic boxes
    ! and for meshes that are not aligned with one dominant direction.
    diag = norm2(xmax - xmin)
    default_tolerance = max(1.0e-10_rp, 1.0e-8_rp * max(1.0_rp, diag))
  end function default_tolerance

end program create_periodic_zones
