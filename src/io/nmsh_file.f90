! Copyright (c) 2019-2021, The Neko Authors
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
!> Neko binary mesh data
module nmsh_file
  use generic_file, only: generic_file_t
  use comm, only: NEKO_COMM, pe_rank, pe_size
  use num_types, only: rp, dp, i4, i8
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBLS
  use utils, only: neko_error
  use point, only: point_t
  use tuple, only: tuple4_i4_t
  use nmsh, only: nmsh_hex_t, nmsh_quad_t, nmsh_zone_t, nmsh_curve_el_t
  use element, only: element_t
  use datadist, only: linear_dist_t
  use neko_mpi_types, only: MPI_NMSH_HEX, MPI_NMSH_QUAD, MPI_NMSH_ZONE, &
       MPI_NMSH_CURVE, MPI_INTEGER_SIZE
  use mpi_f08, only: MPI_Wtime, MPI_Status, MPI_File, MPI_OFFSET_KIND, &
       MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_MODE_RDONLY, MPI_INFO_NULL, &
       MPI_File_open, MPI_File_close, MPI_File_read_all, MPI_File_write_all, &
       MPI_File_write_at_all, MPI_File_read_at_all, MPI_INTEGER, MPI_SUM, &
       MPI_Exscan, MPI_Barrier, MPI_Type_size, MPI_Allreduce, MPI_File_sync
  use logger, only: neko_log, LOG_SIZE
  implicit none

  private
  !> Specifices the maximum number of elements any rank is allowed to write (for nmsh).
  !! Needed in order to generate large meshes where an individual write might exceed 2GB.
  integer, parameter :: max_write_nel = 8000000
  !> Interface for Neko nmsh files
  type, public, extends(generic_file_t) :: nmsh_file_t
   contains
     procedure :: read => nmsh_file_read
     procedure :: write => nmsh_file_write
  end type nmsh_file_t

contains

  !> Load a mesh from a binary Neko nmsh file
  subroutine nmsh_file_read(this, data)
    class(nmsh_file_t) :: this
    class(*), target, intent(inout) :: data
    type(nmsh_hex_t), allocatable :: nmsh_hex(:)
    type(nmsh_quad_t), allocatable :: nmsh_quad(:)
    type(nmsh_zone_t), allocatable :: nmsh_zone(:)
    type(nmsh_curve_el_t), allocatable :: nmsh_curve(:)
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, mpi_el_offset
    integer :: i, j, ierr, element_offset
    integer :: nmsh_quad_size, nmsh_hex_size, nmsh_zone_size
    integer :: nelv, gdim, nzones, ncurves
    integer :: el_idx, el_idx_glb
    type(point_t), target :: p(8)
    type(linear_dist_t) :: dist
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: t_start, t_end

    call this%check_exists()

    select type (data)
    type is (mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select


    call neko_log%section("Mesh")
    call neko_log%message('Reading a binary Neko file ' // this%get_fname())

    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)
    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)
    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)

    call MPI_File_open(NEKO_COMM, trim(this%get_fname()), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

    if (ierr > 0) then
       call neko_error('Could not open the mesh file ' // this%get_fname() // &
            'for reading!')
    end if
    call MPI_File_read_all(fh, nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)

    write(log_buf,1) gdim, nelv
1   format('gdim = ', i1, ', nelements = ', i9)
    call neko_log%message(log_buf)

    if (gdim .eq. 2) then
       call MPI_File_close(fh, ierr)
       call nmsh_file_read_2d(this, msh)
    else

       dist = linear_dist_t(nelv, pe_rank, pe_size, NEKO_COMM)
       nelv = dist%num_local()
       element_offset = dist%start_idx()

       call msh%init(gdim, nelv)

       call neko_log%message('Reading elements')

       if (msh%gdim .eq. 2) then
          allocate(nmsh_quad(msh%nelv))
          mpi_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
               int(element_offset, i8) * int(nmsh_quad_size, i8)
          call MPI_File_read_at_all(fh, mpi_offset, &
               nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
          do i = 1, nelv
             do j = 1, 4
                call p(j)%init(nmsh_quad(i)%v(j)%v_xyz, nmsh_quad(i)%v(j)%v_idx)
             end do
             ! swap vertices to keep symmetric vertex numbering in neko
             call msh%add_element(i, nmsh_quad(i)%el_idx, &
                  p(1), p(2), p(4), p(3))
          end do
          deallocate(nmsh_quad)
          mpi_el_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
               int(dist%num_global(), i8) * int(nmsh_quad_size, i8)
       else if (msh%gdim .eq. 3) then
          allocate(nmsh_hex(msh%nelv))
          mpi_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
               int(element_offset, i8) * int(nmsh_hex_size, i8)
          call MPI_File_read_at_all(fh, mpi_offset, &
               nmsh_hex, msh%nelv, MPI_NMSH_HEX, status, ierr)
          do i = 1, nelv
             do j = 1, 8
                call p(j)%init(nmsh_hex(i)%v(j)%v_xyz, nmsh_hex(i)%v(j)%v_idx)
             end do
             ! swap vertices to keep symmetric vertex numbering in neko
             call msh%add_element(i, nmsh_hex(i)%el_idx, &
                  p(1), p(2), p(4), p(3), p(5), p(6), p(8), p(7))
          end do
          deallocate(nmsh_hex)
          mpi_el_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
               int(dist%num_global(), i8) * int(nmsh_hex_size, i8)
       else
          if (pe_rank .eq. 0) call neko_error('Invalid dimension of mesh')
       end if
       call neko_log%message('Reading BC/zone data')

       mpi_offset = mpi_el_offset
       call MPI_File_read_at_all(fh, mpi_offset, &
            nzones, 1, MPI_INTEGER, status, ierr)
       if (nzones .gt. 0) then
          allocate(nmsh_zone(nzones))

          !>
          !!@todo Fix the parallel reading in this part, let each rank read
          !!a piece and pass the pieces around, filtering out matching zones
          !!in the local mesh.
          !!
          mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8)
          call MPI_File_read_at_all(fh, mpi_offset, &
               nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)

          do i = 1, nzones
             el_idx_glb = nmsh_zone(i)%e
             if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
                select case (nmsh_zone(i)%type)
                case (5)
                   call msh%mark_periodic_facet(nmsh_zone(i)%f, el_idx, &
                        nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, &
                        nmsh_zone(i)%glb_pt_ids)
                case (7)
                   call msh%mark_labeled_facet(nmsh_zone(i)%f, el_idx, &
                        nmsh_zone(i)%p_f)
                end select
             end if
          end do
          !Apply facets, important that marking is finished
          do i = 1, nzones
             el_idx_glb = nmsh_zone(i)%e
             if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
                select case (nmsh_zone(i)%type)
                case (5)
                   call msh%apply_periodic_facet(nmsh_zone(i)%f, el_idx, &
                        nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, &
                        nmsh_zone(i)%glb_pt_ids)
                end select
             end if
          end do

          deallocate(nmsh_zone)
       end if
       call neko_log%message('Reading deformation data')

       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8) + &
            int(nzones, i8)*int(nmsh_zone_size, i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            ncurves, 1, MPI_INTEGER, status, ierr)

       if (ncurves .gt. 0) then

          allocate(nmsh_curve(ncurves))
          mpi_offset = mpi_el_offset + int(2 * MPI_INTEGER_SIZE, i8) + &
               int(nzones, i8)*int(nmsh_zone_size, i8)
          call MPI_File_read_at_all(fh, mpi_offset, &
               nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)

          do i = 1, ncurves
             el_idx_glb = nmsh_curve(i)%e
             if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
                call msh%mark_curve_element(el_idx, &
                     nmsh_curve(i)%curve_data, nmsh_curve(i)%type)
             end if

          end do

          deallocate(nmsh_curve)
       end if

       call MPI_File_close(fh, ierr)
       call neko_log%message('Mesh read, setting up connectivity')

       t_start = MPI_WTIME()
       call msh%finalize()
       call MPI_Barrier(NEKO_COMM, ierr)
       t_end = MPI_WTIME()
       write(log_buf, '(A)') 'Done setting up mesh and connectivity'
       call neko_log%message(log_buf)
       write(log_buf, '(A,F9.6)') &
            'Mesh and connectivity setup (excluding read) time (s): ', &
            t_end - t_start
       call neko_log%message(log_buf)

       call neko_log%end_section()
    end if

  end subroutine nmsh_file_read

  !> Load a mesh from a binary Neko nmsh file
  subroutine nmsh_file_read_2d(this, msh)
    class(nmsh_file_t) :: this
    type(mesh_t), pointer, intent(inout) :: msh
    type(nmsh_hex_t), allocatable :: nmsh_hex(:)
    type(nmsh_quad_t), allocatable :: nmsh_quad(:)
    type(nmsh_zone_t), allocatable :: nmsh_zone(:)
    type(nmsh_curve_el_t), allocatable :: nmsh_curve(:)
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, mpi_el_offset
    integer :: i, j, ierr, element_offset, id
    integer :: nmsh_quad_size, nmsh_hex_size, nmsh_zone_size
    integer :: nelv, gdim, nzones, ncurves, ids(4)
    integer :: el_idx_glb, el_idx
    type(point_t) :: p(8)
    type(linear_dist_t) :: dist
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: depth = 1d0
    real(kind=dp) :: coord(3)
    type(tuple4_i4_t) :: glb_pt_ids


    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)
    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)
    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)

    call MPI_File_open(NEKO_COMM, trim(this%get_fname()), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_read_all(fh, nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)

    write(log_buf,2) gdim
2   format('gdim = ', i1, ', no full 2d support, creating thin slab')
    call neko_log%message(log_buf)
    gdim = 3

    dist = linear_dist_t(nelv, pe_rank, pe_size, NEKO_COMM)
    nelv = dist%num_local()
    element_offset = dist%start_idx()

    call msh%init(gdim, nelv)

    allocate(nmsh_quad(msh%nelv))
    mpi_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
         int(element_offset, i8) * int(nmsh_quad_size, i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
    do i = 1, nelv
       do j = 1, 4
          coord = nmsh_quad(i)%v(j)%v_xyz
          coord(3) = 0_rp
          call p(j)%init(coord, nmsh_quad(i)%v(j)%v_idx)
       end do
       do j = 1, 4
          coord = nmsh_quad(i)%v(j)%v_xyz
          coord(3) = depth
          id = nmsh_quad(i)%v(j)%v_idx+msh%glb_nelv*8
          call p(j+4)%init(coord, id)
       end do
       ! swap vertices to keep symmetric vertex numbering in neko
       call msh%add_element(i, nmsh_quad(i)%el_idx, &
            p(1), p(2), p(4), p(3), p(5), p(6), p(8), p(7))
    end do
    deallocate(nmsh_quad)
    mpi_el_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
         int(dist%num_global(), i8) * int(nmsh_quad_size, i8)

    mpi_offset = mpi_el_offset
    call MPI_File_read_at_all(fh, mpi_offset, &
         nzones, 1, MPI_INTEGER, status, ierr)
    if (nzones .gt. 0) then
       allocate(nmsh_zone(nzones))

       !>
       !!@todo Fix the parallel reading in this part, let each rank read
       !!a piece and pass the pieces around, filtering out matching zones
       !!in the local mesh.
       !!
       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)

       do i = 1, nzones
          el_idx_glb = nmsh_zone(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             select case (nmsh_zone(i)%type)
             case (5)
                nmsh_zone(i)%glb_pt_ids(3) = nmsh_zone(i)%glb_pt_ids(1) + &
                     msh%glb_nelv * 8
                nmsh_zone(i)%glb_pt_ids(4) = nmsh_zone(i)%glb_pt_ids(2) + &
                     msh%glb_nelv * 8
                if (nmsh_zone(i)%f .eq. 1 .or. nmsh_zone(i)%f .eq. 2) then
                   ids(1) = nmsh_zone(i)%glb_pt_ids(1)
                   ids(2) = nmsh_zone(i)%glb_pt_ids(3)
                   ids(3) = nmsh_zone(i)%glb_pt_ids(4)
                   ids(4) = nmsh_zone(i)%glb_pt_ids(2)
                else
                   ids(1) = nmsh_zone(i)%glb_pt_ids(1)
                   ids(2) = nmsh_zone(i)%glb_pt_ids(2)
                   ids(3) = nmsh_zone(i)%glb_pt_ids(4)
                   ids(4) = nmsh_zone(i)%glb_pt_ids(3)
                end if
                nmsh_zone(i)%glb_pt_ids = ids
                call msh%mark_periodic_facet(nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, ids)
             case (7)
                call msh%mark_labeled_facet(nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f)
             end select
          end if
       end do
       !Apply facets, important that marking is finished
       do i = 1, nzones
          el_idx_glb = nmsh_zone(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             select case (nmsh_zone(i)%type)
             case (5)
                call msh%apply_periodic_facet(nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, &
                     nmsh_zone(i)%glb_pt_ids)
             end select
          end if
       end do
       !Do the same for extruded 3d points
       do el_idx = 1, nelv
          call msh%elements(el_idx)%e%facet_order(glb_pt_ids,5)
          call msh%mark_periodic_facet(6, el_idx, &
               5, el_idx, glb_pt_ids%x)
          call msh%elements(el_idx)%e%facet_order(glb_pt_ids,5)
          call msh%mark_periodic_facet(5, el_idx, &
               6, el_idx, glb_pt_ids%x)
       end do
       do el_idx = 1, nelv
          call msh%elements(el_idx)%e%facet_order(glb_pt_ids,5)
          call msh%apply_periodic_facet(6, el_idx, &
               5, el_idx, glb_pt_ids%x)
          call msh%elements(el_idx)%e%facet_order(glb_pt_ids,5)
          call msh%apply_periodic_facet(5, el_idx, &
               6, el_idx, glb_pt_ids%x)
       end do

       deallocate(nmsh_zone)
    end if

    mpi_offset = mpi_el_offset + &
         int(MPI_INTEGER_SIZE, i8) + int(nzones, i8)*int(nmsh_zone_size, i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         ncurves, 1, MPI_INTEGER, status, ierr)

    if (ncurves .gt. 0) then

       allocate(nmsh_curve(ncurves))
       mpi_offset = mpi_el_offset + &
            int(2*MPI_INTEGER_SIZE, i8) + &
            int(nzones, i8)*int(nmsh_zone_size, i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)

       do i = 1, ncurves
          el_idx_glb = nmsh_curve(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             call msh%mark_curve_element(el_idx, &
                  nmsh_curve(i)%curve_data, nmsh_curve(i)%type)
          end if
       end do

       deallocate(nmsh_curve)
    end if

    call MPI_File_close(fh, ierr)

    call msh%finalize()

    call neko_log%end_section()

  end subroutine nmsh_file_read_2d


  !> Write a mesh from to a binary Neko nmsh file
  subroutine nmsh_file_write(this, data, t)
    class(nmsh_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(nmsh_quad_t), allocatable :: nmsh_quad(:)
    type(nmsh_hex_t), allocatable :: nmsh_hex(:)
    type(nmsh_zone_t), allocatable :: nmsh_zone(:)
    type(nmsh_curve_el_t), allocatable :: nmsh_curve(:)
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, mpi_el_offset
    integer :: i, j, ierr, k
    integer :: nmsh_quad_size, nmsh_hex_size, nmsh_zone_size, nmsh_curve_size
    integer :: nzones, nzones_glb, nzones_offset
    integer :: ncurves, ncurves_glb, ncurves_offset
    integer :: el_idx, el_idx_glb
    class(element_t), pointer :: ep
    integer(i4), dimension(8), parameter :: vcyc_to_sym = &
         [1, 2, 4, 3, 5, 6, 8, 7] ! cyclic to symmetric vertex mapping

    select type (data)
    type is (mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select

    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)
    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)
    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)
    call MPI_Type_size(MPI_NMSH_CURVE, nmsh_curve_size, ierr)

    call neko_log%message('Writing data as a binary Neko file ' // this%get_fname())

    call MPI_File_open(NEKO_COMM, trim(this%get_fname()), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

    call MPI_File_write_all(fh, msh%glb_nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, msh%gdim, 1, MPI_INTEGER, status, ierr)

    call msh%reset_periodic_ids()

    if (msh%gdim .eq. 2) then
       allocate(nmsh_quad(msh%nelv))
       do i = 1, msh%nelv
          ep => msh%elements(i)%e
          nmsh_quad(i)%el_idx = ep%id()
          do j = 1, 4
             nmsh_quad(i)%v(j)%v_idx = ep%pts(vcyc_to_sym(j))%p%id()
             nmsh_quad(i)%v(j)%v_xyz = ep%pts(vcyc_to_sym(j))%p%x
          end do
       end do
       mpi_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
            int(msh%offset_el, i8) * int(nmsh_quad_size, i8)
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
       deallocate(nmsh_quad)
       mpi_el_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
            int(msh%glb_nelv, i8) * int(nmsh_quad_size, i8)
    else if (msh%gdim .eq. 3) then
       allocate(nmsh_hex(msh%nelv))
       do i = 1, msh%nelv
          ep => msh%elements(i)%e
          nmsh_hex(i)%el_idx = ep%id()
          do j = 1, 8
             nmsh_hex(i)%v(j)%v_idx = ep%pts(vcyc_to_sym(j))%p%id()
             nmsh_hex(i)%v(j)%v_xyz = ep%pts(vcyc_to_sym(j))%p%x
          end do
       end do
       mpi_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
            int(msh%offset_el, i8) * int(nmsh_hex_size, i8)
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_HEX, min(msh%nelv, max_write_nel), MPI_NMSH_HEX, status, ierr)
       do i = 1, msh%nelv/max_write_nel
          mpi_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
               int(msh%offset_el+i*max_write_nel, i8) * int(nmsh_hex_size, i8)
          call MPI_File_write_at_all(fh, mpi_offset, &
               nmsh_HEX(i*max_write_nel+1), &
               min(msh%nelv-i*max_write_nel, max_write_nel), &
               MPI_NMSH_HEX, status, ierr)
       end do
       deallocate(nmsh_hex)
       mpi_el_offset = int(2 * MPI_INTEGER_SIZE, i8) + &
            int(msh%glb_nelv, i8) * int(nmsh_hex_size, i8)
    else
       call neko_error('Invalid dimension of mesh')
    end if

    nzones = msh%periodic%size
    do i = 1, NEKO_MSH_MAX_ZLBLS
       nzones = nzones + msh%labeled_zones(i)%size
    end do

    call MPI_Allreduce(nzones, nzones_glb, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    nzones_offset = 0
    call MPI_Exscan(nzones, nzones_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    mpi_offset = mpi_el_offset
    call MPI_File_write_at_all(fh, mpi_offset, &
         nzones_glb, 1, MPI_INTEGER, status, ierr)

    if (nzones_glb .gt. 0) then
       allocate(nmsh_zone(nzones))

       if (nzones .gt. 0) then
          nmsh_zone(:)%type = 0

          j = 1
          do i = 1, msh%periodic%size
             nmsh_zone(j)%e = msh%elements(msh%periodic%facet_el(i)%x(2))%e%id()
             nmsh_zone(j)%f = msh%periodic%facet_el(i)%x(1)
             nmsh_zone(j)%p_e = msh%periodic%p_facet_el(i)%x(2)
             nmsh_zone(j)%p_f = msh%periodic%p_facet_el(i)%x(1)
             nmsh_zone(j)%glb_pt_ids = msh%periodic%p_ids(i)%x
             nmsh_zone(j)%type = 5
             j = j + 1
          end do

          do k = 1, NEKO_MSH_MAX_ZLBLS
             do i = 1, msh%labeled_zones(k)%size
                nmsh_zone(j)%e = &
                     msh%elements(msh%labeled_zones(k)%facet_el(i)%x(2))%e%id()
                nmsh_zone(j)%f = msh%labeled_zones(k)%facet_el(i)%x(1)
                nmsh_zone(j)%p_f = k
                nmsh_zone(j)%type = 7
                j = j + 1
             end do
          end do
       end if

       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8) + &
            int(nzones_offset, i8) * int(nmsh_zone_size, i8)
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)

       deallocate(nmsh_zone)
    end if

    ncurves = msh%curve%size


    call MPI_Allreduce(ncurves, ncurves_glb, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    ncurves_offset = 0
    call MPI_Exscan(ncurves, ncurves_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8) + &
         int(nzones_glb, i8)*int(nmsh_zone_size, i8)

    call MPI_File_write_at_all(fh, mpi_offset, &
         ncurves_glb, 1, MPI_INTEGER, status, ierr)

    if (ncurves_glb .gt. 0) then

       allocate(nmsh_curve(ncurves))

       do i = 1, ncurves
          nmsh_curve(i)%type = 0
       end do

       do i = 1, ncurves
          nmsh_curve(i)%e = msh%elements(msh%curve%curve_el(i)%el_idx)%e%id()
          nmsh_curve(i)%curve_data = msh%curve%curve_el(i)%curve_data
          nmsh_curve(i)%type = msh%curve%curve_el(i)%curve_type
       end do

       mpi_offset = mpi_el_offset + int(2*MPI_INTEGER_SIZE, i8) + &
            int(nzones_glb, i8) * int(nmsh_zone_size, i8) + &
            int(ncurves_offset, i8) * int(nmsh_curve_size, i8)

       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)
       deallocate(nmsh_curve)
    end if

    call MPI_File_sync(fh, ierr)
    call MPI_File_close(fh, ierr)
    call neko_log%message('Done')

    !
    ! Re-apply periodic facets
    ! (necessary if the mesh is going to be used after I/O)
    !
    do i = 1, msh%periodic%size
       el_idx_glb = msh%elements(msh%periodic%facet_el(i)%x(2))%e%id()
       if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
          call msh%apply_periodic_facet(msh%periodic%facet_el(i)%x(1), el_idx, &
               msh%periodic%p_facet_el(i)%x(1), msh%periodic%p_facet_el(i)%x(2), &
               msh%periodic%p_ids(i)%x)
       end if
    end do

  end subroutine nmsh_file_write

end module nmsh_file

