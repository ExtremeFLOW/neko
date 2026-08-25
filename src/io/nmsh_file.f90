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
  use generic_file, only : generic_file_t
  use comm, only : NEKO_COMM, pe_rank, pe_size
  use num_types, only : rp, dp, i4, i8
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS
  use utils, only : neko_error
  use point, only : point_t
  use tuple, only : tuple4_i4_t
  use nmsh, only : nmsh_hex_t, nmsh_quad_t, nmsh_zone_t, nmsh_curve_el_t
  use element, only : element_t
  use datadist, only : linear_dist_t
  use neko_mpi_types, only : MPI_NMSH_HEX, MPI_NMSH_QUAD, MPI_NMSH_ZONE, &
       MPI_NMSH_CURVE, MPI_INTEGER_SIZE
  use mpi_f08, only : MPI_Wtime, MPI_Status, MPI_File, MPI_OFFSET_KIND, &
       MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_MODE_RDONLY, MPI_INFO_NULL, &
       MPI_File_open, MPI_File_close, MPI_File_read_all, MPI_File_write_all, &
       MPI_File_write_at_all, MPI_File_read_at_all, MPI_INTEGER, MPI_SUM, &
       MPI_MAX, MPI_Exscan, MPI_Barrier, MPI_Type_size, MPI_Allreduce, &
       MPI_File_sync, MPI_Sendrecv, MPI_Get_count
  use logger, only : neko_log, LOG_SIZE
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
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, mpi_el_offset
    integer :: i, j, ierr, element_offset
    integer :: nmsh_quad_size, nmsh_hex_size, nmsh_zone_size
    integer :: nelv, gdim, nzones, ncurves
    type(point_t), target :: p(8)
    type(linear_dist_t) :: dist
    character(len=LOG_SIZE) :: log_buf
    real(kind=dp) :: t_start, t_end

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
          mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8)
          call nmsh_file_read_zones(fh, mpi_offset, nzones, msh)
       end if
       call neko_log%message('Reading deformation data')

       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8) + &
            int(nzones, i8)*int(nmsh_zone_size, i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            ncurves, 1, MPI_INTEGER, status, ierr)

       if (ncurves .gt. 0) then
          mpi_offset = mpi_el_offset + int(2 * MPI_INTEGER_SIZE, i8) + &
               int(nzones, i8)*int(nmsh_zone_size, i8)
          call nmsh_file_read_curves(fh, mpi_offset, ncurves, msh)
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
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, mpi_el_offset
    integer :: i, j, ierr, element_offset, id
    integer :: nmsh_quad_size, nmsh_hex_size, nmsh_zone_size
    integer :: nelv, gdim, nzones, ncurves
    integer :: el_idx
    type(point_t) :: p(8)
    type(linear_dist_t) :: dist
    character(len=LOG_SIZE) :: log_buf
    real(kind=dp) :: depth = 1.0_dp
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
       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE, i8)
       call nmsh_file_read_zones_2d(fh, mpi_offset, nzones, msh)

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
    end if

    mpi_offset = mpi_el_offset + &
         int(MPI_INTEGER_SIZE, i8) + int(nzones, i8)*int(nmsh_zone_size, i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         ncurves, 1, MPI_INTEGER, status, ierr)

    if (ncurves .gt. 0) then
       mpi_offset = mpi_el_offset + &
            int(2*MPI_INTEGER_SIZE, i8) + &
            int(nzones, i8)*int(nmsh_zone_size, i8)
       call nmsh_file_read_curves(fh, mpi_offset, ncurves, msh, depth)
    end if

    call MPI_File_close(fh, ierr)

    call msh%finalize()

    call neko_log%end_section()

  end subroutine nmsh_file_read_2d

  !> Read zone data in chunks and ring-pass between ranks,
  !! avoiding storage of the full zone list on every rank.
  !! Uses two ring traversals: first marks all matching facets on every
  !! rank, then applies periodic facets (apply must follow mark globally
  !! since apply mutates point ids that mark reads via get_facet_ids).
  subroutine nmsh_file_read_zones(fh, base_offset, nzones, msh)
    type(MPI_File), intent(inout) :: fh
    integer(kind=MPI_OFFSET_KIND), intent(in) :: base_offset
    integer, intent(in) :: nzones
    type(mesh_t), intent(inout) :: msh
    type(nmsh_zone_t), allocatable :: zone_send(:), zone_recv(:)
    type(linear_dist_t) :: dist
    type(MPI_Status) :: status
    integer(kind=MPI_OFFSET_KIND) :: mpi_offset
    integer :: nmsh_zone_size, nlocal, max_recv, n_recv
    integer :: i, ierr, src, dst, step, el_idx, el_idx_glb

    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)

    dist = linear_dist_t(nzones, pe_rank, pe_size, NEKO_COMM)
    nlocal = dist%num_local()

    allocate(zone_send(nlocal))

    mpi_offset = base_offset + &
         int(dist%start_idx(), i8) * int(nmsh_zone_size, i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         zone_send, nlocal, MPI_NMSH_ZONE, status, ierr)

    call MPI_Allreduce(nlocal, max_recv, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)
    allocate(zone_recv(max_recv))

    ! Mark pass: process local chunk
    do i = 1, nlocal
       el_idx_glb = zone_send(i)%e
       if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
          select case (zone_send(i)%type)
          case (5)
             call msh%mark_periodic_facet(zone_send(i)%f, el_idx, &
                  zone_send(i)%p_f, zone_send(i)%p_e, &
                  zone_send(i)%glb_pt_ids)
          case (7)
             call msh%mark_labeled_facet(zone_send(i)%f, el_idx, &
                  zone_send(i)%p_f)
          end select
       end if
    end do

    ! Mark pass: ring-shift through all other ranks
    do step = 1, pe_size - 1
       src = modulo(pe_rank - step + pe_size, pe_size)
       dst = modulo(pe_rank + step, pe_size)
       call MPI_Sendrecv(zone_send, nlocal, MPI_NMSH_ZONE, dst, 0, &
            zone_recv, max_recv, MPI_NMSH_ZONE, src, 0, &
            NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_NMSH_ZONE, n_recv, ierr)
       do i = 1, n_recv
          el_idx_glb = zone_recv(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             select case (zone_recv(i)%type)
             case (5)
                call msh%mark_periodic_facet(zone_recv(i)%f, el_idx, &
                     zone_recv(i)%p_f, zone_recv(i)%p_e, &
                     zone_recv(i)%glb_pt_ids)
             case (7)
                call msh%mark_labeled_facet(zone_recv(i)%f, el_idx, &
                     zone_recv(i)%p_f)
             end select
          end if
       end do
    end do

    ! Apply pass (periodic only): process local chunk first
    do i = 1, nlocal
       el_idx_glb = zone_send(i)%e
       if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
          if (zone_send(i)%type .eq. 5) then
             call msh%apply_periodic_facet(zone_send(i)%f, el_idx, &
                  zone_send(i)%p_f, zone_send(i)%p_e, &
                  zone_send(i)%glb_pt_ids)
          end if
       end if
    end do

    ! Apply pass: ring-shift through all other ranks
    do step = 1, pe_size - 1
       src = modulo(pe_rank - step + pe_size, pe_size)
       dst = modulo(pe_rank + step, pe_size)
       call MPI_Sendrecv(zone_send, nlocal, MPI_NMSH_ZONE, dst, 1, &
            zone_recv, max_recv, MPI_NMSH_ZONE, src, 1, &
            NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_NMSH_ZONE, n_recv, ierr)
       do i = 1, n_recv
          el_idx_glb = zone_recv(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             if (zone_recv(i)%type .eq. 5) then
                call msh%apply_periodic_facet(zone_recv(i)%f, el_idx, &
                     zone_recv(i)%p_f, zone_recv(i)%p_e, &
                     zone_recv(i)%glb_pt_ids)
             end if
          end if
       end do
    end do

    deallocate(zone_send)
    deallocate(zone_recv)

  end subroutine nmsh_file_read_zones

  !> Compute the extruded periodic point ids for the 2D thin-slab path.
  !! Pure function of the zone fields; safe to invoke independently from
  !! the mark and apply ring passes.
  pure subroutine nmsh_file_extrude_periodic_ids(f, glb_pt_ids, glb_nelv, ids)
    integer, intent(in) :: f
    integer, intent(in) :: glb_pt_ids(4)
    integer, intent(in) :: glb_nelv
    integer, intent(out) :: ids(4)
    integer :: pt(4)

    pt(1) = glb_pt_ids(1)
    pt(2) = glb_pt_ids(2)
    pt(3) = glb_pt_ids(1) + glb_nelv * 8
    pt(4) = glb_pt_ids(2) + glb_nelv * 8

    if (f .eq. 1 .or. f .eq. 2) then
       ids(1) = pt(1)
       ids(2) = pt(3)
       ids(3) = pt(4)
       ids(4) = pt(2)
    else
       ids(1) = pt(1)
       ids(2) = pt(2)
       ids(3) = pt(4)
       ids(4) = pt(3)
    end if
  end subroutine nmsh_file_extrude_periodic_ids

  !> Read zone data in chunks and ring-pass between ranks for the
  !! 2D thin-slab path. Two ring traversals as in nmsh_file_read_zones,
  !! but periodic entries also need their glb_pt_ids reordered to
  !! account for the extruded layer.
  subroutine nmsh_file_read_zones_2d(fh, base_offset, nzones, msh)
    type(MPI_File), intent(inout) :: fh
    integer(kind=MPI_OFFSET_KIND), intent(in) :: base_offset
    integer, intent(in) :: nzones
    type(mesh_t), intent(inout) :: msh
    type(nmsh_zone_t), allocatable :: zone_send(:), zone_recv(:)
    type(linear_dist_t) :: dist
    type(MPI_Status) :: status
    integer(kind=MPI_OFFSET_KIND) :: mpi_offset
    integer :: nmsh_zone_size, nlocal, max_recv, n_recv
    integer :: i, ierr, src, dst, step, el_idx, el_idx_glb
    integer :: ids(4)

    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)

    dist = linear_dist_t(nzones, pe_rank, pe_size, NEKO_COMM)
    nlocal = dist%num_local()

    allocate(zone_send(nlocal))

    mpi_offset = base_offset + &
         int(dist%start_idx(), i8) * int(nmsh_zone_size, i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         zone_send, nlocal, MPI_NMSH_ZONE, status, ierr)

    call MPI_Allreduce(nlocal, max_recv, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)
    allocate(zone_recv(max_recv))

    ! Mark pass: process local chunk. We recompute the extruded ids on
    ! demand rather than storing them back into zone_send, so the buffer
    ! we ring-shift always carries the original on-disk values.
    do i = 1, nlocal
       el_idx_glb = zone_send(i)%e
       if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
          select case (zone_send(i)%type)
          case (5)
             call nmsh_file_extrude_periodic_ids(zone_send(i)%f, &
                  zone_send(i)%glb_pt_ids, msh%glb_nelv, ids)
             call msh%mark_periodic_facet(zone_send(i)%f, el_idx, &
                  zone_send(i)%p_f, zone_send(i)%p_e, ids)
          case (7)
             call msh%mark_labeled_facet(zone_send(i)%f, el_idx, &
                  zone_send(i)%p_f)
          end select
       end if
    end do

    ! Mark pass: ring-shift through all other ranks
    do step = 1, pe_size - 1
       src = modulo(pe_rank - step + pe_size, pe_size)
       dst = modulo(pe_rank + step, pe_size)
       call MPI_Sendrecv(zone_send, nlocal, MPI_NMSH_ZONE, dst, 0, &
            zone_recv, max_recv, MPI_NMSH_ZONE, src, 0, &
            NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_NMSH_ZONE, n_recv, ierr)
       do i = 1, n_recv
          el_idx_glb = zone_recv(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             select case (zone_recv(i)%type)
             case (5)
                call nmsh_file_extrude_periodic_ids(zone_recv(i)%f, &
                     zone_recv(i)%glb_pt_ids, msh%glb_nelv, ids)
                call msh%mark_periodic_facet(zone_recv(i)%f, el_idx, &
                     zone_recv(i)%p_f, zone_recv(i)%p_e, ids)
             case (7)
                call msh%mark_labeled_facet(zone_recv(i)%f, el_idx, &
                     zone_recv(i)%p_f)
             end select
          end if
       end do
    end do

    ! Apply pass (periodic only): local chunk
    do i = 1, nlocal
       el_idx_glb = zone_send(i)%e
       if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
          if (zone_send(i)%type .eq. 5) then
             call nmsh_file_extrude_periodic_ids(zone_send(i)%f, &
                  zone_send(i)%glb_pt_ids, msh%glb_nelv, ids)
             call msh%apply_periodic_facet(zone_send(i)%f, el_idx, &
                  zone_send(i)%p_f, zone_send(i)%p_e, ids)
          end if
       end if
    end do

    ! Apply pass: ring-shift through all other ranks
    do step = 1, pe_size - 1
       src = modulo(pe_rank - step + pe_size, pe_size)
       dst = modulo(pe_rank + step, pe_size)
       call MPI_Sendrecv(zone_send, nlocal, MPI_NMSH_ZONE, dst, 1, &
            zone_recv, max_recv, MPI_NMSH_ZONE, src, 1, &
            NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_NMSH_ZONE, n_recv, ierr)
       do i = 1, n_recv
          el_idx_glb = zone_recv(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             if (zone_recv(i)%type .eq. 5) then
                call nmsh_file_extrude_periodic_ids(zone_recv(i)%f, &
                     zone_recv(i)%glb_pt_ids, msh%glb_nelv, ids)
                call msh%apply_periodic_facet(zone_recv(i)%f, el_idx, &
                     zone_recv(i)%p_f, zone_recv(i)%p_e, ids)
             end if
          end if
       end do
    end do

    deallocate(zone_send)
    deallocate(zone_recv)

  end subroutine nmsh_file_read_zones_2d

  !> Read curve element data in chunks and ring-pass between ranks,
  !! avoiding storage of the full curve list on every rank.
  !! If @a extrusion_depth is present, copy the four 2D curve edges to the
  !! corresponding upper edges of the extruded slab.
  subroutine nmsh_file_read_curves(fh, base_offset, ncurves, msh, &
       extrusion_depth)
    type(MPI_File), intent(inout) :: fh
    integer(kind=MPI_OFFSET_KIND), intent(in) :: base_offset
    integer, intent(in) :: ncurves
    type(mesh_t), intent(inout) :: msh
    real(kind=dp), intent(in), optional :: extrusion_depth
    type(nmsh_curve_el_t), allocatable :: curve_send(:), curve_recv(:)
    type(linear_dist_t) :: dist
    type(MPI_Status) :: status
    integer(kind=MPI_OFFSET_KIND) :: mpi_offset
    integer :: nmsh_curve_size, nlocal, max_recv, n_recv
    integer :: i, j, ierr, src, dst, step, el_idx, el_idx_glb

    call MPI_Type_size(MPI_NMSH_CURVE, nmsh_curve_size, ierr)

    dist = linear_dist_t(ncurves, pe_rank, pe_size, NEKO_COMM)
    nlocal = dist%num_local()

    allocate(curve_send(nlocal))

    mpi_offset = base_offset + &
         int(dist%start_idx(), i8) * int(nmsh_curve_size, i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         curve_send, nlocal, MPI_NMSH_CURVE, status, ierr)

    if (present(extrusion_depth)) then
       do i = 1, nlocal
          curve_send(i)%curve_data(:,5:8) = &
               curve_send(i)%curve_data(:,1:4)
          curve_send(i)%type(5:8) = curve_send(i)%type(1:4)
          do j = 5, 8
             if (curve_send(i)%type(j) .eq. 4) then
                curve_send(i)%curve_data(3,j) = extrusion_depth
             end if
          end do
       end do
    end if

    call MPI_Allreduce(nlocal, max_recv, 1, MPI_INTEGER, MPI_MAX, &
         NEKO_COMM, ierr)
    allocate(curve_recv(max_recv))

    do i = 1, nlocal
       el_idx_glb = curve_send(i)%e
       if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
          call msh%mark_curve_element(el_idx, &
               curve_send(i)%curve_data, curve_send(i)%type)
       end if
    end do

    do step = 1, pe_size - 1
       src = modulo(pe_rank - step + pe_size, pe_size)
       dst = modulo(pe_rank + step, pe_size)
       call MPI_Sendrecv(curve_send, nlocal, MPI_NMSH_CURVE, dst, 0, &
            curve_recv, max_recv, MPI_NMSH_CURVE, src, 0, &
            NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_NMSH_CURVE, n_recv, ierr)
       do i = 1, n_recv
          el_idx_glb = curve_recv(i)%e
          if (msh%htel%get(el_idx_glb, el_idx) .eq. 0) then
             call msh%mark_curve_element(el_idx, &
                  curve_recv(i)%curve_data, curve_recv(i)%type)
          end if
       end do
    end do

    deallocate(curve_send)
    deallocate(curve_recv)

  end subroutine nmsh_file_read_curves


  !> Write a mesh from to a binary Neko nmsh file
  subroutine nmsh_file_write(this, data, t)
    class(nmsh_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=dp), intent(in), optional :: t
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

    call neko_log%message('Writing data as a binary Neko file ' // &
         this%get_fname())

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
               msh%periodic%p_facet_el(i)%x(1), &
               msh%periodic%p_facet_el(i)%x(2), msh%periodic%p_ids(i)%x)
       end if
    end do

  end subroutine nmsh_file_write

end module nmsh_file
