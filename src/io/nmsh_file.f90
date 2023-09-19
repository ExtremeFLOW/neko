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
  use generic_file
  use comm
  use mesh
  use utils
  use point
  use tuple
  use nmsh
  use element
  use datadist
  use mpi_types
  use mpi_f08
  use logger
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
    integer :: el_idx
    type(point_t) :: p(8)
    type(linear_dist_t) :: dist
    character(len=LOG_SIZE) :: log_buf

    select type(data)
    type is(mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select


    call neko_log%section("Mesh")
    call neko_log%message('Reading a binary Neko file ' // this%fname)

    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)
    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)
    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)

    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
      
    if (ierr > 0) then
       call neko_error('Could not open the mesh file ' // this%fname  // &
                       'for reading!')
    end if 
    call MPI_File_read_all(fh, nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)

    write(log_buf,1) gdim, nelv
1      format('gdim = ', i1, ', nelements =', i9)
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
       mpi_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(element_offset,i8) * int(nmsh_quad_size,i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
       do i = 1, nelv
          do j = 1, 4
             p(j) = point_t(nmsh_quad(i)%v(j)%v_xyz, nmsh_quad(i)%v(j)%v_idx)
          end do
          ! swap vertices to keep symmetric vertex numbering in neko
          call msh%add_element(i, p(1), p(2), p(4), p(3))
       end do
       deallocate(nmsh_quad)
       mpi_el_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(dist%num_global(),i8) * int(nmsh_quad_size,i8)
    else if (msh%gdim .eq. 3) then
       allocate(nmsh_hex(msh%nelv))
       mpi_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(element_offset,i8) * int(nmsh_hex_size,i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_hex, msh%nelv, MPI_NMSH_HEX, status, ierr)
       do i = 1, nelv
          do j = 1, 8
             p(j) = point_t(nmsh_hex(i)%v(j)%v_xyz, nmsh_hex(i)%v(j)%v_idx)
          end do
          ! swap vertices to keep symmetric vertex numbering in neko
          call msh%add_element(i, &
               p(1), p(2), p(4), p(3), p(5), p(6), p(8), p(7))
       end do
       deallocate(nmsh_hex)
       mpi_el_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(dist%num_global(),i8) * int(nmsh_hex_size,i8)
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
       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE,i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)
       
       do i = 1, nzones
          el_idx = nmsh_zone(i)%e
          if (el_idx .gt. msh%offset_el .and. &
               el_idx .le. msh%offset_el + msh%nelv) then             
             el_idx = el_idx - msh%offset_el
             select case(nmsh_zone(i)%type)
             case(1)
                call msh%mark_wall_facet(nmsh_zone(i)%f, el_idx)
             case(2)
                call msh%mark_inlet_facet(nmsh_zone(i)%f, el_idx)
             case(3)
                call msh%mark_outlet_facet(nmsh_zone(i)%f, el_idx)
             case(4)
                call msh%mark_sympln_facet(nmsh_zone(i)%f, el_idx)
             case(5)
                call msh%mark_periodic_facet(nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, nmsh_zone(i)%glb_pt_ids)
             case(6)
                call msh%mark_outlet_normal_facet(nmsh_zone(i)%f, el_idx)
             case(7)
                call msh%mark_labeled_facet(nmsh_zone(i)%f, el_idx,nmsh_zone(i)%p_f)
             end select
          end if
       end do
       !Apply facets, important that marking is finished
       do i = 1, nzones
          el_idx = nmsh_zone(i)%e
          if (el_idx .gt. msh%offset_el .and. &
               el_idx .le. msh%offset_el + msh%nelv) then             
             el_idx = el_idx - msh%offset_el
             select case(nmsh_zone(i)%type)
             case(5)
                call msh%apply_periodic_facet(nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, nmsh_zone(i)%glb_pt_ids)
             end select
          end if
       end do

       deallocate(nmsh_zone)
    end if
    call neko_log%message('Reading deformation data')

    mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE,i8) + int(nzones,i8)*int(nmsh_zone_size,i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         ncurves, 1, MPI_INTEGER, status, ierr)

    if (ncurves .gt. 0) then
       
       allocate(nmsh_curve(ncurves))
       mpi_offset = mpi_el_offset + int(2 * MPI_INTEGER_SIZE,i8) + int(nzones,i8)*int(nmsh_zone_size,i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)
       
       do i = 1, ncurves 
          el_idx = nmsh_curve(i)%e - msh%offset_el
          if (el_idx .gt. 0 .and. &
              el_idx .le. msh%nelv) then             
             call msh%mark_curve_element(el_idx, &
                  nmsh_curve(i)%curve_data, nmsh_curve(i)%type)
          end if
             
       end do
       
       deallocate(nmsh_curve)
    end if

    call MPI_File_close(fh, ierr)
    call neko_log%message('Mesh read, setting up connectivity')

    call msh%finalize()
    call neko_log%message('Done setting up mesh and connectivity')
    
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
    integer :: el_idx
    type(point_t) :: p(8)
    type(linear_dist_t) :: dist
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: depth = 1d0
    real(kind=dp) :: coord(3)
    type(tuple4_i4_t) :: glb_pt_ids


    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)
    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)
    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)

    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_read_all(fh, nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)

    write(log_buf,2) gdim
2      format('gdim = ', i1, ', no full 2d support, creating thin slab')
    call neko_log%message(log_buf)
    gdim = 3
    
    dist = linear_dist_t(nelv, pe_rank, pe_size, NEKO_COMM)
    nelv = dist%num_local()
    element_offset = dist%start_idx()
    
    call msh%init(gdim, nelv)
   
    allocate(nmsh_quad(msh%nelv))
    mpi_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(element_offset,i8) * int(nmsh_quad_size,i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
    do i = 1, nelv
       do j = 1, 4
          coord = nmsh_quad(i)%v(j)%v_xyz
          coord(3) = 0_rp
          p(j) = point_t(coord, nmsh_quad(i)%v(j)%v_idx)
       end do
       do j = 1, 4
          coord = nmsh_quad(i)%v(j)%v_xyz
          coord(3) = depth
          id = nmsh_quad(i)%v(j)%v_idx+msh%glb_nelv*8
          p(j+4) = point_t(coord, id)
       end do
       ! swap vertices to keep symmetric vertex numbering in neko
       call msh%add_element(i, &
            p(1), p(2), p(4), p(3), p(5), p(6), p(8), p(7))
    end do
    deallocate(nmsh_quad)
    mpi_el_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(dist%num_global(),i8) * int(nmsh_quad_size,i8)

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
       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE,i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)
       
       do i = 1, nzones
          el_idx = nmsh_zone(i)%e
          if (el_idx .gt. msh%offset_el .and. &
               el_idx .le. msh%offset_el + msh%nelv) then             
             el_idx = el_idx - msh%offset_el
             select case(nmsh_zone(i)%type)
             case(1)
                call msh%mark_wall_facet(nmsh_zone(i)%f, el_idx)
             case(2)
                call msh%mark_inlet_facet(nmsh_zone(i)%f, el_idx)
             case(3)
                call msh%mark_outlet_facet(nmsh_zone(i)%f, el_idx)
             case(4)
                call msh%mark_sympln_facet(nmsh_zone(i)%f, el_idx)
             case(5)
                nmsh_zone(i)%glb_pt_ids(3) = nmsh_zone(i)%glb_pt_ids(1)+msh%glb_nelv*8
                nmsh_zone(i)%glb_pt_ids(4) = nmsh_zone(i)%glb_pt_ids(2)+msh%glb_nelv*8
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
             case(6)
                call msh%mark_outlet_normal_facet(nmsh_zone(i)%f, el_idx)
             case(7)
                call msh%mark_labeled_facet(nmsh_zone(i)%f, el_idx,nmsh_zone(i)%p_f)
             end select
          end if
       end do
       !Apply facets, important that marking is finished
       do i = 1, nzones
          el_idx = nmsh_zone(i)%e
          if (el_idx .gt. msh%offset_el .and. &
               el_idx .le. msh%offset_el + msh%nelv) then             
             el_idx = el_idx - msh%offset_el
             select case(nmsh_zone(i)%type)
             case(5)
                call msh%apply_periodic_facet(nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, nmsh_zone(i)%glb_pt_ids)
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
         int(MPI_INTEGER_SIZE,i8) + int(nzones,i8)*int(nmsh_zone_size,i8)
    call MPI_File_read_at_all(fh, mpi_offset, &
         ncurves, 1, MPI_INTEGER, status, ierr)

    if (ncurves .gt. 0) then
       
       allocate(nmsh_curve(ncurves))
       mpi_offset = mpi_el_offset + &
            int(2*MPI_INTEGER_SIZE,i8) + int(nzones,i8)*int(nmsh_zone_size,i8)
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)
       
       do i = 1, ncurves 
          el_idx = nmsh_curve(i)%e - msh%offset_el
          if (el_idx .gt. 0 .and. &
              el_idx .le. msh%nelv) then             
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
    integer :: i, j, ierr, nelgv, element_offset, k
    integer :: nmsh_quad_size, nmsh_hex_size, nmsh_zone_size, nmsh_curve_size
    integer :: nzones, ncurves 
    class(element_t), pointer :: ep
    integer(i4),  dimension(8), parameter :: vcyc_to_sym = (/1, 2, 4, 3, 5, &
         & 6, 8, 7/) ! cyclic to symmetric vertex mapping

    select type(data)
    type is (mesh_t)
       msh => data
    class default
       call neko_error('Invalid output data')
    end select

    call MPI_Type_size(MPI_NMSH_QUAD, nmsh_quad_size, ierr)
    call MPI_Type_size(MPI_NMSH_HEX, nmsh_hex_size, ierr)
    call MPI_Type_size(MPI_NMSH_ZONE, nmsh_zone_size, ierr)
    call MPI_Type_size(MPI_NMSH_CURVE, nmsh_curve_size, ierr)

    call MPI_Reduce(msh%nelv, nelgv, 1, MPI_INTEGER, &
         MPI_SUM, 0, NEKO_COMM, ierr)
    element_offset = 0
    call MPI_Exscan(msh%nelv, element_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    call neko_log%message('Writing data as a binary Neko file ' // this%fname)

    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

    call MPI_File_write_all(fh, nelgv, 1, MPI_INTEGER, status, ierr)
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
       mpi_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(element_offset,i8) * int(nmsh_quad_size,i8)
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
       deallocate(nmsh_quad)
       mpi_el_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(nelgv,i8) * int(nmsh_quad_size,i8)
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
       mpi_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(element_offset,i8) * int(nmsh_hex_size,i8)
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_HEX, min(msh%nelv,max_write_nel), MPI_NMSH_HEX, status, ierr)
       do i = 1, msh%nelv/max_write_nel
          mpi_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(element_offset+i*max_write_nel,i8) * int(nmsh_hex_size,i8)
          call MPI_File_write_at_all(fh, mpi_offset, &
               nmsh_HEX(i*max_write_nel+1), min(msh%nelv-i*max_write_nel,max_write_nel), MPI_NMSH_HEX, status, ierr)
       end do
       deallocate(nmsh_hex)       
       mpi_el_offset = int(2 * MPI_INTEGER_SIZE,i8) + int(nelgv,i8) * int(nmsh_hex_size,i8)
    else 
       call neko_error('Invalid dimension of mesh')
    end if

    nzones = msh%wall%size + msh%inlet%size + msh%outlet%size + &
         msh%sympln%size + msh%periodic%size + msh%outlet_normal%size 

    do i = 1, NEKO_MSH_MAX_ZLBLS
       nzones = nzones + msh%labeled_zones(i)%size
    end do
    mpi_offset = mpi_el_offset
    call MPI_File_write_at_all(fh, mpi_offset, &
         nzones, 1, MPI_INTEGER, status, ierr)

    if (nzones .gt. 0) then
       allocate(nmsh_zone(nzones))
       
       nmsh_zone(:)%type = 0
       
       j = 1
       do i = 1, msh%wall%size
          nmsh_zone(j)%e = msh%wall%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone(j)%f = msh%wall%facet_el(i)%x(1)
          nmsh_zone(j)%type = 1
          j = j + 1
       end do

       do i = 1, msh%inlet%size
          nmsh_zone(j)%e = msh%inlet%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone(j)%f = msh%inlet%facet_el(i)%x(1)
          nmsh_zone(j)%type = 2
          j = j + 1
       end do

       do i = 1, msh%outlet%size
          nmsh_zone(j)%e = msh%outlet%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone(j)%f = msh%outlet%facet_el(i)%x(1)
          nmsh_zone(j)%type = 3
          j = j + 1
       end do

       do i = 1, msh%sympln%size
          nmsh_zone(j)%e = msh%sympln%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone(j)%f = msh%sympln%facet_el(i)%x(1)
          nmsh_zone(j)%type = 4
          j = j + 1
       end do

       do i = 1, msh%periodic%size
          nmsh_zone(j)%e = msh%periodic%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone(j)%f = msh%periodic%facet_el(i)%x(1)
          nmsh_zone(j)%p_e = msh%periodic%p_facet_el(i)%x(2)
          nmsh_zone(j)%p_f = msh%periodic%p_facet_el(i)%x(1)
          nmsh_zone(j)%glb_pt_ids = msh%periodic%p_ids(i)%x
          nmsh_zone(j)%type = 5
          j = j + 1
       end do
       do i = 1, msh%outlet_normal%size
          nmsh_zone(j)%e = msh%outlet_normal%facet_el(i)%x(2) + msh%offset_el
          nmsh_zone(j)%f = msh%outlet_normal%facet_el(i)%x(1)
          nmsh_zone(j)%type = 6
          j = j + 1
       end do
       do k = 1, NEKO_MSH_MAX_ZLBLS
          do i = 1, msh%labeled_zones(k)%size
             nmsh_zone(j)%e = msh%labeled_zones(k)%facet_el(i)%x(2) + msh%offset_el
             nmsh_zone(j)%f = msh%labeled_zones(k)%facet_el(i)%x(1)
             nmsh_zone(j)%p_f = k
             nmsh_zone(j)%type = 7
             j = j + 1
          end do
       end do
  
       
       mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE,i8)
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)
       
       deallocate(nmsh_zone)
    end if
 
    ncurves = msh%curve%size 
    mpi_offset = mpi_el_offset + int(MPI_INTEGER_SIZE,i8) + int(nzones,i8)*int(nmsh_zone_size,i8)

    call MPI_File_write_at_all(fh, mpi_offset, &
         ncurves, 1, MPI_INTEGER, status, ierr)

    if (ncurves .gt. 0) then
       allocate(nmsh_curve(ncurves))
       do i = 1, ncurves
          nmsh_curve(i)%type = 0
       end do
       
       do i = 1, ncurves
          nmsh_curve(i)%e = msh%curve%curve_el(i)%el_idx + msh%offset_el
          nmsh_curve(i)%curve_data = msh%curve%curve_el(i)%curve_data
          nmsh_curve(i)%type = msh%curve%curve_el(i)%curve_type
       end do
       mpi_offset = mpi_el_offset + int(2*MPI_INTEGER_SIZE,i8) + int(nzones,i8)*int(nmsh_zone_size,i8)
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)
       
       deallocate(nmsh_curve)
    end if
   
    call MPI_File_sync(fh, ierr)
    call MPI_File_close(fh, ierr)
    call neko_log%message('Done')

  end subroutine nmsh_file_write
  
end module nmsh_file
  
