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
  use nmsh
  use datadist
  use mpi_types
  use mpi_f08
  use logger
  implicit none
  
  private

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
    integer :: i, j, ierr, nelgv, element_offset
    integer :: nmsh_quad_size, nmsh_hex_size, nmsh_zone_size
    class(element_t), pointer :: ep
    integer :: nelv, gdim, nread, nzones, ncurves
    integer :: el_idx, ids(4), bcs(7), thing
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
    call MPI_File_read_all(fh, nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)

    write(log_buf,1) gdim, nelv
1      format('gdim = ', i1, ', nelements =', i7)
    call neko_log%message(log_buf)
    
    dist = linear_dist_t(nelv, pe_rank, pe_size, NEKO_COMM)
    nelv = dist%num_local()
    element_offset = dist%start_idx()
    
    call mesh_init(msh, gdim, nelv)
   

    if (msh%gdim .eq. 2) then
       allocate(nmsh_quad(msh%nelv))
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_quad_size
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
       do i = 1, nelv
          do j = 1, 4
             p(j) = point_t(nmsh_quad(i)%v(j)%v_xyz, nmsh_quad(i)%v(j)%v_idx)
          end do
          call mesh_add_element(msh, i, p(1), p(2), p(3), p(4))
       end do
       deallocate(nmsh_quad)
       mpi_el_offset = 2 * MPI_INTEGER_SIZE + dist%num_global() * nmsh_quad_size
    else if (msh%gdim .eq. 3) then
       allocate(nmsh_hex(msh%nelv))
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_hex_size
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_hex, msh%nelv, MPI_NMSH_HEX, status, ierr)
       do i = 1, nelv
          do j = 1, 8
             p(j) = point_t(nmsh_hex(i)%v(j)%v_xyz, nmsh_hex(i)%v(j)%v_idx)
          end do
          call mesh_add_element(msh, i, &
               p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))
       end do
       deallocate(nmsh_hex)
       mpi_el_offset = 2 * MPI_INTEGER_SIZE + dist%num_global() * nmsh_hex_size
    else        
       if (pe_rank .eq. 0) call neko_error('Invalid dimension of mesh')
    end if

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
       mpi_offset = mpi_el_offset + MPI_INTEGER_SIZE
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)
       
       do i = 1, nzones
          el_idx = nmsh_zone(i)%e
          if (el_idx .gt. msh%offset_el .and. &
               el_idx .le. msh%offset_el + msh%nelv) then             
             el_idx = el_idx - msh%offset_el
             select case(nmsh_zone(i)%type)
             case(1)
                call mesh_mark_wall_facet(msh, nmsh_zone(i)%f, el_idx)
             case(2)
                call mesh_mark_inlet_facet(msh, nmsh_zone(i)%f, el_idx)
             case(3)
                call mesh_mark_outlet_facet(msh, nmsh_zone(i)%f, el_idx)
             case(4)
                call mesh_mark_sympln_facet(msh, nmsh_zone(i)%f, el_idx)
             case(5)
                call mesh_mark_periodic_facet(msh, nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, nmsh_zone(i)%glb_pt_ids)
             case(6)
                call mesh_mark_outlet_normal_facet(msh, nmsh_zone(i)%f, el_idx)
             case(7)
                call mesh_mark_labeled_facet(msh, nmsh_zone(i)%f, el_idx,nmsh_zone(i)%p_f)
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
                call mesh_apply_periodic_facet(msh, nmsh_zone(i)%f, el_idx, &
                     nmsh_zone(i)%p_f, nmsh_zone(i)%p_e, nmsh_zone(i)%glb_pt_ids)
             end select
          end if
       end do


       deallocate(nmsh_zone)
    end if

    mpi_offset = mpi_el_offset + MPI_INTEGER_SIZE + nzones*nmsh_zone_size
    call MPI_File_read_at_all(fh, mpi_offset, &
         ncurves, 1, MPI_INTEGER, status, ierr)

    if (ncurves .gt. 0) then
       
       allocate(nmsh_curve(ncurves))
       mpi_offset = mpi_el_offset + 2*MPI_INTEGER_SIZE + nzones*nmsh_zone_size
       call MPI_File_read_at_all(fh, mpi_offset, &
            nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)
       
       do i = 1, ncurves 
          el_idx = nmsh_curve(i)%e - msh%offset_el
          if (el_idx .gt. 0 .and. &
              el_idx .le. msh%nelv) then             
             call mesh_mark_curve_element(msh, el_idx, nmsh_curve(i)%curve_data, nmsh_curve(i)%type)
          end if
             
       end do
       
       deallocate(nmsh_curve)
    end if

    call MPI_File_close(fh, ierr)

    call mesh_finalize(msh)
    
    call neko_log%end_section()
       
  end subroutine nmsh_file_read

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

    call MPI_File_write_all(fh, msh%nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, msh%gdim, 1, MPI_INTEGER, status, ierr)

    call mesh_reset_periodic_ids(msh)

    if (msh%gdim .eq. 2) then
       allocate(nmsh_quad(msh%nelv))       
       do i = 1, msh%nelv
          ep => msh%elements(i)%e
          nmsh_quad(i)%el_idx = ep%id()
          do j = 1, 4
             nmsh_quad(i)%v(j)%v_idx = ep%pts(j)%p%id()
             nmsh_quad(i)%v(j)%v_xyz = ep%pts(j)%p%x
          end do
       end do
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_quad_size
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_quad, msh%nelv, MPI_NMSH_QUAD, status, ierr)
       deallocate(nmsh_quad)
       mpi_el_offset = 2 * MPI_INTEGER_SIZE + msh%nelv * nmsh_quad_size
    else if (msh%gdim .eq. 3) then
       allocate(nmsh_hex(msh%nelv))       
       do i = 1, msh%nelv
          ep => msh%elements(i)%e
          nmsh_hex(i)%el_idx = ep%id()
          do j = 1, 8
             nmsh_hex(i)%v(j)%v_idx = ep%pts(j)%p%id()
             nmsh_hex(i)%v(j)%v_xyz = ep%pts(j)%p%x
          end do
       end do
       mpi_offset = 2 * MPI_INTEGER_SIZE + element_offset * nmsh_hex_size
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_HEX, msh%nelv, MPI_NMSH_HEX, status, ierr)
       deallocate(nmsh_hex)       
       mpi_el_offset = 2 * MPI_INTEGER_SIZE + msh%nelv * nmsh_hex_size
    else 
       call neko_error('Invalid dimension of mesh')
    end if

    nzones = msh%wall%size + msh%inlet%size + msh%outlet%size + &
         msh%sympln%size + msh%periodic%size + msh%outlet_normal%size 

    do i = 1, msh%max_labels
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
       do k = 1, msh%max_labels
          do i = 1, msh%labeled_zones(k)%size
             nmsh_zone(j)%e = msh%labeled_zones(k)%facet_el(i)%x(2) + msh%offset_el
             nmsh_zone(j)%f = msh%labeled_zones(k)%facet_el(i)%x(1)
             nmsh_zone(j)%p_f = k
             nmsh_zone(j)%type = 7
             j = j + 1
          end do
       end do
  
       
       mpi_offset = mpi_el_offset + MPI_INTEGER_SIZE
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_zone, nzones, MPI_NMSH_ZONE, status, ierr)
       
       deallocate(nmsh_zone)
    end if
 
    ncurves = msh%curve%size 
    mpi_offset = mpi_el_offset + MPI_INTEGER_SIZE + nzones*nmsh_zone_size

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
       mpi_offset = mpi_el_offset + 2*MPI_INTEGER_SIZE + nzones*nmsh_zone_size
       call MPI_File_write_at_all(fh, mpi_offset, &
            nmsh_curve, ncurves, MPI_NMSH_CURVE, status, ierr)
       
       deallocate(nmsh_curve)
    end if
   
    call MPI_File_close(fh, ierr)
    call neko_log%message('Done')

  end subroutine nmsh_file_write
  
end module nmsh_file
  
