! Copyright (c) 2019-2022, The Neko Authors
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
!> Legacy VTK file format
!! @details This module defines interface to read/write legacy VTK file
module vtk_file
  use num_types
  use generic_file
  use utils
  use mesh
  use field
  use dofmap
  use mesh_field
  use tet_mesh
  use tri_mesh
  use logger
  use comm
  implicit none
  private

  !> Interface for legacy VTK files
  type, public, extends(generic_file_t) :: vtk_file_t
   contains
     procedure :: read => vtk_file_read
     procedure :: write => vtk_file_write
  end type vtk_file_t

contains

  !> Write data in legacy VTK
  subroutine vtk_file_write(this, data, t)
    class(vtk_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(mesh_t), pointer :: msh => null()
    type(field_t), pointer :: fld => null()
    type(mesh_fld_t), pointer :: mfld => null()
    type(dofmap_t), pointer :: dm => null()
    type(tet_mesh_t), pointer :: tet_msh => null()
    type(tri_mesh_t), pointer :: tri_msh => null()
    character(len=10) :: id_str
    integer:: suffix_pos

    select type(data)
    type is (mesh_t)
       msh => data
    type is(field_t)
       msh => data%msh
       fld => data
    type is(mesh_fld_t)
       msh => data%msh
       mfld => data
    type is (dofmap_t)
       dm => data
    type is (tet_mesh_t)
       tet_msh => data
    type is (tri_mesh_t)
       tri_msh => data
    class default
       call neko_log%error('Invalid data')
    end select

    if (pe_size .gt. 1) then
       write(id_str,'(i10.10)') pe_rank
       suffix_pos = filename_suffix_pos(this%fname)
       open(unit=9, file=trim(this%fname(1:suffix_pos-1))//id_str//'.vtk')
    else
       open(unit=9, file=trim(this%fname))
    end if

    ! Write legacy header
    write(9, fmt='(A)') '# vtk DataFile Version 2.0'
    write(9, fmt='(A)') 'Neko'
    write(9, fmt='(A)') 'ASCII'

    if (associated(msh)) then
       write(9, fmt='(A)') 'DATASET UNSTRUCTURED_GRID'

       call vtk_file_write_mesh(9, msh)

       if (associated(mfld)) then
          call vtk_file_write_cell_data(9, mfld)
       else if (associated(fld)) then
          call vtk_file_write_point_data(9, fld)
       end if
    else if (associated(dm)) then
       write(9, fmt='(A)') 'DATASET POLYDATA'

       call vtk_file_write_dofmap_coordinates(9, dm)

       call vtk_file_write_dofmap_data(9, dm)
    else if (associated(tet_msh)) then
       write(9, fmt='(A)') 'DATASET UNSTRUCTURED_GRID'
       call vtk_file_write_tet_mesh(9, tet_msh)
    else if (associated(tri_msh)) then
       write(9, fmt='(A)') 'DATASET UNSTRUCTURED_GRID'
       call vtk_file_write_tri_mesh(9, tri_msh)
    else
       call neko_error('Invalid data')
    end if

    close(9)
  end subroutine vtk_file_write

  subroutine vtk_file_read(this, data)
    class(vtk_file_t) :: this
    class(*), target, intent(inout) :: data

    call neko_error('VTK file read not implemented')
  end subroutine vtk_file_read

  !> Write a mesh in legacy VTK format
  subroutine vtk_file_write_mesh(unit, msh)
    integer :: unit
    type(mesh_t), intent(inout) :: msh
    integer :: i, j, vtk_type
    integer,  dimension(8), parameter :: vcyc_to_sym = (/1, 2, 4, 3, &
                                                         5, 6, 8, 7/)
    ! Dump coordinates
    write(unit, fmt='(A,I8,A)') 'POINTS', msh%mpts,' double'
    do i = 1, msh%mpts
       write(unit, fmt='(F15.8,F15.8,F15.8)') real(msh%points(i)%x,dp)
    end do

    ! Dump cells
    write(unit, fmt='(A,I8,I8)')  'CELLS', msh%nelv, msh%nelv*(msh%npts+1)
    j = 0
    do i = 1, msh%nelv
       write(unit, fmt='(I8,8I8)') msh%npts, &
            (msh%get_local(msh%elements(i)%e%pts(vcyc_to_sym(j))%p) - 1, &
            j=1, msh%npts)
    end do

    ! Dump cell type for each element
    write(unit, fmt='(A,I8)') 'CELL_TYPES', msh%nelv
    vtk_type = 9
    if (msh%gdim .eq. 3) vtk_type = 12
    do i = 1, msh%nelv
       write(unit, fmt='(I2)') vtk_type
    end do

  end subroutine vtk_file_write_mesh

  !> Write a mesh field @a mfld as cell data
  subroutine vtk_file_write_cell_data(unit, mfld)
    integer :: unit
    type(mesh_fld_t), intent(in) :: mfld
    integer :: i

    write(unit, fmt='(A,I8)') 'CELL_DATA', mfld%msh%nelv
    write(unit, fmt='(A,A,A,I8)') 'SCALARS ', trim(mfld%name), ' int', 1
    write(unit, fmt='(A)') 'LOOKUP_TABLE default'

    do i = 1, mfld%msh%nelv
       write(unit, fmt='(I8)') mfld%data(i)
    end do

  end subroutine vtk_file_write_cell_data

  !> Write a field @a fld as point data
  !! @note High-order fields will be interpolated down
  !! to the low-order mesh
  subroutine vtk_file_write_point_data(unit, fld)
    integer :: unit
    type(field_t), intent(inout) :: fld
    real(kind=dp), allocatable :: point_data(:)
    integer :: i, j, lx, ly, lz, id(8)

    if ( (fld%Xh%lx - 1 .gt. 1) .or. &
         (fld%Xh%ly - 1 .gt. 1) .or. &
         (fld%Xh%lz - 1 .gt. 1)) then
       call neko_log%warning("Interpolate high-order data onto a low-order mesh")
    end if

    write(unit, fmt='(A,I8)') 'POINT_DATA', fld%msh%mpts
    write(unit, fmt='(A,A,A,I8)') 'SCALARS ', trim(fld%name), ' double', 1
    write(unit, fmt='(A)') 'LOOKUP_TABLE default'

    lx = fld%Xh%lx
    ly = fld%Xh%ly
    lz = fld%Xh%lz
    allocate(point_data(fld%msh%mpts))

    do i = 1, fld%msh%nelv
       do j = 1, fld%msh%npts
          id(j) = fld%msh%get_local(fld%msh%elements(i)%e%pts(j)%p)
       end do

       point_data(id(1)) = real(fld%x(1,1,1,i),dp)
       point_data(id(2)) = real(fld%x(lx,1,1,i),dp)
       point_data(id(3)) = real(fld%x(1,ly,1,i),dp)
       point_data(id(4)) = real(fld%x(lx,ly,1,i),dp)

       if (fld%msh%gdim .eq. 3) then
          point_data(id(5)) = real(fld%x(1,1,lz,i),dp)
          point_data(id(6)) = real(fld%x(lx,1,lz,i),dp)
          point_data(id(7)) = real(fld%x(1,ly,lz,i),dp)
          point_data(id(8)) = real(fld%x(lx,ly,lz,i),dp)
       end if

    end do

    write(unit, *) point_data

    deallocate(point_data)

  end subroutine vtk_file_write_point_data

  !> Write xyz-coordinates of a dofmap @a dm as points
  subroutine vtk_file_write_dofmap_coordinates(unit, dm)
    integer :: unit
    type(dofmap_t), intent(inout) :: dm
    integer :: i,j,k,l

    write(unit, fmt='(A,I8,A)') 'POINTS', size(dm%x),' double'

    do i = 1, dm%msh%nelv
       do l = 1, dm%Xh%lz
          do k = 1, dm%Xh%ly
             do j = 1, dm%Xh%lx
                write(unit, fmt='(F15.8,F15.8,F15.8)') &
                     real(dm%x(j,k,l,i),dp),&
                     real(dm%y(j,k,l,i),dp),&
                     real(dm%z(j,k,l,i),dp)
             end do
          end do
       end do
    end do

    write(unit, fmt='(A,I8,I8)') 'VERTICES', size(dm%x), 2*size(dm%x)
    do i = 1, size(dm%x)
       write(unit, fmt='(I8,I8)') 1,i-1
    end do


  end subroutine vtk_file_write_dofmap_coordinates

  !> Write a dofmap @a dm data as point data
  subroutine vtk_file_write_dofmap_data(unit, dm)
    integer :: unit
    type(dofmap_t), intent(inout) :: dm
    integer :: i, j, k, l

    write(unit, fmt='(A,I8)') 'POINT_DATA', size(dm%dof)
    write(unit, fmt='(A,A,A,I8)') 'SCALARS ', 'dof_id', ' integer', 1
    write(unit, fmt='(A)') 'LOOKUP_TABLE default'

    do i = 1, dm%msh%nelv
       do l = 1, dm%Xh%lz
          do k = 1, dm%Xh%ly
             do j = 1, dm%Xh%lx
                write(unit, fmt='(I8)') real(dm%dof(j,k,l,i),dp)
             end do
          end do
       end do
    end do

    write(unit, fmt='(A,A,A,I8)') 'SCALARS ', 'shared_dof', ' integer', 1
    write(unit, fmt='(A)') 'LOOKUP_TABLE default'

    do i = 1, dm%msh%nelv
       do l = 1, dm%Xh%lz
          do k = 1, dm%Xh%ly
             do j = 1, dm%Xh%lx
                if (dm%shared_dof(j,k,l,i)) then
                   write(unit, fmt='(I8)') 1
                else
                   write(unit, fmt='(I8)') 0
                end if
             end do
          end do
       end do
    end do

  end subroutine vtk_file_write_dofmap_data

  !> Write a tetrahedral mesh in legacy VTK format
  subroutine vtk_file_write_tet_mesh(unit, tet_msh)
    integer :: unit
    type(tet_mesh_t), intent(inout) :: tet_msh
    integer, parameter :: npts = 4
    integer :: i, j, vtk_type

    ! Dump coordinates
    write(unit, fmt='(A,I8,A)') 'POINTS', tet_msh%msh%mpts,' double'
    do i = 1, tet_msh%msh%mpts
       write(unit, fmt='(F15.8,F15.8,F15.8)') real(tet_msh%msh%points(i)%x,dp)
    end do

    ! Dump cells
    write(unit, fmt='(A,I8,I8)')  'CELLS', tet_msh%nelv, tet_msh%nelv*(npts+1)
    j = 0
    do i = 1, tet_msh%nelv
       write(unit, fmt='(I8,8I8)') npts, &
            (tet_msh%msh%get_local(tet_msh%el(i)%pts(j)%p) - 1, &
            j=1, npts)
    end do

    ! Dump cell type for each element
    write(unit, fmt='(A,I8)') 'CELL_TYPES', tet_msh%nelv
    vtk_type = 10
    do i = 1, tet_msh%nelv
       write(unit, fmt='(I2)') vtk_type
    end do

  end subroutine vtk_file_write_tet_mesh

  !> Write a triangular mesh in legacy VTK format
  subroutine vtk_file_write_tri_mesh(unit, tri_msh)
    integer :: unit
    type(tri_mesh_t), intent(inout) :: tri_msh
    integer, parameter :: npts = 3
    integer :: i, j, vtk_type

    ! Dump coordinates
    write(unit, fmt='(A,I8,A)') 'POINTS', tri_msh%mpts,' double'
    do i = 1, tri_msh%mpts
       write(unit, fmt='(F15.8,F15.8,F15.8)') real(tri_msh%points(i)%x,dp)
    end do

    ! Dump cells
    write(unit, fmt='(A,I8,I8)')  'CELLS', tri_msh%nelv, tri_msh%nelv*(npts+1)
    j = 0
    do i = 1, tri_msh%nelv
       write(unit, fmt='(I8,8I8)') npts, &
            (tri_msh%el(i)%pts(j)%p%id() - 1, j=1, npts)
    end do

    ! Dump cell type for each element
    write(unit, fmt='(A,I8)') 'CELL_TYPES', tri_msh%nelv
    vtk_type = 5
    do i = 1, tri_msh%nelv
       write(unit, fmt='(I2)') vtk_type
    end do

  end subroutine vtk_file_write_tri_mesh

end module vtk_file
