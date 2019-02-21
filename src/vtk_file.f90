!> Legacy VTK file format
!! @details This modules defines interface to read/write legacy VTK file
!
module vtk_file
  use num_types
  use generic_file
  use mesh
  implicit none
  private
  
  !> Interface for legacy VTK files
  type, extends(generic_file_t) :: vtk_file_t
   contains
     procedure :: read => vtk_file_read
     procedure :: write => vtk_file_write
  end type vtk_file_t

  public :: vtk_file_t

contains

  !> Write data in legacy VTK
  subroutine vtk_file_write(this, data)
    class(vtk_file_t) :: this
    class(*), target, intent(in) :: data
    type(mesh_t), pointer :: msh
    integer :: i, j, npts, vtk_type

    select type(data)
    type is (mesh_t)
       msh => data
    class default
       write(*,*) 'Invalid data'
       return
    end select

    open(unit=9, file=trim(this%fname))

    ! Write legacy header
    write(9, fmt='(A)') '# vtk DataFile Version 2.0'
    write(9, fmt='(A)') 'Neko'
    write(9, fmt='(A)') 'ASCII'
    write(9, fmt='(A)') 'DATASET UNSTRUCTURED_GRID'

    npts = 4
    if (msh%dim .eq. 3) npts = 8

    ! Dump coordinates
    write(9, fmt='(A,I8,A)') 'POINTS', msh%lelv*npts,' double'
    if (msh%dim .eq. 2) then
       do i = 1, msh%lelv
          do j = 1, npts
             write(9, fmt='(F15.8,F15.8,F15.8)') msh%xc(j,i), msh%yc(j,i), 0d0
          end do
       end do
    else
       do i = 1, msh%lelv
          do j = 1, npts
             write(9, fmt='(F15.8,F15.8,F15.8)') &
                  msh%xc(j,i), msh%yc(j,i), msh%zc(j,i)
          end do
       end do
    end if

    ! Dump cells
    write(9, fmt='(A,I8,I8)')  'CELLS', msh%lelv, msh%lelv*(npts+1)
    j = 0
    do i = 1, msh%lelv
       write(9, *) npts,(j, j=(i-1)*npts,(i-1)*npts + (npts-1))
    end do

    ! Dump cell type for each element
    write(9, fmt='(A,I8)') 'CELL_TYPES', msh%lelv
    vtk_type = 9
    if (msh%dim .eq. 3) vtk_type = 12
    do i = 1, msh%lelv
       write(9, fmt='(I2)') vtk_type
    end do
    close(9)
  end subroutine vtk_file_write

  subroutine vtk_file_read(this, data)
    class(vtk_file_t) :: this
    class(*), target, intent(inout) :: data
  end subroutine vtk_file_read

end module vtk_file
