!> Legacy VTK file format
!! @details This modules defines interface to read/write legacy VTK file
!
module vtk_file
  use num_types
  use generic_file
  use utils
  use mesh
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
  subroutine vtk_file_write(this, data)
    class(vtk_file_t), intent(in) :: this
    class(*), target, intent(in) :: data
    type(mesh_t), pointer :: msh
    integer :: i, j, vtk_type

    select type(data)
    type is (mesh_t)
       msh => data
    class default
       call neko_error('Invalid data')
    end select

    open(unit=9, file=trim(this%fname))

    ! Write legacy header
    write(9, fmt='(A)') '# vtk DataFile Version 2.0'
    write(9, fmt='(A)') 'Neko'
    write(9, fmt='(A)') 'ASCII'
    write(9, fmt='(A)') 'DATASET UNSTRUCTURED_GRID'

    ! Dump coordinates (yes we're keeping duplicates)
    write(9, fmt='(A,I8,A)') 'POINTS', msh%nelv*msh%npts,' double'
    do i = 1, msh%nelv
       do j = 1, msh%npts
          write(9, fmt='(F15.8,F15.8,F15.8)') msh%elements(i)%e%pts(j)%p%x
       end do
    end do

    ! Dump cells
    write(9, fmt='(A,I8,I8)')  'CELLS', msh%nelv, msh%nelv*(msh%npts+1)
    j = 0
    do i = 1, msh%nelv
       write(9, *) msh%npts,(msh%elements(i)%e%pts(j)%p%id() - 1, j=1, msh%npts)
    end do

    ! Dump cell type for each element
    write(9, fmt='(A,I8)') 'CELL_TYPES', msh%nelv
    vtk_type = 9
    if (msh%gdim .eq. 3) vtk_type = 12
    do i = 1, msh%nelv
       write(9, fmt='(I2)') vtk_type
    end do
    close(9)
  end subroutine vtk_file_write

  subroutine vtk_file_read(this, data)
    class(vtk_file_t) :: this
    class(*), target, intent(inout) :: data
  end subroutine vtk_file_read

end module vtk_file
