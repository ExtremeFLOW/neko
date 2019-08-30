!> Legacy VTK file format
!! @details This modules defines interface to read/write legacy VTK file
!
module vtk_file
  use num_types
  use generic_file
  use utils
  use mesh
  use field
  use mesh_field
  use mpi
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
    type(field_t), pointer :: fld
    type(mesh_fld_t), pointer :: mfld
    integer :: i, j, vtk_type
    integer :: nid, np, ierr
    character(len=80) :: suffix,fname
    character(len=10) :: id_str
    integer:: suffix_pos

    call MPI_Comm_rank(MPI_COMM_WORLD, nid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
    
    select type(data)
    type is (mesh_t)
       msh => data
       nullify(fld)
       nullify(mfld)
    type is(field_t)
       msh => data%msh
       fld => data
       nullify(mfld)
    type is(mesh_fld_t)
       msh => data%msh
       mfld => data
       nullify(fld)
    class default
       call neko_error('Invalid data')
    end select

    if (np .gt. 1) then
       write(id_str,'(i10.10)') nid
       suffix_pos = scan(trim(this%fname), '.', back=.true.)
       open(unit=9, file=trim(this%fname(1:suffix_pos-1))//id_str//'.vtk')
    else
       open(unit=9, file=trim(this%fname))
    end if

    ! Write legacy header
    write(9, fmt='(A)') '# vtk DataFile Version 2.0'
    write(9, fmt='(A)') 'Neko'
    write(9, fmt='(A)') 'ASCII'
    write(9, fmt='(A)') 'DATASET UNSTRUCTURED_GRID'

    ! Dump coordinates
    write(9, fmt='(A,I8,A)') 'POINTS', msh%mpts,' double'
    do i = 1, msh%mpts
       write(9, fmt='(F15.8,F15.8,F15.8)') msh%points(i)%x
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

    if (associated(mfld)) then
       write(9, fmt='(A,I8)') 'CELL_DATA', msh%nelv
       write(9, fmt='(A,A,A,I8)') 'SCALARS ', trim(mfld%name), ' int', 1
       write(9, fmt='(A)') 'LOOKUP_TABLE default'
       do i = 1, msh%nelv
          write(9, fmt='(I8)') mfld%data(i)
       end do
    else if (associated(fld)) then
       !> @todo dump field data (scalar/vector etc)
    end if
    
    close(9)
  end subroutine vtk_file_write

  subroutine vtk_file_read(this, data)
    class(vtk_file_t) :: this
    class(*), target, intent(inout) :: data
  end subroutine vtk_file_read

end module vtk_file
