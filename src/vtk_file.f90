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
  subroutine vtk_file_write(this, data)
    class(vtk_file_t), intent(in) :: this
    class(*), target, intent(in) :: data
    type(mesh_t), pointer :: msh
    type(field_t), pointer :: fld
    type(mesh_fld_t), pointer :: mfld
    character(len=80) :: suffix,fname
    character(len=10) :: id_str
    integer:: suffix_pos

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

    if (pe_size .gt. 1) then
       write(id_str,'(i10.10)') pe_rank
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

    call vtk_file_write_mesh(9, msh)

    if (associated(mfld)) then
       call vtk_file_write_cell_data(9, mfld)
    else if (associated(fld)) then
       call vtk_file_write_point_data(9, fld)
    end if
    
    close(9)
  end subroutine vtk_file_write

  subroutine vtk_file_read(this, data)
    class(vtk_file_t) :: this
    class(*), target, intent(inout) :: data
  end subroutine vtk_file_read

  !> Write a mesh in legacy VTK format
  subroutine vtk_file_write_mesh(unit, msh)
    integer :: unit
    type(mesh_t), intent(in) :: msh
    integer :: i, j, vtk_type

    ! Dump coordinates
    write(unit, fmt='(A,I8,A)') 'POINTS', msh%mpts,' double'
    do i = 1, msh%mpts
       write(unit, fmt='(F15.8,F15.8,F15.8)') msh%points(i)%x
    end do

    ! Dump cells
    write(unit, fmt='(A,I8,I8)')  'CELLS', msh%nelv, msh%nelv*(msh%npts+1)
    j = 0
    do i = 1, msh%nelv
       write(unit, *) msh%npts, &
            (msh%elements(i)%e%pts(j)%p%id() - 1, j=1, msh%npts)
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
    type(field_t), intent(in) :: fld
    type(point_t), target :: p1, p2, p3, p4
    real(kind=dp), allocatable :: point_data(:)
    integer :: i, j, lx, ly, lz, id(8)

    if ( (fld%Xh%lx - 1 .gt. 1) .or. &
         (fld%Xh%ly - 1 .gt. 1) .or. &
         (fld%Xh%lz - 1 .gt. 1)) then
       if (pe_rank .eq. 0) then
          call neko_warning("Interpolate high-order data onto a low-order mesh")
       end if
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
          id = fld%msh%elements(i)%e%pts(j)%p%id()
       end do

       point_data(id(1)) = fld%x(1,1,1,i,1)
       point_data(id(2)) = fld%x(lx,1,1,i,1)
       point_data(id(3)) = fld%x(1,ly,1,i,1)
       point_data(id(4)) = fld%x(lx,ly,1,i,1)       

       if (fld%msh%gdim .eq. 3) then
          point_data(id(5)) = fld%x(1,1,lz,i,1)
          point_data(id(6)) = fld%x(lx,1,lz,i,1)
          point_data(id(7)) = fld%x(1,ly,lz,i,1)
          point_data(id(8)) = fld%x(lx,ly,lz,i,1)       
       end if

    end do

    write(unit, *) point_data
    
    deallocate(point_data)       
    
  end subroutine vtk_file_write_point_data

end module vtk_file
