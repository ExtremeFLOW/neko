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
    real(kind=dp), intent(in), optional :: t
    type(mesh_t), pointer :: msh => null()
    type(field_t), pointer :: fld => null()
    type(mesh_fld_t), pointer :: mfld => null()
    type(dofmap_t), pointer :: dm => null()
    character(len=80) :: suffix,fname
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
    class default
       call neko_error('Invalid data')
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
       
    else
       call neko_error('Invalid data')
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
    type(mesh_t), intent(inout) :: msh
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
       write(unit, fmt='(I8,8I8)') msh%npts, &
            (mesh_get_local_point(msh, msh%elements(i)%e%pts(j)%p) - 1, &
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
          id(j) = mesh_get_local(fld%msh, fld%msh%elements(i)%e%pts(j)%p)
       end do

       point_data(id(1)) = fld%x(1,1,1,i)
       point_data(id(2)) = fld%x(lx,1,1,i)
       point_data(id(4)) = fld%x(1,ly,1,i)
       point_data(id(3)) = fld%x(lx,ly,1,i)       

       if (fld%msh%gdim .eq. 3) then
          point_data(id(5)) = fld%x(1,1,lz,i)
          point_data(id(6)) = fld%x(lx,1,lz,i)
          point_data(id(8)) = fld%x(1,ly,lz,i)
          point_data(id(7)) = fld%x(lx,ly,lz,i)       
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
                     dm%x(j,k,l,i), dm%y(j,k,l,i), dm%z(j,k,l,i)
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
                 write(unit, fmt='(I8)') dm%dof(j,k,l,i)
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

end module vtk_file
