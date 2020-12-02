!> NEKTON fld file format
!! @details this module defines interface to write NEKTON's fld fields
module fld_file
  use generic_file
  use field
  use dofmap
  use mesh
  use utils
  use comm
  use mpi_types
  implicit none
  private


  !> Interface for NEKTON fld files
  type, public, extends(generic_file_t) :: fld_file_t
   contains
     procedure :: read => fld_file_read
     procedure :: write => fld_file_write
  end type fld_file_t

contains

  !> Write a field from a NEKTON fld file
  !! @note currently limited to one scalar field (field_t) containing
  !! double precision data
  subroutine fld_file_write(this, data)
    class(fld_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    type(field_t), pointer :: f
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    type(dofmap_t), pointer :: dof
    integer :: status(MPI_STATUS_SIZE)
    character(len=132) :: hdr
    character :: rdcode(10)
    character(len=6) :: id_str
    character(len=80) :: fname
    integer :: i, ierr, fh, n, j,k,l,el, suffix_pos
    integer, allocatable :: idx(:)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    real(kind=dp), allocatable :: tmp(:)
    real(kind=sp), parameter :: test_pattern = 6.54321
    logical :: write_mesh
   
    select type(data)
    type is (field_t)
       f => data
       msh => f%msh
       Xh => f%Xh
       dof => f%dof
    class default
       call neko_error('Invalid data')
    end select

    
    !
    ! Create fld header for NEKTON's multifile output
    !

    write_mesh = (this%counter .eq. 0)
    
    ! Build rdcode note that for field_t, we only support scalar
    ! fields at the moment
    rdcode = ' '
    i = 1
    if (write_mesh) then
       rdcode(i) = 'X'
       i = i + 1
    end if
    rdcode(i) = 'P'


    !> @todo fix support for single precision output?
    write(hdr, 1) 8,Xh%lx, Xh%ly, Xh%lz,msh%glb_nelv,msh%glb_nelv,&
         0d0,1,1,1,(rdcode(i),i=1,10)
1   format('#std',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,&
         1x,i9,1x,i6,1x,i6,1x,10a)

    ! Store global idx of each element
    allocate(idx(msh%nelv))
    do i = 1, msh%nelv
       idx(i) = msh%elements(i)%e%id()
    end do

    ! Change to NEKTON's fld file format
    suffix_pos = filename_suffix_pos(this%fname)
    write(id_str, '(a,i5.5)') 'f', this%counter
    fname = trim(this%fname(1:suffix_pos-1))//'0.'//id_str
    
    call MPI_File_open(NEKO_COMM, trim(fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

    call MPI_File_write_all(fh, hdr, 132, MPI_CHARACTER,  status, ierr)
    mpi_offset = 132 * MPI_CHARACTER_SIZE

    call MPI_File_write_all(fh, test_pattern, 1, MPI_REAL, status, ierr)
    mpi_offset = mpi_offset  + MPI_REAL_SIZE

    byte_offset = mpi_offset + &
         int(msh%offset_el, 8) * int(MPI_INTEGER_SIZE, 8)
    call MPI_File_write_at_all(fh, byte_offset, idx, msh%nelv, &
         MPI_INTEGER, status, ierr)
    mpi_offset = mpi_offset + int(msh%glb_nelv, 8) * int(MPI_INTEGER_SIZE, 8)

    deallocate(idx)
    
    n = 3*(Xh%lx * Xh%ly * Xh%lz * msh%nelv)
    allocate(tmp(n))
       
    if (write_mesh) then
       i = 1
       do el = 1, msh%nelv
          do l = 1, Xh%lz
             do k = 1, Xh%ly
                do j = 1, Xh%lx
                   tmp(i) = dof%x(j,k,l,el)
                   i = i +1
                end do
             end do
          end do
          do l = 1, Xh%lz
             do k = 1, Xh%ly
                do j = 1, Xh%lx
                   tmp(i) = dof%y(j,k,l,el)
                   i = i +1
                end do
             end do
          end do
          do l = 1, Xh%lz
             do k = 1, Xh%ly
                do j = 1, Xh%lx
                   tmp(i) = dof%z(j,k,l,el)
                   i = i +1
                end do
             end do
          end do
       end do
       
       byte_offset = mpi_offset + int(msh%offset_el, 8) * &
            (int(3 * (Xh%lx * Xh%ly * Xh%lz), 8) * &
            int(MPI_DOUBLE_PRECISION_SIZE, 8))
       call MPI_File_write_at_all(fh, byte_offset, tmp, n, &
            MPI_DOUBLE_PRECISION, status, ierr)
       mpi_offset = mpi_offset + int(msh%glb_nelv, 8) * &
            (int(3 * (Xh%lx * Xh%ly * Xh%lz), 8) * &
            int(MPI_DOUBLE_PRECISION_SIZE, 8))
    end if
 
    i = 1
    do el = 1, msh%nelv
       do l = 1, Xh%lz
          do k = 1, Xh%ly
             do j = 1, Xh%lx
                tmp(i) = real(f%x(j,k,l,el))
                i = i + 1
             end do
          end do
       end do
    end do

    byte_offset = mpi_offset + int(msh%offset_el, 8) * &
         (int((Xh%lx * Xh%ly * Xh%lz), 8) * &
         int(MPI_DOUBLE_PRECISION_SIZE, 8))
    call MPI_File_write_at_all(fh, byte_offset, f%x, n/3, &
         MPI_DOUBLE_PRECISION, status, ierr)

    deallocate(tmp)
    
    call MPI_File_close(fh, ierr)

    ! Write metadata file 
    if (pe_rank .eq. 0) then
       open(unit=9, file=trim(this%fname(1:suffix_pos-1))//'.nek5000', &
            status='replace')
       write(9, fmt='(A,A,A)') 'filetemplate:         ', &
            this%fname(1:suffix_pos-1),'%01d.f%05d'
       write(9, fmt='(A)') 'firsttimestep: 0'
       write(9, fmt='(A,i5)') 'numtimesteps: ', this%counter + 1
       write(9, fmt='(A)') 'type: binary'
       close(9)
    end if

    this%counter = this%counter + 1
    
  end subroutine fld_file_write
  
  !> Load a field from a NEKTON fld file
  subroutine fld_file_read(this, data)
    class(fld_file_t) :: this
    class(*), target, intent(inout) :: data

    call neko_error('Not implemented yet!')
    
  end subroutine fld_file_read


  
end module fld_file
