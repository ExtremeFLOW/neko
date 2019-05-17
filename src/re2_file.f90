!> NEKTON mesh data in re2 format
!! @details This module is used to read/write binary NEKTION mesh data
module re2_file
  use generic_file
  use num_types
  use utils
  use mesh
  use point
  use mpi
  use mpi_types
  use re2
  implicit none
  private
  

  !> Interface for NEKTON re2 files
  type, public, extends(generic_file_t) :: re2_file_t
   contains
     procedure :: read => re2_file_read
     procedure :: write => re2_file_write
  end type re2_file_t

contains

  !> Load a binary NEKTON mesh from a re2 file
  subroutine re2_file_read(this, data)
    class(re2_file_t) :: this
    class(*), target, intent(inout) :: data
    type(re2_xy_t), allocatable :: re2_data_xy(:)
    type(re2_xyz_t), allocatable :: re2_data_xyz(:)
    type(mesh_t), pointer :: msh
    character(len=5) :: hdr_ver
    character(len=54) :: hdr_str
    integer :: i, j, k, fh, nel, ndim, nelv, ierr, el_idx, pt_idx
    integer :: rsize, csize
    integer :: status(MPI_STATUS_SIZE)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset
    real(kind=sp), allocatable :: xyz(:)
    real(kind=sp) :: test
    type(point_t) :: p(8)

    select type(data)
    type is (mesh_t)
       msh => data
    end select

    
    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    write(*, '(A,A)') " Reading binary NEKTON file ", this%fname
    read(9, '(a5,i9,i3,i9,a54)') hdr_ver, nel, ndim, nelv, hdr_str
    write(*,1) ndim, nelv
1   format(1x,'ndim = ', i1, ', nelements =', i7)
    close(9)


    call MPI_Type_size(MPI_CHARACTER, csize, ierr)
    call MPI_Type_size(MPI_REAL, rsize, ierr)

    call MPI_File_open(MPI_COMM_WORLD, trim(this%fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    
    if (ierr .ne. 0) then
       call neko_error("Can't open binary NEKTON file ")
    end if
    

    call mesh_init(msh, ndim, nelv)
   
    ! Set offset (header)
    mpi_offset = RE2_HDR_SIZE * csize

    call MPI_File_read_at(fh, mpi_offset, test, 1, MPI_REAL, status, ierr)
    mpi_offset = mpi_offset + rsize
    
    pt_idx = 1
    el_idx = 1
    if (ndim .eq. 2) then
       allocate(re2_data_xy(nelv))
       call MPI_File_read_at(fh, mpi_offset, &
            re2_data_xy, nelv, MPI_RE2_DATA_XY, status, ierr)
       do i = 1, nelv
          do j = 1, 8             
             p(j) = point_t(dble(re2_data_xy(i)%x(j)), &
                  dble(re2_data_xy(i)%y(j)), 0d0, pt_idx)
             pt_idx = pt_idx + 1
          end do

          call mesh_add_element(msh, el_idx, p(1), p(2), p(3), p(4))
          el_idx = el_idx + 1
       end do
       deallocate(re2_data_xy)
    else if (ndim .eq. 3) then
       allocate(re2_data_xyz(nelv))
       call MPI_File_read_at(fh, mpi_offset, &
            re2_data_xyz, nelv, MPI_RE2_DATA_XYZ, status, ierr)
       do i = 1, nelv
          do j = 1, 8             
             p(j) = point_t(dble(re2_data_xyz(i)%x(j)), &
                  dble(re2_data_xyz(i)%y(j)),&
                  dble(re2_data_xyz(i)%z(j)), pt_idx)
             pt_idx = pt_idx + 1
          end do

          call mesh_add_element(msh, el_idx, &
               p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))          
          el_idx = el_idx + 1
       end do
       deallocate(re2_data_xyz)
    end if
    write(*,*) 'Done'

    !> @todo Add support for curved side data

    close(9)

    
  end subroutine re2_file_read

  subroutine re2_file_write(this, data)
    class(re2_file_t), intent(in) :: this
    class(*), target, intent(in) :: data
  end subroutine re2_file_write

end module re2_file
