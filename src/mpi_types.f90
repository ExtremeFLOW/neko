!> MPI dervied types
module mpi_types
  use mpi
  use re2
  implicit none
  private

  integer :: MPI_RE2_DATA_XYZ !< mpi dervied type for 3D NEKTON re2 data
  integer :: MPI_RE2_DATA_XY !< MPI dervied type for 2D NEKTON re2 data

  ! Public dervied typese
  public :: MPI_RE2_DATA_XYZ, MPI_RE2_DATA_XY

  ! Public subroutines
  public :: mpi_types_init, mpi_types_free

contains

  !> Define all MPI dervied types
  subroutine mpi_types_init
    call mpi_type_re2_xyz_init
  end subroutine mpi_types_init

  !> Define a MPI derived type for a 3d re2 data
  subroutine mpi_type_re2_xyz_init
    type(re2_xyz_t) :: re2_data
    integer(kind=MPI_ADDRESS_KIND) :: disp(4), base    
    integer :: type(4), len(4), ierr

    call MPI_Get_address(re2_data%rgroup, disp(1), ierr)
    call MPI_Get_address(re2_data%x, disp(2), ierr)
    call MPI_Get_address(re2_data%y, disp(3), ierr)
    call MPI_Get_address(re2_data%z, disp(4), ierr)

    base = disp(1)
    disp(1) = disp(1) - base
    disp(2) = disp(2) - base
    disp(3) = disp(3) - base
    disp(4) = disp(4) - base

    len = 8
    len(1) = 1
    type = MPI_REAL

    call MPI_Type_create_struct(4, len, disp, type, MPI_RE2_DATA_XYZ, ierr)
    call MPI_Type_commit(MPI_RE2_DATA_XYZ, ierr)

  end subroutine mpi_type_re2_xyz_init

    !> Define a MPI derived type for a 2d re2 data
  subroutine mpi_type_re2_xy_init
    type(re2_xy_t) :: re2_data
    integer(kind=MPI_ADDRESS_KIND) :: disp(3), base    
    integer :: type(3), len(3), ierr

    call MPI_Get_address(re2_data%rgroup, disp(1), ierr)
    call MPI_Get_address(re2_data%x, disp(2), ierr)
    call MPI_Get_address(re2_data%y, disp(3), ierr)

    base = disp(1)
    disp(1) = disp(1) - base
    disp(2) = disp(2) - base
    disp(3) = disp(3) - base

    len = 8
    len(1) = 1
    type = MPI_REAL

    call MPI_Type_create_struct(3, len, disp, type, MPI_RE2_DATA_XY, ierr)
    call MPI_Type_commit(MPI_RE2_DATA_XY, ierr)

  end subroutine mpi_type_re2_xy_init

  !> Deallocate all dervied MPI types
  subroutine mpi_types_free
    call mpi_type_re2_xyz_free
    call mpi_type_re2_xy_free
  end subroutine mpi_types_free

  !> Deallocate re2 xyz dervied MPI type
  subroutine mpi_type_re2_xyz_free
    integer ierr
    call MPI_Type_free(MPI_RE2_DATA_XYZ, ierr)
  end subroutine mpi_type_re2_xyz_free

    !> Deallocate re2 xyz dervied MPI type
  subroutine mpi_type_re2_xy_free
    integer ierr
    call MPI_Type_free(MPI_RE2_DATA_XY, ierr)
  end subroutine mpi_type_re2_xy_free
  
end module mpi_types
