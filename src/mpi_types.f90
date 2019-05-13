!> MPI dervied types
module mpi_types
  use mpi
  use re2
  implicit none
  private

  integer :: MPI_RE2_HDR            !< MPI dervied type for a NEKTON re2 header

  ! Public dervied types
  public :: MPI_RE2_HDR

  ! Public subroutines
  public :: mpi_types_init, mpi_types_free

contains

  !> Define all MPI dervied types
  subroutine mpi_types_init
    call mpi_type_re2_init
  end subroutine mpi_types_init

  !> Define a MPI derived type for a re2 hdr
  subroutine mpi_type_re2_init
    type(re2_hdr_t) :: re2_hdr
    integer(kind=MPI_ADDRESS_KIND) :: disp(6), base    
    integer :: type(6), len(6), ierr

    call MPI_Get_address(re2_hdr%hdr_ver, disp(1), ierr)
    call MPI_Get_address(re2_hdr%nel, disp(2), ierr)
    call MPI_Get_address(re2_hdr%ndim, disp(3), ierr)
    call MPI_Get_address(re2_hdr%nelv, disp(4), ierr)
    call MPI_Get_address(re2_hdr%hdr_str, disp(5), ierr)
    call MPI_Get_address(re2_hdr%endian_test, disp(6), ierr)

    base = disp(1)
    disp(1) = disp(1) - base
    disp(2) = disp(2) - base
    disp(3) = disp(3) - base
    disp(4) = disp(4) - base
    disp(5) = disp(5) - base
    disp(6) = disp(6) - base

    len = 1
    len(1) = 5
    len(5) = 54

    type(1) = MPI_CHARACTER
    type(2:4) = MPI_INTEGER
    type(5) = MPI_CHARACTER
    type(6) = MPI_REAL

    call MPI_Type_create_struct(6, len, disp, type, MPI_RE2_HDR, ierr)
    call MPI_Type_commit(MPI_RE2_HDR, ierr)

  end subroutine mpi_type_re2_init

  !> Deallocate all dervied MPI types
  subroutine mpi_types_free
    call mpi_type_re2_free
  end subroutine mpi_types_free

  !> Deallocate re2 hdr dervied MPI type
  subroutine mpi_type_re2_free
    integer ierr
    call MPI_Type_free(MPI_RE2_HDR, ierr)
  end subroutine mpi_type_re2_free
  
end module mpi_types
