!> MPI dervied types
module mpi_types
  use comm
  use mpi
  use re2
  use nmsh
  use parameters
  implicit none
  private

  integer :: MPI_NMSH_HEX    !< MPI dervied type for 3D Neko nmsh data
  integer :: MPI_NMSH_QUAD   !< MPI dervied type for 2D Neko nmsh data

  integer :: MPI_RE2_DATA_XYZ !< MPI dervied type for 3D NEKTON re2 data
  integer :: MPI_RE2_DATA_XY  !< MPI dervied type for 2D NEKTON re2 data

  integer :: MPI_NEKO_PARAMS    !< MPI dervied type for parameters

  integer :: MPI_REAL_SIZE             !< Size of MPI type real
  integer :: MPI_DOUBLE_PRECISION_SIZE !< Size of MPI type double precision
  integer :: MPI_CHARACTER_SIZE        !< Size of MPI type character
  integer :: MPI_INTEGER_SIZE          !< size of MPI type integer

  ! Public dervied types and size definitions
  public :: MPI_NMSH_HEX, MPI_NMSH_QUAD, &
       MPI_RE2_DATA_XYZ, MPI_RE2_DATA_XY, &
       MPI_REAL_SIZE, MPI_DOUBLE_PRECISION_SIZE, &
       MPI_CHARACTER_SIZE, MPI_INTEGER_SIZE, &
       MPI_NEKO_PARAMS

  ! Public subroutines
  public :: mpi_types_init, mpi_types_free

contains

  !> Define all MPI dervied types
  subroutine mpi_types_init
    integer :: ierr

    ! Define derived types
    call mpi_type_nmsh_hex_init
    call mpi_type_nmsh_quad_init

    call mpi_type_re2_xyz_init
    call mpi_type_re2_xy_init

    call mpi_type_neko_params_init

    ! Check sizes of MPI types
    call MPI_Type_size(MPI_REAL, MPI_REAL_SIZE, ierr)
    call MPI_Type_size(MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION_SIZE, ierr)
    call MPI_Type_size(MPI_CHARACTER, MPI_CHARACTER_SIZE, ierr)
    call MPI_Type_size(MPI_INTEGER, MPI_INTEGER_SIZE, ierr)

    call MPI_Barrier(NEKO_COMM, ierr)

  end subroutine mpi_types_init

  !> Define a MPI dervied type for a 3d nmsh hex
  subroutine mpi_type_nmsh_hex_init
    type(nmsh_hex_t) :: nmsh_hex
    integer(kind=MPI_ADDRESS_KIND) :: disp(17), base
    integer :: type(17), len(17), i, ierr

    call MPI_Get_address(nmsh_hex%el_idx, disp(1), ierr)
    call MPI_Get_address(nmsh_hex%v(1)%v_idx, disp(2), ierr)
    call MPI_Get_address(nmsh_hex%v(1)%v_xyz, disp(3), ierr)
    call MPI_Get_address(nmsh_hex%v(2)%v_idx, disp(4), ierr)
    call MPI_Get_address(nmsh_hex%v(2)%v_xyz, disp(5), ierr)
    call MPI_Get_address(nmsh_hex%v(3)%v_idx, disp(6), ierr)
    call MPI_Get_address(nmsh_hex%v(3)%v_xyz, disp(7), ierr)
    call MPI_Get_address(nmsh_hex%v(4)%v_idx, disp(8), ierr)
    call MPI_Get_address(nmsh_hex%v(4)%v_xyz, disp(9), ierr)
    call MPI_Get_address(nmsh_hex%v(5)%v_idx, disp(10), ierr)
    call MPI_Get_address(nmsh_hex%v(5)%v_xyz, disp(11), ierr)
    call MPI_Get_address(nmsh_hex%v(6)%v_idx, disp(12), ierr)
    call MPI_Get_address(nmsh_hex%v(6)%v_xyz, disp(13), ierr)
    call MPI_Get_address(nmsh_hex%v(7)%v_idx, disp(14), ierr)
    call MPI_Get_address(nmsh_hex%v(7)%v_xyz, disp(15), ierr)
    call MPI_Get_address(nmsh_hex%v(8)%v_idx, disp(16), ierr)
    call MPI_Get_address(nmsh_hex%v(8)%v_xyz, disp(17), ierr)
    

    base = disp(1)
    do i = 1, 17
       disp(i) = disp(i) - base
    end do

    len(1) = 1
    len(2:16:2) = 1
    len(3:17:2) = 3

    type(1) = MPI_INTEGER
    type(2:16:2) = MPI_INTEGER
    type(3:17:2) = MPI_DOUBLE_PRECISION
    call MPI_Type_create_struct(17, len, disp, type, MPI_NMSH_HEX, ierr)
    call MPI_Type_commit(MPI_NMSH_HEX, ierr)    
  end subroutine mpi_type_nmsh_hex_init

    !> Define a MPI dervied type for a 2d nmsh quad
  subroutine mpi_type_nmsh_quad_init
    type(nmsh_hex_t) :: nmsh_quad
    integer(kind=MPI_ADDRESS_KIND) :: disp(9), base
    integer :: type(9), len(9), i, ierr

    call MPI_Get_address(nmsh_quad%el_idx, disp(1), ierr)
    call MPI_Get_address(nmsh_quad%v(1)%v_idx, disp(2), ierr)
    call MPI_Get_address(nmsh_quad%v(1)%v_xyz, disp(3), ierr)
    call MPI_Get_address(nmsh_quad%v(2)%v_idx, disp(4), ierr)
    call MPI_Get_address(nmsh_quad%v(2)%v_xyz, disp(5), ierr)
    call MPI_Get_address(nmsh_quad%v(3)%v_idx, disp(6), ierr)
    call MPI_Get_address(nmsh_quad%v(3)%v_xyz, disp(7), ierr)
    call MPI_Get_address(nmsh_quad%v(4)%v_idx, disp(8), ierr)
    call MPI_Get_address(nmsh_quad%v(4)%v_xyz, disp(9), ierr)
    

    base = disp(1)
    do i = 1, 9
       disp(i) = disp(i) - base
    end do

    len(1) = 1
    len(2:8:2) = 1
    len(3:9:2) = 3

    type(1) = MPI_INTEGER
    type(2:8:2) = MPI_INTEGER
    type(3:9:2) = MPI_DOUBLE_PRECISION
    call MPI_Type_create_struct(9, len, disp, type, MPI_NMSH_QUAD, ierr)
    call MPI_Type_commit(MPI_NMSH_QUAD, ierr)    
  end subroutine mpi_type_nmsh_quad_init
  
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

    len = 4
    len(1) = 1
    type = MPI_REAL

    call MPI_Type_create_struct(3, len, disp, type, MPI_RE2_DATA_XY, ierr)
    call MPI_Type_commit(MPI_RE2_DATA_XY, ierr)

  end subroutine mpi_type_re2_xy_init

  !> Define a MPI derived type for parameters
  subroutine mpi_type_neko_params_init
    type(param_t) :: param_data
    integer(kind=MPI_ADDRESS_KIND) :: disp(7), base    
    integer :: type(7), len(7), ierr

    call MPI_Get_address(param_data%dt, disp(1), ierr)
    call MPI_Get_address(param_data%nsteps, disp(2), ierr)
    call MPI_Get_address(param_data%rho, disp(3), ierr)
    call MPI_Get_address(param_data%mu, disp(4), ierr)
    call MPI_Get_address(param_data%Re, disp(5), ierr)
    call MPI_Get_address(param_data%uinf, disp(6), ierr)
    call MPI_Get_address(param_data%output_bdry, disp(7), ierr)

    base = disp(1)
    disp(1) = disp(1) - base
    disp(2) = disp(2) - base
    disp(3) = disp(3) - base
    disp(4) = disp(4) - base
    disp(5) = disp(5) - base
    disp(6) = disp(6) - base
    disp(7) = disp(7) - base

    len(1:4) = 1
    len(6) = 3
    len(7) = 1
    
    type(1) = MPI_DOUBLE_PRECISION
    type(2) = MPI_INTEGER
    type(3) = MPI_DOUBLE_PRECISION
    type(4) = MPI_DOUBLE_PRECISION
    type(5) = MPI_DOUBLE_PRECISION
    type(6) = MPI_DOUBLE_PRECISION
    type(7) = MPI_LOGICAL
    
    call MPI_Type_create_struct(6, len, disp, type, MPI_NEKO_PARAMS, ierr)
    call MPI_Type_commit(MPI_NEKO_PARAMS, ierr)

  end subroutine mpi_type_neko_params_init

  !> Deallocate all dervied MPI types
  subroutine mpi_types_free
    call mpi_type_nmsh_hex_free
    call mpi_Type_nmsh_quad_free
    call mpi_type_re2_xyz_free
    call mpi_type_re2_xy_free
    call mpi_type_neko_params_free
  end subroutine mpi_types_free

  !> Deallocate nmsh hex derived MPI type
  subroutine mpi_type_nmsh_hex_free
    integer ierr
    call MPI_Type_free(MPI_NMSH_HEX, ierr)
  end subroutine mpi_type_nmsh_hex_free

    !> Deallocate nmsh quad derived MPI type
  subroutine mpi_type_nmsh_quad_free
    integer ierr
    call MPI_Type_free(MPI_NMSH_QUAD, ierr)
  end subroutine mpi_type_nmsh_quad_free
  
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

  !> Deallocate parameters dervied MPI type
  subroutine mpi_type_neko_params_free
    integer ierr
    call MPI_Type_free(MPI_NEKO_PARAMS, ierr)
  end subroutine mpi_type_neko_params_free
  
end module mpi_types
