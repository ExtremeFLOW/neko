! Copyright (c) 2019-2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> MPI derived types
module mpi_types
  use comm
  use re2
  use nmsh
  use stl
  implicit none
  private

  type(MPI_Datatype) :: MPI_NMSH_HEX    !< MPI derived type for 3D Neko nmsh data
  type(MPI_Datatype) :: MPI_NMSH_QUAD   !< MPI derived type for 2D Neko nmsh data
  type(MPI_Datatype) :: MPI_NMSH_ZONE   !< MPI derived type for Neko nmsh zone data
  type(MPI_Datatype) :: MPI_NMSH_CURVE   !< MPI derived type for Neko nmsh curved elements

  type(MPI_Datatype) :: MPI_RE2V1_DATA_XYZ !< MPI derived type for 3D NEKTON re2 data
  type(MPI_Datatype) :: MPI_RE2V1_DATA_XY  !< MPI derived type for 2D NEKTON re2 data
  type(MPI_Datatype) :: MPI_RE2V1_DATA_CV  !< MPI derived type for NEKTON re2 cv data
  type(MPI_Datatype) :: MPI_RE2V1_DATA_BC  !< MPI derived type for NEKTON re2 bc data

  type(MPI_Datatype) :: MPI_RE2V2_DATA_XYZ !< MPI derived type for 3D NEKTON re2 data
  type(MPI_Datatype) :: MPI_RE2V2_DATA_XY  !< MPI derived type for 2D NEKTON re2 data
  type(MPI_Datatype) :: MPI_RE2V2_DATA_CV  !< MPI derived type for NEKTON re2 cv data
  type(MPI_Datatype) :: MPI_RE2V2_DATA_BC  !< MPI derived type for NEKTON re2 bc data

  type(MPI_Datatype) :: MPI_STL_HEADER     !< MPI Derived type for a STL header
  type(MPI_Datatype) :: MPI_STL_TRIANGLE   !< MPI derived type for a STL triangle

  integer :: MPI_REAL_SIZE             !< Size of MPI type real
  integer :: MPI_DOUBLE_PRECISION_SIZE !< Size of MPI type double precision
  integer :: MPI_CHARACTER_SIZE        !< Size of MPI type character
  integer :: MPI_INTEGER_SIZE          !< Size of MPI type integer
  integer :: MPI_LOGICAL_SIZE          !< Size of MPI type logical
  integer :: MPI_REAL_PREC_SIZE        !< Size of working precision REAL types

  ! Public derived types and size definitions
  public :: MPI_NMSH_HEX, MPI_NMSH_QUAD, MPI_NMSH_ZONE, &
       MPI_NMSH_CURVE, &
       MPI_RE2V1_DATA_XYZ, MPI_RE2V1_DATA_XY, &
       MPI_RE2V1_DATA_CV, MPI_RE2V1_DATA_BC, &
       MPI_RE2V2_DATA_XYZ, MPI_RE2V2_DATA_XY, &
       MPI_RE2V2_DATA_CV, MPI_RE2V2_DATA_BC, &
       MPI_REAL_SIZE, MPI_DOUBLE_PRECISION_SIZE, &
       MPI_CHARACTER_SIZE, MPI_INTEGER_SIZE, &
       MPI_LOGICAL_SIZE, MPI_REAL_PREC_SIZE, &
       MPI_STL_HEADER, MPI_STL_TRIANGLE

  ! Public subroutines
  public :: mpi_types_init, mpi_types_free

contains

  !> Define all MPI derived types
  subroutine mpi_types_init
    integer :: ierr

    ! Define derived types
    call mpi_type_nmsh_hex_init
    call mpi_type_nmsh_quad_init
    call mpi_type_nmsh_zone_init
    call mpi_type_nmsh_curve_init

    call mpi_type_re2_xyz_init
    call mpi_type_re2_xy_init
    call mpi_type_re2_cv_init
    call mpi_type_re2_bc_init

    call mpi_type_stl_header_init
    call mpi_type_stl_triangle_init

    ! Check sizes of MPI types
    call MPI_Type_size(MPI_REAL, MPI_REAL_SIZE, ierr)
    call MPI_Type_size(MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION_SIZE, ierr)
    call MPI_Type_size(MPI_CHARACTER, MPI_CHARACTER_SIZE, ierr)
    call MPI_Type_size(MPI_INTEGER, MPI_INTEGER_SIZE, ierr)
    call MPI_Type_size(MPI_LOGICAL, MPI_LOGICAL_SIZE, ierr)
    call MPI_Type_size(MPI_REAL_PRECISION, MPI_REAL_PREC_SIZE, ierr)

    call MPI_Barrier(NEKO_COMM, ierr)

  end subroutine mpi_types_init

  !> Define a MPI derived type for a 3d nmsh hex
  subroutine mpi_type_nmsh_hex_init
    type(nmsh_hex_t) :: nmsh_hex
    type(MPI_Datatype) :: type(17)
    integer(kind=MPI_ADDRESS_KIND) :: disp(17), base
    integer :: len(17), i, ierr

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
       disp(i) = MPI_Aint_diff(disp(i), base)
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

  !> Define a MPI derived type for a 2d nmsh quad
  subroutine mpi_type_nmsh_quad_init
    type(nmsh_quad_t) :: nmsh_quad
    type(MPI_Datatype) :: type(9)
    integer(kind=MPI_ADDRESS_KIND) :: disp(9), base
    integer :: len(9), i, ierr

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
       disp(i) = MPI_Aint_diff(disp(i), base)
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

  !> Define a MPI derived type for a nmsh zone
  subroutine mpi_type_nmsh_zone_init
    type(nmsh_zone_t) :: nmsh_zone
    type(MPI_Datatype) :: type(6)
    integer(kind=MPI_ADDRESS_KIND) :: disp(6), base
    integer :: len(6), i, ierr

    call MPI_Get_address(nmsh_zone%e, disp(1), ierr)
    call MPI_Get_address(nmsh_zone%f, disp(2), ierr)
    call MPI_Get_address(nmsh_zone%p_e, disp(3), ierr)
    call MPI_Get_address(nmsh_zone%p_f, disp(4), ierr)
    call MPI_Get_address(nmsh_zone%glb_pt_ids, disp(5), ierr)
    call MPI_Get_address(nmsh_zone%type, disp(6), ierr)

    base = disp(1)
    do i = 1, 6
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1:4) = 1
    len(5) = 4
    len(6) = 1
    type = MPI_INTEGER

    call MPI_Type_create_struct(6, len, disp, type, MPI_NMSH_ZONE, ierr)
    call MPI_Type_commit(MPI_NMSH_ZONE, ierr)
    
  end subroutine mpi_type_nmsh_zone_init
 
  !> Define a MPI derived type for a nmsh curved element
  subroutine mpi_type_nmsh_curve_init
    type(nmsh_curve_el_t) :: nmsh_curve_el
    type(MPI_Datatype) :: type(3)
    integer(kind=MPI_ADDRESS_KIND) :: disp(3), base
    integer :: len(3), i, ierr

    call MPI_Get_address(nmsh_curve_el%e, disp(1), ierr)
    call MPI_Get_address(nmsh_curve_el%curve_data, disp(2), ierr)
    call MPI_Get_address(nmsh_curve_el%type, disp(3), ierr)

    base = disp(1)
    do i = 1, 3
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1) = 1
    len(2) = 5*12
    len(3) = 12
    type(1) = MPI_INTEGER
    type(2) = MPI_DOUBLE_PRECISION
    type(3) = MPI_INTEGER

    call MPI_Type_create_struct(3, len, disp, type, MPI_NMSH_CURVE, ierr)
    call MPI_Type_commit(MPI_NMSH_CURVE, ierr)
    
  end subroutine mpi_type_nmsh_curve_init
  
 
  !> Define a MPI derived type for a 3d re2 data
  subroutine mpi_type_re2_xyz_init
    type(re2v1_xyz_t) :: re2v1_data
    type(re2v2_xyz_t) :: re2v2_data
    type(MPI_Datatype) :: type(4)
    integer(kind=MPI_ADDRESS_KIND) :: disp(4), base    
    integer :: len(4), ierr, i

    !
    ! Setup version 1
    !
    
    call MPI_Get_address(re2v1_data%rgroup, disp(1), ierr)
    call MPI_Get_address(re2v1_data%x, disp(2), ierr)
    call MPI_Get_address(re2v1_data%y, disp(3), ierr)
    call MPI_Get_address(re2v1_data%z, disp(4), ierr)

    base = disp(1)
    do i = 1, 4
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len = 8
    len(1) = 1
    type = MPI_REAL

    call MPI_Type_create_struct(4, len, disp, type, MPI_RE2V1_DATA_XYZ, ierr)
    call MPI_Type_commit(MPI_RE2V1_DATA_XYZ, ierr)

    !
    ! Setup version 2
    !
    
    call MPI_Get_address(re2v2_data%rgroup, disp(1), ierr)
    call MPI_Get_address(re2v2_data%x, disp(2), ierr)
    call MPI_Get_address(re2v2_data%y, disp(3), ierr)
    call MPI_Get_address(re2v2_data%z, disp(4), ierr)

    base = disp(1)
    do i = 1, 4
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len = 8
    len(1) = 1
    type = MPI_DOUBLE_PRECISION

    call MPI_Type_create_struct(4, len, disp, type, MPI_RE2V2_DATA_XYZ, ierr)
    call MPI_Type_commit(MPI_RE2V2_DATA_XYZ, ierr)

  end subroutine mpi_type_re2_xyz_init

  !> Define a MPI derived type for a 2d re2 data
  subroutine mpi_type_re2_xy_init
    type(re2v1_xy_t) :: re2v1_data
    type(re2v2_xy_t) :: re2v2_data
    type(MPI_Datatype) :: type(3)
    integer(kind=MPI_ADDRESS_KIND) :: disp(3), base    
    integer :: len(3), ierr, i

    !
    ! Setup version 1
    !

    call MPI_Get_address(re2v1_data%rgroup, disp(1), ierr)
    call MPI_Get_address(re2v1_data%x, disp(2), ierr)
    call MPI_Get_address(re2v1_data%y, disp(3), ierr)

    base = disp(1)
    do i = 1, 3
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len = 4
    len(1) = 1
    type = MPI_REAL

    call MPI_Type_create_struct(3, len, disp, type, MPI_RE2V1_DATA_XY, ierr)
    call MPI_Type_commit(MPI_RE2V1_DATA_XY, ierr)
    
    !
    ! Setup version 2
    !

    call MPI_Get_address(re2v2_data%rgroup, disp(1), ierr)
    call MPI_Get_address(re2v2_data%x, disp(2), ierr)
    call MPI_Get_address(re2v2_data%y, disp(3), ierr)

    base = disp(1)
    do i = 1, 3
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len = 4
    len(1) = 1
    type = MPI_DOUBLE_PRECISION

    call MPI_Type_create_struct(3, len, disp, type, MPI_RE2V2_DATA_XY, ierr)
    call MPI_Type_commit(MPI_RE2V2_DATA_XY, ierr)

  end subroutine mpi_type_re2_xy_init

  !> Define a MPI derived type for re2 cv data
  subroutine mpi_type_re2_cv_init
    type(re2v1_curve_t) :: re2v1_data
    type(re2v2_curve_t) :: re2v2_data
    type(MPI_Datatype) :: type(4)
    integer(kind=MPI_ADDRESS_KIND) :: disp(4), base
    integer :: len(4), ierr, i

    !
    ! Setup version 1
    !
    
    call MPI_Get_address(re2v1_data%elem, disp(1), ierr)
    call MPI_Get_address(re2v1_data%zone, disp(2), ierr)
    call MPI_Get_address(re2v1_data%point, disp(3), ierr)
    call MPI_Get_address(re2v1_data%type, disp(4), ierr)

    base = disp(1)
    do i = 1, 4
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1:2) = 1
    len(3) = 5
    len(4) = 4
    type(1:2) = MPI_INTEGER
    type(3) = MPI_REAL
    type(4) = MPI_CHARACTER

    call MPI_Type_create_struct(4, len, disp, type, MPI_RE2V1_DATA_CV, ierr)
    call MPI_Type_commit(MPI_RE2V1_DATA_CV, ierr)

    !
    ! Setup version 2
    !
    
    call MPI_Get_address(re2v2_data%elem, disp(1), ierr)
    call MPI_Get_address(re2v2_data%zone, disp(2), ierr)
    call MPI_Get_address(re2v2_data%point, disp(3), ierr)
    call MPI_Get_address(re2v2_data%type, disp(4), ierr)

    base = disp(1)
    do i = 1, 4
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1:2) = 1
    len(3) = 5
    len(4) = 8
    type(1:2) = MPI_DOUBLE_PRECISION
    type(3) = MPI_DOUBLE_PRECISION
    type(4) = MPI_CHARACTER

    call MPI_Type_create_struct(4, len, disp, type, MPI_RE2V2_DATA_CV, ierr)
    call MPI_Type_commit(MPI_RE2V2_DATA_CV, ierr)
    
  end subroutine mpi_type_re2_cv_init
  
  !> Define a MPI derived type for re2 bc data
  subroutine mpi_type_re2_bc_init
    type(re2v1_bc_t) :: re2v1_data
    type(re2v2_bc_t) :: re2v2_data
    type(MPI_Datatype) :: type(4)
    integer(kind=MPI_ADDRESS_KIND) :: disp(4), base
    integer :: len(4), ierr, i

    !
    ! Setup version 1
    !

    call MPI_Get_address(re2v1_data%elem, disp(1), ierr)
    call MPI_Get_address(re2v1_data%face, disp(2), ierr)
    call MPI_Get_address(re2v1_data%bc_data, disp(3), ierr)
    call MPI_Get_address(re2v1_data%type, disp(4), ierr)

    base = disp(1)
    do i = 1, 4
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1:2) = 1
    len(3) = 5
    len(4) = 4
    type(1:2) = MPI_INTEGER
    type(3) = MPI_REAL
    type(4) = MPI_CHARACTER

    call MPI_Type_create_struct(4, len, disp, type, MPI_RE2V1_DATA_BC, ierr)
    call MPI_Type_commit(MPI_RE2V1_DATA_BC, ierr)

    !
    ! Setup version 2
    !

    call MPI_Get_address(re2v2_data%elem, disp(1), ierr)
    call MPI_Get_address(re2v2_data%face, disp(2), ierr)
    call MPI_Get_address(re2v2_data%bc_data, disp(3), ierr)
    call MPI_Get_address(re2v2_data%type, disp(4), ierr)

    base = disp(1)
    do i = 1, 4
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1:2) = 1
    len(3) = 5
    len(4) = 8
    type(1:2) = MPI_DOUBLE_PRECISION
    type(3) = MPI_DOUBLE_PRECISIOn
    type(4) = MPI_CHARACTER

    call MPI_Type_create_struct(4, len, disp, type, MPI_RE2V2_DATA_BC, ierr)
    call MPI_Type_commit(MPI_RE2V2_DATA_BC, ierr)
    
  end subroutine mpi_type_re2_bc_init

  !> Define a MPI dervied type for a STL header
  subroutine mpi_type_stl_header_init
    type(stl_hdr_t) :: stl_hdr
    type(MPI_Datatype) :: type(2)
    integer(kind=MPI_ADDRESS_KIND) :: disp(2), base
    integer :: len(2), ierr, i

    call MPI_Get_address(stl_hdr%hdr, disp(1), ierr)
    call MPI_Get_address(stl_hdr%ntri, disp(2), ierr)

    base = disp(1)
    do i = 1, 2
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1) = 80
    len(2) = 1

    type(1) = MPI_CHARACTER
    type(2) = MPI_INTEGER

    call MPI_Type_create_struct(2, len, disp, type, MPI_STL_HEADER, ierr)
    call MPI_Type_commit(MPI_STL_HEADER, ierr)
      
  end subroutine mpi_type_stl_header_init

  !> Define a MPI derived type for a STL triangle
  subroutine mpi_type_stl_triangle_init
    type(stl_triangle_t) :: tri
    type(MPI_Datatype) :: type(5)
    integer(kind=MPI_ADDRESS_KIND) :: disp(5), base
    integer :: len(5), i, ierr

    call MPI_Get_address(tri%n, disp(1), ierr)
    call MPI_Get_address(tri%v1, disp(2), ierr)
    call MPI_Get_address(tri%v2, disp(3), ierr)
    call MPI_Get_address(tri%v3, disp(4), ierr)
    call MPI_Get_address(tri%attrib, disp(5), ierr)

    base = disp(1)
    do i = 1, 5
       disp(i) = MPI_Aint_diff(disp(i), base)
    end do

    len(1:4) = 3
    len(5) = 1

    type(1:4) = MPI_REAL
    type(5) = MPI_INTEGER2

    call MPI_Type_create_struct(5, len, disp, type, MPI_STL_TRIANGLE, ierr)
    call MPI_Type_commit(MPI_STL_TRIANGLE, ierr)
        
  end subroutine mpi_type_stl_triangle_init

  !> Deallocate all derived MPI types
  subroutine mpi_types_free
    call mpi_type_nmsh_hex_free
    call mpi_Type_nmsh_quad_free
    call mpi_Type_nmsh_zone_free
    call mpi_Type_nmsh_curve_free
    call mpi_type_re2_xyz_free
    call mpi_type_re2_xy_free
    call mpi_type_re2_bc_free
    call mpi_type_stl_header_free
    call mpi_type_stl_triangle_free
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

  !> Deallocate nmsh zone derived MPI type
  subroutine mpi_type_nmsh_zone_free
    integer ierr
    call MPI_Type_free(MPI_NMSH_ZONE, ierr)
  end subroutine mpi_type_nmsh_zone_free
  
  !> Deallocate nmsh curve derived MPI type
  subroutine mpi_type_nmsh_curve_free
    integer ierr
    call MPI_Type_free(MPI_NMSH_CURVE, ierr)
  end subroutine mpi_type_nmsh_curve_free
  
  !> Deallocate re2 xyz derived MPI type
  subroutine mpi_type_re2_xyz_free
    integer ierr
    call MPI_Type_free(MPI_RE2V1_DATA_XYZ, ierr)
    call MPI_Type_free(MPI_RE2V2_DATA_XYZ, ierr)
  end subroutine mpi_type_re2_xyz_free

  !> Deallocate re2 xyz derived MPI type
  subroutine mpi_type_re2_xy_free
    integer ierr
    call MPI_Type_free(MPI_RE2V1_DATA_XY, ierr)
    call MPI_Type_free(MPI_RE2V2_DATA_XY, ierr)
  end subroutine mpi_type_re2_xy_free

  !> Deallocate re2 cv derived MPI type
  subroutine mpi_type_re2_cv_free
    integer ierr
    call MPI_Type_free(MPI_RE2V1_DATA_CV, ierr)
    call MPI_Type_free(MPI_RE2V2_DATA_CV, ierr)
  end subroutine mpi_type_re2_cv_free
  
  !> Deallocate re2 bc derived MPI type
  subroutine mpi_type_re2_bc_free
    integer ierr
    call MPI_Type_free(MPI_RE2V1_DATA_BC, ierr)
    call MPI_Type_free(MPI_RE2V2_DATA_BC, ierr)
  end subroutine mpi_type_re2_bc_free

  !> Deallocate STL header dervied MPI type
  subroutine mpi_type_stl_header_free
    integer ierr
    call MPI_Type_free(MPI_STL_HEADER, ierr)
  end subroutine mpi_type_stl_header_free
  
  !> Deallocate STL triangle derived MPI type
  subroutine mpi_type_stl_triangle_free
    integer ierr
    call MPI_Type_free(MPI_STL_TRIANGLE, ierr)
  end subroutine mpi_type_stl_triangle_free
end module mpi_types
