!> NEKTON re2 format
module re2
  use num_types
  implicit none

  !> NEKTON re2 header size
  integer, parameter :: RE2_HDR_SIZE = 80

  !> NEKTION re2 endian test
  real(kind=sp), parameter :: RE2_ENDIAN_TEST = 6.54321

  !> NEKTON re2 element data (version 1)
  type :: re2v1_t
     real(kind=sp) :: rgroup
  end type re2v1_t

  !> NEKTON re2 element data (3d) (version 1)
  type, public, extends(re2v1_t) :: re2v1_xyz_t
     real(kind=sp), dimension(8) :: x
     real(kind=sp), dimension(8) :: y
     real(kind=sp), dimension(8) :: z
  end type re2v1_xyz_t

  !> NEKTON re2 element data (2d) (version 1)
  type, public, extends(re2v1_t) :: re2v1_xy_t
     real(kind=sp), dimension(4) :: x
     real(kind=sp), dimension(4) :: y
  end type re2v1_xy_t

  !> NEKTON re2 curve data (version 1)
  type, public :: re2v1_curve_t
     integer :: elem
     integer :: zone
     real(kind=sp), dimension(5) :: point
     character(len=4) :: type
  end type re2v1_curve_t

  !> NEKTON re2 bc data (version 1)
  type, public :: re2v1_bc_t
     integer :: elem
     integer :: face
     real(kind=sp), dimension(5) :: bc_data
     character(len=4) :: type
  end type re2v1_bc_t

  !> NEKTON re2 element data (version 2)
  type :: re2v2_t
     real(kind=dp) :: rgroup
  end type re2v2_t

  !> NEKTON re2 element data (3d) (version 2)
  type, public, extends(re2v2_t) :: re2v2_xyz_t
     real(kind=dp), dimension(8) :: x
     real(kind=dp), dimension(8) :: y
     real(kind=dp), dimension(8) :: z
  end type re2v2_xyz_t

  !> NEKTON re2 element data (2d) (version 2)
  type, public, extends(re2v1_t) :: re2v2_xy_t
     real(kind=dp), dimension(4) :: x
     real(kind=dp), dimension(4) :: y
  end type re2v2_xy_t

  !> NEKTON re2 curve data (version 2)
  type, public :: re2v2_curve_t
     real(kind=dp) :: elem
     real(kind=dp) :: zone
     real(kind=dp), dimension(5) :: point
     character(len=8) :: type
  end type re2v2_curve_t

  !> NEKTON re2 bc data (version 2)
  type, public :: re2v2_bc_t
     real(kind=dp) :: elem
     real(kind=dp) :: face
     real(kind=dp), dimension(5) :: bc_data
     character(len=8) :: type
  end type re2v2_bc_t

end module re2
