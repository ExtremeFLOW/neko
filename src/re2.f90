!> NEKTON re2 format
module re2
  use num_types
  implicit none

  !> NEKTON re2 header size
  integer, parameter :: RE2_HDR_SIZE = 80

  !> NEKTION re2 endian test
  real(kind=sp), parameter :: RE2_ENDIAN_TEST = 6.54321

  !> NEKTON re2 element data
  type :: re2_t
     real(kind=sp) :: rgroup
  end type re2_t

  !> NEKTON re2 element data (3d)
  type, public, extends(re2_t) :: re2_xyz_t    
     real(kind=sp), dimension(8) :: x
     real(kind=sp), dimension(8) :: y
     real(kind=sp), dimension(8) :: z
  end type re2_xyz_t

  !> NEKTON re2 element data (2d)
  type, public, extends(re2_t) :: re2_xy_t    
     real(kind=sp), dimension(4) :: x
     real(kind=sp), dimension(4) :: y
  end type re2_xy_t

end module re2
