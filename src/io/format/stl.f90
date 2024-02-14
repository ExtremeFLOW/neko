!> Stereolithography format
module stl
  use num_types, only : sp, i2
  implicit none
  private

  !> Defines a STL hdr
  type, public :: stl_hdr_t
     character(len=80) :: hdr
     integer :: ntri
  end type stl_hdr_t

  !> Defines a STL triangle
  type, public :: stl_triangle_t
     real(kind=sp) :: n(3)
     real(kind=sp) :: v1(3)
     real(kind=sp) :: v2(3)
     real(kind=sp) :: v3(3)
     integer(kind=i2) :: attrib
  end type stl_triangle_t

end module stl
