!> NEKTON re2 format
module re2
  use num_types
  implicit none

  !> NEKTON re2 header
  type, public :: re2_hdr_t
     character(len=5) :: hdr_ver
     integer :: nel
     integer :: ndim
     integer :: nelv
     character(len=54) :: hdr_str
     real(kind=sp) :: endian_test
  end type re2_hdr_t

end module re2
