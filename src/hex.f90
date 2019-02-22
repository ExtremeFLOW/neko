!> Defines a hexahedron element
!
module hex
  use num_types
  use element
  implicit none
  private

  type, public, extends(element_t) :: hex_t
   contains
     procedure, pass(this) :: hex_init
     procedure, pass(this) :: diameter => hex_diameter
     generic, public :: init => hex_init
  end type hex_t

contains
  
  subroutine hex_init(this, dim, npts)
    class(hex_t), intent(inout) :: this
    integer, intent(inout) :: dim
    integer, intent(in) :: npts

    call this%init(dim)

  end subroutine hex_init

  function hex_diameter(this) result(res)
    class(hex_t), intent(in) :: this
    integer :: res
  end function hex_diameter
  
end module hex
