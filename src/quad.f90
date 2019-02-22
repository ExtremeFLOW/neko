!> Defines a quadrilateral element
!
module quad
  use num_types
  use element
  implicit none
  private

  type, public, extends(element_t) :: quad_t
   contains
     procedure, pass(this) :: quad_init
     procedure, pass(this) :: diameter => quad_diameter
     generic, public :: init => quad_init
  end type quad_t

contains
  
  subroutine quad_init(this, dim, npts)
    class(quad_t), intent(inout) :: this
    integer, intent(inout) :: dim
    integer, intent(in) :: npts

    call this%init(dim)

  end subroutine quad_init

  function quad_diameter(this) result(res)
    class(quad_t), intent(in) :: this
    integer :: res
  end function quad_diameter
  
end module quad
