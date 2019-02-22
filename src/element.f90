module element
  use num_types
  implicit none
  private

  !> Base type for an element
  type, public, abstract :: element_t
     integer, private :: gdim            !< Geometric dimension
     integer, private :: npts
   contains
     procedure, pass(this) :: element_init
     procedure(element_diameter), pass(this), deferred :: diameter
     procedure :: n_points => element_npts
     generic, public :: init => element_init
  end type element_t

  abstract interface
     function element_diameter(this) result(res)
       import :: element_t
       class(element_t), intent(in) :: this
       integer :: res
     end function element_diameter
  end interface

contains

  subroutine element_init(this, dim)
    class(element_t), intent(inout)  :: this
    integer, intent(in) :: dim
    this%gdim = dim
    write(*,*) 'element'
  end subroutine element_init
  
  function element_npts(this) result(npts)
    class(element_t), intent(in) :: this
    integer :: npts
    npts = this%npts
  end function element_npts

end module element
