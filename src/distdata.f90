!> Distributed mesh data
module distdata
  use stack
  use tuple
  implicit none
  
  type, public ::  distdata_t
     type(stack_i4t2_t) :: shared_facet !< Elemenets with shared facets
  end type distdata_t

contains

  !> Initialise a distdata type
  subroutine distdata_init(distdata)
    type(distdata_t), intent(inout) :: distdata

    call distdata%shared_facet%init()
    
  end subroutine distdata_init

  !> Free a distdata type
  subroutine distdata_free(distdata)
    type(distdata_t), intent(inout) :: distdata

    call distdata%shared_facet%free()
    
  end subroutine distdata_free

  !> Mark an element's facet as shared
  subroutine distdata_set_shared(distdata, element, side)
    type(distdata_t), intent(inout) :: distdata
    integer, intent(in), value :: element !< Element index (local numbering)
    integer, intent(in), value :: side    !< Facet index
    type(tuple_i4_t) :: t

    t = (/ element, side /)
    call distdata%shared_facet%push(t)
    
  end subroutine distdata_set_shared
  
end module distdata
