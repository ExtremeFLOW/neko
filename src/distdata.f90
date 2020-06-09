!> Distributed mesh data
module distdata
  use stack
  use tuple
  use uset
  implicit none
  
  type, public :: distdata_t
     type(stack_i4t2_t) :: shared_facet !< Elemenets with shared facets
     type(stack_i4t2_t) :: shared_edge  !< Elements with shared edges
     type(uset_i4_t) :: shared_point    !< List of shared points

  end type distdata_t

contains

  !> Initialise a distdata type
  subroutine distdata_init(distdata)
    type(distdata_t), intent(inout) :: distdata

    call distdata%shared_facet%init()
    call distdata%shared_edge%init()   
    call distdata%shared_point%init()
    
  end subroutine distdata_init

  !> Free a distdata type
  subroutine distdata_free(distdata)
    type(distdata_t), intent(inout) :: distdata

    call distdata%shared_facet%free()
    call distdata%shared_edge%free()
    call distdata%shared_point%free()
    
  end subroutine distdata_free

  !> Mark an element's facet as shared
  subroutine distdata_set_shared_facet(distdata, element, side)
    type(distdata_t), intent(inout) :: distdata
    integer, intent(in), value :: element !< Element index (local numbering)
    integer, intent(in), value :: side    !< Facet index
    type(tuple_i4_t) :: t

    t = (/ element, side /)
    call distdata%shared_facet%push(t)
    
  end subroutine distdata_set_shared_facet

  !> Mark an element's edge as shared
  !! @attention only defined for elements where facet .ne. edges
  subroutine distdata_set_shared_edge(distdata, element, edge)
    type(distdata_t), intent(inout) :: distdata
    integer, intent(in), value :: element !< Element index (local numbering)
    integer, intent(in), value :: edge    !< Edge index
    type(tuple_i4_t) :: t

    t = (/ element, edge /)
    call distdata%shared_edge%push(t)
    
  end subroutine distdata_set_shared_edge

  !> Mark a point as shared
  subroutine distdata_set_shared_point(distdata, point)
    type(distdata_t), intent(inout) :: distdata
    integer, value :: point !< Point index (local numbering)

    call distdata%shared_point%add(point)
    
  end subroutine distdata_set_shared_point
  
end module distdata
