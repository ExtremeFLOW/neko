!> Distributed mesh data
module distdata
  use stack
  use tuple
  use uset
  implicit none
  
  type, public :: distdata_t
     type(stack_i4t2_t) :: shared_el_facet !< Elemenets with shared facets
     
     type(uset_i4_t) :: shared_facet    !< List of shared facets
     type(uset_i4_t) :: shared_edge     !< List of shared edges
     type(uset_i4_t) :: shared_point    !< List of shared points
     
     integer, allocatable :: local_to_global_facet(:)!< Local to global (facets)
     integer, allocatable :: local_to_global_edge(:) !< Local to global (edges)
     
  end type distdata_t

contains

  !> Initialise a distdata type
  subroutine distdata_init(distdata)
    type(distdata_t), intent(inout) :: distdata

    call distdata%shared_el_facet%init()

    call distdata%shared_facet%init()
    call distdata%shared_edge%init()   
    call distdata%shared_point%init()
    
  end subroutine distdata_init

  !> Free a distdata type
  subroutine distdata_free(distdata)
    type(distdata_t), intent(inout) :: distdata

    call distdata%shared_el_facet%free()
    
    call distdata%shared_facet%free()
    call distdata%shared_edge%free()
    call distdata%shared_point%free()

    if (allocated(distdata%local_to_global_facet)) then
       deallocate(distdata%local_to_global_facet)
    end if
    
    if (allocated(distdata%local_to_global_edge)) then
       deallocate(distdata%local_to_global_edge)
    end if
    
  end subroutine distdata_free

  !> Mark an element's facet as shared
  subroutine distdata_set_shared_el_facet(distdata, element, side)
    type(distdata_t), intent(inout) :: distdata
    integer, intent(in), value :: element !< Element index (local numbering)
    integer, intent(in), value :: side    !< Facet index
    type(tuple_i4_t) :: t

    t = (/ element, side /)
    call distdata%shared_el_facet%push(t)
    
  end subroutine distdata_set_shared_el_facet

  !> Mark an element's edge as shared
  !! @attention only defined for elements where facet .ne. edges
  subroutine distdata_set_shared_edge(distdata, edge)
    type(distdata_t), intent(inout) :: distdata
    integer, value :: edge      !< Edge index (local numbering) 

    call distdata%shared_edge%add(edge)
    
  end subroutine distdata_set_shared_edge

  !> Mark a point as shared
  subroutine distdata_set_shared_point(distdata, point)
    type(distdata_t), intent(inout) :: distdata
    integer, value :: point !< Point index (local numbering)

    call distdata%shared_point%add(point)
    
  end subroutine distdata_set_shared_point

  !> Set local to global mapping (facets)
  subroutine distdata_set_local_to_global_facet(distdata, local, global)
    type(distdata_t), intent(inout) :: distdata
    integer, intent(in), value :: local  !< Local facet index
    integer, intent(in), value :: global !< Global facet index

    distdata%local_to_global_facet(local) = global
    
  end subroutine distdata_set_local_to_global_facet
  !> Set local to global mapping (edges)
  subroutine distdata_set_local_to_global_edge(distdata, local, global)
    type(distdata_t), intent(inout) :: distdata
    integer, intent(in) , value :: local  !< Local edge index
    integer, intent(in) , value :: global !< Global edge index

    distdata%local_to_global_edge(local) = global
    
  end subroutine distdata_set_local_to_global_edge
  
end module distdata
