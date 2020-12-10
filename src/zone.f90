!> Defines a zone as a subset of facets in a mesh
module zone
  use tuple
  use stack
  use utils
  implicit none
  
  type :: zone_t
     type(tuple_i4_t), allocatable :: facet_el(:)
     integer :: size = 0
     logical, private :: finalized = .false.
     type(stack_i4t2_t), private :: scratch
  end type zone_t

contains

  !> Initialize a zone
  subroutine zone_init(z, size)
    type(zone_t), intent(inout) :: z
    integer, optional :: size

    call zone_free(z)

    if (present(size)) then
       call z%scratch%init(size)
    else
       call z%scratch%init()
    end if
    
  end subroutine zone_init

  !> Deallocate a zone
  subroutine zone_free(z)
    type(zone_t), intent(inout) :: z
    if (allocated(z%facet_el)) then
       deallocate(z%facet_el)
    end if

    z%finalized = .false.
    z%size = 0

    call z%scratch%free()
    
  end subroutine zone_free

  !> Finalize a zone list
  !! @details Create a static list of (facet,el) tuples
  subroutine zone_finalize(z)
    type(zone_t), intent(inout) :: z
    type(tuple_i4_t), pointer :: tp(:)
    integer :: i
    
    if (.not. z%finalized) then

       allocate(z%facet_el(z%scratch%size()))
       
       tp => z%scratch%array()
       do i = 1, z%scratch%size()
          z%facet_el(i) = tp(i)
       end do

       z%size = z%scratch%size()

       call z%scratch%clear()

       z%finalized = .true.
       
    end if
    
  end subroutine zone_finalize

  !> Add a (facet, el) tuple to an unfinalized zone
  subroutine zone_add_facet(z, facet, el)
    type(zone_t), intent(inout) :: z
    integer, intent(in) :: facet
    integer, intent(in) :: el
    type(tuple_i4_t) :: t

    if (z%finalized) then
       call neko_error('Zone already finalized')
    end if

    t = (/ facet, el /)
    call z%scratch%push(t)
    
  end subroutine zone_add_facet
  
end module zone
