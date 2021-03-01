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
   contains
     procedure, pass(z) :: init => zone_init
     procedure, pass(z) :: free => zone_free
     procedure, pass(z) :: finalize => zone_finalize
     procedure, pass(z) :: add_facet => zone_add_facet
  end type zone_t

  type, extends(zone_t) :: zone_periodic_t
     type(tuple_i4_t), allocatable :: p_facet_el(:)
     type(stack_i4t2_t), private :: p_scratch
     type(tuple4_i4_t), allocatable :: p_ids(:)
     type(stack_i4t4_t), private :: p_id_scratch
   contains
     procedure, pass(z) :: init => zone_periodic_init
     procedure, pass(z) :: free => zone_periodic_free
     procedure, pass(z) :: finalize => zone_periodic_finalize
     procedure, pass(z) :: add_periodic_facet => zone_periodic_add_facet
  end type zone_periodic_t
  
contains

  !> Initialize a zone
  subroutine zone_init(z, size)
    class(zone_t), intent(inout) :: z
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
    class(zone_t), intent(inout) :: z
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
    class(zone_t), intent(inout) :: z
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
    class(zone_t), intent(inout) :: z
    integer, intent(in) :: facet   !< Facet in the zone
    integer, intent(in) :: el      !< Element  in the zone
    type(tuple_i4_t) :: t

    if (z%finalized) then
       call neko_error('Zone already finalized')
    end if

    t = (/ facet, el /)
    call z%scratch%push(t)
    
  end subroutine zone_add_facet

    !> Initialize a periodic zone
  subroutine zone_periodic_init(z, size)
    class(zone_periodic_t), intent(inout) :: z
    integer, optional :: size

    call z%free()

    if (present(size)) then
       call zone_init(z, size)
       call z%p_scratch%init(size)
       call z%p_id_scratch%init(size)
    else
       call zone_init(z)
       call z%p_scratch%init()
       call z%p_id_scratch%init()
    end if
    
  end subroutine zone_periodic_init

  !> Deallocate a zone
  subroutine zone_periodic_free(z)
    class(zone_periodic_t), intent(inout) :: z

    call zone_free(z)

    if (allocated(z%p_facet_el)) then
       deallocate(z%p_facet_el)
    end if

    if (allocated(z%p_ids)) then
       deallocate(z%p_ids)
    end if

    call z%p_scratch%free()
    call z%p_id_scratch%free()
    
  end subroutine zone_periodic_free

  !> Finalize a periodic zone list
  !! @details Create a static list of (facet,el) tuples
  subroutine zone_periodic_finalize(z)
    class(zone_periodic_t), intent(inout) :: z
    type(tuple_i4_t), pointer :: tp(:)
    type(tuple4_i4_t), pointer :: tp2(:)
    integer :: i
    
    if (.not. z%finalized) then

       call zone_finalize(z)

       if (z%size .ne. z%p_scratch%size()) then
          call neko_error('Zone size mismatch')
       end if

       allocate(z%p_facet_el(z%size))
       allocate(z%p_ids(z%size))
       
       tp => z%p_scratch%array()
       do i = 1, z%size
          z%p_facet_el(i) = tp(i)
       end do
       tp2 => z%p_id_scratch%array()
       do i = 1, z%size
          z%p_ids(i) = tp2(i)
       end do

       call z%p_scratch%clear()
       call z%p_id_scratch%clear()

    end if
    
  end subroutine zone_periodic_finalize

  !> Add a (facet, el) tuple to an unfinalized zone
  subroutine zone_periodic_add_facet(z, facet, el, p_facet, p_el, pids)
    class(zone_periodic_t), intent(inout) :: z
    integer, intent(in) :: facet   !< Facet in the zone
    integer, intent(in) :: el      !< Element  in the zone
    integer, intent(in) :: p_facet !< Facet at periodic length
    integer, intent(in) :: p_el    !< Element at periodic length
    integer, intent(in) :: pids(4) !< Element at periodic length ids
    type(tuple_i4_t) :: t
    type(tuple4_i4_t) :: t2

    if (z%finalized) then
       call neko_error('Zone already finalized')
    end if

    call z%add_facet(facet, el)

    t = (/ p_facet, p_el /)
    call z%p_scratch%push(t)
    t2 = pids
    call z%p_id_scratch%push(t2)
    
  end subroutine zone_periodic_add_facet
  
end module zone
