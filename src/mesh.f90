!> Defines a mesh
module mesh
  use num_types
  use point
  use element
  use hex
  use quad
  use utils
  implicit none

  type, private :: mesh_element_t
     class(element_t), allocatable :: e
  end type mesh_element_t

  type mesh_t

     integer :: nelv            !< Number of elements
     integer :: npts            !< Number of points per element
     integer :: gdim            !< Geometric dimension

     type(point_t), allocatable :: points(:) !< List of points
     type(mesh_element_t), allocatable :: elements(:) !< List of elements
  end type mesh_t

  !> Add an element to the mesh
  interface mesh_add_element
     module procedure mesh_add_quad, mesh_add_hex
  end interface mesh_add_element
  

contains 

  !> Initialise a mesh @a m
  subroutine mesh_init(m, gdim, nelv)
    type(mesh_t), intent(inout) :: m
    integer, intent(in) :: gdim
    integer, intent(in) :: nelv
    integer :: i, npts
    type(point_t) :: p(4)
    
    call mesh_free(m)

    m%nelv = nelv
    m%gdim = gdim

    allocate(m%elements(nelv))
    if (gdim .eq. 3) then
       do i = 1, nelv
          allocate(hex_t::m%elements(i)%e)
       end do
       m%npts = NEKO_HEX_NPTS
    else if (gdim .eq. 2) then
       do i = 1, nelv
          allocate(quad_t::m%elements(i)%e)
       end do
       m%npts = NEKO_QUAD_NPTS
    else
       call neko_error("Invalid dimension")
    end if

    allocate(m%points(m%npts*m%nelv))
    
  end subroutine mesh_init
  
  subroutine mesh_free(m)
    type(mesh_t), intent(inout) :: m
    integer :: i
    
    if (allocated(m%points)) then
       deallocate(m%points)
    end if

    if (allocated(m%elements)) then
       do i = 1, m%nelv
          deallocate(m%elements(i)%e)
       end do
    end if

  end subroutine mesh_free

  !> Add a quadrilateral element to the mesh @a m
  subroutine mesh_add_quad(m, el, p1, p2, p3, p4)
    type(mesh_t), target, intent(inout) :: m
    integer, intent(inout) :: el
    type(point_t), intent(inout) :: p1, p2, p3, p4
    class(element_t), pointer :: ep
    integer :: pt_offset

    pt_offset = m%npts * (el - 1)

    m%points(pt_offset + 1) = p1
    m%points(pt_offset + 2) = p2
    m%points(pt_offset + 3) = p3
    m%points(pt_offset + 4) = p4
    
    ep => m%elements(el)%e
    select type(ep)
    type is (quad_t)
       call ep%init(el, m%points(pt_offset + 1), m%points(pt_offset + 2), &
            m%points(pt_offset + 3), m%points(pt_offset + 4))
    class default
       call neko_error('Invalid element type')
    end select
    
  end subroutine mesh_add_quad

  !> Add a hexahedral element to the mesh @a m
  subroutine mesh_add_hex(m, el, p1, p2, p3, p4, p5, p6, p7, p8)
    type(mesh_t), target, intent(inout) :: m
    integer, intent(inout) :: el
    type(point_t), intent(inout) :: p1, p2, p3, p4, p5, p6, p7, p8
    class(element_t), pointer :: ep
    integer :: pt_offset

    pt_offset = m%npts * (el - 1)

    m%points(pt_offset + 1) = p1
    m%points(pt_offset + 2) = p2
    m%points(pt_offset + 3) = p3
    m%points(pt_offset + 4) = p4
    m%points(pt_offset + 5) = p5
    m%points(pt_offset + 6) = p6
    m%points(pt_offset + 7) = p7
    m%points(pt_offset + 8) = p8

    ep => m%elements(el)%e
    select type(ep)
    type is (hex_t)
       call ep%init(el, m%points(pt_offset + 1), m%points(pt_offset + 2), &
            m%points(pt_offset + 3), m%points(pt_offset + 4), &
            m%points(pt_offset + 5), m%points(pt_offset + 6), &
            m%points(pt_offset + 7), m%points(pt_offset + 8))
    class default
       call neko_error('Invalid element type')
    end select
    
  end subroutine mesh_add_hex


end module mesh
