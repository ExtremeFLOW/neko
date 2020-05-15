!> Defines a mesh
module mesh
  use num_types
  use point
  use element
  use hex
  use quad
  use utils
  use htable
  implicit none

  type, private :: mesh_element_t
     class(element_t), allocatable :: e
  end type mesh_element_t

  type mesh_t

     integer :: nelv            !< Number of elements
     integer :: npts            !< Number of points per element
     integer :: gdim            !< Geometric dimension
     integer :: mpts            !< Number of (unique) points in the mesh

     type(point_t), allocatable :: points(:) !< list of points
     type(mesh_element_t), allocatable :: elements(:) !< List of elements
     
     !> @todo flush this table once mesh is finalized
     type(htable_i4_t) :: htp   !< Table of unique points

     logical :: lconn           !< Valid connectivity
     
  end type mesh_t

  !> Add an element to the mesh
  interface mesh_add_element
     module procedure mesh_add_quad, mesh_add_hex
  end interface mesh_add_element

  private :: mesh_add_quad, mesh_add_hex

contains 

  !> Initialise a mesh @a m
  subroutine mesh_init(m, gdim, nelv)
    type(mesh_t), intent(inout) :: m
    integer, intent(in) :: gdim
    integer, intent(in) :: nelv
    integer :: i
    
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

    call m%htp%init(m%npts*m%nelv, i)
    m%mpts = 0

    m%lconn = .false.
    
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
    integer, value :: el
    type(point_t), intent(inout) :: p1, p2, p3, p4
    class(element_t), pointer :: ep
    integer :: p(4)

    ! Connectivity invalidated if a new element is added        
    m%lconn = .false.           

    call mesh_add_point(m, p1, p(1))
    call mesh_add_point(m, p2, p(2))
    call mesh_add_point(m, p3, p(3))
    call mesh_add_point(m, p4, p(4))

    ep => m%elements(el)%e
    select type(ep)
    type is (quad_t)
       call ep%init(el, m%points(p(1)), m%points(p(2)), &
            m%points(p(3)), m%points(p(4)))
    class default
       call neko_error('Invalid element type')
    end select
        
  end subroutine mesh_add_quad

  !> Add a hexahedral element to the mesh @a m
  subroutine mesh_add_hex(m, el, p1, p2, p3, p4, p5, p6, p7, p8)
    type(mesh_t), target, intent(inout) :: m
    integer, value :: el
    type(point_t), intent(inout) :: p1, p2, p3, p4, p5, p6, p7, p8
    class(element_t), pointer :: ep
    integer :: p(8)

    ! Connectivity invalidated if a new element is added        
    m%lconn = .false.
    
    call mesh_add_point(m, p1, p(1))
    call mesh_add_point(m, p2, p(2))
    call mesh_add_point(m, p3, p(3))
    call mesh_add_point(m, p4, p(4))
    call mesh_add_point(m, p5, p(5))
    call mesh_add_point(m, p6, p(6))
    call mesh_add_point(m, p7, p(7))
    call mesh_add_point(m, p8, p(8))

    ep => m%elements(el)%e
    select type(ep)
    type is (hex_t)
       call ep%init(el, m%points(p(1)), m%points(p(2)), &
            m%points(p(3)), m%points(p(4)), &
            m%points(p(5)), m%points(p(6)), &
            m%points(p(7)), m%points(p(8)))
    class default
       call neko_error('Invalid element type')
    end select

  end subroutine mesh_add_hex

  !> Add a unique point to the mesh
  subroutine mesh_add_point(m, p, idx)
    type(mesh_t), intent(inout) :: m
    type(point_t), intent(inout) :: p
    integer, intent(inout) :: idx
    integer :: tmp
   
    tmp = p%id()
    
    if (tmp .le. 0) then
       call neko_error("Invalid point id")
    end if

    if (m%htp%get(tmp, idx) .gt. 0) then
       m%mpts = m%mpts + 1
       call m%htp%set(tmp, m%mpts)
       m%points(m%mpts) = p
       idx = m%mpts
    end if
    
  end subroutine mesh_add_point


end module mesh
