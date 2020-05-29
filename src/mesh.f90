!> Defines a mesh
module mesh
  use num_types
  use point
  use element
  use hex
  use quad
  use utils
  use htable
  use mpi
  use comm
  use datadist
  implicit none

  type, private :: mesh_element_t
     class(element_t), allocatable :: e
  end type mesh_element_t

  type mesh_t

     integer :: nelv            !< Number of elements
     integer :: npts            !< Number of points per element
     integer :: gdim            !< Geometric dimension
     integer :: mpts            !< Number of (unique) points in the mesh

     integer :: glb_nelv        !< Global number of elements
     integer :: offset_el       !< Element offset
     
     type(point_t), allocatable :: points(:) !< list of points
     type(mesh_element_t), allocatable :: elements(:) !< List of elements
     
     !> @todo flush this table once mesh is finalized
     type(htable_i4_t) :: htp   !< Table of unique points

     logical :: lconn = .false.     !< valid connectivity
     logical :: finalized = .false. !< Valid mesh
  end type mesh_t

  !> Initialise a mesh
  interface mesh_init
     module procedure mesh_init_nelv, mesh_init_dist
  end interface mesh_init
  
  !> Add an element to the mesh
  interface mesh_add_element
     module procedure mesh_add_quad, mesh_add_hex
  end interface mesh_add_element

  private :: mesh_init_common, mesh_add_quad, mesh_add_hex

contains 

  !> Initialise a mesh @a m with @a nelv elements
  subroutine mesh_init_nelv(m, gdim, nelv)
    type(mesh_t), intent(inout) :: m !< Mesh
    integer, intent(in) :: gdim      !< Geometric dimension
    integer, intent(in) :: nelv      !< Local number of elements
    integer :: ierr
    
    call mesh_free(m)

    m%nelv = nelv
    m%gdim = gdim

    call MPI_Allreduce(m%nelv, m%glb_nelv, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    m%offset_el = 0
    call MPI_Exscan(m%nelv, m%offset_el, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    call mesh_init_common(m)
    
  end subroutine mesh_init_nelv

  !> Initialise a mesh @a m based on a distribution @a dist
  subroutine mesh_init_dist(m, gdim, dist)
    type(mesh_t), intent(inout) :: m        !< Mesh
    integer, intent(in) :: gdim             !< Geometric dimension
    type(linear_dist_t), intent(in) :: dist !< Data distribution

    call mesh_free(m)
    
    m%nelv = dist%num_local()
    m%glb_nelv = dist%num_global()
    m%offset_el = dist%start_idx()
    m%gdim = gdim

    call mesh_init_common(m)
    
  end subroutine mesh_init_dist

  subroutine mesh_init_common(m)
    type(mesh_t), intent(inout) :: m
    integer :: i

    allocate(m%elements(m%nelv))
    if (m%gdim .eq. 3) then
       do i = 1, m%nelv
          allocate(hex_t::m%elements(i)%e)
       end do
       m%npts = NEKO_HEX_NPTS
    else if (m%gdim .eq. 2) then
       do i = 1, m%nelv
          allocate(quad_t::m%elements(i)%e)
       end do
       m%npts = NEKO_QUAD_NPTS
    else
       call neko_error("Invalid dimension")
    end if

    allocate(m%points(m%npts*m%nelv))

    call m%htp%init(m%npts*m%nelv, i)
    m%mpts = 0    

  end subroutine mesh_init_common
  
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
    integer :: p(4), el_glb_idx

    ! Connectivity invalidated if a new element is added        
    m%lconn = .false.           

    call mesh_add_point(m, p1, p(1))
    call mesh_add_point(m, p2, p(2))
    call mesh_add_point(m, p3, p(3))
    call mesh_add_point(m, p4, p(4))

    ep => m%elements(el)%e
    el_glb_idx = el + m%offset_el
    select type(ep)
    type is (quad_t)
       call ep%init(el_glb_idx, &
            m%points(p(1)), m%points(p(2)), &
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
    integer :: p(8), el_glb_idx

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
    el_glb_idx = el + m%offset_el
    select type(ep)
    type is (hex_t)
       call ep%init(el_glb_idx, &
            m%points(p(1)), m%points(p(2)), &
            m%points(p(3)), m%points(p(4)), &
            m%points(p(5)), m%points(p(6)), &
            m%points(p(7)), m%points(p(8)))
    class default
       call neko_error('Invalid element type')
    end select

  end subroutine mesh_add_hex

  !> Add a unique point to the mesh
  !! @todo remove hash table is only necessary for legacy formats
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
