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
  use stack
  use tuple
  use htable
  use datadist
  use distdata
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
     
     type(htable_i4_t) :: htp   !< Table of unique points (global->local)

     integer, allocatable :: facet_neigh(:,:)  !< Facet to neigh. element table
     class(htable_t), allocatable :: facet_map !< Facet to element's id tuple 
                                               !! \f$ t=(odd, even) \f$

     type(stack_i4_t), allocatable :: point_neigh(:) !< Point to neigh. table

     type(distdata_t) :: distdata              !< Mesh distributed data

     logical :: lconn = .false.                !< valid connectivity
     logical :: ldist = .false.                !< valid distributed data

  end type mesh_t

  !> Initialise a mesh
  interface mesh_init
     module procedure mesh_init_nelv, mesh_init_dist
  end interface mesh_init
  
  !> Add an element to the mesh
  interface mesh_add_element
     module procedure mesh_add_quad, mesh_add_hex
  end interface mesh_add_element

  !> Get local id for a mesh entity
  !! @todo Add similar mappings for element ids
  interface mesh_get_local
     module procedure mesh_get_local_point
  end interface mesh_get_local

  private :: mesh_init_common, mesh_add_quad, mesh_add_hex, &
       mesh_generate_external_facet_conn, mesh_generate_external_point_conn
  

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
    type(tuple_i4_t) :: facet_data

    allocate(m%elements(m%nelv))
    if (m%gdim .eq. 3) then
       do i = 1, m%nelv
          allocate(hex_t::m%elements(i)%e)
       end do
       m%npts = NEKO_HEX_NPTS

       allocate(htable_i4t4_t::m%facet_map)
       select type (fmp => m%facet_map)
       type is(htable_i4t4_t)
          call fmp%init(m%nelv, facet_data)
       end select

       allocate(m%facet_neigh(6, m%nelv))
    else if (m%gdim .eq. 2) then
       do i = 1, m%nelv
          allocate(quad_t::m%elements(i)%e)
       end do
       m%npts = NEKO_QUAD_NPTS

       allocate(htable_i4t2_t::m%facet_map)       
       select type (fmp => m%facet_map)
       type is(htable_i4t2_t)
          call fmp%init(m%nelv, facet_data)
       end select

       allocate(m%facet_neigh(4, m%nelv))
    else
       call neko_error("Invalid dimension")
    end if

    !> @todo resize onces final size is known
    allocate(m%points(m%npts*m%nelv))

    !> @todo resize onces final size is known
    allocate(m%point_neigh(m%npts*m%nelv))
    do i = 1, m%npts*m%nelv
       call m%point_neigh(i)%init()
    end do
    
    call m%htp%init(m%npts*m%nelv, i)
   
    call distdata_init(m%distdata)
    
    m%mpts = 0    

  end subroutine mesh_init_common
  
  !> Deallocate a mesh %a m
  subroutine mesh_free(m)
    type(mesh_t), intent(inout) :: m
    integer :: i
    
    call m%htp%free()

    if (allocated(m%points)) then
       deallocate(m%points)
    end if

    if (allocated(m%elements)) then
       do i = 1, m%nelv
          deallocate(m%elements(i)%e)
       end do
    end if

    if (allocated(m%facet_map)) then
       select type (fmp => m%facet_map)
       type is(htable_i4t2_t)
          call fmp%free()
       type is(htable_i4t4_t)
          call fmp%free()
       end select
       deallocate(m%facet_map)
    end if

    if (allocated(m%facet_neigh)) then
       deallocate(m%facet_neigh)
    end if

    if (allocated(m%point_neigh)) then
       do i = 1, m%npts * m%nelv
          call m%point_neigh(i)%free()
       end do
       deallocate(m%point_neigh)
    end if

  end subroutine mesh_free

  !> Generate element-to-element connectivity
  subroutine mesh_generate_conn(m)
    type(mesh_t), intent(inout) :: m
    type(tuple_i4_t) :: edge, facet_key
    type(tuple4_i4_t) :: face
    integer :: i, j, k, el_glb_idx, n_sides, n_nodes

    if (m%lconn) return

    if (m%gdim .eq. 2) then
       n_sides = 4
       n_nodes = 2
    else
       n_sides = 6
       n_nodes = 4
    end if

    !
    ! Find all (local) boundaries
    !
    
    !> @note We have to sweep through the facet map twice to make sure
    !! that both odd and even sides are marked
    !! @todo These loop nests needs a lot of love...
    select type (fmp => m%facet_map)
    type is(htable_i4t2_t)
       do k = 1, 2              
          do i = 1, m%nelv
             el_glb_idx = i + m%offset_el
             do j = 1, n_sides
                call m%elements(i)%e%facet_id(edge, j)
                
                ! Assume that all facets are on the exterior
                facet_key = (/ 0, 0 /)
                
                if (fmp%get(edge, facet_key) .gt. 0) then
                   if (mod(j, 2) .gt. 0) then
                      facet_key%x(1) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(2)
                   else
                      facet_key%x(2) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(1)
                   end if
                   call fmp%set(edge, facet_key)
                else
                   if (mod(j, 2) .gt. 0) then
                      facet_key%x(1) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(2)
                   else
                      facet_key%x(2) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(1)
                   end if
                   call fmp%set(edge, facet_key)
                end if
             end do
          end do
       end do
    type is(htable_i4t4_t)
       do k = 1, 2
          do i = 1, m%nelv
             el_glb_idx = i + m%offset_el
             do j = 1, n_sides
                call m%elements(i)%e%facet_id(face, j)
                
                ! Assume that all facets are on the exterior
                facet_key = (/ 0, 0 /)
                
                if (fmp%get(face, facet_key) .gt. 0) then
                   if(mod(j, 2) .gt. 0) then
                      facet_key%x(1) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(2)
                   else 
                      facet_key%x(2) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(1)

                   end if
                   call fmp%set(face, facet_key)                             
                else
                   if(mod(j, 2) .gt. 0) then
                      facet_key%x(1) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(2)
                   else
                      facet_key%x(2) = el_glb_idx
                      m%facet_neigh(j, i) = facet_key%x(1)
                   end if
                   call fmp%set(face, facet_key)               
                end if
          end do
       end do
    end do
    class default
       call neko_error('Invalid facet map')
    end select

    !
    ! Find all external (between PEs) boundaries
    !
    if (pe_size .gt. 1) then
       call mesh_generate_external_facet_conn(m)

       call mesh_generate_external_point_conn(m)
    end if
    
    m%lconn = .true.
    
  end subroutine mesh_generate_conn
 
  !> Generate element-element connectivity via facets between PEs
  subroutine mesh_generate_external_facet_conn(m)
    type(mesh_t), intent(inout) :: m
    type(tuple_i4_t) :: edge, facet_key
    type(tuple4_i4_t) :: face
    type(stack_i4_t) :: buffer
    integer, allocatable :: recv_buffer(:)
    integer :: i, j, k, el_glb_idx, n_sides, n_nodes, facet, element
    integer :: max_recv, ierr, src, dst, n_recv, recv_side, neigh_el
    integer :: status(MPI_STATUS_SIZE)

    if (m%gdim .eq. 2) then
       n_sides = 4
       n_nodes = 2
    else
       n_sides = 6
       n_nodes = 4
    end if

    call buffer%init()
        
    ! Build send buffers containing
    ! [el_glb_idx, side number, facet data (global ids of points)]
    do i = 1, m%nelv
       el_glb_idx = i + m%offset_el
       do j = 1, n_sides
          facet = j             ! Adhere to standards...
          if (m%facet_neigh(j, i) .eq. 0) then
             if (n_nodes .eq. 2) then
                call m%elements(i)%e%facet_id(edge, j)                
                call buffer%push(el_glb_idx)
                call buffer%push(facet)
                do k = 1, n_nodes
                   call buffer%push(edge%x(k))
                end do
             else
                call m%elements(i)%e%facet_id(face, j)
                call buffer%push(el_glb_idx)
                call buffer%push(facet)
                do k = 1, n_nodes
                   call buffer%push(face%x(k))
                end do
             end if
          end if
       end do
    end do


    call MPI_Allreduce(buffer%size(), max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    allocate(recv_buffer(max_recv))
    
    do i = 1, pe_size - 1
       src = modulo(pe_rank - i + pe_size, pe_size)
       dst = modulo(pe_rank + i, pe_size)

       call MPI_Sendrecv(buffer%array(), buffer%size(), MPI_INTEGER, dst, 0, &
            recv_buffer, max_recv, MPI_INTEGER, src, 0, NEKO_COMM, status, ierr)

       call MPI_Get_count(status, MPI_INTEGER, n_recv, ierr)

       select type (fmp => m%facet_map)
       type is(htable_i4t2_t)
          do j = 1, n_recv, n_nodes + 2
             neigh_el = recv_buffer(j)
             recv_side = recv_buffer(j+1)

             edge = (/ recv_buffer(j+2), recv_buffer(j+3) /)
             
             facet_key = (/ 0 , 0 /)
             !Check if the face is present on this PE
             if (fmp%get(edge, facet_key) .eq. 0) then
                ! Determine opposite side and update neighbor
                if (mod(recv_side, 2) .eq. 1) then
                   m%facet_neigh(recv_side + 1, &
                        facet_key%x(2) - m%offset_el) = -neigh_el
                else  if (mod(recv_side, 2) .eq. 0) then
                   m%facet_neigh(recv_side - 1, &
                        facet_key%x(1) - m%offset_el) = -neigh_el
                end if
             end if
             
          end do
       type is(htable_i4t4_t)
          do j = 1, n_recv, n_nodes + 2
             neigh_el = recv_buffer(j)
             recv_side = recv_buffer(j+1)

             face = (/ recv_buffer(j+2), recv_buffer(j+3), &
                  recv_buffer(j+4), recv_buffer(j+5) /)
             
             facet_key = (/ 0 , 0 /)
             !Check if the face is present on this PE
             if (fmp%get(face, facet_key) .eq. 0) then
                ! Determine opposite side and update neighbor
                if (mod(recv_side, 2) .eq. 1) then
                   element = facet_key%x(2) - m%offset_el
                   facet = recv_side + 1
                   m%facet_neigh(facet, element) = -neigh_el
                else  if (mod(recv_side, 2) .eq. 0) then
                   element = facet_key%x(1) - m%offset_el
                   facet  = recv_side - 1
                   m%facet_neigh(facet, element) = -neigh_el
                end if
                call distdata_set_shared(m%distdata, element, facet)       
             end if
             
          end do
       end select

    end do

    deallocate(recv_buffer)

    call buffer%free()

  end subroutine mesh_generate_external_facet_conn

  !> Generate element-element connectivity via points between PEs
  subroutine mesh_generate_external_point_conn(m)
    type(mesh_t), intent(inout) :: m
    type(tuple_i4_t) :: edge, facet_key
    type(tuple4_i4_t) :: face
    type(stack_i4_t) :: send_buffer
    integer, allocatable :: recv_buffer(:)
    integer :: i, j, k, el_glb_idx, n_sides, n_nodes, facet, element
    integer :: max_recv, ierr, src, dst, n_recv, recv_side, neigh_el
    integer :: pt_glb_idx, pt_loc_idx, num_neigh
    integer, pointer :: neighs(:)
    integer :: status(MPI_STATUS_SIZE)

    
    call send_buffer%init()
    
    ! Build send buffers containing
    ! [pt_glb_idx, #neigh, neigh id_1 ....neigh_id_n] 
    do i = 1, m%mpts
       pt_glb_idx = m%points(i)%id() ! Adhere to standards...
       num_neigh = m%point_neigh(i)%size()
       call send_buffer%push(pt_glb_idx)
       call send_buffer%push(num_neigh)

       neighs => m%point_neigh(i)%array()
       do j = 1, m%point_neigh(i)%size()
          call send_buffer%push(neighs(j))
       end do
    end do

    call MPI_Allreduce(send_buffer%size(), max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)
    allocate(recv_buffer(max_recv))
       
    do i = 1, pe_size - 1
       src = modulo(pe_rank - i + pe_size, pe_size)
       dst = modulo(pe_rank + i, pe_size)

       call MPI_Sendrecv(send_buffer%array(), send_buffer%size(), &
            MPI_INTEGER, dst, 0, recv_buffer, max_recv, MPI_INTEGER, src, 0, &
            NEKO_COMM, status, ierr)

       call MPI_Get_count(status, MPI_INTEGER, n_recv, ierr)

       j = 1
       do while (j .le. n_recv)
          pt_glb_idx = recv_buffer(j)
          num_neigh = recv_buffer(j + 1)
          ! Check if the point is present on this PE
          pt_loc_idx = mesh_have_point_glb_idx(m, pt_glb_idx)
          if (pt_loc_idx .gt. 0) then
             do k = 1, num_neigh
                neigh_el = -recv_buffer(j + 2 + k)
                !> @todo Store shared points in distdata
                call m%point_neigh(pt_loc_idx)%push(neigh_el)
             end do
          end if
          j = j + (2 + num_neigh)          
       end do
       
    end do

    deallocate(recv_buffer)
    call send_buffer%free()
    
  end subroutine mesh_generate_external_point_conn
  
  
  !> Add a quadrilateral element to the mesh @a m
  subroutine mesh_add_quad(m, el, p1, p2, p3, p4)
    type(mesh_t), target, intent(inout) :: m
    integer, value :: el
    type(point_t), intent(inout) :: p1, p2, p3, p4
    class(element_t), pointer :: ep
    integer :: p(4), el_glb_idx, i, p_local_idx

    ! Connectivity invalidated if a new element is added        
    m%lconn = .false.           

    call mesh_add_point(m, p1, p(1))
    call mesh_add_point(m, p2, p(2))
    call mesh_add_point(m, p3, p(3))
    call mesh_add_point(m, p4, p(4))

    ep => m%elements(el)%e
    el_glb_idx = el + m%offset_el

    do i = 1,4
       p_local_idx = mesh_get_local(m, m%points(p(i)))
       call m%point_neigh(p_local_idx)%push(el_glb_idx)
    end do
    
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
    integer :: p(8), el_glb_idx, i, p_local_idx

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

    do i = 1,8
       p_local_idx = mesh_get_local(m, m%points(p(i)))
       call m%point_neigh(p_local_idx)%push(el_glb_idx)
    end do
    
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

  !> Return the local id of a point @a p
  function mesh_get_local_point(m, p) result(local_id)
    type(mesh_t), intent(inout) :: m
    type(point_t), intent(inout) :: p
    integer :: local_id
    integer :: tmp

    !> @todo why do we still need to do this?
    tmp = p%id()

    if (m%htp%get(tmp, local_id) .gt. 0) then
       call neko_error('Invalid global id')
    end if
    
  end function mesh_get_local_point

  !> Check if the mesh has a point given its global index
  !! @return The local id of the point (if present) otherwise -1
  !! @todo Consider moving this to distdata
  function mesh_have_point_glb_idx(m, index) result(local_id)
    type(mesh_t), intent(inout) :: m 
    integer, intent(inout) :: index  !< Global index
    integer :: local_id

    if (m%htp%get(index, local_id) .eq. 1) then
       local_id = -1
    end if
        
  end function mesh_have_point_glb_idx

end module mesh
