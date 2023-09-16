! Copyright (c) 2018-2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Defines a mesh
module mesh
  use num_types, only : rp, dp, i8
 ! use mpi_f08
  use point, only : point_t
  use element, only : element_t
  use hex
  use quad
  use utils, only : neko_error
  use stack, only : stack_i4_t, stack_i8_t, stack_i4t4_t, stack_i4t2_t
  use tuple, only : tuple_i4_t, tuple4_i4_t
  use htable
  use datadist
  use distdata
  use comm    
  use zone, only : zone_t, zone_periodic_t
  use math
  use uset, only : uset_i8_t
  use curve, only : curve_t
  implicit none
  private
  
  integer, public, parameter :: NEKO_MSH_MAX_ZLBLS = 20 !< Max num. zone labels

  type, private :: mesh_element_t
     class(element_t), allocatable :: e
  end type mesh_element_t

  type, public ::  mesh_t
     integer :: nelv            !< Number of elements
     integer :: npts            !< Number of points per element
     integer :: gdim            !< Geometric dimension
     integer :: mpts            !< Number of (unique) points in the mesh
     integer :: mfcs            !< Number of (unique) faces in the mesh
     integer :: meds            !< Number of (unique) edges in the mesh

     integer :: glb_nelv        !< Global number of elements
     integer :: glb_mpts        !< Global number of unique points
     integer :: glb_mfcs        !< Global number of unique faces
     integer :: glb_meds        !< Global number of unique edges
     
     integer :: offset_el       !< Element offset
     integer :: max_pts_id      !< Max local point id
     
     type(point_t), allocatable :: points(:) !< list of points
     type(mesh_element_t), allocatable :: elements(:) !< List of elements
     logical, allocatable :: dfrmd_el(:) !< List of elements
     
     type(htable_i4_t) :: htp   !< Table of unique points (global->local)
     type(htable_i4t4_t) :: htf !< Table of unique faces (facet->local id)
     type(htable_i4t2_t) :: hte !< Table of unique edges (edge->local id)

     integer, allocatable :: facet_neigh(:,:)  !< Facet to neigh. element table

     !> Facet to element's id tuple and the mapping of the
     !! points between lower id element and higher
     !! \f$ t=(low_id element, element with higher global id) \f$
     class(htable_t), allocatable :: facet_map 
     type(stack_i4_t), allocatable :: point_neigh(:) !< Point to neigh. table

     type(distdata_t) :: ddata            !< Mesh distributed data
     logical, allocatable :: neigh(:)     !< Neighbouring ranks
     integer, allocatable :: neigh_order(:) !< Neighbour order

     integer(2), allocatable :: facet_type(:,:) !< Facet type     
     
     type(zone_t) :: wall                 !< Zone of wall facets
     type(zone_t) :: inlet                !< Zone of inlet facets
     type(zone_t) :: outlet               !< Zone of outlet facets
     type(zone_t) :: outlet_normal        !< Zone of outlet normal facets
     type(zone_t) :: sympln               !< Zone of symmetry plane facets
     type(zone_t), allocatable :: labeled_zones(:) !< Zones with labeled facets
     type(zone_periodic_t) :: periodic             !< Zones with periodic facets
     type(curve_t) :: curve                        !< Set of curved elements

     logical :: lconn = .false.                !< valid connectivity
     logical :: ldist = .false.                !< valid distributed data
     logical :: lnumr = .false.                !< valid numbering
     logical :: lgenc = .true.                 !< generate connectivity

     !> enables user to specify a deformation
     !! that is applied to all x,y,z coordinates generated with this mesh
     procedure(mesh_deform), pass(msh), pointer  :: apply_deform => null()
   contains
     procedure, private, pass(m) :: init_nelv => mesh_init_nelv
     procedure, private, pass(m) :: init_dist => mesh_init_dist
     procedure, private, pass(m) :: add_quad => mesh_add_quad
     procedure, private, pass(m) :: add_hex => mesh_add_hex
     procedure, private, pass(m) :: get_local_point => mesh_get_local_point
     procedure, private, pass(m) :: get_local_edge => mesh_get_local_edge
     procedure, private, pass(m) :: get_local_facet => mesh_get_local_facet
     procedure, private, pass(m) :: get_global_edge => mesh_get_global_edge
     procedure, private, pass(m) :: get_global_facet => mesh_get_global_facet
     procedure, private, pass(m) :: is_shared_point => mesh_is_shared_point
     procedure, private, pass(m) :: is_shared_edge => mesh_is_shared_edge
     procedure, private, pass(m) :: is_shared_facet => mesh_is_shared_facet
     procedure, pass(m) :: free => mesh_free
     procedure, pass(m) :: finalize => mesh_finalize
     procedure, pass(m) :: mark_wall_facet => mesh_mark_wall_facet
     procedure, pass(m) :: mark_inlet_facet => mesh_mark_inlet_facet
     procedure, pass(m) :: mark_outlet_facet => mesh_mark_outlet_facet
     procedure, pass(m) :: mark_sympln_facet => mesh_mark_sympln_facet
     procedure, pass(m) :: mark_periodic_facet => mesh_mark_periodic_facet
     procedure, pass(m) :: mark_outlet_normal_facet => &
          mesh_mark_outlet_normal_facet
     procedure, pass(m) :: mark_labeled_facet => mesh_mark_labeled_facet
     procedure, pass(m) :: mark_curve_element => mesh_mark_curve_element
     procedure, pass(m) :: apply_periodic_facet => mesh_apply_periodic_facet
     procedure, pass(m) :: all_deformed => mesh_all_deformed
     procedure, pass(m) :: get_facet_ids => mesh_get_facet_ids
     procedure, pass(m) :: reset_periodic_ids => mesh_reset_periodic_ids
     procedure, pass(m) :: create_periodic_ids => mesh_create_periodic_ids
     procedure, pass(m) :: generate_conn => mesh_generate_conn
     !> Initialise a mesh
     generic :: init => init_nelv, init_dist
     !> Add an element to the mesh
     generic :: add_element => add_quad, add_hex
     !> Get local id for a mesh entity
     !! @todo Add similar mappings for element ids
     generic :: get_local => get_local_point, get_local_edge, get_local_facet
     !> Get global id for a mesh entity
     !! @todo Add similar mappings for element ids
     generic :: get_global => get_global_edge, get_global_facet
     !> Check if a mesh entity is shared
     generic :: is_shared => is_shared_point, is_shared_edge, is_shared_facet
  end type mesh_t

  abstract interface
     subroutine mesh_deform(msh, x, y, z, lx, ly, lz)
       import mesh_t       
       import rp
       class(mesh_t) :: msh
       integer, intent(in) :: lx, ly, lz
       real(kind=rp), intent(inout) :: x(lx, ly, lz, msh%nelv)
       real(kind=rp), intent(inout) :: y(lx, ly, lz, msh%nelv)
       real(kind=rp), intent(inout) :: z(lx, ly, lz, msh%nelv)
     end subroutine mesh_deform
  end interface

  public :: mesh_deform
  
contains 

  !> Initialise a mesh @a m with @a nelv elements
  subroutine mesh_init_nelv(m, gdim, nelv)
    class(mesh_t), intent(inout) :: m !< Mesh
    integer, intent(in) :: gdim       !< Geometric dimension
    integer, intent(in) :: nelv       !< Local number of elements
    integer :: ierr
    
    call m%free()

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
    class(mesh_t), intent(inout) :: m       !< Mesh
    integer, intent(in) :: gdim             !< Geometric dimension
    type(linear_dist_t), intent(in) :: dist !< Data distribution

    call m%free()
    
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

    m%max_pts_id = 0
    
    allocate(m%elements(m%nelv))
    allocate(m%dfrmd_el(m%nelv))
    if (m%gdim .eq. 3) then
       do i = 1, m%nelv
          allocate(hex_t::m%elements(i)%e)
       end do
       m%npts = NEKO_HEX_NPTS
       !> Only intialize if we generate connectivity
       if (m%lgenc) then
          allocate(htable_i4t4_t::m%facet_map)
          select type (fmp => m%facet_map)
          type is(htable_i4t4_t)
             call fmp%init(m%nelv, facet_data)
          end select

          allocate(m%facet_neigh(NEKO_HEX_NFCS, m%nelv))

          call m%htf%init(m%nelv * NEKO_HEX_NFCS, i)
          call m%hte%init(m%nelv * NEKO_HEX_NEDS, i)
       end if
    else if (m%gdim .eq. 2) then
       do i = 1, m%nelv
          allocate(quad_t::m%elements(i)%e)
       end do
       m%npts = NEKO_QUAD_NPTS
       if (m%lgenc) then
          allocate(htable_i4t2_t::m%facet_map)       
          select type (fmp => m%facet_map)
          type is(htable_i4t2_t)
             call fmp%init(m%nelv, facet_data)
          end select

          allocate(m%facet_neigh(NEKO_QUAD_NEDS, m%nelv))

          call m%hte%init(m%nelv * NEKO_QUAD_NEDS, i)
       end if
    else
       call neko_error("Invalid dimension")
    end if

    !> @todo resize onces final size is known
    allocate(m%points(m%npts*m%nelv))

    !> @todo resize onces final size is known
    !! Only init if we generate connectivity
    if (m%lgenc) then
       allocate(m%point_neigh(m%gdim*m%npts*m%nelv))
       do i = 1, m%gdim*m%npts*m%nelv
          call m%point_neigh(i)%init()
       end do
    end if

    allocate(m%facet_type(2 * m%gdim, m%nelv))
    m%facet_type = 0
    
    call m%htp%init(m%npts*m%nelv, i)

    call m%wall%init(m%nelv)
    call m%inlet%init(m%nelv)
    call m%outlet%init(m%nelv)
    call m%outlet_normal%init(m%nelv)
    call m%sympln%init(m%nelv)
    call m%periodic%init(m%nelv)

    allocate(m%labeled_zones(NEKO_MSH_MAX_ZLBLS))
    do i = 1, NEKO_MSH_MAX_ZLBLS
       call m%labeled_zones(i)%init(m%nelv)
    end do

    call m%curve%init(m%nelv)
   
    call distdata_init(m%ddata)

    allocate(m%neigh(0:pe_size-1))
    m%neigh = .false.
    
    m%mpts = 0
    m%mfcs = 0
    m%meds = 0

  end subroutine mesh_init_common
  
  !> Deallocate a mesh @a m
  subroutine mesh_free(m)
    class(mesh_t), intent(inout) :: m
    integer :: i
    
    call m%htp%free()
    call m%htf%free()
    call m%hte%free()
    call distdata_free(m%ddata)
    
    if (allocated(m%points)) then
       deallocate(m%points)
    end if
    if (allocated(m%dfrmd_el)) then
       deallocate(m%dfrmd_el)
    end if

    if (allocated(m%elements)) then
       do i = 1, m%nelv
          call m%elements(i)%e%free()
          deallocate(m%elements(i)%e)
       end do
       deallocate(m%elements)
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
       do i = 1, m%gdim * m%npts * m%nelv
          call m%point_neigh(i)%free()
       end do
       deallocate(m%point_neigh)
    end if

    if (allocated(m%facet_type)) then
       deallocate(m%facet_type)
    end if
    if (allocated(m%labeled_zones)) then
       do i = 1, NEKO_MSH_MAX_ZLBLS
          call m%labeled_zones(i)%free()
       end do
       deallocate(m%labeled_zones)
    end if

    if (allocated(m%neigh)) then
       deallocate(m%neigh)
    end if

    if (allocated(m%neigh_order)) then
       deallocate(m%neigh_order)
    end if

    call m%wall%free()
    call m%inlet%free()
    call m%outlet%free()
    call m%outlet_normal%free()
    call m%sympln%free()
    call m%periodic%free()
    
  end subroutine mesh_free

  subroutine mesh_finalize(m)
    class(mesh_t), target, intent(inout) :: m
    integer :: i


    call mesh_generate_flags(m)
    call mesh_generate_conn(m)

    call m%wall%finalize()
    call m%inlet%finalize()
    call m%outlet%finalize()
    call m%outlet_normal%finalize()
    call m%sympln%finalize()
    call m%periodic%finalize()
    do i = 1, NEKO_MSH_MAX_ZLBLS
       call m%labeled_zones(i)%finalize()
    end do
    call m%curve%finalize()

  end subroutine mesh_finalize

  subroutine mesh_generate_flags(m)
    type(mesh_t), intent(inout) :: m
    real(kind=dp) :: u(3), v(3), w(3), temp
    integer :: e

    do e = 1, m%nelv
       if (m%gdim .eq. 2) then
          m%dfrmd_el(e) = .false.
          u = m%elements(e)%e%pts(2)%p%x - m%elements(e)%e%pts(1)%p%x
          v = m%elements(e)%e%pts(3)%p%x - m%elements(e)%e%pts(1)%p%x
          temp = u(1)*v(1) + u(2)*v(2)
          if(.not. abscmp(temp, 0d0)) m%dfrmd_el(e) = .true.
       else
          m%dfrmd_el(e) = .false.
          u = m%elements(e)%e%pts(2)%p%x - m%elements(e)%e%pts(1)%p%x
          v = m%elements(e)%e%pts(3)%p%x - m%elements(e)%e%pts(1)%p%x
          w = m%elements(e)%e%pts(5)%p%x - m%elements(e)%e%pts(1)%p%x
          temp = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
          if(.not. abscmp(temp, 0d0)) m%dfrmd_el(e) = .true.
          temp = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
          if(.not. abscmp(temp, 0d0)) m%dfrmd_el(e) = .true.
          u = m%elements(e)%e%pts(7)%p%x - m%elements(e)%e%pts(8)%p%x
          v = m%elements(e)%e%pts(6)%p%x - m%elements(e)%e%pts(8)%p%x
          w = m%elements(e)%e%pts(4)%p%x - m%elements(e)%e%pts(8)%p%x
          temp = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
          if(.not. abscmp(temp, 0d0)) m%dfrmd_el(e) = .true.
          temp = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
          if(.not. abscmp(temp, 0d0)) m%dfrmd_el(e) = .true.
       end if
    end do
  end subroutine mesh_generate_flags

  !> Set all elements as if they are deformed
  subroutine mesh_all_deformed(m)
    class(mesh_t), intent(inout) :: m
    m%dfrmd_el = .true.
  end subroutine mesh_all_deformed
  
  !> Generate element-to-element connectivity
  subroutine mesh_generate_conn(m)
    class(mesh_t), target, intent(inout) :: m
    type(tuple_i4_t) :: edge
    type(tuple4_i4_t) :: face, face_comp
    type(tuple_i4_t) :: facet_data
    type(stack_i4_t) :: neigh_order

    integer :: i, j, k, ierr, el_glb_idx, n_sides, n_nodes, src, dst

    if (m%lconn) return

    if (.not. m%lgenc) return

    if (m%gdim .eq. 2) then
       n_sides = 4
       n_nodes = 2
    else
       n_sides = 6
       n_nodes = 4
    end if

    ! Compute global number of unique points
    call MPI_Allreduce(m%max_pts_id, m%glb_mpts, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

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
                facet_data%x = (/  0, 0/)
                
                !check it this face has shown up earlier
                if (fmp%get(edge, facet_data) .eq. 0) then
                  !if element is already recognized on face
                  if (facet_data%x(1) .eq. el_glb_idx ) then
                     m%facet_neigh(j, i) = facet_data%x(2)
                  else if( facet_data%x(2) .eq. el_glb_idx) then
                     m%facet_neigh(j, i) = facet_data%x(1)
                  !if this is the second element, arrange so low id is first
                  else if(facet_data%x(1) .gt. el_glb_idx) then 
                    facet_data%x(2) = facet_data%x(1)
                    facet_data%x(1) = el_glb_idx
                    m%facet_neigh(j, i) = facet_data%x(2)
                    call fmp%set(edge, facet_data)                             
                  else if(facet_data%x(1) .lt. el_glb_idx) then 
                    facet_data%x(2) = el_glb_idx
                    m%facet_neigh(j, i) = facet_data%x(1)
                    call fmp%set(edge, facet_data)
                 end if
                else
                   facet_data%x(1) = el_glb_idx
                   m%facet_neigh(j, i) = facet_data%x(2)
                   call fmp%set(edge, facet_data)               
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
               
                facet_data%x = (/ 0, 0/)

                !check it this face has shown up earlier
                if (fmp%get(face, facet_data) .eq. 0) then
                  !if element is already recognized on face
                  if (facet_data%x(1) .eq. el_glb_idx ) then
                     m%facet_neigh(j, i) = facet_data%x(2)
                     call m%elements(i)%e%facet_id(face_comp, j+(2*mod(j,2)-1)) 
                     if (face_comp .eq. face) then
                       facet_data%x(2) = el_glb_idx
                       m%facet_neigh(j, i) = facet_data%x(1)
                       call fmp%set(face, facet_data)
                    end if
                  else if( facet_data%x(2) .eq. el_glb_idx) then
                    m%facet_neigh(j, i) = facet_data%x(1)
                  !if this is the second element, arrange so low id is first
                  else if(facet_data%x(1) .gt. el_glb_idx) then 
                    facet_data%x(2) = facet_data%x(1)
                    facet_data%x(1) = el_glb_idx
                    m%facet_neigh(j, i) = facet_data%x(2)
                    call fmp%set(face, facet_data)                             
                  else if(facet_data%x(1) .lt. el_glb_idx) then 
                    facet_data%x(2) = el_glb_idx
                    m%facet_neigh(j, i) = facet_data%x(1)
                    call fmp%set(face, facet_data)
                 end if
                else
                   facet_data%x(1) = el_glb_idx
                   m%facet_neigh(j, i) = 0
                   call fmp%set(face, facet_data)               
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
       
       call mesh_generate_external_point_conn(m)

       !
       ! Generate neighbour exchange order
       !
       call neigh_order%init(pe_size)
              
       do i = 1, pe_size - 1
          src = modulo(pe_rank - i + pe_size, pe_size)
          dst = modulo(pe_rank + i, pe_size)
          if (m%neigh(src) .or. m%neigh(dst)) then
             j = i ! adhere to standards...
             call neigh_order%push(j)      
          end if
       end do

       allocate(m%neigh_order(neigh_order%size()))
       select type(order => neigh_order%data)
       type is (integer)
          do i = 1, neigh_order%size()
             m%neigh_order(i) = order(i)
          end do
       end select
       call neigh_order%free()
       
       call mesh_generate_external_facet_conn(m)
    else
       allocate(m%neigh_order(1))
       m%neigh_order = 1
    end if

    !
    ! Find all internal/extenral edge connections
    ! (Note it needs to be called after external point connections has
    ! been established)
    !
    if (m%gdim .eq. 3) then
       call mesh_generate_edge_conn(m)
    end if
    

    call mesh_generate_facet_numbering(m)
    
    m%lconn = .true.
    
  end subroutine mesh_generate_conn
 
  !> Generate element-element connectivity via facets between PEs
  subroutine mesh_generate_external_facet_conn(m)
    type(mesh_t), intent(inout) :: m
    type(tuple_i4_t) :: edge, edge2
    type(tuple4_i4_t) :: face,  face2
    type(tuple_i4_t) :: facet_data
    type(stack_i4_t) :: buffer
    type(MPI_Status) :: status
    type(MPI_Request) :: send_req, recv_req
    integer, allocatable :: recv_buffer(:)
    integer :: i, j, k, el_glb_idx, n_sides, n_nodes, facet, element, l
    integer :: max_recv, ierr, src, dst, n_recv, recv_side, neigh_el


    if (m%gdim .eq. 2) then
       n_sides = 4
       n_nodes = 2
    else
       n_sides = 6
       n_nodes = 4
    end if

    call buffer%init()
        
    ! Build send buffers containing
    ! [el_glb_idx, side number, facet_id (global ids of points)]
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
    
    do i = 1, size(m%neigh_order)
       src = modulo(pe_rank - m%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + m%neigh_order(i), pe_size)

       if (m%neigh(src)) then
          call MPI_Irecv(recv_buffer, max_recv, MPI_INTEGER, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (m%neigh(dst)) then
          call MPI_Isend(buffer%array(), buffer%size(), MPI_INTEGER, &
               dst, 0, NEKO_COMM, send_req, ierr)
       end if

       if (m%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER, n_recv, ierr)

          select type (fmp => m%facet_map)
          type is(htable_i4t2_t)
             do j = 1, n_recv, n_nodes + 2
                neigh_el = recv_buffer(j)
                recv_side = recv_buffer(j+1)

                edge = (/ recv_buffer(j+2), recv_buffer(j+3) /)

                facet_data = (/ 0, 0 /)
                !Check if the face is present on this PE
                if (fmp%get(edge, facet_data) .eq. 0) then
                   element = facet_data%x(1) - m%offset_el
                   !Check which side is connected
                   do l = 1, n_sides
                      call m%elements(element)%e%facet_id(edge2, l)
                      if(edge2 .eq. edge) then
                         facet = l
                         exit
                      end if
                   end do
                   m%facet_neigh(facet, element) = -neigh_el
                   facet_data%x(2) = -neigh_el

                   !  Update facet map
                   call fmp%set(edge, facet_data)

                   call distdata_set_shared_el_facet(m%ddata, element, facet)

                   if (m%hte%get(edge, facet) .eq. 0) then
                      call distdata_set_shared_facet(m%ddata, facet)
                   else
                      call neko_error("Invalid shared edge")
                   end if

                end if

             end do
          type is(htable_i4t4_t)
             do j = 1, n_recv, n_nodes + 2
                neigh_el = recv_buffer(j)
                recv_side = recv_buffer(j+1)

                face%x = (/ recv_buffer(j+2), recv_buffer(j+3), &
                     recv_buffer(j+4), recv_buffer(j+5) /)


                facet_data%x = (/ 0, 0 /)

                !Check if the face is present on this PE
                if (fmp%get(face, facet_data) .eq. 0) then
                   ! Determine opposite side and update neighbor
                   element = facet_data%x(1) - m%offset_el
                   do l = 1, 6
                      call m%elements(element)%e%facet_id(face2, l)
                      if(face2 .eq. face) then
                         facet = l
                         exit
                      end if
                   end do
                   m%facet_neigh(facet, element) = -neigh_el
                   facet_data%x(2) = -neigh_el                   

                   ! Update facet map
                   call fmp%set(face, facet_data)

                   call distdata_set_shared_el_facet(m%ddata, element, facet)

                   if (m%htf%get(face, facet) .eq. 0) then
                      call distdata_set_shared_facet(m%ddata, facet)
                   else
                      call neko_error("Invalid shared face")
                   end if


                end if

             end do
          end select
       end if

       if (m%neigh(dst)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
       end if
       
    end do

       
    deallocate(recv_buffer)

    call buffer%free()

  end subroutine mesh_generate_external_facet_conn

  !> Generate element-element connectivity via points between PEs
  subroutine mesh_generate_external_point_conn(m)
    type(mesh_t), intent(inout) :: m
    type(stack_i4_t) :: send_buffer
    type(MPI_Status) :: status
    integer, allocatable :: recv_buffer(:)
    integer :: i, j, k
    integer :: max_recv, ierr, src, dst, n_recv, neigh_el
    integer :: pt_glb_idx, pt_loc_idx, num_neigh
    integer, pointer :: neighs(:)

    
    call send_buffer%init(m%mpts * 2)
    
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
             m%neigh(src) = .true.
             do k = 1, num_neigh
                neigh_el = -recv_buffer(j + 1 + k)
                call m%point_neigh(pt_loc_idx)%push(neigh_el)
                call distdata_set_shared_point(m%ddata, pt_loc_idx)
             end do
          end if
          j = j + (2 + num_neigh)          
       end do
       
    end do

    deallocate(recv_buffer)
    call send_buffer%free()
    
  end subroutine mesh_generate_external_point_conn

  !> Generate element-element connectivity via edges
  !! both between internal and between PEs
  !! @attention only for elements where facet .ne. edges
  subroutine mesh_generate_edge_conn(m)
    type(mesh_t), target, intent(inout) :: m
    type(htable_iter_i4t2_t) :: it
    type(tuple_i4_t), pointer :: edge
    type(uset_i8_t), target :: edge_idx, ghost, owner
    type(stack_i8_t), target :: send_buff
    type(htable_i8_t) :: glb_to_loc
    type(MPI_Status) :: status
    type(MPI_Request) :: send_req, recv_req
    integer, pointer :: p1(:), p2(:), ns_id(:)
    integer :: i, j, id, ierr, num_edge_glb, edge_offset, num_edge_loc
    integer :: k, l , shared_offset, glb_nshared, n_glb_id
    integer(kind=i8) :: C, glb_max, glb_id
    integer(kind=i8), pointer :: glb_ptr
    integer(kind=i8), allocatable :: recv_buff(:)
    logical :: shared_edge
    type(stack_i4_t), target :: non_shared_edges
    integer :: max_recv, src, dst, n_recv


    !>@todo move this into distdata
    allocate(m%ddata%local_to_global_edge(m%meds))

    call edge_idx%init(m%hte%num_entries())
    call send_buff%init(m%hte%num_entries())
    call owner%init(m%hte%num_entries())

    call glb_to_loc%init(32, i)

    !
    ! Determine/ constants used to generate unique global edge numbers
    ! for shared edges 
    !
    C = int(m%glb_nelv, i8) * int(NEKO_HEX_NEDS, i8)

    num_edge_glb = 2* m%meds
    call MPI_Allreduce(MPI_IN_PLACE, num_edge_glb, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM,  ierr)

    glb_max = int(num_edge_glb, i8)

    call non_shared_edges%init(m%hte%num_entries())

    call it%init(m%hte)
    do while(it%next())       
       edge => it%key()
       call it%data(id)

       k = mesh_have_point_glb_idx(m, edge%x(1))
       l = mesh_have_point_glb_idx(m, edge%x(2))
       p1 => m%point_neigh(k)%array()
       p2 => m%point_neigh(l)%array()

       shared_edge = .false.
       
       ! Find edge neighbor from point neighbors 
       do i = 1, m%point_neigh(k)%size()
          do j = 1, m%point_neigh(l)%size()
             if ((p1(i) .eq. p2(j)) .and. &
                  (p1(i) .lt. 0) .and. (p2(j) .lt. 0)) then
                call distdata_set_shared_edge(m%ddata, id)
                shared_edge = .true.
             end if
          end do
       end do

       ! Generate a unique id for the shared edge as,
       ! ((e1 * C) + e2 )) + glb_max if e1 > e2
       ! ((e2 * C) + e1 )) + glb_max if e2 > e1     
       if (shared_edge) then
          glb_id = ((int(edge%x(1), i8)) + int(edge%x(2), i8)*C) + glb_max
          call glb_to_loc%set(glb_id, id)
          call edge_idx%add(glb_id)
          call owner%add(glb_id) ! Always assume the PE is the owner
          call send_buff%push(glb_id)
       else
          call non_shared_edges%push(id)
       end if
    end do

    ! Determine start offset for global numbering of locally owned edges
    edge_offset = 0
    num_edge_loc = non_shared_edges%size()
    call MPI_Exscan(num_edge_loc, edge_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    edge_offset = edge_offset + 1

    ! Construct global numbering of locally owned edges
    ns_id => non_shared_edges%array()
    do i = 1, non_shared_edges%size()
       call distdata_set_local_to_global_edge(m%ddata, ns_id(i), edge_offset)
       edge_offset = edge_offset + 1          
    end do
    nullify(ns_id)
    
    !
    ! Renumber shared edges into integer range
    !
    
    call MPI_Allreduce(send_buff%size(), max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    call ghost%init(send_buff%size())

    allocate(recv_buff(max_recv))

    do i = 1, size(m%neigh_order)
       src = modulo(pe_rank - m%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + m%neigh_order(i), pe_size)

       if (m%neigh(src)) then
          call MPI_Irecv(recv_buff, max_recv, MPI_INTEGER8, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (m%neigh(dst)) then
          ! We should use the %array() procedure, which works great for
          ! GNU, Intel and NEC, but it breaks horribly on Cray when using
          ! certain data types
          select type(sbarray=>send_buff%data)
          type is (integer(i8))
             call MPI_Isend(sbarray, send_buff%size(), MPI_INTEGER8, &
                  dst, 0, NEKO_COMM, send_req, ierr)
          end select
       end if

       if (m%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER8, n_recv, ierr)
          
          do j = 1, n_recv
             if ((edge_idx%element(recv_buff(j))) .and. (src .lt. pe_rank)) then
                call ghost%add(recv_buff(j))
                call owner%remove(recv_buff(j))
             end if
          end do
       end if

       if (m%neigh(dst)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
       end if
    end do

   
    ! Determine start offset for global numbering of shared edges
    glb_nshared = num_edge_loc
    call MPI_Allreduce(MPI_IN_PLACE, glb_nshared, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    shared_offset = 0
    call MPI_Exscan(owner%size(), shared_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    shared_offset = shared_offset + glb_nshared + 1
    
    ! Renumber locally owned set of shared edges
    call send_buff%clear()
    call owner%iter_init()
    do while (owner%iter_next())
       glb_ptr => owner%iter_value()
       if (glb_to_loc%get(glb_ptr, id) .eq. 0) then
          call distdata_set_local_to_global_edge(m%ddata, id, shared_offset)

          ! Add new number to send data as [old_glb_id new_glb_id] for each edge
          call send_buff%push(glb_ptr)    ! Old glb_id integer*8
          glb_id = int(shared_offset, i8) ! Waste some space here...
          call send_buff%push(glb_id)     ! New glb_id integer*4

          shared_offset = shared_offset + 1
       else
          call neko_error('Invalid edge id')
       end if
    end do
    nullify(glb_ptr)
       
    ! Determine total number of unique edges in the mesh
    ! (This can probably be done in a clever way...)
    m%glb_meds = shared_offset -1 
    call MPI_Allreduce(MPI_IN_PLACE, m%glb_meds, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, IERR)    

    !
    ! Update ghosted edges with new global id
    !

    call MPI_Allreduce(send_buff%size(), max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    deallocate(recv_buff)
    allocate(recv_buff(max_recv))


    do i = 1, size(m%neigh_order)
       src = modulo(pe_rank - m%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + m%neigh_order(i), pe_size)

       if (m%neigh(src)) then
          call MPI_Irecv(recv_buff, max_recv, MPI_INTEGER8, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (m%neigh(dst)) then
          ! We should use the %array() procedure, which works great for
          ! GNU, Intel and NEC, but it breaks horribly on Cray when using
          ! certain data types
          select type(sbarray=>send_buff%data)
          type is (integer(i8))
             call MPI_Isend(sbarray, send_buff%size(), MPI_INTEGER8, &
                  dst, 0, NEKO_COMM, send_req, ierr)
          end select
       end if
       
       if (m%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER8, n_recv, ierr)
          
          do j = 1, n_recv, 2
             if (ghost%element(recv_buff(j))) then
                if (glb_to_loc%get(recv_buff(j), id) .eq. 0) then
                   n_glb_id = int(recv_buff(j + 1 ), 4)
                   call distdata_set_local_to_global_edge(m%ddata, id, n_glb_id)
                else
                   call neko_error('Invalid edge id')
                end if
             end if
          end do
       end if

       if (m%neigh(dst)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
       end if
    end do

    deallocate(recv_buff)
    call glb_to_loc%free()
    call send_buff%free()
    call edge_idx%free()
    call non_shared_edges%free()
    call ghost%free()
    call owner%free()

  end subroutine mesh_generate_edge_conn

  !> Generate a unique facet numbering
  subroutine mesh_generate_facet_numbering(m)
    type(mesh_t), target, intent(inout) :: m
    type(htable_iter_i4t4_t), target :: face_it
    type(htable_iter_i4t2_t), target :: edge_it
    type(tuple4_i4_t), pointer :: face, fd(:)
    type(tuple_i4_t), pointer :: edge, ed(:)
    type(tuple_i4_t) :: facet_data
    type(tuple4_i4_t) :: recv_face
    type(tuple_i4_t) :: recv_edge
    type(stack_i4t4_t) :: face_owner
    type(htable_i4t4_t) :: face_ghost
    type(stack_i4t2_t) :: edge_owner
    type(htable_i4t2_t) :: edge_ghost
    type(stack_i4_t) :: send_buff
    type(MPI_Status) :: status
    type(MPI_Request) :: send_req, recv_req
    integer, allocatable :: recv_buff(:)
    integer :: non_shared_facets, shared_facets, facet_offset    
    integer :: id, glb_nshared, shared_offset, owned_facets
    integer :: i, j, ierr, max_recv, src, dst, n_recv

    shared_facets = m%ddata%shared_facet%size()

    !>@todo move this into distdata
    if (m%gdim .eq. 2) then
       allocate(m%ddata%local_to_global_facet(m%meds))
       call edge_owner%init(m%meds)
       call edge_ghost%init(64, i)       
       non_shared_facets = m%hte%num_entries() - shared_facets
    else
       allocate(m%ddata%local_to_global_facet(m%mfcs))
       call face_owner%init(m%mfcs)
       call face_ghost%init(64, i)       
       non_shared_facets = m%htf%num_entries() - shared_facets
    end if
    
    !> @todo Move this into distdata as a method...
    
    facet_offset = 0
    call MPI_Exscan(non_shared_facets, facet_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    facet_offset = facet_offset + 1
    
    ! Determine ownership of shared facets
    if (m%gdim .eq. 2) then
       call edge_it%init(m%hte)
       do while (edge_it%next())
          call edge_it%data(id)
          edge => edge_it%key()
          if (.not. m%ddata%shared_facet%element(id)) then       
             call distdata_set_local_to_global_facet(m%ddata, &
                  id, facet_offset)
             facet_offset = facet_offset + 1
          else
             select type(fmp => m%facet_map)
             type is(htable_i4t2_t)
                if (fmp%get(edge, facet_data) .eq. 0) then
                   if (facet_data%x(2) .lt. 0) then
                      if (abs(facet_data%x(2)) .lt. (m%offset_el + 1)) then
                         call edge_ghost%set(edge, id)
                      else
                         call edge_owner%push(edge)
                      end if
                   else
                      call neko_error("Invalid edge neigh.")
                   end if
                end if
             end select
          end if
       end do
       owned_facets = edge_owner%size()
    else
       call face_it%init(m%htf)
       do while (face_it%next())
          call face_it%data(id)
          face => face_it%key()
          if (.not. m%ddata%shared_facet%element(id)) then       
             call distdata_set_local_to_global_facet(m%ddata, &
                  id, facet_offset)
             facet_offset = facet_offset + 1
          else
             select type(fmp => m%facet_map)
             type is(htable_i4t4_t)
                if (fmp%get(face, facet_data) .eq. 0) then
                   if (facet_data%x(2) .lt. 0) then
                      if (abs(facet_data%x(2)) .lt. (m%offset_el + 1)) then
                         call face_ghost%set(face, id)
                      else
                         call face_owner%push(face)
                      end if
                   else
                      call neko_error("Invalid face neigh.")
                   end if
                end if
             end select
          end if
       end do
       owned_facets = face_owner%size()
    end if
    
    ! Determine start offset for global numbering of shared facets
    glb_nshared = non_shared_facets
    call MPI_Allreduce(MPI_IN_PLACE, glb_nshared, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    
    shared_offset = 0
    call MPI_Exscan(owned_facets, shared_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    shared_offset = shared_offset + glb_nshared + 1

    if (m%gdim .eq. 2) then
       
       if (owned_facets .gt. 32)  then
          call send_buff%init(owned_facets)
       else
          call send_buff%init()
       end if
       
       ed => edge_owner%array()
       do i = 1, edge_owner%size()
          if (m%hte%get(ed(i), id) .eq. 0) then
             call distdata_set_local_to_global_facet(m%ddata, id, shared_offset)
             
             ! Add new number to send buffer
             ! [edge id1 ... edge idn new_glb_id]
             do j = 1, 2
                call send_buff%push(ed(i)%x(j))
             end do
             call send_buff%push(shared_offset)
             
             shared_offset = shared_offset + 1
          end if
       end do
       
    else
       
       if (owned_facets .gt. 32)  then
          call send_buff%init(owned_facets)
       else
          call send_buff%init()
       end if
           
       fd => face_owner%array()
       do i = 1, face_owner%size()
          if (m%htf%get(fd(i), id) .eq. 0) then
             call distdata_set_local_to_global_facet(m%ddata, id, shared_offset)

             ! Add new number to send buffer
             ! [face id1 ... face idn new_glb_id]
             do j = 1, 4
                call send_buff%push(fd(i)%x(j))
             end do
             call send_buff%push(shared_offset)
             
             shared_offset = shared_offset + 1
          end if
       end do
       nullify(fd)
       
    end if

    ! Determine total number of unique facets in the mesh
    ! (This can probably be done in a clever way...)
    m%glb_mfcs = shared_offset - 1
    call MPI_Allreduce(MPI_IN_PLACE, m%glb_mfcs, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, IERR)    

    !
    ! Update ghosted facets with new global id
    !
    
    call MPI_Allreduce(send_buff%size(), max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    allocate(recv_buff(max_recv))    

    !> @todo Since we now the neigh. we can actually do p2p here...
    do i = 1, size(m%neigh_order)
       src = modulo(pe_rank - m%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + m%neigh_order(i), pe_size)

       if (m%neigh(src)) then
          call MPI_Irecv(recv_buff, max_recv, MPI_INTEGER, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (m%neigh(dst)) then
          call MPI_Isend(send_buff%array(), send_buff%size(), MPI_INTEGER, &
               dst, 0, NEKO_COMM, send_req, ierr)
       end if
       
       if (m%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER, n_recv, ierr)
          
          if (m%gdim .eq. 2) then
             do j = 1, n_recv, 3
             
                recv_edge = (/recv_buff(j), recv_buff(j+1)/)
             
                ! Check if the PE has the shared edge
                if (edge_ghost%get(recv_edge, id) .eq. 0) then
                   call distdata_set_local_to_global_facet(m%ddata, &
                        id, recv_buff(j+2))
                end if
             end do
          else             
             do j = 1, n_recv, 5
                
                recv_face = (/recv_buff(j), recv_buff(j+1), &
                     recv_buff(j+2), recv_buff(j+3) /)
                
                ! Check if the PE has the shared face
                if (face_ghost%get(recv_face, id) .eq. 0) then
                   call distdata_set_local_to_global_facet(m%ddata, &
                        id, recv_buff(j+4))
                end if
             end do
          end if
       end if

       if (m%neigh(dst)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
       end if
       
    end do

    if (m%gdim .eq. 2) then
       call edge_owner%free()
       call edge_ghost%free()
    else
       call face_owner%free()
       call face_ghost%free()
    end if
    
    call send_buff%free()
    deallocate(recv_buff)
       
  end subroutine mesh_generate_facet_numbering
  
  
  !> Add a quadrilateral element to the mesh @a m
  subroutine mesh_add_quad(m, el, p1, p2, p3, p4)
    class(mesh_t), target, intent(inout) :: m
    integer, value :: el
    type(point_t), intent(inout) :: p1, p2, p3, p4
    class(element_t), pointer :: ep
    integer :: p(4), el_glb_idx, i, p_local_idx
    type(tuple_i4_t) :: e

    ! Connectivity invalidated if a new element is added        
    m%lconn = .false.

    ! Numbering invalidated if a new element is added
    m%lnumr = .false.
    
    call mesh_add_point(m, p1, p(1))
    call mesh_add_point(m, p2, p(2))
    call mesh_add_point(m, p3, p(3))
    call mesh_add_point(m, p4, p(4))

    ep => m%elements(el)%e
    el_glb_idx = el + m%offset_el

    do i = 1, NEKO_QUAD_NPTS
       p_local_idx = m%get_local(m%points(p(i)))
       call m%point_neigh(p_local_idx)%push(el_glb_idx)
    end do
    
    select type(ep)
    type is (quad_t)
       call ep%init(el_glb_idx, &
            m%points(p(1)), m%points(p(2)), &
            m%points(p(3)), m%points(p(4)))

       do i = 1, NEKO_QUAD_NEDS
          call ep%facet_id(e, i)
          call mesh_add_edge(m, e)
       end do

    class default
       call neko_error('Invalid element type')
    end select
        
  end subroutine mesh_add_quad

  !> Add a hexahedral element to the mesh @a m
  subroutine mesh_add_hex(m, el, p1, p2, p3, p4, p5, p6, p7, p8)
    class(mesh_t), target, intent(inout) :: m
    integer, value :: el
    type(point_t), intent(inout) :: p1, p2, p3, p4, p5, p6, p7, p8
    class(element_t), pointer :: ep
    integer :: p(8), el_glb_idx, i, p_local_idx
    type(tuple4_i4_t) :: f
    type(tuple_i4_t) :: e

    ! Connectivity invalidated if a new element is added        
    m%lconn = .false.

    ! Numbering invalidated if a new element is added
    m%lnumr = .false.
    
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
    if (m%lgenc) then
       do i = 1, NEKO_HEX_NPTS
          p_local_idx = m%get_local(m%points(p(i)))
          call m%point_neigh(p_local_idx)%push(el_glb_idx)
       end do
    end if
    select type(ep)
    type is (hex_t)
       call ep%init(el_glb_idx, &
            m%points(p(1)), m%points(p(2)), &
            m%points(p(3)), m%points(p(4)), &
            m%points(p(5)), m%points(p(6)), &
            m%points(p(7)), m%points(p(8)))
       
       if (m%lgenc) then
          do i = 1, NEKO_HEX_NFCS
             call ep%facet_id(f, i)
             call mesh_add_face(m, f)
          end do

          do i = 1, NEKO_HEX_NEDS
             call ep%edge_id(e, i)
             call mesh_add_edge(m, e)
          end do
       end if
       
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

    m%max_pts_id = max(m%max_pts_id, tmp)
    
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

  !> Add a unique face represented as a 4-tuple to the mesh
  subroutine mesh_add_face(m, f)
    type(mesh_t), intent(inout) :: m
    type(tuple4_i4_t), intent(inout) :: f
    integer :: idx

    if (m%htf%get(f, idx) .gt. 0) then
       m%mfcs = m%mfcs + 1
       call m%htf%set(f, m%mfcs)
    end if
    
  end subroutine mesh_add_face
  
  !> Add a unique edge represented as a 2-tuple to the mesh
  subroutine mesh_add_edge(m, e)
    type(mesh_t), intent(inout) :: m
    type(tuple_i4_t), intent(inout) :: e
    integer :: idx

    if (m%hte%get(e, idx) .gt. 0) then
       m%meds = m%meds + 1
       call m%hte%set(e, m%meds)
    end if
    
  end subroutine mesh_add_edge

  !> Mark facet @a f in element @a e as a wall
  subroutine mesh_mark_wall_facet(m, f, e)
    class(mesh_t), intent(inout) :: m
    integer, intent(inout) :: f
    integer, intent(inout) :: e

    if (e .gt. m%nelv) then
       call neko_error('Invalid element index')
    end if

    if ((m%gdim .eq. 2 .and. f .gt. 4) .or. &
         (m%gdim .eq. 3 .and. f .gt. 6)) then
       call neko_error('Invalid facet index')
    end if
    m%facet_type(f, e) = 2
    call m%wall%add_facet(f, e)
    
  end subroutine mesh_mark_wall_facet

  !> Mark element @a e as a curve element
  subroutine mesh_mark_curve_element(m, e, curve_data, curve_type)
    class(mesh_t), intent(inout) :: m
    integer, intent(in) :: e
    real(kind=dp), dimension(5,12), intent(in) :: curve_data 
    integer, dimension(12), intent(in) :: curve_type 

    if (e .gt. m%nelv) then
       call neko_error('Invalid element index')
    end if
    if ((m%gdim .eq. 2 .and. sum(curve_type(5:8)) .gt. 0) ) then
       call neko_error('Invalid curve element')
    end if
    call m%curve%add_element(e, curve_data, curve_type)
    
  end subroutine mesh_mark_curve_element


  !> Mark facet @a f in element @a e as an inlet
  subroutine mesh_mark_inlet_facet(m, f, e)
    class(mesh_t), intent(inout) :: m
    integer, intent(in) :: f
    integer, intent(in) :: e

    if (e .gt. m%nelv) then
       call neko_error('Invalid element index')
    end if

    if ((m%gdim .eq. 2 .and. f .gt. 4) .or. &
         (m%gdim .eq. 3 .and. f .gt. 6)) then
       call neko_error('Invalid facet index')
    end if
    m%facet_type(f, e) = 2
    call m%inlet%add_facet(f, e)
    
  end subroutine mesh_mark_inlet_facet
 
  !> Mark facet @a f in element @a e with label 
  subroutine mesh_mark_labeled_facet(m, f, e, label)
    class(mesh_t), intent(inout) :: m
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: label

    if (e .gt. m%nelv) then
       call neko_error('Invalid element index')
    end if

    if ((m%gdim .eq. 2 .and. f .gt. 4) .or. &
         (m%gdim .eq. 3 .and. f .gt. 6)) then
       call neko_error('Invalid facet index')
    end if
    call m%labeled_zones(label)%add_facet(f, e)
    m%facet_type(f,e) = -label
    
  end subroutine mesh_mark_labeled_facet

 
  !> Mark facet @a f in element @a e as an outlet normal
  subroutine mesh_mark_outlet_normal_facet(m, f, e)
    class(mesh_t), intent(inout) :: m
    integer, intent(inout) :: f
    integer, intent(inout) :: e

    if (e .gt. m%nelv) then
       call neko_error('Invalid element index')
    end if

    if ((m%gdim .eq. 2 .and. f .gt. 4) .or. &
         (m%gdim .eq. 3 .and. f .gt. 6)) then
       call neko_error('Invalid facet index')
    end if
    m%facet_type(f, e) = 1
    call m%outlet_normal%add_facet(f, e)
    
  end subroutine mesh_mark_outlet_normal_facet


  !> Mark facet @a f in element @a e as an outlet
  subroutine mesh_mark_outlet_facet(m, f, e)
    class(mesh_t), intent(inout) :: m
    integer, intent(inout) :: f
    integer, intent(inout) :: e

    if (e .gt. m%nelv) then
       call neko_error('Invalid element index')
    end if

    if ((m%gdim .eq. 2 .and. f .gt. 4) .or. &
         (m%gdim .eq. 3 .and. f .gt. 6)) then
       call neko_error('Invalid facet index')
    end if
    m%facet_type(f, e) = 1
    call m%outlet%add_facet(f, e)
    
  end subroutine mesh_mark_outlet_facet

  !> Mark facet @a f in element @a e as a symmetry plane
  subroutine mesh_mark_sympln_facet(m, f, e)
    class(mesh_t), intent(inout) :: m
    integer, intent(inout) :: f
    integer, intent(inout) :: e

    if (e .gt. m%nelv) then
       call neko_error('Invalid element index')
    end if

    if ((m%gdim .eq. 2 .and. f .gt. 4) .or. &
         (m%gdim .eq. 3 .and. f .gt. 6)) then
       call neko_error('Invalid facet index')
    end if
    m%facet_type(f, e) = 2
    call m%sympln%add_facet(f, e)
    
  end subroutine mesh_mark_sympln_facet

  !> Mark facet @a f in element @a e as periodic with (@a pf, @a pe)
  subroutine mesh_mark_periodic_facet(m, f, e, pf, pe, pids)
    class(mesh_t), intent(inout) :: m
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: pf
    integer, intent(in) :: pe
    integer, intent(inout) :: pids(4)
    integer, dimension(4) :: org_ids
    
    call mesh_get_facet_ids(m, f, e, org_ids)
    call m%periodic%add_periodic_facet(f, e, pf, pe, pids, org_ids)
  end subroutine mesh_mark_periodic_facet

  !> Get original ids of periodic points
  subroutine mesh_get_facet_ids(m, f, e, pids)
    class(mesh_t), intent(inout) :: m
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(inout) :: pids(4)
    type(point_t), pointer :: pi
    type(tuple4_i4_t) :: t
    type(tuple_i4_t) :: t2
    integer :: i

    select type(ele => m%elements(e)%e)
    type is(hex_t)
       call ele%facet_order(t,f)
       pids = t%x
    type is(quad_t)
       call ele%facet_order(t2,f)
       pids(1) = t2%x(1)
       pids(2) = t2%x(2)
       pids(3) = 0
       pids(4) = 0
    end select
  end subroutine mesh_get_facet_ids
  
  !> Reset ids of periodic points to their original ids
  subroutine mesh_reset_periodic_ids(m)
    class(mesh_t), intent(inout) :: m
    integer :: i,j
    integer :: f
    integer :: e
    integer :: pf
    integer :: pe
    integer :: org_ids(4), pids(4)
    type(point_t), pointer :: pi
    integer, dimension(4, 6) :: face_nodes = reshape((/1,5,7,3,&
                                                       2,6,8,4,&
                                                       1,2,6,5,&
                                                       3,4,8,7,&
                                                       1,2,4,3,&
                                                       5,6,8,7/),&
                                                       (/4,6/))
    integer, dimension(2, 4) :: edge_nodes = reshape((/1,3,&
                                                                2,4,&
                                                                1,2,&
                                                                3,4 /),&
                                                                (/2,4/))
 
    do i = 1, m%periodic%size
       e = m%periodic%facet_el(i)%x(2) 
       f = m%periodic%facet_el(i)%x(1)
       pe = m%periodic%p_facet_el(i)%x(2)
       pf = m%periodic%p_facet_el(i)%x(1)
       pids = m%periodic%p_ids(i)%x
       call mesh_get_facet_ids(m, f, e, pids)
       m%periodic%p_ids(i)%x = pids
    end do
    do i = 1, m%periodic%size
       e = m%periodic%facet_el(i)%x(2) 
       f = m%periodic%facet_el(i)%x(1)
       org_ids = m%periodic%org_ids(i)%x
       select type(ele => m%elements(e)%e)
       type is(hex_t)
       do j = 1, 4
          pi => ele%pts(face_nodes(j,f))%p
          call pi%set_id(org_ids(j))
       end do
       type is(quad_t)
       do j = 1, 2
          pi => ele%pts(edge_nodes(j,f))%p
          call pi%set_id(org_ids(j))
       end do
       end select
    end do
  end subroutine mesh_reset_periodic_ids
  
  !> Creates common ids for matching periodic points.
  subroutine mesh_create_periodic_ids(m, f, e, pf, pe)
    class(mesh_t), intent(inout) :: m
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: pf
    integer, intent(in) :: pe
    type(point_t), pointer :: pi, pj
    real(kind=dp) :: L(3)
    integer :: i, j, id, p_local_idx
    type(tuple4_i4_t) :: ft
    type(tuple_i4_t) :: et
    integer, dimension(4, 6) :: face_nodes = reshape((/1,5,7,3,&
                                                       2,6,8,4,&
                                                       1,2,6,5,&
                                                       3,4,8,7,&
                                                       1,2,4,3,&
                                                       5,6,8,7/),&
                                                       (/4,6/))
    integer, dimension(2, 4) :: edge_nodes = reshape((/1,3,&
                                                       2,4,&
                                                       1,2,&
                                                       3,4 /),&
                                                      (/2,4/))
  
    select type(ele => m%elements(e)%e)
    type is(hex_t)
    select type(elp => m%elements(pe)%e)
    type is(hex_t)
       L = 0d0
       do i = 1, 4
          L = L + ele%pts(face_nodes(i,f))%p%x(1:3) - &
               elp%pts(face_nodes(i,pf))%p%x(1:3)
       end do
       L = L/4
       do i = 1, 4
          pi => ele%pts(face_nodes(i,f))%p
          do j = 1, 4
             pj => elp%pts(face_nodes(j,pf))%p
             if (norm2(pi%x(1:3) - pj%x(1:3) - L) .lt. 1d-7) then
                id = min(pi%id(), pj%id())
                call pi%set_id(id)
                call pj%set_id(id)
                p_local_idx = m%get_local(m%points(id))
                if (m%lgenc) then
                   id = ele%id()
                   call m%point_neigh(p_local_idx)%push(id)
                   id = elp%id()
                   call m%point_neigh(p_local_idx)%push(id)
                end if
             end if
          end do
       end do

       if (m%lgenc) then
          do i = 1, NEKO_HEX_NFCS
             call ele%facet_id(ft, i)
             call mesh_add_face(m, ft)
             call elp%facet_id(ft, i)
             call mesh_add_face(m, ft)
          end do

          do i = 1, NEKO_HEX_NEDS
             call ele%edge_id(et, i)
             call mesh_add_edge(m, et)
             call elp%edge_id(et, i)
             call mesh_add_edge(m, et)
          end do
       end if
    end select
    type is(quad_t)
    select type(elp => m%elements(pe)%e)
    type is(quad_t)
       L = 0d0
       do i = 1, 2
          L = L + ele%pts(edge_nodes(i,f))%p%x(1:3) - &
               elp%pts(edge_nodes(i,pf))%p%x(1:3)
       end do
       L = L/2
       do i = 1, 2
          pi => ele%pts(edge_nodes(i,f))%p
          do j = 1, 2
             pj => elp%pts(edge_nodes(j,pf))%p
             !whatabout thie tolerance?
             if (norm2(pi%x(1:3) - pj%x(1:3) - L) .lt. 1d-7) then
                id = min(pi%id(), pj%id())
                call pi%set_id(id)
                call pj%set_id(id)
                p_local_idx = m%get_local(m%points(id))
                if (m%lgenc) then
                   id = ele%id()
                   call m%point_neigh(p_local_idx)%push(id)
                   id = elp%id()
                   call m%point_neigh(p_local_idx)%push(id)
                end if
             end if
          end do
       end do
       if (m%lgenc) then
          do i = 1, NEKO_QUAD_NEDS
             call ele%facet_id(et, i)
             call mesh_add_edge(m, et)
             call elp%facet_id(et, i)
             call mesh_add_edge(m, et)
          end do
       end if
    end select
    end select
  end subroutine mesh_create_periodic_ids

  !> Replaces the periodic point's id with a common id for matching
  !! periodic points
  subroutine mesh_apply_periodic_facet(m, f, e, pf, pe, pids)
    class(mesh_t), intent(inout) :: m
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: pf
    integer, intent(in) :: pe
    integer, intent(inout) :: pids(4)
    type(point_t), pointer :: pi
    integer :: i, id, p_local_idx
    type(tuple4_i4_t) :: ft
    type(tuple_i4_t) :: et
    integer, dimension(4, 6) :: face_nodes = reshape((/1,5,7,3,&
                                                       2,6,8,4,&
                                                       1,2,6,5,&
                                                       3,4,8,7,&
                                                       1,2,4,3,&
                                                       5,6,8,7/),&
                                                       (/4,6/))
  
    select type(ele => m%elements(e)%e)
    type is(hex_t)
       do i = 1, 4
          pi => ele%pts(face_nodes(i,f))%p
          call pi%set_id(pids(i))
          call mesh_add_point(m,pi,id)
          p_local_idx = m%get_local(m%points(id))
          id = ele%id()
          if (m%lgenc) then
             call m%point_neigh(p_local_idx)%push(id)
          end if
       end do
       if (m%lgenc) then
          do i = 1, NEKO_HEX_NFCS
             call ele%facet_id(ft, i)
             call mesh_add_face(m, ft)
          end do

          do i = 1, NEKO_HEX_NEDS
             call ele%edge_id(et, i)
             call mesh_add_edge(m, et)
          end do
       end if
    end select

  end subroutine mesh_apply_periodic_facet

  !> Return the local id of a point @a p
  function mesh_get_local_point(m, p) result(local_id)
    class(mesh_t), intent(inout) :: m
    type(point_t), intent(inout) :: p
    integer :: local_id
    integer :: tmp

    !> @todo why do we still need to do this?
    tmp = p%id()

    if (m%htp%get(tmp, local_id) .gt. 0) then
       call neko_error('Invalid global id (local point)')
    end if
    
  end function mesh_get_local_point

  !> Return the local id of an edge @a e
  !! @attention only defined for gdim .ne. 2
  function mesh_get_local_edge(m, e) result(local_id)
    class(mesh_t), intent(inout) :: m
    type(tuple_i4_t), intent(inout) :: e
    integer :: local_id

    if (m%hte%get(e, local_id) .gt. 0) then
       call neko_error('Invalid global id (local edge)')
    end if
    
  end function mesh_get_local_edge

  !> Return the local id of a face @a f
  function mesh_get_local_facet(m, f) result(local_id)
    class(mesh_t), intent(inout) :: m
    type(tuple4_i4_t), intent(inout) :: f
    integer :: local_id

    if (m%htf%get(f, local_id) .gt. 0) then
       call neko_error('Invalid global id (local facet)')
    end if
    
  end function mesh_get_local_facet

  !> Return the global id of an edge @a e
  function mesh_get_global_edge(m, e) result(global_id)
    class(mesh_t), intent(inout) :: m
    type(tuple_i4_t), intent(inout) :: e
    integer :: global_id

    global_id = m%get_local(e)

    if (m%gdim .eq. 2) then
       if (pe_size .gt. 1) then
          global_id = m%ddata%local_to_global_facet(global_id)
       end if
    else 
       if (pe_size .gt. 1) then
          global_id = m%ddata%local_to_global_edge(global_id)
       end if
    end if

  end function mesh_get_global_edge

  !> Return the local id of a face @a f
  function mesh_get_global_facet(m, f) result(global_id)
    class(mesh_t), intent(inout) :: m
    type(tuple4_i4_t), intent(inout) :: f
    integer :: global_id

    global_id = mesh_get_local_facet(m, f)
    
    if (pe_size .gt. 1) then
       global_id = m%ddata%local_to_global_facet(global_id)
    end if
    
  end function mesh_get_global_facet

  
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


  !> Check if a point is shared
  function mesh_is_shared_point(m, p) result(shared)
    class(mesh_t), intent(inout) :: m
    type(point_t), intent(inout) :: p
    integer :: local_index
    logical shared

    local_index = m%get_local(p)
    shared = m%ddata%shared_point%element(local_index)
    
  end function mesh_is_shared_point
  

  !> Check if an edge is shared
  !! @attention only defined for gdim .ne. 2
  function mesh_is_shared_edge(m, e) result(shared)
    class(mesh_t), intent(inout) :: m
    type(tuple_i4_t), intent(inout) :: e
    integer :: local_index
    logical shared
    local_index = m%get_local(e)
    if (m%gdim .eq. 2) then
       shared = m%ddata%shared_facet%element(local_index)
    else 
       shared = m%ddata%shared_edge%element(local_index)
    end if
  end function mesh_is_shared_edge

  !> Check if a facet is shared
  function mesh_is_shared_facet(m, f) result(shared)
    class(mesh_t), intent(inout) :: m
    type(tuple4_i4_t), intent(inout) :: f
    integer :: local_index
    logical shared

    local_index = m%get_local(f)
    shared = m%ddata%shared_facet%element(local_index)
    
  end function mesh_is_shared_facet

end module mesh
