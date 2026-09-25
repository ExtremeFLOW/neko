! Copyright (c) 2018-2026, The Neko Authors
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
  use point, only : point_t
  use element, only : element_t
  use hex, only : hex_t, NEKO_HEX_NEDS, NEKO_HEX_NFCS, &
       NEKO_HEX_NPTS
  use quad, only : quad_t, NEKO_QUAD_NEDS, NEKO_QUAD_NPTS
  use utils, only : neko_error, neko_warning, nonlinear_index
  use mask, only : mask_t
  use stack, only : stack_i4_t, stack_i8_t
  use tuple, only : tuple_i4_t, tuple4_i4_t
  use htable, only : htable_t, htable_i8_t, htable_i4_t, htable_i4t4_t,&
       htable_i4t2_t
  use datadist, only : linear_dist_t
  use distdata, only : distdata_t
  use comm, only : pe_size, pe_rank, NEKO_COMM
  use facet_zone, only : facet_zone_t, facet_zone_periodic_t
  use math, only : abscmp, sort
  use crystal_router, only : crystal_router_transfer, crystal_router_pack
  use mpi_f08, only : MPI_INTEGER, MPI_MAX, MPI_SUM, MPI_IN_PLACE, &
       MPI_Allreduce, MPI_Exscan, MPI_Request, MPI_Status, MPI_Wait, &
       MPI_Issend, MPI_Irecv, MPI_STATUS_IGNORE, MPI_Integer8, &
       MPI_Get_count
  use uset, only : uset_i8_t
  use curve, only : curve_t
  use logger, only : LOG_SIZE
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  private

  !> Max num. zone labels
  integer, public, parameter :: NEKO_MSH_MAX_ZLBLS = 20
  !> Max length of a zone label
  integer, public, parameter :: NEKO_MSH_MAX_ZLBL_LEN = 40

  type, private :: mesh_element_t
     class(element_t), allocatable :: e
  end type mesh_element_t

  type, public :: mesh_t
     integer :: nelv !< Number of elements
     integer :: npts !< Number of points per element
     integer :: gdim !< Geometric dimension
     integer :: mpts !< Number of (unique) points in the mesh
     integer :: mfcs !< Number of (unique) faces in the mesh
     integer :: meds !< Number of (unique) edges in the mesh

     integer :: glb_nelv !< Global number of elements
     integer :: glb_mpts !< Global number of unique points
     integer :: glb_mfcs !< Global number of unique faces
     integer :: glb_meds !< Global number of unique edges

     integer :: offset_el !< Element offset
     integer :: max_pts_id !< Max local point id

     type(point_t), allocatable :: points(:) !< list of points
     type(mesh_element_t), allocatable :: elements(:) !< List of elements
     logical, allocatable :: dfrmd_el(:) !< List of elements

     !> Table of unique points (global->local)
     !! @details Construction-only scratch: it is what deduplicates points as
     !! elements are added, and is released at the end of @a generate_conn,
     !! by which point every local id still needed lives in @a pt_lid.
     !! @a allocated(htp) is therefore the test for "the mesh is still
     !! being built".
     type(htable_i4_t), private, allocatable :: htp
     type(htable_i4_t) :: htel !< Table of unique elements (global->local)

     !> Local point id of each element's points \f$ (point, element) \f$
     !! @details Lets the point->local id hash table be released once the
     !! connectivity has been generated, see @a edge_lid
     integer, allocatable :: pt_lid(:,:)

     !> Local edge id of each element's edges \f$ (edge, element) \f$
     integer, allocatable :: edge_lid(:,:)

     !> Local facet id of each element's facets \f$ (facet, element) \f$
     !! @attention only allocated for gdim .eq. 3, in two dimensions the
     !! facets of an element are its edges
     integer, allocatable :: face_lid(:,:)

     !> Endpoints (global point ids) of each unique local edge \f$ (2, meds) \f$
     !! @details Only valid while the connectivity is generated, released at
     !! the end of @a generate_conn.
     integer, private, allocatable :: edge_pts(:,:)

     !> Points (global point ids) of each unique local face \f$ (4, mfcs) \f$
     !! @details Only valid while the connectivity is generated, released at
     !! the end of @a generate_conn.
     integer, private, allocatable :: face_pts(:,:)


     integer, allocatable :: facet_neigh(:,:) !< Facet to neigh. element table

     !> Facet to element's id tuple and the mapping of the
     !! points between lower id element and higher
     !! \f$ t=(low_id element, element with higher global id) \f$
     class(htable_t), allocatable :: facet_map
     type(stack_i4_t), allocatable :: point_neigh(:) !< Point to neigh. table

     type(distdata_t) :: ddata !< Mesh distributed data
     logical, allocatable :: neigh(:) !< Neighbouring ranks
     integer, allocatable :: neigh_order(:) !< Neighbour order

     integer(2), allocatable :: facet_type(:,:) !< Facet type

     type(facet_zone_t), allocatable :: labeled_zones(:) !< Zones with labeled facets
     type(facet_zone_periodic_t) :: periodic !< Zones with periodic facets
     type(curve_t) :: curve !< Set of curved elements

     logical :: lconn = .false. !< valid connectivity
     logical :: ldist = .false. !< valid distributed data
     logical :: lnumr = .false. !< valid numbering
     logical :: lgenc = .true. !< generate connectivity

     logical :: is_submesh = .false. !< is this mesh a subset of another mesh?

     !> enables user to specify a deformation
     !! that is applied to all x,y,z coordinates generated with this mesh
     procedure(mesh_deform), pass(msh), pointer :: apply_deform => null()
   contains
     procedure, private, pass(this) :: init_nelv => mesh_init_nelv
     procedure, private, pass(this) :: init_dist => mesh_init_dist
     procedure, private, pass(this) :: add_quad => mesh_add_quad
     procedure, private, pass(this) :: add_hex => mesh_add_hex
     procedure, private, pass(this) :: add_point => mesh_add_point
     procedure, pass(this) :: get_global_edge => mesh_get_global_edge
     procedure, pass(this) :: get_global_facet => mesh_get_global_facet
     procedure, pass(this) :: is_shared_point => mesh_is_shared_point
     procedure, pass(this) :: is_shared_edge => mesh_is_shared_edge
     procedure, pass(this) :: is_shared_facet => mesh_is_shared_facet
     procedure, pass(this) :: free => mesh_free
     procedure, pass(this) :: finalize => mesh_finalize
     procedure, pass(this) :: mark_periodic_facet => mesh_mark_periodic_facet
     procedure, pass(this) :: mark_labeled_facet => mesh_mark_labeled_facet
     procedure, pass(this) :: mark_curve_element => mesh_mark_curve_element
     procedure, pass(this) :: apply_periodic_facet => mesh_apply_periodic_facet
     procedure, pass(this) :: all_deformed => mesh_all_deformed
     procedure, pass(this) :: get_facet_ids => mesh_get_facet_ids
     procedure, pass(this) :: reset_periodic_ids => mesh_reset_periodic_ids
     procedure, pass(this) :: create_periodic_ids => mesh_create_periodic_ids
     procedure, pass(this) :: generate_conn => mesh_generate_conn
     procedure, pass(this) :: have_point_glb_idx => mesh_have_point_glb_idx
     procedure, pass(this) :: subset_by_mask => mesh_subset_by_mask

     !> Check the correct orientation of the rst coordindates.
     procedure, pass(this) :: check_right_handedness => &
          mesh_check_right_handedness
     !> Initialise a mesh
     generic :: init => init_nelv, init_dist
     !> Add an element to the mesh
     generic :: add_element => add_quad, add_hex
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

  public :: mesh_deform, parallelepiped_signed_volume


contains

  !> Initialise a mesh @a this with @a nelv elements
  subroutine mesh_init_nelv(this, gdim, nelv)
    class(mesh_t), intent(inout) :: this !< Mesh
    integer, intent(in) :: gdim !< Geometric dimension
    integer, intent(in) :: nelv !< Local number of elements
    integer :: ierr
    logical :: lgenc
    character(len=LOG_SIZE) :: log_buf

    ! Preserve the caller's connectivity-generation choice across free(), which
    ! resets lgenc to its .true. default (see mesh_free). Without this, setting
    ! `msh%lgenc = .false.` before init has no effect, and connectivity is always
    ! generated.
    lgenc = this%lgenc
    call this%free()
    this%lgenc = lgenc

    this%nelv = nelv
    this%gdim = gdim

    if (this%nelv < 1) then
       write(log_buf, '(A,I0,A)') 'MPI rank ', pe_rank, ' has zero elements'
       call neko_warning(log_buf)
    end if

    call MPI_Allreduce(this%nelv, this%glb_nelv, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    this%offset_el = 0
    call MPI_Exscan(this%nelv, this%offset_el, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    call mesh_init_common(this)

  end subroutine mesh_init_nelv

  !> Initialise a mesh @a this based on a distribution @a dist
  subroutine mesh_init_dist(this, gdim, dist)
    class(mesh_t), intent(inout) :: this !< Mesh
    integer, intent(in) :: gdim !< Geometric dimension
    type(linear_dist_t), intent(in) :: dist !< Data distribution
    logical :: lgenc
    character(len=LOG_SIZE) :: log_buf

    ! Preserve the caller's connectivity-generation choice across free()
    ! (see mesh_init_nelv / mesh_free).
    lgenc = this%lgenc
    call this%free()
    this%lgenc = lgenc

    this%nelv = dist%num_local()
    if (this%nelv < 1) then
       write(log_buf, '(A,I0,A)') 'MPI rank ', pe_rank, ' has zero elements'
       call neko_warning(log_buf)
    end if
    this%glb_nelv = dist%num_global()
    this%offset_el = dist%start_idx()
    this%gdim = gdim

    call mesh_init_common(this)

  end subroutine mesh_init_dist

  subroutine mesh_init_common(this)
    type(mesh_t), intent(inout) :: this
    integer :: i
    type(tuple_i4_t) :: facet_data

    this%max_pts_id = 0

    allocate(this%elements(this%nelv))
    allocate(this%dfrmd_el(this%nelv))
    if (this%gdim .eq. 3) then
       do i = 1, this%nelv
          allocate(hex_t::this%elements(i)%e)
       end do
       this%npts = NEKO_HEX_NPTS
       !> Only intialize if we generate connectivity
       if (this%lgenc) then
          allocate(htable_i4t4_t::this%facet_map)
          select type (fmp => this%facet_map)
          type is (htable_i4t4_t)
             call fmp%init(this%nelv, facet_data)
          end select

          allocate(this%facet_neigh(NEKO_HEX_NFCS, this%nelv))

          allocate(this%edge_lid(NEKO_HEX_NEDS, this%nelv))
          allocate(this%face_lid(NEKO_HEX_NFCS, this%nelv))
       end if
    else if (this%gdim .eq. 2) then
       do i = 1, this%nelv
          allocate(quad_t::this%elements(i)%e)
       end do
       this%npts = NEKO_QUAD_NPTS
       if (this%lgenc) then
          allocate(htable_i4t2_t::this%facet_map)
          select type (fmp => this%facet_map)
          type is (htable_i4t2_t)
             call fmp%init(this%nelv, facet_data)
          end select

          allocate(this%facet_neigh(NEKO_QUAD_NEDS, this%nelv))

          allocate(this%edge_lid(NEKO_QUAD_NEDS, this%nelv))
       end if
    else
       call neko_error("Invalid dimension")
    end if

    !> @todo resize onces final size is known
    allocate(this%points(this%npts*this%nelv))

    ! Only init if we generate connectivity; point_neigh is sized and
    ! allocated by generate_conn, once the number of unique points is known
    if (this%lgenc) then
       allocate(this%pt_lid(this%npts, this%nelv))
    end if

    allocate(this%facet_type(2 * this%gdim, this%nelv))
    this%facet_type = 0

    allocate(this%htp)
    call this%htp%init(this%npts*this%nelv, i)
    call this%htel%init(this%nelv, i)

    call this%periodic%init(this%nelv)

    allocate(this%labeled_zones(NEKO_MSH_MAX_ZLBLS))
    do i = 1, NEKO_MSH_MAX_ZLBLS
       call this%labeled_zones(i)%init(this%nelv)
    end do

    call this%curve%init(this%nelv)

    call this%ddata%init()

    allocate(this%neigh(0:pe_size-1))
    this%neigh = .false.

    this%mpts = 0
    this%mfcs = 0
    this%meds = 0

  end subroutine mesh_init_common

  !> Deallocate a mesh @a this
  subroutine mesh_free(this)
    class(mesh_t), intent(inout) :: this
    integer :: i

    if (allocated(this%htp)) then
       call this%htp%free()
       deallocate(this%htp)
    end if
    call this%htel%free()
    call this%ddata%free()
    call this%curve%free()

    if (allocated(this%pt_lid)) then
       deallocate(this%pt_lid)
    end if

    if (allocated(this%edge_lid)) then
       deallocate(this%edge_lid)
    end if

    if (allocated(this%face_lid)) then
       deallocate(this%face_lid)
    end if

    if (allocated(this%edge_pts)) then
       deallocate(this%edge_pts)
    end if

    if (allocated(this%face_pts)) then
       deallocate(this%face_pts)
    end if

    if (allocated(this%dfrmd_el)) then
       deallocate(this%dfrmd_el)
    end if

    if (allocated(this%elements)) then
       do i = 1, this%nelv
          call this%elements(i)%e%free()
          deallocate(this%elements(i)%e)
       end do
       deallocate(this%elements)
    end if

    if (allocated(this%facet_map)) then
       select type (fmp => this%facet_map)
       type is (htable_i4t2_t)
          call fmp%free()
       type is (htable_i4t4_t)
          call fmp%free()
       end select
       deallocate(this%facet_map)
    end if

    if (allocated(this%facet_neigh)) then
       deallocate(this%facet_neigh)
    end if

    if (allocated(this%point_neigh)) then
       do i = 1, size(this%point_neigh)
          call this%point_neigh(i)%free()
       end do
       ! This causes Cray Fortran to take a long vacation
       !deallocate(this%point_neigh)
    end if

    if (allocated(this%facet_type)) then
       deallocate(this%facet_type)
    end if
    if (allocated(this%labeled_zones)) then
       do i = 1, NEKO_MSH_MAX_ZLBLS
          call this%labeled_zones(i)%free()
       end do
       deallocate(this%labeled_zones)
    end if

    if (allocated(this%neigh)) then
       deallocate(this%neigh)
    end if

    if (allocated(this%neigh_order)) then
       deallocate(this%neigh_order)
    end if

    if (allocated(this%points)) then
       deallocate(this%points)
    end if

    call this%periodic%free()
    this%lconn = .false.
    this%lnumr = .false.
    this%ldist = .false.
    this%lgenc = .true.

  end subroutine mesh_free

  subroutine mesh_finalize(this)
    class(mesh_t), target, intent(inout) :: this
    integer :: i

    call mesh_generate_flags(this)
    call mesh_generate_conn(this)

    call this%periodic%finalize()
    do i = 1, NEKO_MSH_MAX_ZLBLS
       call this%labeled_zones(i)%finalize()
    end do
    call this%curve%finalize()

    ! Due to a bug, right handedness check disabled for the time being.
    !call this%check_right_handedness()

  end subroutine mesh_finalize

  subroutine mesh_generate_flags(this)
    type(mesh_t), intent(inout) :: this
    real(kind=dp) :: u(3), v(3), w(3), temp
    integer :: e

    do e = 1, this%nelv
       if (this%gdim .eq. 2) then
          this%dfrmd_el(e) = .false.
          u = this%elements(e)%e%pts(2)%p%x - this%elements(e)%e%pts(1)%p%x
          v = this%elements(e)%e%pts(3)%p%x - this%elements(e)%e%pts(1)%p%x
          temp = u(1)*v(1) + u(2)*v(2)
          if(.not. abscmp(temp, 0d0)) this%dfrmd_el(e) = .true.
       else
          this%dfrmd_el(e) = .false.
          u = this%elements(e)%e%pts(2)%p%x - this%elements(e)%e%pts(1)%p%x
          v = this%elements(e)%e%pts(3)%p%x - this%elements(e)%e%pts(1)%p%x
          w = this%elements(e)%e%pts(5)%p%x - this%elements(e)%e%pts(1)%p%x
          temp = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
          if(.not. abscmp(temp, 0d0)) this%dfrmd_el(e) = .true.
          temp = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
          if(.not. abscmp(temp, 0d0)) this%dfrmd_el(e) = .true.
          u = this%elements(e)%e%pts(7)%p%x - this%elements(e)%e%pts(8)%p%x
          v = this%elements(e)%e%pts(6)%p%x - this%elements(e)%e%pts(8)%p%x
          w = this%elements(e)%e%pts(4)%p%x - this%elements(e)%e%pts(8)%p%x
          temp = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
          if(.not. abscmp(temp, 0d0)) this%dfrmd_el(e) = .true.
          temp = u(1)*w(1) + u(2)*w(2) + u(3)*w(3)
          if(.not. abscmp(temp, 0d0)) this%dfrmd_el(e) = .true.
       end if
    end do
  end subroutine mesh_generate_flags

  !> Set all elements as if they are deformed
  subroutine mesh_all_deformed(this)
    class(mesh_t), intent(inout) :: this
    this%dfrmd_el = .true.
  end subroutine mesh_all_deformed

  !> Generate element-to-element connectivity
  subroutine mesh_generate_conn(this)
    class(mesh_t), target, intent(inout) :: this
    type(tuple_i4_t) :: edge
    type(tuple4_i4_t) :: face, face_comp
    type(tuple_i4_t) :: facet_data
    type(stack_i4_t) :: neigh_order
    class(element_t), pointer :: ep
    integer :: p_local_idx
    integer :: el, id
    integer :: i, j, k, ierr, el_glb_idx, n_sides, n_nodes, src, dst

    if (this%lconn) return

    if (.not. this%lgenc) return

    !If we generate connectivity, we do that here.
    !
    ! Register every point first; add_point hands back the local id directly,
    ! and once all of them are in, mpts is the exact size of point_neigh
    do el = 1, this%nelv
       ep => this%elements(el)%e
       do i = 1, this%npts
          call this%add_point(ep%pts(i)%p, p_local_idx)
          this%pt_lid(i, el) = p_local_idx
       end do
    end do

    ! Temporary workaround to avoid long vacations with Cray Fortran
    if (allocated(this%point_neigh)) then
       deallocate(this%point_neigh)
    end if
    allocate(this%point_neigh(this%mpts))
    do i = 1, this%mpts
       call this%point_neigh(i)%init(size = 4)
    end do

    do el = 1, this%nelv
       !should stack have inout on what we push? would be neat with in
       id = this%elements(el)%e%id()
       do i = 1, this%npts
          call this%point_neigh(this%pt_lid(i, el))%push(id)
       end do
    end do

    !
    ! Enumerate the unique edges and faces of the local mesh
    ! (Note it needs to be called after all points have been added)
    !
    call mesh_generate_edge_lid(this)

    if (this%gdim .eq. 3) then
       call mesh_generate_face_lid(this)
    end if


    if (this%gdim .eq. 2) then
       n_sides = 4
       n_nodes = 2
    else
       n_sides = 6
       n_nodes = 4
    end if

    ! Compute global number of unique points
    call MPI_Allreduce(this%max_pts_id, this%glb_mpts, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    !
    ! Find all (local) boundaries
    !

    !> @note We have to sweep through the facet map twice to make sure
    !! that both odd and even sides are marked
    !! @todo These loop nests needs a lot of love...
    select type (fmp => this%facet_map)
    type is(htable_i4t2_t)
       do k = 1, 2
          do i = 1, this%nelv
             el_glb_idx = i + this%offset_el
             do j = 1, n_sides
                call this%elements(i)%e%facet_id(edge, j)

                ! Assume that all facets are on the exterior
                facet_data%x = [0, 0]

                !check it this face has shown up earlier
                if (fmp%get(edge, facet_data) .eq. 0) then
                   !if element is already recognized on face
                   if (facet_data%x(1) .eq. el_glb_idx ) then
                      this%facet_neigh(j, i) = facet_data%x(2)
                   else if( facet_data%x(2) .eq. el_glb_idx) then
                      this%facet_neigh(j, i) = facet_data%x(1)
                      !if this is the second element, arrange so low id is first
                   else if(facet_data%x(1) .gt. el_glb_idx) then
                      facet_data%x(2) = facet_data%x(1)
                      facet_data%x(1) = el_glb_idx
                      this%facet_neigh(j, i) = facet_data%x(2)
                      call fmp%set(edge, facet_data)
                   else if(facet_data%x(1) .lt. el_glb_idx) then
                      facet_data%x(2) = el_glb_idx
                      this%facet_neigh(j, i) = facet_data%x(1)
                      call fmp%set(edge, facet_data)
                   end if
                else
                   facet_data%x(1) = el_glb_idx
                   this%facet_neigh(j, i) = facet_data%x(2)
                   call fmp%set(edge, facet_data)
                end if
             end do
          end do
       end do
    type is(htable_i4t4_t)

       do k = 1, 2
          do i = 1, this%nelv
             el_glb_idx = i + this%offset_el
             do j = 1, n_sides
                call this%elements(i)%e%facet_id(face, j)

                facet_data%x = (/ 0, 0/)

                !check it this face has shown up earlier
                if (fmp%get(face, facet_data) .eq. 0) then
                   !if element is already recognized on face
                   if (facet_data%x(1) .eq. el_glb_idx ) then
                      this%facet_neigh(j, i) = facet_data%x(2)
                      call this%elements(i)%e%facet_id(face_comp, &
                           j + (2*mod(j, 2) - 1))
                      if (face_comp .eq. face) then
                         facet_data%x(2) = el_glb_idx
                         this%facet_neigh(j, i) = facet_data%x(1)
                         call fmp%set(face, facet_data)
                      end if
                   else if( facet_data%x(2) .eq. el_glb_idx) then
                      this%facet_neigh(j, i) = facet_data%x(1)
                      !if this is the second element, arrange so low id is first
                   else if(facet_data%x(1) .gt. el_glb_idx) then
                      facet_data%x(2) = facet_data%x(1)
                      facet_data%x(1) = el_glb_idx
                      this%facet_neigh(j, i) = facet_data%x(2)
                      call fmp%set(face, facet_data)
                   else if(facet_data%x(1) .lt. el_glb_idx) then
                      facet_data%x(2) = el_glb_idx
                      this%facet_neigh(j, i) = facet_data%x(1)
                      call fmp%set(face, facet_data)
                   end if
                else
                   facet_data%x(1) = el_glb_idx
                   this%facet_neigh(j, i) = 0
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

       call mesh_generate_external_point_conn(this)

       !
       ! Generate neighbour exchange order
       !
       call neigh_order%init(pe_size)

       do i = 1, pe_size - 1
          src = modulo(pe_rank - i + pe_size, pe_size)
          dst = modulo(pe_rank + i, pe_size)
          if (this%neigh(src) .or. this%neigh(dst)) then
             j = i ! adhere to standards...
             call neigh_order%push(j)
          end if
       end do

       allocate(this%neigh_order(neigh_order%size()))
       select type(order => neigh_order%data)
       type is (integer)
          do i = 1, neigh_order%size()
             this%neigh_order(i) = order(i)
          end do
       end select
       call neigh_order%free()

       call mesh_generate_external_facet_conn(this)
    else
       allocate(this%neigh_order(1))
       this%neigh_order = 1
    end if

    !
    ! Find all internal/extenral edge connections
    ! (Note it needs to be called after external point connections has
    ! been established)
    !
    if (this%gdim .eq. 3) then
       call mesh_generate_edge_conn(this)
    end if


    call mesh_generate_facet_numbering(this)

    ! The endpoints are only needed while the numbering is generated, from
    ! here on an edge or a face is addressed by its element and local number
    deallocate(this%edge_pts)
    if (allocated(this%face_pts)) then
       deallocate(this%face_pts)
    end if

    ! The global->local point table has served its purpose; every local id
    ! needed from here on is in pt_lid, edge_lid and face_lid
    call this%htp%free()
    deallocate(this%htp)

    this%lconn = .true.

  end subroutine mesh_generate_conn

  !> Enumerate the unique edges of the local mesh
  !! @details Assigns a local id to each distinct edge, stored per element
  !! and local edge number in @a edge_lid, and collects the endpoints of
  !! every unique edge in @a edge_pts. Edges are deduplicated by chaining
  !! them on the local id of their lowest numbered endpoint; the chains hold
  !! only the edges meeting at a point, so a lookup scans a handful of
  !! entries. Ids are handed out in first seen order, matching the numbering
  !! the previously used hash table produced.
  !! @attention Requires all points to have been added to the mesh
  subroutine mesh_generate_edge_lid(this)
    type(mesh_t), target, intent(inout) :: this
    type(tuple_i4_t) :: e
    class(element_t), pointer :: ep
    integer, allocatable :: chain(:), head(:)
    integer :: el, i, id, n_eds

    if (this%gdim .eq. 3) then
       n_eds = NEKO_HEX_NEDS
    else
       n_eds = NEKO_QUAD_NEDS
    end if

    allocate(this%edge_pts(2, n_eds * this%nelv))
    allocate(chain(n_eds * this%nelv))

    ! Chains are indexed by local point id, bounded by the number of points
    allocate(head(this%npts * this%nelv))
    head = 0

    this%meds = 0
    do el = 1, this%nelv
       ep => this%elements(el)%e
       select type (ep)
       type is (hex_t)
          do i = 1, NEKO_HEX_NEDS
             call ep%edge_id(e, i)
             call mesh_add_edge(this, e, head, chain, id)
             this%edge_lid(i, el) = id
          end do
       type is (quad_t)
          do i = 1, NEKO_QUAD_NEDS
             call ep%facet_id(e, i)
             call mesh_add_edge(this, e, head, chain, id)
             this%edge_lid(i, el) = id
          end do
       end select
    end do

    deallocate(chain)
    deallocate(head)

  end subroutine mesh_generate_edge_lid

  !> Enumerate the unique faces of the local mesh
  !! @details Assigns a local id to each distinct face, stored per element
  !! and local facet number in @a face_lid, and collects the points of every
  !! unique face in @a face_pts. Faces are deduplicated with the same point
  !! chaining as @a mesh_generate_edge_lid, and ids are handed out in first
  !! seen order, matching the numbering the previously used hash table
  !! produced.
  !! @attention Requires all points to have been added to the mesh
  subroutine mesh_generate_face_lid(this)
    type(mesh_t), target, intent(inout) :: this
    type(tuple4_i4_t) :: f
    class(element_t), pointer :: ep
    integer, allocatable :: chain(:), head(:)
    integer :: el, i, id

    allocate(this%face_pts(4, NEKO_HEX_NFCS * this%nelv))
    allocate(chain(NEKO_HEX_NFCS * this%nelv))

    ! Chains are indexed by local point id, bounded by the number of points
    allocate(head(this%npts * this%nelv))
    head = 0

    this%mfcs = 0
    do el = 1, this%nelv
       ep => this%elements(el)%e
       select type (ep)
       type is (hex_t)
          do i = 1, NEKO_HEX_NFCS
             call ep%facet_id(f, i)
             call mesh_add_face(this, f, head, chain, id)
             this%face_lid(i, el) = id
          end do
       end select
    end do

    deallocate(chain)
    deallocate(head)

  end subroutine mesh_generate_face_lid

  !> Generate element-element connectivity via facets between PEs
  subroutine mesh_generate_external_facet_conn(this)
    type(mesh_t), intent(inout) :: this
    type(tuple_i4_t) :: edge, edge2
    type(tuple4_i4_t) :: face, face2
    type(tuple_i4_t) :: facet_data
    type(stack_i4_t) :: buffer
    type(MPI_Status) :: status
    type(MPI_Request) :: send_req, recv_req
    integer, allocatable :: recv_buffer(:)
    integer :: i, j, k, el_glb_idx, n_sides, n_nodes, facet, element, l
    integer :: max_recv, ierr, src, dst, n_recv, recv_side, neigh_el


    if (this%gdim .eq. 2) then
       n_sides = 4
       n_nodes = 2
    else
       n_sides = 6
       n_nodes = 4
    end if

    call buffer%init()

    ! Build send buffers containing
    ! [el_glb_idx, side number, facet_id (global ids of points)]
    do i = 1, this%nelv
       el_glb_idx = i + this%offset_el
       do j = 1, n_sides
          facet = j ! Adhere to standards...
          if (this%facet_neigh(j, i) .eq. 0) then
             if (n_nodes .eq. 2) then
                call this%elements(i)%e%facet_id(edge, j)
                call buffer%push(el_glb_idx)
                call buffer%push(facet)
                do k = 1, n_nodes
                   call buffer%push(edge%x(k))
                end do
             else
                call this%elements(i)%e%facet_id(face, j)
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

    do i = 1, size(this%neigh_order)
       src = modulo(pe_rank - this%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + this%neigh_order(i), pe_size)

       if (this%neigh(src)) then
          call MPI_Irecv(recv_buffer, max_recv, MPI_INTEGER, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (this%neigh(dst)) then
          ! Synchronous send so the buffer is never eager-buffered as an
          ! unexpected message, which exhausts the MPI internal buffer pool
          ! (SIGBUS) under flat MPI at high rank counts. Deadlock-safe: the
          ! matching recv is pre-posted at the top of the peer's iteration.
          call MPI_Issend(buffer%array(), buffer%size(), MPI_INTEGER, &
               dst, 0, NEKO_COMM, send_req, ierr)
       end if

       if (this%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER, n_recv, ierr)

          select type (fmp => this%facet_map)
          type is(htable_i4t2_t)
             do j = 1, n_recv, n_nodes + 2
                neigh_el = recv_buffer(j)
                recv_side = recv_buffer(j+1)

                edge = (/ recv_buffer(j+2), recv_buffer(j+3) /)

                facet_data = (/ 0, 0 /)
                !Check if the face is present on this PE
                if (fmp%get(edge, facet_data) .eq. 0) then
                   element = facet_data%x(1) - this%offset_el
                   !Check which side is connected
                   do l = 1, n_sides
                      call this%elements(element)%e%facet_id(edge2, l)
                      if(edge2 .eq. edge) then
                         facet = l
                         exit
                      end if
                   end do
                   this%facet_neigh(facet, element) = -neigh_el
                   facet_data%x(2) = -neigh_el

                   !  Update facet map
                   call fmp%set(edge, facet_data)

                   call this%ddata%set_shared_el_facet(element, facet)

                   call this%ddata%set_shared_facet( &
                        this%edge_lid(facet, element))

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
                   element = facet_data%x(1) - this%offset_el
                   do l = 1, 6
                      call this%elements(element)%e%facet_id(face2, l)
                      if(face2 .eq. face) then
                         facet = l
                         exit
                      end if
                   end do
                   this%facet_neigh(facet, element) = -neigh_el
                   facet_data%x(2) = -neigh_el

                   ! Update facet map
                   call fmp%set(face, facet_data)

                   call this%ddata%set_shared_el_facet(element, facet)

                   call this%ddata%set_shared_facet( &
                        this%face_lid(facet, element))

                end if

             end do
          end select
       end if

       if (this%neigh(dst)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
       end if

    end do


    deallocate(recv_buffer)

    call buffer%free()

  end subroutine mesh_generate_external_facet_conn

  !> Generate element-element connectivity via points between PEs
  !! @details Uses a canonical-owner rendezvous routed through the crystal
  !! router instead of a dense O(P) all-to-all. Each local point is hashed to
  !! an owner rank, `mod(glb_idx, P)`; the owner gathers every rank holding
  !! that point and reflects, back to each holder, the other holders' element
  !! lists. Shared points (held by more than one rank) are thus discovered in
  !! O(log P) communication stages, the neigh array becomes symmetric by
  !! construction, and no O(P) buffers are ever allocated.
  subroutine mesh_generate_external_point_conn(this)
    type(mesh_t), intent(inout) :: this
    type(stack_i8_t) :: cr_buf
    integer(i8), allocatable :: buf(:), body(:)
    integer(i8), pointer :: cr_data(:)
    integer, allocatable :: gkey(:), gperm(:), rpos(:)
    integer, contiguous, pointer :: neighs(:)
    integer :: i, j, k, n, p, owner, num_neigh, nrec, rlen
    integer :: pt_glb_idx, pt_loc_idx, src_rank, neigh_el, rk, rp

    !
    ! Phase 1: route every local point's element list to its canonical owner.
    !   record payload = [glb_idx, origin, elems...]
    !
    call cr_buf%init(this%mpts * 4)
    allocate(body(8))
    do i = 1, this%mpts
       pt_glb_idx = this%points(i)%id() ! Adhere to standards...
       num_neigh = this%point_neigh(i)%size()
       if (2 + num_neigh .gt. size(body)) then
          deallocate(body)
          allocate(body(2 + num_neigh))
       end if
       body(1) = int(pt_glb_idx, i8) ! glb_idx
       body(2) = int(pe_rank, i8) ! origin
       neighs => this%point_neigh(i)%array()
       do j = 1, num_neigh
          body(2 + j) = int(neighs(j), i8) ! element ids
       end do
       owner = modulo(pt_glb_idx, pe_size)
       call crystal_router_pack(cr_buf, owner, body(1:2 + num_neigh))
    end do
    deallocate(body)

    n = cr_buf%size()
    allocate(buf(max(n, 1)))
    if (n .gt. 0) then
       cr_data => cr_buf%array()
       buf(1:n) = cr_data(1:n)
    end if
    call cr_buf%free()

    call crystal_router_transfer(buf, n)

    !
    ! Phase 2: at the owner, group received records by glb_idx and reflect,
    !   to each holder, the element lists of every *other* holder.
    !   reply = [dest=holder, len=2+num_neigh, glb_idx, src_rank, elems...]
    !
    ! Index the received records and sort their keys (glb_idx) so equal keys
    ! form contiguous runs; per-point holder counts are small, so the
    ! all-pairs reflection within a run is cheap.
    nrec = 0
    p = 1
    do while (p .le. n)
       nrec = nrec + 1
       p = p + 2 + int(buf(p + 1))
    end do

    allocate(gkey(max(nrec, 1)), gperm(max(nrec, 1)), rpos(max(nrec, 1)))
    nrec = 0
    p = 1
    do while (p .le. n)
       nrec = nrec + 1
       rpos(nrec) = p ! record start offset in buf
       gkey(nrec) = int(buf(p + 2)) ! glb_idx
       p = p + 2 + int(buf(p + 1))
    end do
    if (nrec .gt. 0) call sort(gkey, gperm, nrec)

    call cr_buf%init(max(n, 1))
    i = 1
    do while (i .le. nrec)
       ! [i, j) is the run of records sharing the same glb_idx
       j = i
       do while (j .le. nrec)
          if (gkey(j) .ne. gkey(i)) exit
          j = j + 1
       end do
       ! All-pairs reflection within the run (skip singletons = unshared):
       ! send source holder p's record body (glb_idx, origin, elems) to
       ! recipient holder k, addressed to k's origin rank.
       if (j - i .gt. 1) then
          do k = i, j - 1 ! recipient holder
             rk = rpos(gperm(k))
             do p = i, j - 1 ! source holder
                if (p .eq. k) cycle
                rp = rpos(gperm(p))
                call crystal_router_pack(cr_buf, int(buf(rk + 3)), &
                     buf(rp + 2 : rp + 1 + int(buf(rp + 1))))
             end do
          end do
       end if
       i = j
    end do
    deallocate(gkey, gperm, rpos)

    n = cr_buf%size()
    if (allocated(buf)) deallocate(buf)
    allocate(buf(max(n, 1)))
    if (n .gt. 0) then
       cr_data => cr_buf%array()
       buf(1:n) = cr_data(1:n)
    end if
    call cr_buf%free()

    call crystal_router_transfer(buf, n)

    !
    ! Phase 3: finalise locally. Each reply names a remote holder of one of our
    !   points; mark it as a neighbour and absorb its (remote) element list.
    !
    p = 1
    do while (p .le. n)
       rlen = int(buf(p + 1))
       pt_glb_idx = int(buf(p + 2))
       src_rank = int(buf(p + 3))
       pt_loc_idx = this%have_point_glb_idx(pt_glb_idx)
       if (pt_loc_idx .gt. 0) then
          this%neigh(src_rank) = .true.
          call this%ddata%set_shared_point(pt_loc_idx)
          do k = 1, rlen - 2
             neigh_el = -int(buf(p + 3 + k))
             call this%point_neigh(pt_loc_idx)%push(neigh_el)
          end do
       end if
       p = p + 2 + rlen
    end do

    if (allocated(buf)) deallocate(buf)

  end subroutine mesh_generate_external_point_conn

  !> Generate element-element connectivity via edges
  !! both between internal and between PEs
  !! @attention only for elements where facet .ne. edges
  subroutine mesh_generate_edge_conn(this)
    type(mesh_t), target, intent(inout) :: this
    integer, allocatable :: edge_lp(:,:)
    logical, allocatable :: shared_edges(:)
    type(uset_i8_t), target :: edge_idx, ghost, owner
    type(stack_i8_t), target :: send_buff
    type(htable_i8_t) :: glb_to_loc
    type(MPI_Status) :: status
    type(MPI_Request) :: send_req, recv_req
    integer, contiguous, pointer :: p1(:), p2(:), ns_id(:)
    integer :: i, j, id, lid, ierr, num_edge_glb, edge_offset, num_edge_loc
    integer :: k, l , shared_offset, glb_nshared, n_glb_id
    integer(kind=i8) :: C, glb_max, glb_id
    integer(kind=i8), pointer :: glb_ptr
    integer(kind=i8), allocatable :: recv_buff(:)
    type(stack_i4_t), target :: non_shared_edges
    integer :: max_recv, src, dst, n_recv


    !>@todo move this into distdata
    allocate(this%ddata%local_to_global_edge(this%meds))

    call edge_idx%init(this%meds)
    call send_buff%init(this%meds)
    call owner%init(this%meds)

    call glb_to_loc%init(32, i)

    !
    ! Determine/ constants used to generate unique global edge numbers
    ! for shared edges
    !
    C = int(this%glb_nelv, i8) * int(NEKO_HEX_NEDS, i8)

    num_edge_glb = 2* this%meds
    call MPI_Allreduce(MPI_IN_PLACE, num_edge_glb, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    glb_max = int(num_edge_glb, i8)

    call non_shared_edges%init(this%meds)

    ! Resolve both endpoints of every edge to a local point id up front, so
    ! that the neighbour search below reads nothing but plain arrays
    allocate(edge_lp(2, this%meds))
    do lid = 1, this%meds
       id = this%edge_pts(1, lid)
       edge_lp(1, lid) = this%have_point_glb_idx(id)
       id = this%edge_pts(2, lid)
       edge_lp(2, lid) = this%have_point_glb_idx(id)
    end do

    !
    ! An edge is shared when both of its endpoints see the same remote
    ! element. Every iteration writes only its own flag, so this search,
    ! the expensive part of the numbering, is the part that threads
    !
    allocate(shared_edges(this%meds))
    !$omp parallel do private(lid, k, l, p1, p2, i, j)
    do lid = 1, this%meds
       k = edge_lp(1, lid)
       l = edge_lp(2, lid)
       p1 => this%point_neigh(k)%array()
       p2 => this%point_neigh(l)%array()

       shared_edges(lid) = .false.

       ! Find edge neighbor from point neighbors
       do i = 1, this%point_neigh(k)%size()
          do j = 1, this%point_neigh(l)%size()
             if ((p1(i) .eq. p2(j)) .and. &
                  (p1(i) .lt. 0) .and. (p2(j) .lt. 0)) then
                shared_edges(lid) = .true.
             end if
          end do
       end do
    end do
    !$omp end parallel do
    deallocate(edge_lp)

    ! The bookkeeping stays ordered; the order the ids are pushed in is what
    ! the global numbering below is built from
    do lid = 1, this%meds
       id = lid

       ! Generate a unique id for the shared edge as,
       ! ((e1 * C) + e2 )) + glb_max if e1 > e2
       ! ((e2 * C) + e1 )) + glb_max if e2 > e1
       if (shared_edges(id)) then
          call this%ddata%set_shared_edge(id)
          glb_id = ((int(this%edge_pts(1, id), i8)) + &
               int(this%edge_pts(2, id), i8)*C) + glb_max
          call glb_to_loc%set(glb_id, id)
          call edge_idx%add(glb_id)
          call owner%add(glb_id) ! Always assume the PE is the owner
          call send_buff%push(glb_id)
       else
          call non_shared_edges%push(id)
       end if
    end do
    deallocate(shared_edges)

    ! Determine start offset for global numbering of locally owned edges
    edge_offset = 0
    num_edge_loc = non_shared_edges%size()
    call MPI_Exscan(num_edge_loc, edge_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    edge_offset = edge_offset + 1

    ! Construct global numbering of locally owned edges
    ns_id => non_shared_edges%array()
    do i = 1, non_shared_edges%size()
       call this%ddata%set_local_to_global_edge(ns_id(i), edge_offset)
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

    do i = 1, size(this%neigh_order)
       src = modulo(pe_rank - this%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + this%neigh_order(i), pe_size)

       if (this%neigh(src)) then
          call MPI_Irecv(recv_buff, max_recv, MPI_INTEGER8, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (this%neigh(dst)) then
          ! We should use the %array() procedure, which works great for
          ! GNU, Intel and NEC, but it breaks horribly on Cray when using
          ! certain data types
          select type(sbarray=>send_buff%data)
          type is (integer(i8))
             ! Synchronous send to avoid eager-buffering the key list as an
             ! unexpected message (FJMPI buffer-pool exhaustion / SIGBUS under
             ! flat MPI); matching recv is pre-posted by the peer.
             call MPI_Issend(sbarray, send_buff%size(), MPI_INTEGER8, &
                  dst, 0, NEKO_COMM, send_req, ierr)
          end select
       end if

       if (this%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER8, n_recv, ierr)

          do j = 1, n_recv
             if ((edge_idx%element(recv_buff(j))) .and. (src .lt. pe_rank)) then
                call ghost%add(recv_buff(j))
                call owner%remove(recv_buff(j))
             end if
          end do
       end if

       if (this%neigh(dst)) then
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
          call this%ddata%set_local_to_global_edge(id, shared_offset)

          ! Add new number to send data as [old_glb_id new_glb_id] for each edge
          call send_buff%push(glb_ptr) ! Old glb_id integer*8
          glb_id = int(shared_offset, i8) ! Waste some space here...
          call send_buff%push(glb_id) ! New glb_id integer*4

          shared_offset = shared_offset + 1
       else
          call neko_error('Invalid edge id')
       end if
    end do
    nullify(glb_ptr)

    ! Determine total number of unique edges in the mesh
    ! (This can probably be done in a clever way...)
    this%glb_meds = shared_offset -1
    call MPI_Allreduce(MPI_IN_PLACE, this%glb_meds, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, IERR)

    !
    ! Update ghosted edges with new global id
    !

    call MPI_Allreduce(send_buff%size(), max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    deallocate(recv_buff)
    allocate(recv_buff(max_recv))


    do i = 1, size(this%neigh_order)
       src = modulo(pe_rank - this%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + this%neigh_order(i), pe_size)

       if (this%neigh(src)) then
          call MPI_Irecv(recv_buff, max_recv, MPI_INTEGER8, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (this%neigh(dst)) then
          ! We should use the %array() procedure, which works great for
          ! GNU, Intel and NEC, but it breaks horribly on Cray when using
          ! certain data types
          select type(sbarray=>send_buff%data)
          type is (integer(i8))
             ! Synchronous send to avoid eager-buffering the key list as an
             ! unexpected message (FJMPI buffer-pool exhaustion / SIGBUS under
             ! flat MPI); matching recv is pre-posted by the peer.
             call MPI_Issend(sbarray, send_buff%size(), MPI_INTEGER8, &
                  dst, 0, NEKO_COMM, send_req, ierr)
          end select
       end if

       if (this%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER8, n_recv, ierr)

          do j = 1, n_recv, 2
             if (ghost%element(recv_buff(j))) then
                if (glb_to_loc%get(recv_buff(j), id) .eq. 0) then
                   n_glb_id = int(recv_buff(j + 1 ), 4)
                   call this%ddata%set_local_to_global_edge(id, n_glb_id)
                else
                   call neko_error('Invalid edge id')
                end if
             end if
          end do
       end if

       if (this%neigh(dst)) then
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
  subroutine mesh_generate_facet_numbering(this)
    type(mesh_t), target, intent(inout) :: this
    integer, contiguous, pointer :: fd(:), ed(:)
    type(tuple4_i4_t) :: face
    type(tuple_i4_t) :: edge
    type(tuple_i4_t) :: facet_data
    type(tuple4_i4_t) :: recv_face
    type(tuple_i4_t) :: recv_edge
    type(stack_i4_t) :: face_owner
    type(htable_i4t4_t) :: face_ghost
    type(stack_i4_t) :: edge_owner
    type(htable_i4t2_t) :: edge_ghost
    type(stack_i4_t) :: send_buff
    type(MPI_Status) :: status
    type(MPI_Request) :: send_req, recv_req
    integer, allocatable :: recv_buff(:)
    integer :: non_shared_facets, shared_facets, facet_offset
    integer :: id, lid, glb_nshared, shared_offset, owned_facets
    integer :: i, j, ierr, max_recv, src, dst, n_recv

    shared_facets = this%ddata%shared_facet%size()

    !>@todo move this into distdata
    if (this%gdim .eq. 2) then
       allocate(this%ddata%local_to_global_facet(this%meds))
       call edge_owner%init(this%meds)
       call edge_ghost%init(64, i)
       non_shared_facets = this%meds - shared_facets
    else
       allocate(this%ddata%local_to_global_facet(this%mfcs))
       call face_owner%init(this%mfcs)
       call face_ghost%init(64, i)
       non_shared_facets = this%mfcs - shared_facets
    end if

    !> @todo Move this into distdata as a method...

    facet_offset = 0
    call MPI_Exscan(non_shared_facets, facet_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    facet_offset = facet_offset + 1

    ! Determine ownership of shared facets
    if (this%gdim .eq. 2) then
       do lid = 1, this%meds
          id = lid
          if (.not. this%ddata%shared_facet%element(id)) then
             call this%ddata%set_local_to_global_facet(id, facet_offset)
             facet_offset = facet_offset + 1
          else
             edge%x = this%edge_pts(:, id)
             select type(fmp => this%facet_map)
             type is(htable_i4t2_t)
                if (fmp%get(edge, facet_data) .eq. 0) then
                   if (facet_data%x(2) .lt. 0) then
                      if (abs(facet_data%x(2)) .lt. (this%offset_el + 1)) then
                         call edge_ghost%set(edge, id)
                      else
                         call edge_owner%push(id)
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
       do lid = 1, this%mfcs
          id = lid
          if (.not. this%ddata%shared_facet%element(id)) then
             call this%ddata%set_local_to_global_facet(id, facet_offset)
             facet_offset = facet_offset + 1
          else
             face%x = this%face_pts(:, id)
             select type(fmp => this%facet_map)
             type is(htable_i4t4_t)
                if (fmp%get(face, facet_data) .eq. 0) then
                   if (facet_data%x(2) .lt. 0) then
                      if (abs(facet_data%x(2)) .lt. (this%offset_el + 1)) then
                         call face_ghost%set(face, id)
                      else
                         call face_owner%push(id)
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

    if (this%gdim .eq. 2) then

       if (owned_facets .gt. 32) then
          call send_buff%init(owned_facets)
       else
          call send_buff%init()
       end if

       ed => edge_owner%array()
       do i = 1, edge_owner%size()
          id = ed(i)
          call this%ddata%set_local_to_global_facet(id, shared_offset)

          ! Add new number to send buffer
          ! [edge id1 ... edge idn new_glb_id]
          do j = 1, 2
             call send_buff%push(this%edge_pts(j, id))
          end do
          call send_buff%push(shared_offset)

          shared_offset = shared_offset + 1
       end do
       nullify(ed)

    else

       if (owned_facets .gt. 32) then
          call send_buff%init(owned_facets)
       else
          call send_buff%init()
       end if

       fd => face_owner%array()
       do i = 1, face_owner%size()
          id = fd(i)
          call this%ddata%set_local_to_global_facet(id, shared_offset)

          ! Add new number to send buffer
          ! [face id1 ... face idn new_glb_id]
          do j = 1, 4
             call send_buff%push(this%face_pts(j, id))
          end do
          call send_buff%push(shared_offset)

          shared_offset = shared_offset + 1
       end do
       nullify(fd)

    end if

    ! Determine total number of unique facets in the mesh
    ! (This can probably be done in a clever way...)
    this%glb_mfcs = shared_offset - 1
    call MPI_Allreduce(MPI_IN_PLACE, this%glb_mfcs, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, IERR)

    !
    ! Update ghosted facets with new global id
    !

    call MPI_Allreduce(send_buff%size(), max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    allocate(recv_buff(max_recv))

    !> @todo Since we now the neigh. we can actually do p2p here...
    do i = 1, size(this%neigh_order)
       src = modulo(pe_rank - this%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + this%neigh_order(i), pe_size)

       if (this%neigh(src)) then
          call MPI_Irecv(recv_buff, max_recv, MPI_INTEGER, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (this%neigh(dst)) then
          ! Synchronous send to avoid eager-buffered unexpected messages
          ! (FJMPI buffer-pool exhaustion / SIGBUS under flat MPI); the
          ! matching recv is pre-posted at the top of the peer's iteration.
          call MPI_Issend(send_buff%array(), send_buff%size(), MPI_INTEGER, &
               dst, 0, NEKO_COMM, send_req, ierr)
       end if

       if (this%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER, n_recv, ierr)

          if (this%gdim .eq. 2) then
             do j = 1, n_recv, 3

                recv_edge = (/recv_buff(j), recv_buff(j+1)/)

                ! Check if the PE has the shared edge
                if (edge_ghost%get(recv_edge, id) .eq. 0) then
                   call this%ddata%set_local_to_global_facet(id, recv_buff(j+2))
                end if
             end do
          else
             do j = 1, n_recv, 5

                recv_face = (/recv_buff(j), recv_buff(j+1), &
                     recv_buff(j+2), recv_buff(j+3) /)

                ! Check if the PE has the shared face
                if (face_ghost%get(recv_face, id) .eq. 0) then
                   call this%ddata%set_local_to_global_facet(id, recv_buff(j+4))
                end if
             end do
          end if
       end if

       if (this%neigh(dst)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
       end if

    end do

    if (this%gdim .eq. 2) then
       call edge_owner%free()
       call edge_ghost%free()
    else
       call face_owner%free()
       call face_ghost%free()
    end if

    call send_buff%free()
    deallocate(recv_buff)

  end subroutine mesh_generate_facet_numbering


  !> Add a quadrilateral element to the mesh @a this
  subroutine mesh_add_quad(this, el, el_glb, p1, p2, p3, p4)
    class(mesh_t), target, intent(inout) :: this
    integer, value :: el, el_glb
    type(point_t), target, intent(inout) :: p1, p2, p3, p4
    integer :: p(4)
    type(tuple_i4_t) :: e

    ! Connectivity invalidated if a new element is added
    this%lconn = .false.

    ! Numbering invalidated if a new element is added
    this%lnumr = .false.

    call this%add_point(p1, p(1))
    call this%add_point(p2, p(2))
    call this%add_point(p3, p(3))
    call this%add_point(p4, p(4))

    select type (ep => this%elements(el)%e)
    type is (quad_t)
       call ep%init(el_glb, &
            this%points(p(1)), this%points(p(2)), &
            this%points(p(3)), this%points(p(4)))


    class default
       call neko_error('Invalid element type')
    end select

  end subroutine mesh_add_quad

  !> Add a hexahedral element to the mesh @a this
  subroutine mesh_add_hex(this, el, el_glb, p1, p2, p3, p4, p5, p6, p7, p8)
    class(mesh_t), target, intent(inout) :: this
    integer, value :: el, el_glb
    type(point_t), target, intent(inout) :: p1, p2, p3, p4, p5, p6, p7, p8
    integer :: p(8)
    type(tuple4_i4_t) :: f
    type(tuple_i4_t) :: e

    ! Connectivity invalidated if a new element is added
    this%lconn = .false.

    ! Numbering invalidated if a new element is added
    this%lnumr = .false.

    call this%add_point(p1, p(1))
    call this%add_point(p2, p(2))
    call this%add_point(p3, p(3))
    call this%add_point(p4, p(4))
    call this%add_point(p5, p(5))
    call this%add_point(p6, p(6))
    call this%add_point(p7, p(7))
    call this%add_point(p8, p(8))

    ! Global to local mapping
    call this%htel%set(el_glb, el)

    select type (ep => this%elements(el)%e)
    type is (hex_t)
       call ep%init(el_glb, &
            this%points(p(1)), this%points(p(2)), &
            this%points(p(3)), this%points(p(4)), &
            this%points(p(5)), this%points(p(6)), &
            this%points(p(7)), this%points(p(8)))
    class default
       call neko_error('Invalid element type')
    end select

  end subroutine mesh_add_hex

  !> Add a unique point to the mesh
  subroutine mesh_add_point(this, p, idx)
    class(mesh_t), intent(inout) :: this
    type(point_t), intent(inout) :: p
    integer, intent(inout) :: idx
    integer :: tmp

    tmp = p%id()

    this%max_pts_id = max(this%max_pts_id, tmp)

    if (tmp .le. 0) then
       call neko_error("Invalid point id")
    end if

    if (this%htp%get(tmp, idx) .gt. 0) then
       this%mpts = this%mpts + 1
       call this%htp%set(tmp, this%mpts)
       this%points(this%mpts) = p
       idx = this%mpts
    end if

  end subroutine mesh_add_point

  !> Add a unique face represented as a 4-tuple to the mesh
  !! @details Returns the local id of @a f in @a idx, adding it to the set of
  !! unique faces if it has not been seen before. @a head and @a chain hold
  !! the lookup chains described in @a mesh_generate_face_lid
  subroutine mesh_add_face(this, f, head, chain, idx)
    type(mesh_t), intent(inout) :: this
    type(tuple4_i4_t), intent(inout) :: f
    integer, intent(inout) :: head(:)
    integer, intent(inout) :: chain(:)
    integer, intent(out) :: idx
    integer :: lp

    ! The tuple is ordered, so chaining on the first point leaves any two
    ! faces in the same chain differing in their remaining points
    lp = this%have_point_glb_idx(f%x(1))
    if (lp .lt. 1) then
       call neko_error('Invalid face point')
    end if

    idx = head(lp)
    do while (idx .gt. 0)
       if ((this%face_pts(2, idx) .eq. f%x(2)) .and. &
            (this%face_pts(3, idx) .eq. f%x(3)) .and. &
            (this%face_pts(4, idx) .eq. f%x(4))) return
       idx = chain(idx)
    end do

    this%mfcs = this%mfcs + 1
    idx = this%mfcs
    this%face_pts(1, idx) = f%x(1)
    this%face_pts(2, idx) = f%x(2)
    this%face_pts(3, idx) = f%x(3)
    this%face_pts(4, idx) = f%x(4)
    chain(idx) = head(lp)
    head(lp) = idx

  end subroutine mesh_add_face

  !> Add a unique edge represented as a 2-tuple to the mesh
  !! @details Returns the local id of @a e in @a idx, adding it to the set of
  !! unique edges if it has not been seen before. @a head and @a chain hold
  !! the lookup chains described in @a mesh_generate_edge_lid
  subroutine mesh_add_edge(this, e, head, chain, idx)
    type(mesh_t), intent(inout) :: this
    type(tuple_i4_t), intent(inout) :: e
    integer, intent(inout) :: head(:)
    integer, intent(inout) :: chain(:)
    integer, intent(out) :: idx
    integer :: lp

    ! The tuple is ordered, so chaining on the first point leaves any two
    ! edges in the same chain differing in their second point
    lp = this%have_point_glb_idx(e%x(1))
    if (lp .lt. 1) then
       call neko_error('Invalid edge endpoint')
    end if

    idx = head(lp)
    do while (idx .gt. 0)
       if (this%edge_pts(2, idx) .eq. e%x(2)) return
       idx = chain(idx)
    end do

    this%meds = this%meds + 1
    idx = this%meds
    this%edge_pts(1, idx) = e%x(1)
    this%edge_pts(2, idx) = e%x(2)
    chain(idx) = head(lp)
    head(lp) = idx

  end subroutine mesh_add_edge

  !> Mark element @a e as a curve element
  subroutine mesh_mark_curve_element(this, e, curve_data, curve_type)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: e
    real(kind=dp), dimension(5,12), intent(in) :: curve_data
    integer, dimension(12), intent(in) :: curve_type

    if (e .gt. this%nelv) then
       call neko_error('Invalid element index')
    end if
    if ((this%gdim .eq. 2 .and. sum(curve_type(5:8)) .gt. 0) ) then
       call neko_error('Invalid curve element')
    end if
    call this%curve%add_element(e, curve_data, curve_type)

  end subroutine mesh_mark_curve_element

  !> Mark facet @a f in element @a e with label
  subroutine mesh_mark_labeled_facet(this, f, e, label)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: label

    if (e .gt. this%nelv) then
       call neko_error('Invalid element index')
    end if

    if ((this%gdim .eq. 2 .and. f .gt. 4) .or. &
         (this%gdim .eq. 3 .and. f .gt. 6)) then
       call neko_error('Invalid facet index')
    end if
    call this%labeled_zones(label)%add_facet(f, e)
    this%facet_type(f,e) = -label

  end subroutine mesh_mark_labeled_facet

  !> Mark facet @a f in element @a e as periodic with (@a pf, @a pe)
  subroutine mesh_mark_periodic_facet(this, f, e, pf, pe, pids)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: pf
    integer, intent(in) :: pe
    integer, intent(inout) :: pids(4)
    integer, dimension(4) :: org_ids

    call this%get_facet_ids(f, e, org_ids)
    call this%periodic%add_periodic_facet(f, e, pf, pe, pids, org_ids)
  end subroutine mesh_mark_periodic_facet

  !> Get original ids of periodic points
  subroutine mesh_get_facet_ids(this, f, e, pids)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(inout) :: pids(4)
    type(point_t), pointer :: pi
    type(tuple4_i4_t) :: t
    type(tuple_i4_t) :: t2

    select type(ele => this%elements(e)%e)
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
  subroutine mesh_reset_periodic_ids(this)
    class(mesh_t), intent(inout) :: this
    integer :: i,j
    integer :: f
    integer :: e
    integer :: pf
    integer :: pe
    integer :: org_ids(4), pids(4)
    type(point_t), pointer :: pi
    integer, dimension(4, 6) :: face_nodes = reshape([ &
         1,5,7,3, &
         2,6,8,4, &
         1,2,6,5, &
         3,4,8,7, &
         1,2,4,3, &
         5,6,8,7],&
         [4,6])
    integer, dimension(2, 4) :: edge_nodes = reshape([ &
         1,3, &
         2,4, &
         1,2, &
         3,4],&
         [2,4])

    do i = 1, this%periodic%size
       e = this%periodic%facet_el(i)%x(2)
       f = this%periodic%facet_el(i)%x(1)
       pe = this%periodic%p_facet_el(i)%x(2)
       pf = this%periodic%p_facet_el(i)%x(1)
       pids = this%periodic%p_ids(i)%x
       call this%get_facet_ids(f, e, pids)
       this%periodic%p_ids(i)%x = pids
    end do
    do i = 1, this%periodic%size
       e = this%periodic%facet_el(i)%x(2)
       f = this%periodic%facet_el(i)%x(1)
       org_ids = this%periodic%org_ids(i)%x
       select type(ele => this%elements(e)%e)
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
  subroutine mesh_create_periodic_ids(this, f, e, pf, pe)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: pf
    integer, intent(in) :: pe
    type(point_t), pointer :: pi, pj
    real(kind=dp) :: L(3)
    integer :: i, j, id, match
    type(tuple4_i4_t) :: ft
    type(tuple_i4_t) :: et
    integer :: envvar_len
    character(len=255) :: tol_str
    real(kind=dp) :: tol
    integer, dimension(4, 6) :: face_nodes = reshape([&
         1,5,7,3,&
         2,6,8,4,&
         1,2,6,5,&
         3,4,8,7,&
         1,2,4,3,&
         5,6,8,7],&
         [4,6])
    integer, dimension(2, 4) :: edge_nodes = reshape([&
         1,3,&
         2,4,&
         1,2,&
         3,4 ],&
         [2,4])

    call get_environment_variable("NEKO_PERIODIC_TOL", tol_str, envvar_len)
    if (envvar_len .gt. 0) then
       read(tol_str(1:envvar_len), *) tol
    else
       tol = 1d-7
    end if

    select type(ele => this%elements(e)%e)
    type is(hex_t)
       select type(elp => this%elements(pe)%e)
       type is(hex_t)
          L = 0d0
          do i = 1, 4
             L = L + ele%pts(face_nodes(i,f))%p%x(1:3) - &
                  elp%pts(face_nodes(i,pf))%p%x(1:3)
          end do
          L = L/4
          do i = 1, 4
             pi => ele%pts(face_nodes(i,f))%p
             match = 0
             do j = 1, 4
                pj => elp%pts(face_nodes(j,pf))%p
                if (norm2(pi%x(1:3) - pj%x(1:3) - L) .lt. tol) then
                   id = min(pi%id(), pj%id())
                   call pi%set_id(id)
                   call pj%set_id(id)
                   match = match + 1
                end if
             end do
             if ( match .gt. 1) then
                call neko_error('Multiple matches when creating periodic ids')
             else if (match .eq. 0) then
                call neko_error('Cannot find matching periodic point')
             end if
          end do
       end select
    type is(quad_t)
       select type(elp => this%elements(pe)%e)
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
                if (norm2(pi%x(1:3) - pj%x(1:3) - L) .lt. tol) then
                   id = min(pi%id(), pj%id())
                   call pi%set_id(id)
                   call pj%set_id(id)
                end if
             end do
          end do
       end select
    end select
  end subroutine mesh_create_periodic_ids

  !> Replaces the periodic point's id with a common id for matching
  !! periodic points
  subroutine mesh_apply_periodic_facet(this, f, e, pf, pe, pids)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: f
    integer, intent(in) :: e
    integer, intent(in) :: pf
    integer, intent(in) :: pe
    integer, intent(inout) :: pids(4)
    type(point_t), pointer :: pi
    integer :: i, id
    type(tuple4_i4_t) :: ft
    type(tuple_i4_t) :: et
    integer, dimension(4, 6) :: face_nodes = reshape([&
         1,5,7,3,&
         2,6,8,4,&
         1,2,6,5,&
         3,4,8,7,&
         1,2,4,3,&
         5,6,8,7],&
         [4,6])
    select type(ele => this%elements(e)%e)
    type is(hex_t)
       do i = 1, 4
          pi => ele%pts(face_nodes(i,f))%p
          call pi%set_id(pids(i))
          ! Only register the point while the mesh is being built; once the
          ! connectivity is generated the periodic ids are already known and
          ! the global->local table has been released
          if (allocated(this%htp)) then
             call this%add_point(pi, id)
          end if
       end do
    end select

  end subroutine mesh_apply_periodic_facet

  !> Return the global id of edge @a e in element @a el
  function mesh_get_global_edge(this, el, e) result(global_id)
    class(mesh_t), intent(in) :: this
    integer, intent(in) :: el !< Local element id
    integer, intent(in) :: e !< Local edge number of the element
    integer :: global_id

    global_id = this%edge_lid(e, el)

    if (this%gdim .eq. 2) then
       if (pe_size .gt. 1) then
          global_id = this%ddata%local_to_global_facet(global_id)
       end if
    else
       if (pe_size .gt. 1) then
          global_id = this%ddata%local_to_global_edge(global_id)
       end if
    end if

  end function mesh_get_global_edge

  !> Return the global id of facet @a f in element @a el
  !! @attention only defined for gdim .eq. 3, in two dimensions the facets
  !! of an element are its edges, see @a mesh_get_global_edge
  function mesh_get_global_facet(this, el, f) result(global_id)
    class(mesh_t), intent(in) :: this
    integer, intent(in) :: el !< Local element id
    integer, intent(in) :: f !< Local facet number of the element
    integer :: global_id

    global_id = this%face_lid(f, el)

    if (pe_size .gt. 1) then
       global_id = this%ddata%local_to_global_facet(global_id)
    end if

  end function mesh_get_global_facet


  !> Check if the mesh has a point given its global index
  !! @details The one way to turn a global point id into a local one; a point
  !! carries its global id, so @a p%id() is the key for a caller holding a
  !! point_t. For an element's own corner read @a pt_lid instead.
  !! @return The local id of the point (if present) otherwise -1
  !! @attention Only valid until @a generate_conn releases the global->local
  !! point table
  !! @todo Consider moving this to distdata
  function mesh_have_point_glb_idx(this, index) result(local_id)
    class(mesh_t), intent(inout) :: this
    integer, intent(inout) :: index !< Global index
    integer :: local_id

    if (.not. allocated(this%htp)) then
       call neko_error('have_point_glb_idx is only valid before generate_conn')
    end if

    if (this%htp%get(index, local_id) .eq. 1) then
       local_id = -1
    end if

  end function mesh_have_point_glb_idx


  !> Check if point @a p in element @a el is shared
  function mesh_is_shared_point(this, el, p) result(shared)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: el !< Local element id
    integer, intent(in) :: p !< Local point number of the element
    integer :: local_index
    logical shared

    local_index = this%pt_lid(p, el)
    shared = this%ddata%shared_point%element(local_index)

  end function mesh_is_shared_point


  !> Check if edge @a e in element @a el is shared
  function mesh_is_shared_edge(this, el, e) result(shared)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: el !< Local element id
    integer, intent(in) :: e !< Local edge number of the element
    integer :: local_index
    logical shared
    local_index = this%edge_lid(e, el)
    if (this%gdim .eq. 2) then
       shared = this%ddata%shared_facet%element(local_index)
    else
       shared = this%ddata%shared_edge%element(local_index)
    end if
  end function mesh_is_shared_edge

  !> Check if facet @a f in element @a el is shared
  !! @attention only defined for gdim .eq. 3, in two dimensions the facets
  !! of an element are its edges, see @a mesh_is_shared_edge
  function mesh_is_shared_facet(this, el, f) result(shared)
    class(mesh_t), intent(inout) :: this
    integer, intent(in) :: el !< Local element id
    integer, intent(in) :: f !< Local facet number of the element
    integer :: local_index
    logical shared

    local_index = this%face_lid(f, el)
    shared = this%ddata%shared_facet%element(local_index)

  end function mesh_is_shared_facet

  !> Check the correct orientation of the rst coordindates.
  !! @note Similar algorithm as in Nek5000 `verify` routine.
  subroutine mesh_check_right_handedness(this)
    class(mesh_t), intent(inout) :: this
    integer :: i
    real(kind=rp) :: v(8)
    type(point_t) :: centroid
    logical :: fail

    fail = .false.

    if (this%gdim .eq. 3) then
       do i = 1, this%nelv
          v(1) = parallelepiped_signed_volume( &
               this%elements(i)%e%pts(2)%p%x, &
               this%elements(i)%e%pts(3)%p%x, &
               this%elements(i)%e%pts(5)%p%x, &
               this%elements(i)%e%pts(1)%p%x &
               )

          v(2) = parallelepiped_signed_volume( &
               this%elements(i)%e%pts(4)%p%x, &
               this%elements(i)%e%pts(1)%p%x, &
               this%elements(i)%e%pts(6)%p%x, &
               this%elements(i)%e%pts(2)%p%x &
               )

          v(3) = parallelepiped_signed_volume( &
               this%elements(i)%e%pts(1)%p%x, &
               this%elements(i)%e%pts(4)%p%x, &
               this%elements(i)%e%pts(7)%p%x, &
               this%elements(i)%e%pts(3)%p%x &
               )

          v(4) = parallelepiped_signed_volume( &
               this%elements(i)%e%pts(3)%p%x, &
               this%elements(i)%e%pts(2)%p%x, &
               this%elements(i)%e%pts(8)%p%x, &
               this%elements(i)%e%pts(4)%p%x &
               )

          v(5) = -parallelepiped_signed_volume( &
               this%elements(i)%e%pts(6)%p%x, &
               this%elements(i)%e%pts(7)%p%x, &
               this%elements(i)%e%pts(1)%p%x, &
               this%elements(i)%e%pts(5)%p%x &
               )

          v(6) = -parallelepiped_signed_volume( &
               this%elements(i)%e%pts(8)%p%x, &
               this%elements(i)%e%pts(5)%p%x, &
               this%elements(i)%e%pts(2)%p%x, &
               this%elements(i)%e%pts(6)%p%x &
               )

          v(7) = -parallelepiped_signed_volume( &
               this%elements(i)%e%pts(5)%p%x, &
               this%elements(i)%e%pts(8)%p%x, &
               this%elements(i)%e%pts(3)%p%x, &
               this%elements(i)%e%pts(7)%p%x &
               )

          v(8) = -parallelepiped_signed_volume( &
               this%elements(i)%e%pts(7)%p%x, &
               this%elements(i)%e%pts(6)%p%x, &
               this%elements(i)%e%pts(4)%p%x, &
               this%elements(i)%e%pts(8)%p%x &
               )

          if (v(1) .le. 0.0_rp .or. &
               v(2) .le. 0.0_rp .or. &
               v(3) .le. 0.0_rp .or. &
               v(4) .le. 0.0_rp .or. &
               v(5) .le. 0.0_rp .or. &
               v(6) .le. 0.0_rp .or. &
               v(7) .le. 0.0_rp .or. &
               v(8) .le. 0.0_rp ) then

             centroid = this%elements(i)%e%centroid()

             write(error_unit, '(A, A, I0, A, 3G12.5)') "*** ERROR ***: ", &
                  "Wrong orientation of mesh element ", i, &
                  " with centroid ", centroid%x

             fail = .true.
          end if
       end do
    end if

    if (fail) then
       call neko_error("Some mesh elements are not right-handed")
    end if
  end subroutine mesh_check_right_handedness

  !> Compute a signed volume of a parallelepiped formed by three vectors, in
  !! turn defined via three points, `p1`, `p2`, and `p3` and an `origin`.
  !! @param p1 The first point.
  !! @param p2 The second point.
  !! @param p3 The third point.
  !! @param origin The point defining the origin.
  !! @note Used to check right-handness of the elements: the volumes should be
  !! positive.
  function parallelepiped_signed_volume(p1, p2, p3, origin) result(v)
    real(kind=dp), dimension(3), intent(in) :: p1, p2, p3, origin
    real(kind=dp) :: v
    real(kind=dp) :: vp1(3), vp2(3), vp3(3), cross(3)

    vp1 = p1 - origin
    vp2 = p2 - origin
    vp3 = p3 - origin

    cross(1) = vp1(2)*vp2(3) - vp2(3)*vp1(2)
    cross(2) = vp1(3)*vp2(1) - vp1(1)*vp2(3)
    cross(3) = vp1(1)*vp2(2) - vp1(2)*vp2(1)

    v = cross(1)*vp3(1) + cross(2)*vp3(2) + cross(3)*vp3(3)

  end function parallelepiped_signed_volume

  !> Create a subset of the mesh @a this in @a other based on the provided
  !! mask.
  !! @param this The original mesh.
  !! @param other The subset mesh to be created.
  !! @param mask The mask defining the subset.
  !! @param lx the quadrature degree in x direction.
  !! @param ly the quadrature degree in y direction.
  !! @param lz the quadrature degree in z direction.
  !! @note Partially lifted from nmsh_file.f90.
  subroutine mesh_subset_by_mask(this, other, mask, lx, ly, lz)
    class(mesh_t), intent(in) :: this
    class(mesh_t), intent(inout) :: other
    type(mask_t), intent(in) :: mask
    integer, intent(in) :: lx, ly, lz
    integer :: i, j, k, nelv, lxyz, gdim, e_m, nidx(4), nelv_c, el_c, el, i_m
    type(point_t) :: p(8)
    integer :: p_id = 1

    call other%free()
    lxyz = lx * ly * lz

    ! Initialize
    nelv = mask%size()/lxyz
    call other%init(this%gdim, nelv)

    ! Assign the elements
    if (other%gdim .eq. 2) then
       call neko_error("Subset mesh not implemented for 2d")
    else if (other%gdim .eq. 3) then
       do el = 1, nelv
          i_m = 1 + lxyz * (el - 1)
          nidx = nonlinear_index(mask%get(i_m), lx, ly, lz)
          e_m = nidx(4) ! Actual element from the original mesh
          ! Retrieve the points form the other mesh.
          ! No need to shift points, since original
          ! mesh has done it.
          ! Had to use a new point id to avoid issues at
          ! periodic boundaries
          ! But this means that all points might be incorrectly
          ! marked as unique.
          do j = 1, 8
             call p(j)%init(this%elements(e_m)%e%pts(j)%p%x, p_id)
             p_id = p_id + 1
          end do

          call other%add_element(el, el + other%offset_el, &
               p(1), p(2), p(3), p(4), &
               p(5), p(6), p(7), p(8))
       end do
    else
       if (pe_rank .eq. 0) call neko_error('Invalid dimension of mesh')
    end if

    ! Skip searching for boundaries.

    ! Update the curvature
    nelv_c = this%curve%size
    if (nelv_c .gt. 0) then
       el_c = 1
       el = 1
       ! 2 pointer scan
       do while (el .le. nelv .and. el_c .le. nelv_c)

          i_m = 1 + lxyz * (el - 1)
          nidx = nonlinear_index(mask%get(i_m), lx, ly, lz)
          e_m = nidx(4)

          if (e_m .lt. this%curve%curve_el(el_c)%el_idx) then
             el = el + 1

          else if (e_m .gt. this%curve%curve_el(el_c)%el_idx) then
             el_c = el_c + 1

          else
             call other%mark_curve_element(el, &
                  this%curve%curve_el(el_c)%curve_data, &
                  this%curve%curve_el(el_c)%curve_type)
             el = el + 1
             el_c = el_c + 1
          end if

       end do
    end if

    ! Finalize
    call other%finalize()

    other%is_submesh = .true.

  end subroutine mesh_subset_by_mask

end module mesh
