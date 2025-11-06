! Copyright (c) 2019-2025, The Neko Authors
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
!> Implementation of the `mesh_manager_t` type
module mesh_manager
  use num_types, only : i4, i8, rp, dp
  use json_module, only : json_file

  implicit none
  private

  !> Type for geometrical independent nodes.
  !! @details Geometrical independent nodes are globally unique nodes not
  !! sharing physical coordinates and global id and not located on the
  !! nonconforming interfaces. It is a minimal set of points needed to build
  !! a conforming, consistent element mesh based on verices only (no curvature
  !! information). This type does not include full connectivity information, as
  !! periodicity and nonconforming interfaces are omitted.
  type, public :: mesh_geom_node_ind_t
     !> Number of owned nodes
     integer(i4) :: lown
     !> Number of owned shared nodes
     integer(i4) :: lshr
     !> Local offset of owned nodes
     integer(i4) :: loff
     !> Local number of nodes
     integer(i4) :: lnum
     !> Global indexing of unique nodes
     integer(i8), allocatable, dimension(:) :: gidx
     !> Node owner (MPI rank)
     integer(i4), allocatable, dimension(:) :: ndown
     !> Geometrical mesh dimension
     integer(i4) :: gdim
     !> Physical node coordinates
     real(kind=dp), allocatable, dimension(:,:) :: coord
   contains
     procedure, pass(this) :: init => mesh_geom_node_ind_init
     procedure, pass(this) :: init_type => mesh_geom_node_ind_init_type
     procedure, pass(this) :: free => mesh_geom_node_ind_free
  end type mesh_geom_node_ind_t

  !> Type for geometrical hanging nodes
  !! @details Hanging nodes are not independent ones. They are located at
  !! the centre of the nonconforming face or edge (have unique physical
  !! coordinates) and can be maped to the independent face/edge vertices. This
  !! mapping is not unique and will depend on the element position in a tree.
  !! There are two types of hanging nodes h2 (dependent on 2 independent nodes;
  !! 2D face and 3D edge) and h4 (dependent on 4 independent nodes; 3D face).
  type, public :: mesh_geom_node_hng_t
     !> Local number of nodes
     integer(i4) :: lnum
     !> Global indexing of unique nodes
     integer(i8), allocatable, dimension(:) :: gidx
     !> Number of independent nodes in the dependency list
     integer(i4) :: ndep
     !> Local hanging to independent node mapping
     integer(i4), allocatable, dimension(:,:) :: lmap
     !> Geometrical mesh dimension
     integer(i4) :: gdim
     !> Physical node coordinates
     real(kind=dp), allocatable, dimension(:,:) :: coord
   contains
     procedure, pass(this) :: init => mesh_geom_node_hng_init
     procedure, pass(this) :: init_type => mesh_geom_node_hng_init_type
     procedure, pass(this) :: free => mesh_geom_node_hng_free
  end type mesh_geom_node_hng_t

  !> Type for element geometry information
  !! @details This type collects all the geometrical information including
  !! independent and hanging nodes lists and combines it with element vertices
  !! to node mapping
  type, public :: mesh_geom_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> Geometrical independent nodes
     type(mesh_geom_node_ind_t) :: geom_ind
     !> Geometrical h2-type hanging nodes
     type(mesh_geom_node_hng_t) :: geom_hng_edg
     !> Geometrical h4-type hanging nodes
     type(mesh_geom_node_hng_t) :: geom_hng_fcs
     !> local number of elements
     integer(i4) :: nel
     !> Number of vertices per element
     integer(i4) :: nvrt
     ! Mapping for:
     ! 1<= vnmap(iv,iel) <= nin - independent node
     ! nin < vnmap(iv,iel) <= nin + nhf - face hanging node
     ! nin + nhf < vnmap(iv,iel) <= nin + nhf + nhe - edge hanging node
     ! where
     ! nin = number of local independent nodes
     ! nhf = number of local face hanging nodes
     ! nhe = number of local edge hanging nodes
     ! vnmap uses symmetric vertex notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !             3+--------+4    ^ s
     !             /        /|     |
     !            /        / |     |
     !           /        /  |     |
     !         7+--------+8  +2    +----> r
     !          |        |  /     /
     !          |        | /     /
     !          |        |/     /
     !         5+--------+6    t
     !> Local mapping of element vertices to the nodes
     integer(i4), allocatable, dimension(:,:) :: vmap
   contains
     procedure, pass(this) :: init => mesh_geom_init
     procedure, pass(this) :: init_type => mesh_geom_init_type
     procedure, pass(this) :: free => mesh_geom_free
  end type mesh_geom_t

  !> Type for connectivity information regarding vertices, faces and edges.
  !! @details It contains both object global numbering and communication
  !! information
  type, public :: mesh_conn_obj_t
     !> Number of local objects
     integer(i4) :: lnum
     !> Number of owned objects
     integer(i4) :: lown
     !> Global object offset
     integer(i8) :: goff
     !> Global number of objects
     integer(i8) :: gnum
     !> Number of MPI ranks sharing objects
     integer(i4) :: nrank
     !> Number of shared objects
     integer(i4) :: nshare
     !> Global indexing of unique objects of given type
     integer(i8), allocatable, dimension(:) :: gidx
     !> List of ranks sharing objects
     integer(i4), allocatable, dimension(:) :: rank
     !> List of shared objects
     integer(i4), allocatable, dimension(:) :: share
     !> Offset in the share list
     integer(i4), allocatable, dimension(:) :: off
   contains
     procedure, pass(this) :: init => mesh_conn_obj_init
     procedure, pass(this) :: init_type => mesh_conn_obj_init_type
     procedure, pass(this) :: free => mesh_conn_obj_free
  end type mesh_conn_obj_t

  !> Type for element connectivity information
  type, public :: mesh_conn_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> Connectivity information for vertices, edges and faces
     type(mesh_conn_obj_t) :: conn_vrt, conn_edg, conn_fcs
     !> local number of elements
     integer(i4) :: nel
     !> Number of vertices per element
     integer(i4) :: nvrt
     ! Node info for building communicators; includes global indexing of various
     ! element objects
     ! vertex%lmap uses symmetric vertex notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !             3+--------+4    ^ s
     !             /        /|     |
     !            /        / |     |
     !           /        /  |     |
     !         7+--------+8  +2    +----> r
     !          |        |  /     /
     !          |        | /     /
     !          |        |/     /
     !         5+--------+6    t
     !
     !> Element vertices to object mapping
     integer(i4), allocatable, dimension(:,:) :: vmap

     !> Number of faces per element
     integer(i4) :: nfcs
     ! face%lmap uses symmetric face notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !              +--------+     ^ s
     !             /        /|     |
     !            /    4   / |     |
     !      1--> /        /  |     |
     !          +--------+ 2 +     +----> r
     !          |        |  /     /
     !          |    6   | /     /
     !          |        |/     /
     !          +--------+     t
     !               3
     !> Local element faces to object mapping
     integer(i4), allocatable, dimension(:,:) :: fmap
     ! Face orientation.
     ! In 2D case permutation array takes 2 values
     !          0 => no permutations
     !          1 => row permutations
     ! In 3D case there are 8 possible transformations for faces
     ! Symmetric vertex numbering  (0,1,2,3)
     !                       2--3
     !                       |  |
     !                       0--1
     ! They are numbered according to p4est by 0..7 and stored in falg
     ! { 0, 1, 2, 3 }      0 => identity
     ! { 0, 2, 1, 3 }      1 => T
     ! { 1, 0, 3, 2 }      2 => P_x
     ! { 1, 3, 0, 2 }      3 => P_x T
     ! { 2, 0, 3, 1 }      4 => P_y T
     ! { 2, 3, 0, 1 }      5 => P_y
     ! { 3, 1, 2, 0 }      6 => P_y P_x T
     ! { 3, 2, 1, 0 }      7 => P_y P_x
     ! where   T - transpose
     !               P_x - permutation in x (rows)
     !               P_y - permutation in y (columns)
     !> Face alignment
     integer(i4), allocatable, dimension(:,:) :: falgn

     !> Number of edges per element
     integer(i4) :: nedg
     ! edge%lmap uses symmetric edge notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !              +---2----+     ^ s
     !             /        /|     |
     !           11       12 6     |
     !           /        /  |     |
     !          +---4----+   +     +----> r
     !          |        |  /     /
     !          7        8 10    /
     !          |        |/     /
     !          +---3----+     t
     !
     !> Local element edges to object mapping
     integer(i4), allocatable, dimension(:,:) :: emap
     ! Edge orientation is similar to 2D face
     !> Edge alignment
     integer(i4), allocatable, dimension(:,:) :: ealgn

     ! Topology information: hanging elements, faces, edges and vertices
     ! The element is marked as hanging if it contains at least one hanging
     ! object.
     ! Hanging face
     ! 2D
     !   = -1 if the face is not hanging,
     !   = 0 if the face is the first half,
     !   = 1 if the face is the second half.
     ! 3D
     !   = -1 if the face is not hanging,
     !   = the corner of the full face that it touches:
     ! Hanging edge
     !   = -1 if the edge is not hanging,
     !   =  0 if the edge is the first half of a full edge,
     !        but neither of the two faces touching the
     !        edge is hanging,
     !   =  1 if the edge is the second half of a full edge,
     !        but neither of the two faces touching the
     !        edge is hanging,
     !   =  2 if the edge is the first half of a full edge
     !        and is on the boundary of a full face,
     !   =  3 if the edge is the second half of a full edge
     !        and is on the boundary of a full face,
     !   =  4 if the edge is in the middle of a full face.
     !        See the diagram below for clarification.
     !  o...............o o...............o +---2---+.......o o.......+---3---+
     !  :               : :               : |       |       : :       |       |
     !  :               : :               : 3   2   4       : :       4   3   3
     !  :               : :               : |       |       : :       |       |
     !  +---4---+       : :       +---4---+ +---4---+       : :       +---4---+
     !  |       |       : :       |       | :               : :               :
     !  2   0   4       : :       4   1   2 :               : :               :
     !  |       |       : :       |       | :               : :               :
     !  +---2---+.......o o.......+---3---+ o...............o o...............o
     !
     !                     o                  +-------+
     !                     :                  |\       \
     !                     :                  1 \       \
     !                     :                  |  +-------+
     !                     +-------+          +  |       |
     !                     |\       \         :\ |       |
     !                     0 \       \        : \|       |
     !                     |  +-------+       :  +-------+
     !                     +  |       |       o
     !                      \ |       |
     !                       \|       |
     !                        +-------+
     !> Hanging element list; 1 if at least one hanging face or edge otherwise 0
     integer(i4), allocatable, dimension(:) :: hngel
     !> Hanging face list; position of hanging face; otherwise -1
     integer(i4), allocatable, dimension(:,:) :: hngfc
     !> Hanging edge list;
     integer(i4), allocatable, dimension(:,:) :: hnged

   contains
     procedure, pass(this) :: init => mesh_conn_init
     procedure, pass(this) :: init_type => mesh_conn_init_type
     procedure, pass(this) :: free => mesh_conn_free
  end type mesh_conn_t

  type, public :: mesh_manager_mesh_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> total number of local element (both V- and T-type)
     integer(i4) :: nelt
     !> number of local V-type elements
     integer(i4) :: nelv
     !> global element number
     integer(i8) :: gnelt
     !> global element offset
     integer(i8) :: gnelto
     !> max refinement level across the whole mesh
     integer(i4) :: level_max
     !> global element index
     integer(i8), allocatable, dimension(:) :: gidx
     !> element refinement level
     integer(i4), allocatable, dimension(:) :: level
     !> element group, (not used right now)
     integer(i4), allocatable, dimension(:) :: igrp
     !> Number of faces per element
     integer(i4) :: nfcs
     !> face curvature flag (not used right now)
     integer(i4), allocatable, dimension(:,:) :: crv
     !> face boundary condition: -1- periodic, 0-internal, 0< user specified
     integer(i4), allocatable, dimension(:,:) :: bc
     !> Geometrical information
     type(mesh_geom_t) :: geom
     !> Connectivity information
     type(mesh_conn_t) :: conn
   contains
     !> Constructor for the mesh data
     procedure, pass(this) :: init => mesh_manager_mesh_init
     !> Constructor for the mesh data based on the other mesh type
     procedure, pass(this) :: init_type => mesh_manager_mesh_init_type
     !> Destructor for the mesh data
     procedure, pass(this) :: free => mesh_manager_mesh_free
  end type mesh_manager_mesh_t

  !> Base abstract type for mesh manager.
  type, abstract, public :: mesh_manager_t
     !> Manager type name
     character(len=:), allocatable :: type_name
     !> 3rd-party software activation flag
     logical :: ifstarted
     !> mesh information
     type(mesh_manager_mesh_t) :: mesh
   contains
     !> Constructor for the mesh_manager_t (base) type.
     procedure, pass(this) :: init_base => mesh_manager_init_base
     !> Destructor for the mesh_manager_t (base) type.
     procedure, pass(this) :: free_base => mesh_manager_free_base
     !> Start 3rd-party software (if needed)
     procedure(mesh_manager_start), pass(this), deferred :: start
     !> Stop 3rd-party software (if needed)
     procedure(mesh_manager_free), pass(this), deferred :: stop
     !> The common constructor using a JSON object.
     procedure(mesh_manager_init), pass(this), deferred :: init
     !> Destructor.
     procedure(mesh_manager_free), pass(this), deferred :: free
     !> Import mesh data into current type
     procedure(mesh_manager_free), pass(this), deferred :: import
     !> Import mesh data creating a new variable
     procedure(mesh_manager_import_new), pass(this), deferred :: import_new
  end type mesh_manager_t

  abstract interface
     !> Start 3rd-party software (if needed)
     !! @param[out]  ierr  error flag
     subroutine mesh_manager_start(this, json, ierr)
       import mesh_manager_t, json_file
       class(mesh_manager_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       integer, intent(out) :: ierr
     end subroutine mesh_manager_start

     !> The common constructor using a JSON object.
     !! @param json       The JSON object for the mesh manager.
     !! @param type_name  Manager type name
     subroutine mesh_manager_init(this, json)
       import mesh_manager_t, json_file
       class(mesh_manager_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
     end subroutine mesh_manager_init

     !> Destructor, 3rd-party stopping, importing
     subroutine mesh_manager_free(this)
       import mesh_manager_t
       class(mesh_manager_t), intent(inout) :: this
     end subroutine mesh_manager_free

     !> Import mesh data creating a new variable
     !! @param  mesh_new   new mesh data
     subroutine mesh_manager_import_new(this, mesh_new)
       import mesh_manager_t
       class(mesh_manager_t), intent(inout) :: this
       class(mesh_manager_t), allocatable, intent(inout) :: mesh_new
     end subroutine mesh_manager_import_new
  end interface

  interface
     !> Mesh manager factory. Both constructs and initializes the object.
     !! @param object The object to be initialised.
     !! @param json JSON object initialising the mesh manager.
     module subroutine mesh_manager_factory(object, json)
       class(mesh_manager_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
     end subroutine mesh_manager_factory
  end interface

  interface
     !> Mesh manager allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the type to allocate.
     module subroutine mesh_manager_allocator(object, type_name)
       class(mesh_manager_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
     end subroutine mesh_manager_allocator
  end interface

  public :: mesh_manager_factory

contains

  !> Initialise independent nodes type
  !! @param[in]    lown    number of owned nodes
  !! @param[in]    lshr    number of owned shared nodes
  !! @param[in]    loff    local offset of owned nodes
  !! @param[in]    lnum    local number of nodes
  !! @param[in]    gdim    geometrical dimension
  !! @param[inout] gidx    global node index
  !! @param[inout] ndown   node owner (MPI rank)
  !! @param[inout] coord   node coordinates
  subroutine mesh_geom_node_ind_init(this, lown, lshr, loff, lnum, gdim, &
       gidx, ndown, coord)
    ! argument list
    class(mesh_geom_node_ind_t), intent(inout) :: this
    integer(i4), intent(in) :: lown, lshr, loff, lnum, gdim
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: ndown
    real(kind=dp), allocatable, dimension(:,:), intent(inout) :: coord

    call this%free()

    this%lown = lown
    this%lshr = lshr
    this%loff = loff
    this%lnum = lnum
    this%gdim = gdim

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)
    if (allocated(ndown)) call move_alloc(ndown, this%ndown)
    if (allocated(coord)) call move_alloc(coord, this%coord)

  end subroutine mesh_geom_node_ind_init

  !> Initialise independent nodes type based on another node type
  !! @param[inout] node   node data
  subroutine mesh_geom_node_ind_init_type(this, node)
    ! argument list
    class(mesh_geom_node_ind_t), intent(inout) :: this
    type(mesh_geom_node_ind_t), intent(inout) :: node

    call this%free()

    this%lown = node%lown
    this%lshr = node%lshr
    this%loff = node%loff
    this%lnum = node%lnum
    this%gdim = node%gdim

    if (allocated(node%gidx)) call move_alloc(node%gidx, this%gidx)
    if (allocated(node%ndown)) call move_alloc(node%ndown, this%ndown)
    if (allocated(node%coord)) call move_alloc(node%coord, this%coord)

  end subroutine mesh_geom_node_ind_init_type

  !> Free independent nodes type
  subroutine mesh_geom_node_ind_free(this)
    ! argument list
    class(mesh_geom_node_ind_t), intent(inout) :: this

    this%lown = 0
    this%lshr = 0
    this%loff = 0
    this%lnum = 0
    this%gdim = 0

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%ndown)) deallocate(this%ndown)
    if (allocated(this%coord)) deallocate(this%coord)

  end subroutine mesh_geom_node_ind_free

  !> Initialise hanging nodes type
  !! @param[in]    lnum    local number of nodes
  !! @param[in]    gdim    geometrical dimension
  !! @param[in]    ndep    node dependency
  !! @param[inout] gidx    global node index
  !! @param[inout] lmap    node mapping to independent node
  !! @param[inout] coord   node coordinates
  subroutine mesh_geom_node_hng_init(this, lnum, gdim, ndep, gidx, lmap, coord)
    ! argument list
    class(mesh_geom_node_hng_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, gdim, ndep
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), allocatable, dimension(:, :), intent(inout) :: lmap
    real(kind=dp), allocatable, dimension(:,:), intent(inout) :: coord

    call this%free()

    this%lnum = lnum
    this%ndep = ndep
    this%gdim = gdim

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)
    if (allocated(lmap)) call move_alloc(lmap, this%lmap)
    if (allocated(coord)) call move_alloc(coord, this%coord)

  end subroutine mesh_geom_node_hng_init

  !> Initialise hanging nodes type based on another node type
  !! @param[inout] node   node data
  subroutine mesh_geom_node_hng_init_type(this, node)
    ! argument list
    class(mesh_geom_node_hng_t), intent(inout) :: this
    type(mesh_geom_node_hng_t), intent(inout) :: node

    call this%free()

    this%lnum = node%lnum
    this%ndep = node%ndep
    this%gdim = node%gdim

    if (allocated(node%gidx)) call move_alloc(node%gidx, this%gidx)
    if (allocated(node%lmap)) call move_alloc(node%lmap, this%lmap)
    if (allocated(node%coord)) call move_alloc(node%coord, this%coord)

  end subroutine mesh_geom_node_hng_init_type

  !> Free hanging nodes type
  subroutine mesh_geom_node_hng_free(this)
    ! argument list
    class(mesh_geom_node_hng_t), intent(inout) :: this

    this%lnum = 0
    this%ndep = 0
    this%gdim = 0

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%coord)) deallocate(this%coord)

  end subroutine mesh_geom_node_hng_free

  !> Initialise geometry type
  !! @param[in]    tdim    topological mesh dimension
  !! @param[in]    nel     local element number
  !! @param[inout] vmap    element vertices to node mapping
  subroutine mesh_geom_init(this, tdim, nel, vmap)
    class(mesh_geom_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    integer(i4), allocatable, dimension(:,:), intent(inout) :: vmap

    call this%free()

    this%tdim = tdim
    this%nel = nel
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nvrt = 2**tdim

    if (allocated(vmap)) call move_alloc(vmap, this%vmap)

  end subroutine mesh_geom_init

  !> Initialise geometry type based on another geometry type
  !! @param[inout] geom   geometry data
  subroutine mesh_geom_init_type(this, geom)
    ! argument list
    class(mesh_geom_t), intent(inout) :: this
    type(mesh_geom_t), intent(inout) :: geom

    call this%free()

    call this%geom_ind%init_type(geom%geom_ind)
    call this%geom_hng_edg%init_type(geom%geom_hng_edg)
    call this%geom_hng_fcs%init_type(geom%geom_hng_fcs)

    this%tdim = geom%tdim
    this%nel = geom%nel
    this%nvrt = geom%nvrt

    if (allocated(geom%vmap)) call move_alloc(geom%vmap, this%vmap)

  end subroutine mesh_geom_init_type

  !> Free geometry information
  subroutine mesh_geom_free(this)
    class(mesh_geom_t), intent(inout) :: this

    call this%geom_ind%free()
    call this%geom_hng_edg%free()
    call this%geom_hng_fcs%free()

    this%tdim = 0
    this%nel = 0
    this%nvrt = 0

    if (allocated(this%vmap)) deallocate(this%vmap)

  end subroutine mesh_geom_free

  !> Initialise connectivity object type
  !! @param[in]    lnum    local number of objects
  !! @param[in]    lown    number of owned objects
  !! @param[in]    goff    global object offset
  !! @param[in]    gnum    global number of objects
  !! @param[in]    nrank   number of MPI ranks sharing objects and
  !! @param[in]    nshare  number of shared objects
  !! @param[inout] gidx    global object index
  !! @param[inout] rank    list of ranks sharing object
  !! @param[inout] share   list of shared objects
  !! @param[inout] off     offset in share list
  subroutine mesh_conn_obj_init(this, lnum, lown, goff, gnum, nrank, nshare, &
       gidx, rank, share, off)
    ! argument list
    class(mesh_conn_obj_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, lown, nrank, nshare
    integer(i8), intent(in) :: goff, gnum
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: rank, share, off

    call this%free()

    this%lnum = lnum
    this%lown = lown
    this%goff = goff
    this%gnum = gnum
    this%nrank = nrank
    this%nshare = nshare

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)
    if (allocated(rank)) call move_alloc(rank, this%rank)
    if (allocated(share)) call move_alloc(share, this%share)
    if (allocated(off)) call move_alloc(off, this%off)

  end subroutine mesh_conn_obj_init

  !> Initialise connectivity object type based on another connectivity type
  !! @param[inout] conn   connectivity object data
  subroutine mesh_conn_obj_init_type(this, conn)
    ! argument list
    class(mesh_conn_obj_t), intent(inout) :: this
    type(mesh_conn_obj_t), intent(inout) :: conn

    call this%free()

    this%lnum = conn%lnum
    this%lown = conn%lown
    this%goff = conn%goff
    this%gnum = conn%gnum
    this%nrank = conn%nrank
    this%nshare = conn%nshare

    if (allocated(conn%gidx)) call move_alloc(conn%gidx, this%gidx)
    if (allocated(conn%rank)) call move_alloc(conn%rank, this%rank)
    if (allocated(conn%share)) call move_alloc(conn%share, this%share)
    if (allocated(conn%off)) call move_alloc(conn%off, this%off)

  end subroutine mesh_conn_obj_init_type

  !> Free connectivity object type
  subroutine mesh_conn_obj_free(this)
    ! argument list
    class(mesh_conn_obj_t), intent(inout) :: this

    this%lnum = 0
    this%lown = 0
    this%goff = 0
    this%gnum = 0
    this%nrank = 0
    this%nshare = 0

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%rank)) deallocate(this%rank)
    if (allocated(this%share)) deallocate(this%share)
    if (allocated(this%off)) deallocate(this%off)

  end subroutine mesh_conn_obj_free

  !> Initialise connectivity information
  !! @param[in]    tdim    topological dimension
  !! @param[in]    nel     local element number
  !! @param[inout] vmap    element vertex mapping
  !! @param[inout] fmap    element face mapping
  !! @param[inout] falgn   element face alignment
  !! @param[inout] emap    element edge mapping
  !! @param[inout] ealgn   element edge alignment
  !! @param[inout] hngel   element hanging flag
  !! @param[inout] hngfc   element face hanging flag
  !! @param[inout] hnged   element edge hanging flag
  subroutine mesh_conn_init(this, tdim, nel, vmap, fmap, falgn, emap, ealgn, &
       hngel, hngfc, hnged)
    ! argument list
    class(mesh_conn_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    integer(i4), allocatable, dimension(:), intent(inout) :: hngel
    integer(i4), allocatable, dimension(:,:), intent(inout) :: vmap, fmap, &
         falgn, emap, ealgn, hngfc, hnged

    call this%free()

    this%tdim = tdim
    this%nel = nel
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nvrt = 2**tdim
    this%nfcs = 2 * tdim
    this%nedg = 12 * (tdim -2) !tdim * 2**(tdim - 1)

    if (allocated(vmap)) call move_alloc(vmap, this%vmap)
    if (allocated(fmap)) call move_alloc(fmap, this%fmap)
    if (allocated(falgn)) call move_alloc(falgn, this%falgn)
    if (allocated(emap)) call move_alloc(emap, this%emap)
    if (allocated(ealgn)) call move_alloc(ealgn, this%ealgn)
    if (allocated(hngel)) call move_alloc(hngel, this%hngel)
    if (allocated(hngfc)) call move_alloc(hngfc, this%hngfc)
    if (allocated(hnged)) call move_alloc(hnged, this%hnged)

  end subroutine mesh_conn_init

  !> Initialise connectivity type based on another connectivity type
  !! @param[inout] conn   connectivity data
  subroutine mesh_conn_init_type(this, conn)
    ! argument list
    class(mesh_conn_t), intent(inout) :: this
    type(mesh_conn_t), intent(inout) :: conn

    call this%free()

    call this%conn_vrt%init_type(conn%conn_vrt)
    call this%conn_fcs%init_type(conn%conn_fcs)
    call this%conn_edg%init_type(conn%conn_edg)

    this%tdim = conn%tdim
    this%nel = conn%nel
    this%nvrt = conn%nvrt
    this%nfcs = conn%nfcs
    this%nedg = conn%nedg

    if (allocated(conn%vmap)) call move_alloc(conn%vmap, this%vmap)
    if (allocated(conn%fmap)) call move_alloc(conn%fmap, this%fmap)
    if (allocated(conn%falgn)) call move_alloc(conn%falgn, this%falgn)
    if (allocated(conn%emap)) call move_alloc(conn%emap, this%emap)
    if (allocated(conn%ealgn)) call move_alloc(conn%ealgn, this%ealgn)
    if (allocated(conn%hngel)) call move_alloc(conn%hngel, this%hngel)
    if (allocated(conn%hngfc)) call move_alloc(conn%hngfc, this%hngfc)
    if (allocated(conn%hnged)) call move_alloc(conn%hnged, this%hnged)

  end subroutine mesh_conn_init_type

  !> Free connectivity information
  subroutine mesh_conn_free(this)
    ! argument list
    class(mesh_conn_t), intent(inout) :: this

    call this%conn_vrt%free()
    call this%conn_fcs%free()
    call this%conn_edg%free()

    this%tdim = 0
    this%nel = 0
    this%nvrt = 0
    this%nfcs = 0
    this%nedg = 0

    if (allocated(this%vmap)) deallocate(this%vmap)
    if (allocated(this%fmap)) deallocate(this%fmap)
    if (allocated(this%falgn)) deallocate(this%falgn)
    if (allocated(this%emap)) deallocate(this%emap)
    if (allocated(this%ealgn)) deallocate(this%ealgn)
    if (allocated(this%hngel)) deallocate(this%hngel)
    if (allocated(this%hngfc)) deallocate(this%hngfc)
    if (allocated(this%hnged)) deallocate(this%hnged)

  end subroutine mesh_conn_free

  !>  Initialise mesh data
  !! @param[in]    nelt       total number of local elements
  !! @param[in]    nelv       number of V-type elements
  !! @param[in]    gnelt      global number of all the elements
  !! @param[in]    gnelto     global element offset
  !! @param[in]    level_max  max refinement level across the mesh
  !! @param[in]    tdim       topological mesh dimension
  !! @param[inout] gidx       global element number
  !! @param[inout] level      element refinement level
  !! @param[inout] igrp       element group
  !! @param[inout] crv        face curvature data
  !! @param[inout] bc         face boundary condition
  subroutine mesh_manager_mesh_init(this, nelt, nelv, gnelt, gnelto, &
       level_max, tdim, gidx, level, igrp, crv, bc)
    class(mesh_manager_mesh_t), intent(inout) :: this
    integer(i4), intent(in) :: nelt, nelv, level_max, tdim
    integer(i8), intent(in) :: gnelt, gnelto
    integer(i8), allocatable, dimension(:), intent(inout)  :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: level, igrp
    integer(i4), allocatable, dimension(:,:), intent(inout) :: crv, bc

    call this%free()

    this%tdim = tdim
    this%nelt = nelt
    this%nelv = nelv
    this%gnelt = gnelt
    this%gnelto = gnelto
    this%level_max  = level_max
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nfcs = 2 * tdim

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)
    if (allocated(level)) call move_alloc(level, this%level)
    if (allocated(igrp)) call move_alloc(igrp, this%igrp)
    if (allocated(crv)) call move_alloc(crv, this%crv)
    if (allocated(bc)) call move_alloc(bc, this%bc)

  end subroutine mesh_manager_mesh_init

  !>  Initialise mesh data based on another mesh type
  !! @param[inout] mesh   mesh data
  subroutine mesh_manager_mesh_init_type(this, mesh)
    class(mesh_manager_mesh_t), intent(inout) :: this
    type(mesh_manager_mesh_t), intent(inout) :: mesh

    call this%free()

    call this%geom%init_type(mesh%geom)
    call this%conn%init_type(mesh%conn)

    this%tdim = mesh%tdim
    this%nelt = mesh%nelt
    this%nelv = mesh%nelv
    this%gnelt = mesh%gnelt
    this%gnelto = mesh%gnelto
    this%level_max = mesh%level_max
    this%nfcs = mesh%nfcs

    if (allocated(mesh%gidx)) call move_alloc(mesh%gidx, this%gidx)
    if (allocated(mesh%level)) call move_alloc(mesh%level, this%level)
    if (allocated(mesh%igrp)) call move_alloc(mesh%igrp, this%igrp)
    if (allocated(mesh%crv)) call move_alloc(mesh%crv, this%crv)
    if (allocated(mesh%bc)) call move_alloc(mesh%bc, this%bc)

  end subroutine mesh_manager_mesh_init_type

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine mesh_manager_mesh_free(this)
    class(mesh_manager_mesh_t), intent(inout) :: this

    call this%geom%free()
    call this%conn%free()

    this%tdim = 0
    this%nelt = 0
    this%nelv = 0
    this%gnelt = 0
    this%gnelto = 0
    this%level_max = 0
    this%nfcs = 0

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%igrp)) deallocate(this%igrp)
    if (allocated(this%crv)) deallocate(this%crv)
    if (allocated(this%bc)) deallocate(this%bc)

  end subroutine mesh_manager_mesh_free

  !> Constructor for the `mesh_manager_t` (base) type.
  subroutine mesh_manager_init_base(this, type_name)
    class(mesh_manager_t), intent(inout) :: this
    character(len=*), intent(in) :: type_name

    this%type_name = trim(type_name)
    this%ifstarted = .false.

  end subroutine mesh_manager_init_base

  !> Destructor for the `mesh_manager_t` (base) type.
  subroutine mesh_manager_free_base(this)
    class(mesh_manager_t), intent(inout) :: this

    call this%mesh%free()

    this%ifstarted = .false.

    if (allocated(this%type_name)) deallocate(this%type_name)

  end subroutine mesh_manager_free_base


end module mesh_manager
