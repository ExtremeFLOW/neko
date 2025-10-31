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
     integer(i4) :: lown ! number of owned nodes
     integer(i4) :: lshr ! number of owned shared nodes
     integer(i4) :: loff ! local offset of owned nodes
     integer(i4) :: lnum ! local number of nodes
      ! global indexing of unique nodes
     integer(i8), allocatable, dimension(:) :: gidx
     integer(i4), allocatable, dimension(:) :: ndown ! node owner (mpi rank)
     ! physical node coordinates
     real(kind=dp), allocatable, dimension(:,:) :: coord
   contains
     procedure, public, pass(this) :: init_ind => mesh_geom_node_ind_init
     procedure, public, pass(this) :: free_ind => mesh_geom_node_ind_free
  end type mesh_geom_node_ind_t

  !> Type for geometrical hanging nodes
  !! @details Hanging nodes are not independent ones. They are located at
  !! the centre of the nonconforming face or edge (have unique physical
  !! coordinates) and can be maped to the independent face/edge vertices. This
  !! mapping is not unique and will depend on the element position in a tree.
  !! There are two types of hanging nodes h2 (dependent on 2 independent nodes;
  !! 2D face and 3D edge) and h4 (dependent on 4 independent nodes; 3D face).
  type, public, extends(mesh_geom_node_ind_t) :: mesh_geom_node_hng_t
     ! local hanging to independent node mapping
     integer(i4), allocatable, dimension(:,:) :: lmap
   contains
     procedure, public, pass(this) :: init_hng => mesh_geom_node_hng_init
     procedure, public, pass(this) :: free_hng => mesh_geom_node_hng_free
  end type mesh_geom_node_hng_t

  !> Type for connectivity information regarding vertices, faces and edges.
  !! @details It contains both object global numbering and communication
  !! information
  type, public :: mesh_conn_obj_t
     integer(i4) :: lnum ! number of local objects
     integer(i4) :: lown ! number of owned objects
     integer(i8) :: goff ! global object offset
     integer(i8) :: gnum ! global number of objects
     ! element vertices/faces/edges to object mapping
     integer(i4), allocatable, dimension(:,:) :: lmap
     integer(i4) :: nrank ! number of MPI ranks sharing objects and
     integer(i4) :: nshare ! number of shared objects
     ! global indexing of unique objects of given type
     integer(i8), allocatable, dimension(:) :: lgidx
     ! list of ranks sharing objects
     integer(i4), allocatable, dimension(:) :: lrank
     integer(i4), allocatable, dimension(:) :: lshare ! list of shared objects
     integer(i4), allocatable, dimension(:) :: loff ! offset in the lshare list
   contains
     procedure, public, pass(this) :: init => mesh_conn_obj_init
     procedure, public, pass(this) :: free => mesh_conn_obj_free
  end type mesh_conn_obj_t

  !> Base abstract type for mesh manager.
  type, abstract, public :: mesh_manager_t
     !> Manager type name
     character(len=:), allocatable :: type_name
     !> Geometrical independent nodes
     type(mesh_geom_node_ind_t) :: geom_ind
     !> Geometrical h2-type hanging nodes
     type(mesh_geom_node_hng_t) :: geom_hng_h2
     !> Geometrical h4-type hanging nodes
     type(mesh_geom_node_hng_t) :: geom_hng_h4
     !> Connectivity information of vertices, edges and faces
     type(mesh_conn_obj_t) :: conn_vrt, conn_edg, conn_fcs
   contains
     !> Constructor for the mesh_manager_t (base) type.
     procedure, pass(this) :: init_base => mesh_manager_init_base
     !> Destructor for the mesh_manager_t (base) type.
     procedure, pass(this) :: free_base => mesh_manager_free_base
     !> Destructor for the data in the mesh_manager_t type.
     procedure, pass(this) :: free_base_data => mesh_manager_free_base_data
     !> The common constructor using a JSON object.
     procedure(mesh_manager_init), pass(this), deferred :: init
     !> Destructor.
     procedure(mesh_manager_free), pass(this), deferred :: free
  end type mesh_manager_t

  abstract interface
     !> The common constructor using a JSON object.
     !! @param json       The JSON object for the mesh manager.
     !! @param type_name  Manager type name
     subroutine mesh_manager_init(this, json, type_name)
       import mesh_manager_t, json_file
       class(mesh_manager_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       character(len=*), intent(in) :: type_name
     end subroutine mesh_manager_init
  end interface

  abstract interface
     !> Destructor.
     subroutine mesh_manager_free(this)
       import mesh_manager_t
       class(mesh_manager_t), intent(inout) :: this
     end subroutine mesh_manager_free
  end interface

  interface
     !> Mesh manager factory. Both constructs and initializes the object.
     !! @param object The object to be initialized.
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
  !! @param[in] lown    number of owned nodes
  !! @param[in] lshr    number of owned shared nodes
  !! @param[in] loff    local offset of owned nodes
  !! @param[in] lnum    local number of nodes
  !! @param[in] ndim    geometrical dimension
  subroutine mesh_geom_node_ind_init(this, lown, lshr, loff, lnum, ndim)
    ! argument list
    class(mesh_geom_node_ind_t), intent(inout) :: this
    integer(i4), intent(in) :: lown, lshr, loff, lnum, ndim

    call this%free_ind()

    this%lown = lown
    this%lshr = lshr
    this%loff = loff
    this%lnum = lnum

    allocate(this%gidx(lnum), this%ndown(lnum), this%coord(lnum, ndim))

  end subroutine mesh_geom_node_ind_init

  !> Free independent nodes type
  subroutine mesh_geom_node_ind_free(this)
    ! argument list
    class(mesh_geom_node_ind_t), intent(inout) :: this

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%ndown)) deallocate(this%ndown)
    if (allocated(this%coord)) deallocate(this%coord)

    ! Reset registers
    this%lown = 0
    this%lshr = 0
    this%loff = 0
    this%lnum = 0

  end subroutine mesh_geom_node_ind_free

  !> Initialise hanging nodes type
  !! @param[in] lown    number of owned nodes
  !! @param[in] lshr    number of owned shared nodes
  !! @param[in] loff    local offset of owned nodes
  !! @param[in] lnum    local number of nodes
  !! @param[in] ndim    geometrical dimension
  !! @param[in] ndep    node dependency
  subroutine mesh_geom_node_hng_init(this, lown, lshr, loff, lnum, ndim, ndep)
    ! argument list
    class(mesh_geom_node_hng_t), intent(inout) :: this
    integer(i4), intent(in) :: lown, lshr, loff, lnum, ndim, ndep

    call this%free_hng()

    call this%init_ind(lown, lshr, loff, lnum, ndim)
    allocate(this%lmap(lnum, ndep))

  end subroutine mesh_geom_node_hng_init

  !> Free hanging nodes type
  subroutine mesh_geom_node_hng_free(this)
    ! argument list
    class(mesh_geom_node_hng_t), intent(inout) :: this

    call this%free_ind()
    if (allocated(this%lmap)) deallocate(this%lmap)

  end subroutine mesh_geom_node_hng_free

  !> Initialise connectivity object type
  !! @param[in] lnum    local number of objects
  !! @param[in] lown    number of owned objects
  !! @param[in] goff    global object offset
  !! @param[in] gnum    global number of objects
  !! @param[in] nrank   number of MPI ranks sharing objects and
  !! @param[in] nshare  number of shared objects
  !! @param[in] ndep    node dependency
  subroutine mesh_conn_obj_init(this, lnum, lown, goff, gnum, nrank, nshare, &
       ndep)
    ! argument list
    class(mesh_conn_obj_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, lown, nrank, nshare, ndep
    integer(i8), intent(in) :: goff, gnum

    call this%free()

    this%lnum = lnum
    this%lown = lown
    this%goff = goff
    this%gnum = gnum
    this%nrank = nrank
    this%nshare = nshare

!    allocate(this%lmap(lnum,), this%lgidx(lnum), this%lrank(lnum), &
!         this%lshare(), this%loff())

  end subroutine mesh_conn_obj_init

  !> Free connectivity object type
  subroutine mesh_conn_obj_free(this)
    ! argument list
    class(mesh_conn_obj_t), intent(inout) :: this

    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%lgidx)) deallocate(this%lgidx)
    if (allocated(this%lrank)) deallocate(this%lrank)
    if (allocated(this%lshare)) deallocate(this%lshare)
    if (allocated(this%loff)) deallocate(this%loff)

    ! Reset registers
    this%lnum = 0
    this%lown = 0
    this%goff = 0
    this%gnum = 0
    this%nrank = 0
    this%nshare = 0

  end subroutine mesh_conn_obj_free

  !> Constructor for the `mesh_manager_t` (base) type.
  subroutine mesh_manager_init_base(this, type_name)
    class(mesh_manager_t), intent(inout) :: this
    character(len=*), intent(in) :: type_name

    this%type_name = type_name

  end subroutine mesh_manager_init_base

  !> Destructor for the `mesh_manager_t` (base) type.
  subroutine mesh_manager_free_base(this)
    class(mesh_manager_t), intent(inout) :: this

    call this%free_base_data
    if (allocated(this%type_name)) deallocate(this%type_name)

  end subroutine mesh_manager_free_base

  !> Destructor for the data in `mesh_manager_t` (base) type.
  subroutine mesh_manager_free_base_data(this)
    class(mesh_manager_t), intent(inout) :: this

    call this%geom_ind%free_ind()
    call this%geom_hng_h2%free_hng()
    call this%geom_hng_h4%free_hng()
    call this%conn_vrt%free()
    call this%conn_edg%free()
    call this%conn_fcs%free()

  end subroutine mesh_manager_free_base_data

end module mesh_manager
