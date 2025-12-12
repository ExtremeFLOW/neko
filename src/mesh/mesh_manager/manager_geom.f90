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
!> Implementation of the mesh geometry type for mesh manager
module manager_geom
  use num_types, only : i4, i8, rp, dp

  implicit none
  private

  !> Base abstract type for geometrical independent nodes.
  !! @details Geometrical independent nodes are globally unique nodes not
  !! sharing physical coordinates and global id and not located on the
  !! nonconforming interfaces. It is a minimal set of points needed to build
  !! a conforming, consistent element mesh based on verices only (no curvature
  !! information). This type does not include full connectivity information, as
  !! periodicity and nonconforming interfaces are omitted. This type is used
  !! as well as a base type for the hanging nodes.  Hanging nodes are not
  !! independent ones. They are located at the centre of the nonconforming
  !! face or edge (have unique physical coordinates).
  type, abstract, public :: manager_geom_node_t
     !> Local number of nodes
     integer(i4) :: lnum
     !> Global indexing of unique nodes
     integer(i8), allocatable, dimension(:) :: gidx
     !> Geometrical mesh dimension
     integer(i4) :: gdim
     !> Physical node coordinates
     real(kind=dp), allocatable, dimension(:,:) :: coord
   contains
     procedure, pass(this) :: init_data_base => manager_geom_node_init_data_base
     procedure, pass(this) :: init_type_base => manager_geom_node_init_type_base
     procedure, pass(this) :: free_base => manager_geom_node_free_base
     !> Initialise data from type
     procedure(mesh_node_init_type), pass(this), deferred :: init_type
     !> Free type
     procedure(mesh_node_free), pass(this), deferred :: free
  end type manager_geom_node_t

  abstract interface
     subroutine mesh_node_init_type(this, node)
       import manager_geom_node_t
       class(manager_geom_node_t), intent(inout) :: this
       class(manager_geom_node_t), intent(inout) :: node
     end subroutine mesh_node_init_type

     subroutine mesh_node_free(this)
       import manager_geom_node_t
       class(manager_geom_node_t), intent(inout) :: this
     end subroutine mesh_node_free
  end interface

  !> Type for element geometry information
  !! @details This type collects all the geometrical information including
  !! independent nodes lists and combines it with element vertices for node
  !! mapping.
  type, abstract, public :: manager_geom_t
     !> Did we import simple or complete mesh information set
     logical :: ifcomplete
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> Geometrical independent nodes
     class(manager_geom_node_t), allocatable :: ind
     !> local number of elements
     integer(i4) :: nel
     !> Number of vertices per element
     integer(i4) :: nvrt
     ! Mapping for:
     ! 1<= vmap(iv,iel) <= nin - independent node
     ! vmap uses symmetric vertex notation with (r,s,t) being a local
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
     !> Local mapping of element vertices to the nodes; complete global info
     integer(i4), allocatable, dimension(:, :) :: vmap
     !> Coordinates of element vertices; just local info
     real(dp), allocatable, dimension(:, :, :) :: vcoord
   contains
     procedure, pass(this) :: init_data_base => manager_geom_init_data_base
     procedure, pass(this) :: init_simple_base => manager_geom_init_simple_base
     procedure, pass(this) :: init_type_base => manager_geom_init_type_base
     procedure, pass(this) :: free_data_base => manager_geom_free_data_base
     procedure, pass(this) :: free_base => manager_geom_free_base
     !> Allocate types
     procedure(mesh_geom_init), pass(this), deferred :: init
     !> Initialise data from type
     procedure(mesh_geom_init_type), pass(this), deferred :: init_type
     !> Free type_data
     procedure(mesh_geom_free_data), pass(this), deferred :: free_data
     !> Free type
     procedure(mesh_geom_init), pass(this), deferred :: free
  end type manager_geom_t

  abstract interface
     subroutine mesh_geom_init(this)
       import manager_geom_t
       class(manager_geom_t), intent(inout) :: this
     end subroutine mesh_geom_init

     subroutine mesh_geom_init_type(this, geom)
       import manager_geom_t
       class(manager_geom_t), intent(inout) :: this
       class(manager_geom_t), intent(inout) :: geom
     end subroutine mesh_geom_init_type

     subroutine mesh_geom_free_data(this, ifsave)
       import manager_geom_t
       class(manager_geom_t), intent(inout) :: this
       logical, optional, intent(in) :: ifsave
     end subroutine mesh_geom_free_data
  end interface

contains

  !> Initialise nodes type data
  !! @param[in]    lnum    local number of nodes
  !! @param[in]    gdim    geometrical dimension
  !! @param[inout] gidx    global node index
  !! @param[inout] coord   node coordinates
  subroutine manager_geom_node_init_data_base(this, lnum, gdim, gidx, coord)
    class(manager_geom_node_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, gdim
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    real(kind=dp), allocatable, dimension(:,:), intent(inout) :: coord

    call this%free_base()

    this%lnum = lnum
    this%gdim = gdim

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)
    if (allocated(coord)) call move_alloc(coord, this%coord)

  end subroutine manager_geom_node_init_data_base

  !> Initialise nodes type based on another node type
  !! @param[inout] node   node data
  subroutine manager_geom_node_init_type_base(this, node)
    class(manager_geom_node_t), intent(inout) :: this, node

    call this%free_base()

    this%lnum = node%lnum
    this%gdim = node%gdim

    if (allocated(node%gidx)) call move_alloc(node%gidx, this%gidx)
    if (allocated(node%coord)) call move_alloc(node%coord, this%coord)

  end subroutine manager_geom_node_init_type_base

  !> Free node type
  subroutine manager_geom_node_free_base(this)
    class(manager_geom_node_t), intent(inout) :: this

    this%lnum = 0
    this%gdim = 0

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%coord)) deallocate(this%coord)

  end subroutine manager_geom_node_free_base

  !> Initialise complete geometry type
  !! @param[in]    tdim    topological mesh dimension
  !! @param[in]    nel     local element number
  !! @param[inout] vmap    element vertices to node mapping
  !! @param[in]   ifsave      save component types
  subroutine manager_geom_init_data_base(this, tdim, nel, vmap, ifsave)
    class(manager_geom_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    integer(i4), allocatable, dimension(:, :), intent(inout) :: vmap
    logical, optional, intent(in) :: ifsave

    if (present(ifsave)) then
       call this%free_data_base(ifsave)
    else
       call this%free_data_base()
    end if

    this%ifcomplete = .true.
    this%tdim = tdim
    this%nel = nel
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nvrt = 2**tdim

    if (allocated(vmap)) call move_alloc(vmap, this%vmap)

  end subroutine manager_geom_init_data_base

  !> Initialise simple geometry type
  !! @param[in]    tdim    topological mesh dimension
  !! @param[in]    nel     local element number
  !! @param[inout] vcoord  coordinates of element vertices
  subroutine manager_geom_init_simple_base(this, tdim, nel, vcoord)
    class(manager_geom_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    real(dp), allocatable, dimension(:, :, :), intent(inout) :: vcoord

    call this%free_data_base()

    this%ifcomplete = .false.
    this%tdim = tdim
    this%nel = nel
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nvrt = 2**tdim

    if (allocated(vcoord)) call move_alloc(vcoord, this%vcoord)

  end subroutine manager_geom_init_simple_base

  !> Initialise geometry type based on another geometry type
  !! @param[inout] geom   geometry data
  subroutine manager_geom_init_type_base(this, geom)
    class(manager_geom_t), intent(inout) :: this, geom

    call this%free_data_base()

    if (allocated(this%ind) .and. allocated(geom%ind)) &
            call this%ind%init_type(geom%ind)

    this%ifcomplete = geom%ifcomplete
    this%tdim = geom%tdim
    this%nel = geom%nel
    this%nvrt = geom%nvrt

    if (allocated(geom%vmap)) call move_alloc(geom%vmap, this%vmap)
    if (allocated(geom%vcoord)) call move_alloc(geom%vcoord, this%vcoord)

  end subroutine manager_geom_init_type_base

  !> Free geometry data
  !! @param[in]   ifsave      save component types
  subroutine manager_geom_free_data_base(this, ifsave)
    class(manager_geom_t), intent(inout) :: this
    logical, optional, intent(in) :: ifsave
    logical :: ifsavel

    ifsavel = .false.
    if (present(ifsave)) ifsavel = ifsave

    if (.not.ifsavel) then
       if (allocated(this%ind)) call this%ind%free()
    end if

    this%ifcomplete = .false.
    this%tdim = 0
    this%nel = 0
    this%nvrt = 0

    if (allocated(this%vmap)) deallocate(this%vmap)
    if (allocated(this%vcoord)) deallocate(this%vcoord)

  end subroutine manager_geom_free_data_base

  !> Free geometry type
  subroutine manager_geom_free_base(this)
    class(manager_geom_t), intent(inout) :: this

    call this%free_data_base()
    if (allocated(this%ind)) deallocate(this%ind)

  end subroutine manager_geom_free_base

end module manager_geom
