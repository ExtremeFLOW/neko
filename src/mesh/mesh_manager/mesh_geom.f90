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
module mesh_geom
  use num_types, only : i4, i8, rp, dp

  implicit none
  private

  !> Type base abstract type for geometrical independent nodes.
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
     procedure, pass(this) :: init_data => mesh_geom_node_ind_init_data
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
     procedure, pass(this) :: init_data => mesh_geom_node_hng_init_data
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
     procedure, pass(this) :: init_data => mesh_geom_init_data
     procedure, pass(this) :: init_type => mesh_geom_init_type
     procedure, pass(this) :: free => mesh_geom_free
  end type mesh_geom_t

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
  subroutine mesh_geom_node_ind_init_data(this, lown, lshr, loff, lnum, gdim, &
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

  end subroutine mesh_geom_node_ind_init_data

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
  subroutine mesh_geom_node_hng_init_data(this, lnum, gdim, ndep, gidx, &
       lmap, coord)
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

  end subroutine mesh_geom_node_hng_init_data

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
  subroutine mesh_geom_init_data(this, tdim, nel, vmap)
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

  end subroutine mesh_geom_init_data

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

end module mesh_geom
