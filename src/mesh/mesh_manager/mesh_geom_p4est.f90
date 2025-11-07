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
!> Implementation of the mesh geometry type for p4est mesh manager
module mesh_geom_p4est
  use num_types, only : i4, i8, rp, dp
  use mesh_geom, only : mesh_geom_node_t, mesh_geom_t

  implicit none
  private

  !> p4est type extension for geometrical independent nodes.
  !! @details It adds ownership information
  type, extends(mesh_geom_node_t), public :: mesh_geom_node_ind_p4est_t
     !> Number of owned nodes
     integer(i4) :: lown
     !> Number of owned shared nodes
     integer(i4) :: lshr
     !> Local offset of owned nodes
     integer(i4) :: loff
     !> Node owner (MPI rank)
     integer(i4), allocatable, dimension(:) :: ndown
   contains
     procedure, pass(this) :: init_data => mesh_geom_node_ind_init_data_p4est
     procedure, pass(this) :: init_type => mesh_geom_node_ind_init_type_p4est
     procedure, pass(this) :: free => mesh_geom_node_ind_free_p4est
  end type mesh_geom_node_ind_p4est_t

  !> p4est type extension for geometrical hanging nodes
  !! @details p4est adds mapping to the independent face/edge vertices. This
  !! mapping is not unique and will depend on the element position in a tree.
  !! There are two types of hanging nodes h2 (dependent on 2 independent nodes;
  !! 2D face and 3D edge) and h4 (dependent on 4 independent nodes; 3D face).
  type, extends(mesh_geom_node_t), public :: mesh_geom_node_hng_p4est_t
     !> Number of independent nodes in the dependency list
     integer(i4) :: ndep
     !> Local hanging to independent node mapping
     integer(i4), allocatable, dimension(:,:) :: lmap
   contains
     procedure, pass(this) :: init_data => mesh_geom_node_hng_init_data_p4est
     procedure, pass(this) :: init_type => mesh_geom_node_hng_init_type_p4est
     procedure, pass(this) :: free => mesh_geom_node_hng_free_p4est
  end type mesh_geom_node_hng_p4est_t

  !> Type for element geometry information
  !! @details This type extends a base type with hanging nodes.
  type, extends(mesh_geom_t), public :: mesh_geom_p4est_t
     !> Geometrical h2-type hanging nodes
     type(mesh_geom_node_hng_p4est_t) :: geom_hng_edg
     !> Geometrical h4-type hanging nodes
     type(mesh_geom_node_hng_p4est_t) :: geom_hng_fcs
     ! Mapping including hanging nodes:
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
   contains
     procedure, pass(this) :: init => mesh_geom_init_p4est
     procedure, pass(this) :: init_data => mesh_geom_init_data_p4est
     procedure, pass(this) :: init_type => mesh_geom_init_type_p4est
     procedure, pass(this) :: free_data => mesh_geom_free_data_p4est
     procedure, pass(this) :: free => mesh_geom_free_p4est
  end type mesh_geom_p4est_t

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
  subroutine mesh_geom_node_ind_init_data_p4est(this, lown, lshr, loff, lnum, &
       gdim, gidx, ndown, coord)
    ! argument list
    class(mesh_geom_node_ind_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: lown, lshr, loff, lnum, gdim
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: ndown
    real(kind=dp), allocatable, dimension(:,:), intent(inout) :: coord

    call this%free()
    call this%init_data_base(lnum, gdim, gidx, coord)

    this%lown = lown
    this%lshr = lshr
    this%loff = loff

    if (allocated(ndown)) call move_alloc(ndown, this%ndown)

  end subroutine mesh_geom_node_ind_init_data_p4est

  !> Initialise independent nodes type based on another node type
  !! @param[inout] node   node data
  subroutine mesh_geom_node_ind_init_type_p4est(this, node)
    ! argument list
    class(mesh_geom_node_ind_p4est_t), intent(inout) :: this
    class(mesh_geom_node_t), intent(inout) :: node

    call this%free()
    call this%init_type_base(node)

    select type (node)
    type is(mesh_geom_node_ind_p4est_t)
       this%lown = node%lown
       this%lshr = node%lshr
       this%loff = node%loff

       if (allocated(node%ndown)) call move_alloc(node%ndown, this%ndown)

    end select

  end subroutine mesh_geom_node_ind_init_type_p4est

  !> Free independent nodes type
  subroutine mesh_geom_node_ind_free_p4est(this)
    ! argument list
    class(mesh_geom_node_ind_p4est_t), intent(inout) :: this

    call this%free_base()

    this%lown = 0
    this%lshr = 0
    this%loff = 0

    if (allocated(this%ndown)) deallocate(this%ndown)

  end subroutine mesh_geom_node_ind_free_p4est

  !> Initialise hanging nodes type
  !! @param[in]    lnum    local number of nodes
  !! @param[in]    gdim    geometrical dimension
  !! @param[in]    ndep    node dependency
  !! @param[inout] gidx    global node index
  !! @param[inout] lmap    node mapping to independent node
  !! @param[inout] coord   node coordinates
  subroutine mesh_geom_node_hng_init_data_p4est(this, lnum, gdim, ndep, gidx, &
       lmap, coord)
    ! argument list
    class(mesh_geom_node_hng_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, gdim, ndep
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), allocatable, dimension(:, :), intent(inout) :: lmap
    real(kind=dp), allocatable, dimension(:,:), intent(inout) :: coord

    call this%free()
    call this%init_data_base(lnum, gdim, gidx, coord)

    this%ndep = ndep

    if (allocated(lmap)) call move_alloc(lmap, this%lmap)

  end subroutine mesh_geom_node_hng_init_data_p4est

  !> Initialise hanging nodes type based on another node type
  !! @param[inout] node   node data
  subroutine mesh_geom_node_hng_init_type_p4est(this, node)
    ! argument list
    class(mesh_geom_node_hng_p4est_t), intent(inout) :: this
    class(mesh_geom_node_t), intent(inout) :: node

    call this%free()
    call this%init_type_base(node)

    select type (node)
    type is(mesh_geom_node_hng_p4est_t)
       this%ndep = node%ndep

       if (allocated(node%lmap)) call move_alloc(node%lmap, this%lmap)
    end select

  end subroutine mesh_geom_node_hng_init_type_p4est

  !> Free hanging nodes type
  subroutine mesh_geom_node_hng_free_p4est(this)
    ! argument list
    class(mesh_geom_node_hng_p4est_t), intent(inout) :: this

    call this%free_base()

    this%ndep = 0

    if (allocated(this%lmap)) deallocate(this%lmap)

  end subroutine mesh_geom_node_hng_free_p4est

  !> Allocate types
  subroutine mesh_geom_init_p4est(this)
    class(mesh_geom_p4est_t), intent(inout) :: this

    if (allocated(this%geom_ind))then
       call this%geom_ind%free()
       deallocate(this%geom_ind)
    end if
    allocate(mesh_geom_node_ind_p4est_t::this%geom_ind)

  end subroutine mesh_geom_init_p4est

  !> Initialise geometry type
  !! @param[in]    tdim    topological mesh dimension
  !! @param[in]    nel     local element number
  !! @param[inout] vmap    element vertices to node mapping
  subroutine mesh_geom_init_data_p4est(this, tdim, nel, vmap)
    class(mesh_geom_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    integer(i4), allocatable, dimension(:,:), intent(inout) :: vmap

    call this%free_data()

    call this%init_data_base(tdim, nel, vmap)

  end subroutine mesh_geom_init_data_p4est

  !> Initialise geometry type based on another geometry type
  !! @param[inout] geom   geometry data
  subroutine mesh_geom_init_type_p4est(this, geom)
    ! argument list
    class(mesh_geom_p4est_t), intent(inout) :: this
    class(mesh_geom_t), intent(inout) :: geom

    call this%free_data()

    call this%init_type_base(geom)
    select type (geom)
    type is(mesh_geom_p4est_t)
       call this%geom_hng_edg%init_type(geom%geom_hng_edg)
       call this%geom_hng_fcs%init_type(geom%geom_hng_fcs)
    end select

  end subroutine mesh_geom_init_type_p4est

  !> Free geometry data
  subroutine mesh_geom_free_data_p4est(this)
    class(mesh_geom_p4est_t), intent(inout) :: this

    call this%free_data_base()
    call this%geom_hng_edg%free()
    call this%geom_hng_fcs%free()

  end subroutine mesh_geom_free_data_p4est

  !> Free geometry
  subroutine mesh_geom_free_p4est(this)
    class(mesh_geom_p4est_t), intent(inout) :: this

    call this%free_base()
    call this%geom_hng_edg%free()
    call this%geom_hng_fcs%free()

  end subroutine mesh_geom_free_p4est

end module mesh_geom_p4est
