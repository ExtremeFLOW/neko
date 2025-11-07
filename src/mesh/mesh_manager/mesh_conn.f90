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
!> Implementation of the mesh connectivity type for mesh manager
module mesh_conn
  use num_types, only : i4, i8, rp, dp

  implicit none
  private

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
     procedure, pass(this) :: init_data => mesh_conn_obj_init_data
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
     procedure, pass(this) :: init_data => mesh_conn_init_data
     procedure, pass(this) :: init_type => mesh_conn_init_type
     procedure, pass(this) :: free => mesh_conn_free
  end type mesh_conn_t

contains

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
  subroutine mesh_conn_obj_init_data(this, lnum, lown, goff, gnum, nrank, &
       nshare, gidx, rank, share, off)
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

  end subroutine mesh_conn_obj_init_data

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
  subroutine mesh_conn_init_data(this, tdim, nel, vmap, fmap, falgn, emap, &
       ealgn, hngel, hngfc, hnged)
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

  end subroutine mesh_conn_init_data

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

end module mesh_conn
