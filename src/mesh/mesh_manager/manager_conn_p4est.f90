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
!> Implementation of the mesh connectivity type for p4est mesh manager
module manager_conn_p4est
  use num_types, only : i4, i8, rp, dp
  use utils, only : neko_error
  use manager_conn, only : manager_conn_obj_t, manager_conn_t

  implicit none
  private

  !> Extended type for connectivity information regarding vertices, faces and
  !! edges.
  !! @details It adds sharing information for communication
  type, extends(manager_conn_obj_t), public :: manager_conn_obj_p4est_t
     !> Number of MPI ranks sharing objects
     integer(i4) :: nrank
     !> Number of shared objects
     integer(i4) :: nshare
     !> List of ranks sharing objects
     integer(i4), allocatable, dimension(:) :: rank
     !> List of shared objects
     integer(i4), allocatable, dimension(:) :: share
     !> Offset in the share list
     integer(i4), allocatable, dimension(:) :: off
   contains
     procedure, pass(this) :: init_data => manager_conn_obj_init_data_p4est
     procedure, pass(this) :: init_type => manager_conn_obj_init_type_p4est
     procedure, pass(this) :: free => manager_conn_obj_free_p4est
  end type manager_conn_obj_p4est_t

  !> Type for element connectivity information
  type, extends(manager_conn_t), public :: manager_conn_p4est_t
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
     procedure, pass(this) :: init => manager_conn_init_p4est
     procedure, pass(this) :: init_data => manager_conn_init_data_p4est
     procedure, pass(this) :: init_type => manager_conn_init_type_p4est
     procedure, pass(this) :: free_data => manager_conn_free_data_p4est
     procedure, pass(this) :: free => manager_conn_free_p4est
  end type manager_conn_p4est_t

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
  subroutine manager_conn_obj_init_data_p4est(this, lnum, lown, goff, gnum, &
       nrank, nshare, gidx, rank, share, off)
    class(manager_conn_obj_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, lown, nrank, nshare
    integer(i8), intent(in) :: goff, gnum
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx
    integer(i4), allocatable, dimension(:), intent(inout) :: rank, share, off

    ! sanity check
    if ((lnum .ne. size(gidx)) .or. (nrank .ne. size(rank)) .or. &
         (nshare .ne. size(share)) .or. (nrank + 1 .ne. size(off))) &
         call neko_error('Inconsistent array sizes; p4est%conn_obj')

    call this%free()
    call this%init_data_base(lnum, lown, goff, gnum, gidx)

    this%nrank = nrank
    this%nshare = nshare

    if (allocated(rank)) call move_alloc(rank, this%rank)
    if (allocated(share)) call move_alloc(share, this%share)
    if (allocated(off)) call move_alloc(off, this%off)

  end subroutine manager_conn_obj_init_data_p4est

  !> Initialise connectivity object type based on another connectivity type
  !! @param[inout] conn   connectivity object data
  subroutine manager_conn_obj_init_type_p4est(this, conn)
    class(manager_conn_obj_p4est_t), intent(inout) :: this
    class(manager_conn_obj_t), intent(inout) :: conn

    call this%free()
    call this%init_type_base(conn)

    select type (conn)
    type is(manager_conn_obj_p4est_t)
       this%nrank = conn%nrank
       this%nshare = conn%nshare

       if (allocated(conn%rank)) call move_alloc(conn%rank, this%rank)
       if (allocated(conn%share)) call move_alloc(conn%share, this%share)
       if (allocated(conn%off)) call move_alloc(conn%off, this%off)
    end select

  end subroutine manager_conn_obj_init_type_p4est

  !> Free connectivity object type
  subroutine manager_conn_obj_free_p4est(this)
    class(manager_conn_obj_p4est_t), intent(inout) :: this

    call this%free_base()

    this%nrank = 0
    this%nshare = 0

    if (allocated(this%rank)) deallocate(this%rank)
    if (allocated(this%share)) deallocate(this%share)
    if (allocated(this%off)) deallocate(this%off)

  end subroutine manager_conn_obj_free_p4est

  !> Allocate types
  subroutine manager_conn_init_p4est(this)
    class(manager_conn_p4est_t), intent(inout) :: this

    if (allocated(this%vrt))then
       call this%vrt%free()
       deallocate(this%vrt)
    end if
    allocate(manager_conn_obj_p4est_t::this%vrt)

    if (allocated(this%fcs))then
       call this%fcs%free()
       deallocate(this%fcs)
    end if
    allocate(manager_conn_obj_p4est_t::this%fcs)

    if (allocated(this%edg))then
       call this%edg%free()
       deallocate(this%edg)
    end if
    allocate(manager_conn_obj_p4est_t::this%edg)

  end subroutine manager_conn_init_p4est

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
  !! @param[in]   ifsave      save component types
  subroutine manager_conn_init_data_p4est(this, tdim, nel, vmap, fmap, falgn, &
       emap, ealgn, hngel, hngfc, hnged, ifsave)
    class(manager_conn_p4est_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    integer(i4), allocatable, dimension(:), intent(inout) :: hngel
    integer(i4), allocatable, dimension(:,:), intent(inout) :: vmap, fmap, &
         falgn, emap, ealgn, hngfc, hnged
    logical, optional, intent(in) :: ifsave
    integer(i4) :: nvert, nface, nedge

    nvert = 2**tdim
    nface = 2 * tdim
    nedge = 12 * (tdim - 2)

    ! sanity check
    if ((nvert .ne. size(vmap, 1)) .or. (nel .ne. size(vmap, 2)) .or. &
         (nface .ne. size(fmap, 1)) .or. (nel .ne. size(fmap, 2)) .or. &
         (nface .ne. size(falgn, 1)) .or. (nel .ne. size(falgn, 2)) .or. &
         (nedge .ne. size(emap, 1)) .or. (nel .ne. size(emap, 2)) .or. &
         (nedge .ne. size(ealgn, 1)) .or. (nel .ne. size(ealgn, 2)) .or. &
         (nel .ne. size(hngel)) .or. &
         (nface .ne. size(hngfc, 1)) .or. (nel .ne. size(hngfc, 2)) .or. &
         (nedge .ne. size(hnged, 1)) .or. (nel .ne. size(hnged, 2))) &
         call neko_error('Inconsistent array sizes; p4est%conn')

    if (present(ifsave)) then
       call this%free_data(ifsave)
       call this%init_data_base(tdim, nel, vmap, fmap, emap, ifsave)
    else
       call this%free_data()
       call this%init_data_base(tdim, nel, vmap, fmap, emap)
    end if

    if (allocated(falgn)) call move_alloc(falgn, this%falgn)
    if (allocated(ealgn)) call move_alloc(ealgn, this%ealgn)
    if (allocated(hngel)) call move_alloc(hngel, this%hngel)
    if (allocated(hngfc)) call move_alloc(hngfc, this%hngfc)
    if (allocated(hnged)) call move_alloc(hnged, this%hnged)

  end subroutine manager_conn_init_data_p4est

  !> Initialise connectivity type based on another connectivity type
  !! @param[inout] conn   connectivity data
  subroutine manager_conn_init_type_p4est(this, conn)
    class(manager_conn_p4est_t), intent(inout) :: this
    class(manager_conn_t), intent(inout) :: conn

    call this%free_data()

    call this%init_type_base(conn)

    select type (conn)
    type is(manager_conn_p4est_t)
       if (allocated(conn%falgn)) call move_alloc(conn%falgn, this%falgn)
       if (allocated(conn%ealgn)) call move_alloc(conn%ealgn, this%ealgn)
       if (allocated(conn%hngel)) call move_alloc(conn%hngel, this%hngel)
       if (allocated(conn%hngfc)) call move_alloc(conn%hngfc, this%hngfc)
       if (allocated(conn%hnged)) call move_alloc(conn%hnged, this%hnged)
    end select

  end subroutine manager_conn_init_type_p4est

  !> Free connectivity information
  !! @param[in]   ifsave      save component types
  subroutine manager_conn_free_data_p4est(this, ifsave)
    class(manager_conn_p4est_t), intent(inout) :: this
    logical, optional, intent(in) :: ifsave

    if (present(ifsave)) then
       call this%free_data_base(ifsave)
    else
       call this%free_data_base()
    end if

    if (allocated(this%falgn)) deallocate(this%falgn)
    if (allocated(this%ealgn)) deallocate(this%ealgn)
    if (allocated(this%hngel)) deallocate(this%hngel)
    if (allocated(this%hngfc)) deallocate(this%hngfc)
    if (allocated(this%hnged)) deallocate(this%hnged)

  end subroutine manager_conn_free_data_p4est

  !> Free connectivity information
  subroutine manager_conn_free_p4est(this)
    class(manager_conn_p4est_t), intent(inout) :: this

    call this%free_base()

    if (allocated(this%falgn)) deallocate(this%falgn)
    if (allocated(this%ealgn)) deallocate(this%ealgn)
    if (allocated(this%hngel)) deallocate(this%hngel)
    if (allocated(this%hngfc)) deallocate(this%hngfc)
    if (allocated(this%hnged)) deallocate(this%hnged)

  end subroutine manager_conn_free_p4est

end module manager_conn_p4est
