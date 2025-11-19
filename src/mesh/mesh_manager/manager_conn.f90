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
module manager_conn
  use num_types, only : i4, i8, rp, dp

  implicit none
  private

  !> Abstract base type for connectivity information regarding vertices,
  !! faces and edges.
  !! @details It contains global numbering of objects
  !! information
  type, abstract, public :: manager_conn_obj_t
     !> Number of local objects
     integer(i4) :: lnum
     !> Number of owned objects
     integer(i4) :: lown
     !> Global object offset
     integer(i8) :: goff
     !> Global number of objects
     integer(i8) :: gnum
     !> Global indexing of unique objects of given type
     integer(i8), allocatable, dimension(:) :: gidx
   contains
     procedure, pass(this) :: init_data_base => manager_conn_obj_init_data_base
     procedure, pass(this) :: init_type_base => manager_conn_obj_init_type_base
     procedure, pass(this) :: free_base => manager_conn_obj_free_base
     !> Initialise data from type
     procedure(mesh_obj_init_type), pass(this), deferred :: init_type
     !> Free type
     procedure(mesh_obj_free), pass(this), deferred :: free
  end type manager_conn_obj_t

  abstract interface
     subroutine mesh_obj_init_type(this, conn)
       import manager_conn_obj_t
       class(manager_conn_obj_t), intent(inout) :: this
       class(manager_conn_obj_t), intent(inout) :: conn
     end subroutine mesh_obj_init_type

     subroutine mesh_obj_free(this)
       import manager_conn_obj_t
       class(manager_conn_obj_t), intent(inout) :: this
     end subroutine mesh_obj_free
  end interface

  !> Type for element connectivity information
  type, abstract, public :: manager_conn_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> Connectivity information for vertices, edges and faces
     class(manager_conn_obj_t), allocatable :: conn_vrt, conn_fcs, conn_edg
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

   contains
     procedure, pass(this) :: init_data_base => manager_conn_init_data_base
     procedure, pass(this) :: init_type_base => manager_conn_init_type_base
     procedure, pass(this) :: free_data_base => manager_conn_free_data_base
     procedure, pass(this) :: free_base => manager_conn_free_base
     !> Allocate types
     procedure(mesh_conn_init), pass(this), deferred :: init
     !> Initialise data from type
     procedure(mesh_conn_init_type), pass(this), deferred :: init_type
     !> Free type data
     procedure(mesh_conn_free_data), pass(this), deferred :: free_data
     !> Free type
     procedure(mesh_conn_init), pass(this), deferred :: free
  end type manager_conn_t

  abstract interface
     subroutine mesh_conn_init(this)
       import manager_conn_t
       class(manager_conn_t), intent(inout) :: this
     end subroutine mesh_conn_init

     subroutine mesh_conn_init_type(this, conn)
       import manager_conn_t
       class(manager_conn_t), intent(inout) :: this
       class(manager_conn_t), intent(inout) :: conn
     end subroutine mesh_conn_init_type

     subroutine mesh_conn_free_data(this, ifsave)
       import manager_conn_t
       class(manager_conn_t), intent(inout) :: this
       logical, optional, intent(in) :: ifsave
     end subroutine mesh_conn_free_data
  end interface

contains

  !> Initialise connectivity object type
  !! @param[in]    lnum    local number of objects
  !! @param[in]    lown    number of owned objects
  !! @param[in]    goff    global object offset
  !! @param[in]    gnum    global number of objects
  !! @param[inout] gidx    global object index
  subroutine manager_conn_obj_init_data_base(this, lnum, lown, goff, gnum, gidx)
    class(manager_conn_obj_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, lown
    integer(i8), intent(in) :: goff, gnum
    integer(i8), allocatable, dimension(:), intent(inout) :: gidx

    call this%free_base()

    this%lnum = lnum
    this%lown = lown
    this%goff = goff
    this%gnum = gnum

    if (allocated(gidx)) call move_alloc(gidx, this%gidx)

  end subroutine manager_conn_obj_init_data_base

  !> Initialise connectivity object type based on another connectivity type
  !! @param[inout] conn   connectivity object data
  subroutine manager_conn_obj_init_type_base(this, conn)
    class(manager_conn_obj_t), intent(inout) :: this
    class(manager_conn_obj_t), intent(inout) :: conn

    call this%free_base()

    this%lnum = conn%lnum
    this%lown = conn%lown
    this%goff = conn%goff
    this%gnum = conn%gnum

    if (allocated(conn%gidx)) call move_alloc(conn%gidx, this%gidx)

  end subroutine manager_conn_obj_init_type_base

  !> Free connectivity object type
  subroutine manager_conn_obj_free_base(this)
    class(manager_conn_obj_t), intent(inout) :: this

    this%lnum = 0
    this%lown = 0
    this%goff = 0
    this%gnum = 0

    if (allocated(this%gidx)) deallocate(this%gidx)

  end subroutine manager_conn_obj_free_base

  !> Initialise connectivity information
  !! @param[in]    tdim    topological dimension
  !! @param[in]    nel     local element number
  !! @param[inout] vmap    element vertex mapping
  !! @param[inout] fmap    element face mapping
  !! @param[inout] emap    element edge mapping
  !! @param[in]   ifsave      save component types
  subroutine manager_conn_init_data_base(this, tdim, nel, vmap, fmap, emap, &
       ifsave)
    class(manager_conn_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    integer(i4), allocatable, dimension(:,:), intent(inout) :: vmap, fmap, emap
    logical, optional, intent(in) :: ifsave

    if (present(ifsave)) then
       call this%free_data_base(ifsave)
    else
       call this%free_data_base()
    end if

    this%tdim = tdim
    this%nel = nel
    ! we work with hex/quad only and there is no difference between
    ! topology and geometrical dimensions
    this%nvrt = 2**tdim
    this%nfcs = 2 * tdim
    this%nedg = 12 * (tdim -2) !tdim * 2**(tdim - 1)

    if (allocated(vmap)) call move_alloc(vmap, this%vmap)
    if (allocated(fmap)) call move_alloc(fmap, this%fmap)
    if (allocated(emap)) call move_alloc(emap, this%emap)

  end subroutine manager_conn_init_data_base

  !> Initialise connectivity type based on another connectivity type
  !! @param[inout] conn   connectivity data
  subroutine manager_conn_init_type_base(this, conn)
    class(manager_conn_t), intent(inout) :: this
    class(manager_conn_t), intent(inout) :: conn

    call this%free_data_base()

    if (allocated(this%conn_vrt) .and. allocated(conn%conn_vrt)) &
         call this%conn_vrt%init_type(conn%conn_vrt)
    if (allocated(this%conn_fcs) .and. allocated(conn%conn_fcs)) &
         call this%conn_fcs%init_type(conn%conn_fcs)
    if (allocated(this%conn_edg) .and. allocated(conn%conn_edg)) &
         call this%conn_edg%init_type(conn%conn_edg)

    this%tdim = conn%tdim
    this%nel = conn%nel
    this%nvrt = conn%nvrt
    this%nfcs = conn%nfcs
    this%nedg = conn%nedg

    if (allocated(conn%vmap)) call move_alloc(conn%vmap, this%vmap)
    if (allocated(conn%fmap)) call move_alloc(conn%fmap, this%fmap)
    if (allocated(conn%emap)) call move_alloc(conn%emap, this%emap)

  end subroutine manager_conn_init_type_base

  !> Free connectivity data
  !! @param[in]   ifsave      save component types
  subroutine manager_conn_free_data_base(this, ifsave)
    class(manager_conn_t), intent(inout) :: this
    logical, optional, intent(in) :: ifsave
    logical :: ifsavel

    ifsavel = .false.
    if (present(ifsave)) ifsavel = ifsave

    if (.not.ifsavel) then
       if (allocated(this%conn_vrt)) call this%conn_vrt%free()
       if (allocated(this%conn_fcs)) call this%conn_fcs%free()
       if (allocated(this%conn_edg)) call this%conn_edg%free()
    end if

    this%tdim = 0
    this%nel = 0
    this%nvrt = 0
    this%nfcs = 0
    this%nedg = 0

    if (allocated(this%vmap)) deallocate(this%vmap)
    if (allocated(this%fmap)) deallocate(this%fmap)
    if (allocated(this%emap)) deallocate(this%emap)

  end subroutine manager_conn_free_data_base

  !> Free connectivity type
  subroutine manager_conn_free_base(this)
    class(manager_conn_t), intent(inout) :: this

    call this%free_data_base()
    if (allocated(this%conn_vrt)) deallocate(this%conn_vrt)
    if (allocated(this%conn_fcs)) deallocate(this%conn_fcs)
    if (allocated(this%conn_edg)) deallocate(this%conn_edg)

  end subroutine manager_conn_free_base

end module manager_conn
