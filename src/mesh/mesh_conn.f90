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
!> Implementation of the mesh connectivity type
module mesh_conn
  use num_types, only : i4, i8
  use utils, only : neko_error

  implicit none
  private

  !> Type for connectivity information regarding vertices, faces and edges.
  !! @details It contains global numbering of objects and element mapping
  !! information
  type, public :: mesh_conn_obj_t
     !> Number of local objects
     integer(i4) :: lnum
     !> Global number of objects
     integer(i8) :: gnum
     !> Global indexing of unique objects of given type
     integer(i8), allocatable, dimension(:) :: gidx
     !> Flag indicating sharing with other MPI ranks
     logical, allocatable, dimension(:) :: share
     !> Local number of elements
     integer(i4) :: nel
     !> Number of objects per element
     integer(i4) :: nobj
     !> Element to object mapping
     integer(i4), allocatable, dimension(:,:) :: map
     !> Is object alignment directly specified
     logical :: ifalgn
     !> Edge alignment
     integer(i4), allocatable, dimension(:,:) :: algn
     !> Are hanging object directly specified
     logical :: ifhang
     !> Hanging object; position of hanging object, otherwise -1
     integer(i4), allocatable, dimension(:,:) :: hang
   contains
     procedure, pass(this) :: init => mesh_conn_obj_init
     procedure, pass(this) :: free => mesh_conn_obj_free
  end type mesh_conn_obj_t

  !> Type for element connectivity information
  type, public :: mesh_conn_t
     !> Topological mesh dimension
     integer(i4) :: tdim
     !> Local number of elements
     integer(i4) :: nel
     !> Connectivity information for vertices, edges and faces
     type(mesh_conn_obj_t) :: vrt, fcs, edg
     !> Are elements with hanging object directly specified
     logical :: ifhang
     !> Elements with hanging object
     logical, allocatable, dimension(:) :: hang
   contains
     procedure, pass(this) :: init => mesh_conn_init
     procedure, pass(this) :: free => mesh_conn_free
  end type mesh_conn_t

contains

  !> Initialise object connectivity information
  !! @param[in]    lnum     local number of objects
  !! @param[in]    gnum     global number of objects
  !! @param[in]    nel      local number of elements
  !! @param[in]    nobj     number of objects per element
  !! @param[in]    gidx     object global id
  !! @param[in]    share    sharing with other ranks
  !! @param[in]    map      element to objects mapping
  !! @param[in]    algn     object alignment information
  !! @param[in]    hang     object hanging information
  subroutine mesh_conn_obj_init(this, lnum, gnum, nel, nobj, gidx, share, map, &
       algn, hang)
    class(mesh_conn_obj_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, nel, nobj
    integer(i8), intent(in) :: gnum
    integer(i8), dimension(:), intent(in) :: gidx
    logical, dimension(:), intent(in) :: share
    integer(i4), dimension(:, :), intent(in) :: map
    integer(i4), dimension(:, :), optional, intent(in) :: algn, hang

    call this%free()

    this%lnum = lnum
    this%gnum = gnum
    this%nel = nel
    this%nobj = nobj

    allocate(this%gidx(lnum), this%share(lnum), this%map(nobj, nel))
    ! sanity check
    if (lnum .ne. size(gidx) .or. lnum .ne. size(share) .or. &
         nel .ne. size(map, 2) .or. nobj .ne. size(map, 1)) &
         call neko_error('Inconsistent array sizes; conn_obj')
    this%gidx(:) = gidx(:)
    this%share(:) = share(:)
    this%map(:, :) = map(:, :)

    if (present(algn)) then
       ! sanity check
       if (nel .ne. size(algn, 2) .or. nobj .ne. size(algn, 1)) &
            call neko_error('Inconsistent array sizes; conn_obj%algn')
       this%ifalgn = .true.
       allocate(this%algn(nobj, nel))
       this%algn(:, :) = algn(:, :)
    end if

    if (present(hang)) then
       ! sanity check
       if (nel .ne. size(hang, 2) .or. nobj .ne. size(hang, 1)) &
            call neko_error('Inconsistent array sizes; conn_obj%hang')
       this%ifhang = .true.
       allocate(this%hang(nobj, nel))
       this%hang(:, :) = hang(:, :)
    end if

  end subroutine mesh_conn_obj_init

  !> Free object connectivity information
  subroutine mesh_conn_obj_free(this)
    class(mesh_conn_obj_t), intent(inout) :: this

    this%lnum = 0
    this%gnum = 0
    this%nel = 0
    this%nobj = 0
    this%ifalgn = .false.
    this%ifhang = .false.

    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%share)) deallocate(this%share)
    if (allocated(this%map)) deallocate(this%map)
    if (allocated(this%algn)) deallocate(this%algn)
    if (allocated(this%hang)) deallocate(this%hang)

  end subroutine mesh_conn_obj_free

  !> Initialise mesh connectivity information
  !! @param[in]    tdim     topological dimension
  !! @param[in]    nel      local number of elements
  !! @param[in]    hang     element hanging information
  subroutine mesh_conn_init(this, tdim, nel, hang)
    class(mesh_conn_t), intent(inout) :: this
    integer(i4), intent(in) :: tdim, nel
    logical, dimension(:), optional, intent(in) :: hang

    call this%free()

    this%tdim = tdim
    this%nel = nel

    if (present(hang)) then
       ! sanity check
       if (nel .ne. size(hang)) &
            call neko_error('Inconsistent array sizes; conn%hang')
       this%ifhang = .true.
       allocate(this%hang(nel))
       this%hang(:) = hang(:)
    end if

  end subroutine mesh_conn_init

  !> Free mesh connectivity information
  subroutine mesh_conn_free(this)
    class(mesh_conn_t), intent(inout) :: this

    call this%vrt%free()
    call this%fcs%free()
    call this%edg%free()

    this%tdim = 0
    this%nel = 0
    this%ifhang = .false.

    if (allocated(this%hang)) deallocate(this%hang)

  end subroutine mesh_conn_free



end module mesh_conn
