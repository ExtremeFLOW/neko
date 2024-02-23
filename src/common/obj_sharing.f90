! Copyright (c) 2018-2024, The Neko Authors
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
!> Abstract type for defining object ownership, mapping and sharing
!! @details These abstract types are bases for different type object lists
!! e.g, points, elements or topology objects.
module obj_sharing
  use num_types, only : i4
  use utils, only : neko_error
  use comm, only : NEKO_COMM, MPI_INTEGER, MPI_SUM

  implicit none
  private

  public :: obj_ownership_t, obj_mapping_t, obj_sharing_t

  !> Base type for the lists with defined object ownership
  type, abstract :: obj_ownership_t
     !> Number of local objects
     integer(i4), private :: lnum_ = -1
     !> Global number of objects
     integer(i4), private :: gnum_ = -1
     !> Global object offset based on local owned objects
     integer(i4), private :: goff_ = -1
     !> Number of local owned objects
     integer(i4), private :: lown_ = -1
     !> List of local numbers of local owned objects
     integer(i4), dimension(:), allocatable :: owned_list
     !> Global to local mapping
   contains
     !> Set space for ownership data
     procedure, pass(this) :: init_own => obj_ownership_init_own
     !> Free ownership data
     procedure, pass(this) :: free_own => obj_ownership_free_own
     !> Global to local id mapping
     procedure, pass(this) :: glb_loc_id => obj_ownership_glb_loc_id
     !> Local to global id mapping
     procedure(obj_ownership_loc_glb_id), pass(this), deferred :: loc_glb_id
  end type obj_ownership_t

  !> Base type for the lists with defined object mapping
  type, extends(obj_ownership_t), abstract :: obj_mapping_t
   contains
  end type obj_mapping_t

  !> Base type for the lists with defined object sharing
  type, extends(obj_mapping_t), abstract :: obj_sharing_t
   contains
  end type obj_sharing_t

  !> Returns global object id based on it's local position
  !! @parameter[in]   lid   object local id position
  !! @return gid
  abstract interface
     function obj_ownership_loc_glb_id(this, lid) result(gid)
       import i4
       import obj_ownership_t
       class(obj_ownership_t), intent(in) :: this
       integer(i4), intent(in) :: lid
       integer(i4) :: gid
     end function obj_ownership_loc_glb_id
  end interface

contains

  !> Set space for ownership data
  !! @parameter[in]    lnum     local number of objects
  !! @parameter[inout] gid_list list of object global ids
  !! @parameter[in]    lown     number of local owned
  !! @parameter[inout] gid_list list of owned object
  subroutine obj_ownership_init_own(this, lnum, gid_list, lown, own_list)
    class(obj_ownership_t), intent(inout) :: this
    integer(i4), intent(in) :: lnum, lown
    integer(i4), dimension(:), allocatable, intent(inout) :: gid_list, own_list
    integer(i4) :: ierr

    call this%free_own()

    ! CHANGE IT TO CALCULATE OWNERSHIP INTERNALLY
    ! ADD FLAG TO SKIP OWNERSHIP CALCULATION AND ASSUME ALL OBJECTS ARE OWNED
    if (lnum > 0) then
       this%lnum_ = lnum
       if (allocated(gid_list)) then
          ! place to work on global to local mapping
          !call move_alloc(own_list, this%owned_list)
       else
          call neko_error('Missing object global id list.')
       end if
    else
       this%lnum_ = 0
    end if
    if (lown > 0) then
       if (lown > this%lnum_) &
            & call neko_error('Inconsistent owned object number.')
       this%lown_ = lown
       if (allocated(own_list)) then
          call move_alloc(own_list, this%owned_list)
       else
          call neko_error('Missing owned object list.')
       end if
    else
       this%lown_ = 0
    end if

    call MPI_Allreduce(this%lown_, this%gnum_, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    this%goff_ = 0
    call MPI_Exscan(this%lown_, this%goff_, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  end subroutine obj_ownership_init_own

  !> Free ownership data
  subroutine obj_ownership_free_own(this)
    class(obj_ownership_t), intent(inout) :: this
    this%lnum_ = -1
    this%gnum_ = -1
    this%goff_ = -1
    this%lown_ = -1
    if (allocated(this%owned_list)) deallocate(this%owned_list)
  end subroutine obj_ownership_free_own

  !> Returns local object position based on it's global id
  !! @parameter[in]   gid   object local id position
  !! @return lid
  function obj_ownership_glb_loc_id(this, gid) result(lid)
    class(obj_ownership_t), intent(in) :: this
    integer(i4), intent(in) :: gid
    integer(i4) :: lid
  end function obj_ownership_glb_loc_id

end module obj_sharing
