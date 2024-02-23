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
!>
module point_list
  use num_types, only : i4
  use utils, only : neko_error
  use obj_sharing, only : obj_mapping_t
  use point, only : point_t

  implicit none
  private

  public :: point_list_t

  !> Type for the geometrical points list
  type, extends(obj_mapping_t) :: point_list_t
     !> Geometrical point list
     type(point_t), dimension(:), allocatable :: points
   contains
     !> Free list data
     procedure, pass(this) :: free => point_list_free
     !> Local to global id mapping
     procedure, pass(this) :: loc_glb_id => point_list_loc_glb_id
  end type point_list_t

contains

  !> Free list data
  subroutine point_list_free(this)
    class(point_list_t), intent(inout) :: this
    integer(i4) :: il

    call this%free_own()

    if (allocated(this%points)) then
       deallocate(this%points)
    end if
  end subroutine point_list_free

  !> Returns global point id based on it's local position
  !! @parameter[in]   lid   object local id position
  !! @return gid
  function point_list_loc_glb_id(this, lid) result(gid)
    class(point_list_t), intent(in) :: this
    integer(i4), intent(in) :: lid
    integer(i4) :: gid
    gid = this%points(lid)%id()
  end function point_list_loc_glb_id

end module point_list
