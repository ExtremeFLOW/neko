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
module element_list
  use num_types, only : i4
  use utils, only : neko_error
  use obj_sharing, only : obj_ownership_t
  use element_new, only : mesh_element_t

  implicit none
  private

  public :: element_list_t

  !> Type for the mesh element list
  type, extends(obj_ownership_t) :: element_list_t
     !> Mesh element list
     type(mesh_element_t), dimension(:), allocatable :: elements
   contains
     !> Free list data
     procedure, pass(this) :: free => element_list_free
     !> Local to global id mapping
     procedure, pass(this) :: loc_glb_id => element_list_loc_glb_id
  end type element_list_t

contains

  !> Free list data
  subroutine element_list_free(this)
    class(element_list_t), intent(inout) :: this
    integer(i4) :: il

    call this%free_own()

    if (allocated(this%elements)) then
       do il = 1, size(this%elements)
          if (allocated(this%elements(il)%obj)) then
             call this%elements(il)%obj%free()
             deallocate(this%elements(il)%obj)
          end if
       end do
       deallocate(this%elements)
    end if
  end subroutine element_list_free

  !> Returns global element id based on it's local position
  !! @parameter[in]   lid   object local id position
  !! @return gid
  function element_list_loc_glb_id(this, lid) result(gid)
    class(element_list_t), intent(in) :: this
    integer(i4), intent(in) :: lid
    integer(i4) :: gid
    gid = this%elements(lid)%obj%id()
  end function element_list_loc_glb_id

end module element_list
