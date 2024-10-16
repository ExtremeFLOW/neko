! Copyright (c) 2024, The Neko Authors
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
!> Implements a mesh intersection detector
module intersection_detector
  use num_types, only : dp
  use aabb_tree, only : aabb_tree_t
  use mesh, only : mesh_t
  use hex, only : hex_t
  use point, only : point_t
  use utils, only : neko_error
  use stack, only : stack_i4_t
  use aabb
  implicit none
  private

  !> Intersector detector type
  !!
  !! This type is used to detect intersections between a mesh and a given point
  !! in space. This allow to rapid reduction in the number of elements to be
  !! tested for intersection. The detector is based on a private AABB tree and
  !! is linked to a mesh.
  !!
  !! @note Please note that the current implementation is limited to hexahedral
  !! elements, and we assume that the mesh is static.
  type, public :: intersect_detector_t
     type(mesh_t), private, pointer :: msh => null()
     type(aabb_tree_t), private :: search_tree
   contains
     procedure, pass(this) :: init => intersect_detector_init
     procedure, pass(this) :: free => intersect_detector_free

     generic :: overlap => overlap_stack, overlap_array
     procedure, pass(this) :: overlap_stack => intersect_detector_overlap_stack
     procedure, pass(this) :: overlap_array => intersect_detector_overlap_array
  end type intersect_detector_t

contains

  !> Initialise an intersector detector for a given mesh @a msh
  subroutine intersect_detector_init(this, msh, padding)
    class(intersect_detector_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    real(kind=dp), intent(in), optional :: padding

    call this%free()
    this%msh => msh

    call this%search_tree%init(msh%nelv)
    call this%search_tree%build(msh%elements, padding)

    if (this%search_tree%get_size() .ne. (msh%nelv)) then
       call neko_error("Error building the search tree.")
    end if

  end subroutine intersect_detector_init

  !> Destroy an intersector detector
  subroutine intersect_detector_free(this)
    class(intersect_detector_t), intent(inout) :: this

    if (associated(this%msh)) then
       nullify(this%msh)
    end if

    !> @todo cleanup the aabb tree
    call this%search_tree%free()

  end subroutine intersect_detector_free

  !> Computes the overlap between elements and a given point @a p
  subroutine intersect_detector_overlap_stack(this, p, overlaps)
    class(intersect_detector_t), intent(inout) :: this
    type(point_t), intent(in) :: p
    type(stack_i4_t), intent(inout) :: overlaps

    call this%search_tree%query_overlaps(p, overlaps)

  end subroutine intersect_detector_overlap_stack

  !> Computes the overlap between elements and a given point @a p
  subroutine intersect_detector_overlap_array(this, p, overlaps)
    class(intersect_detector_t), intent(inout) :: this
    type(point_t), intent(in) :: p
    integer, dimension(:), allocatable, intent(inout) :: overlaps

    call this%search_tree%query_overlaps(p, overlaps)

  end subroutine intersect_detector_overlap_array

end module intersection_detector
