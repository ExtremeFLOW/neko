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
 !> Abstract type for abstract aligned polytope class for mesh topology
module polytope_aligned
  use num_types, only : i4
  use utils, only : neko_error
  use polytope, only : polytope_t
  use alignment, only : alignment_t
  implicit none
  private

  public :: polytope_aligned_t, topology_object_t

  !> Base type for an abstract aligned polytope
  !! @details This is an abstract type combining polytope class and alignment
  !! operator. Note that vertices have no alignment. This type corresponds to
  !! building blocs of the higher dimension abstract objects that build mesh
  !! topology.
  type, abstract :: polytope_aligned_t
     !> Polytope pointer
     class(polytope_t), pointer :: polytope => null()
     !> Is the object aligned
     logical :: ifaligned = .false.
     !> Alignment operator
     class(alignment_t), allocatable :: algn_op
   contains
     !> Free polytope and aligned data
     procedure, pass(this) :: free_algn => polytope_free_algn
     !> Return a pointer to the polytope
     procedure, pass(this) :: polyp => polytope_ptr
     !> Return alignment value
     procedure, pass(this) :: algn => polytope_algn_get
  end type polytope_aligned_t

  !> Single topology object allocatable space
  type :: topology_object_t
     class(polytope_aligned_t), allocatable :: obj
  end type topology_object_t

contains

  !> Free polytope and aligned data
  subroutine polytope_free_algn(this)
    class(polytope_aligned_t), intent(inout) :: this
    this%polytope => null()
    this%ifaligned = .false.
    if (allocated(this%algn_op)) deallocate(this%algn_op)
  end subroutine polytope_free_algn

  !> @brief Return pointer to the polytoppe
  !! @parameter[out]  poly   polytope pointer
  subroutine polytope_ptr(this, poly)
    class(polytope_aligned_t), intent(in) :: this
    class(polytope_t), pointer, intent(out) :: poly
    poly => this%polytope
  end subroutine polytope_ptr

  !> @brief Get polytope alignment
  !! @return   algn
  pure function polytope_algn_get(this) result(algn)
    class(polytope_aligned_t), intent(in) :: this
    integer(i4) :: algn
    if (this%ifaligned) then
       algn = this%algn_op%algn()
    else
       algn = -1
    end if
  end function polytope_algn_get

end module polytope_aligned
