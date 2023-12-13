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
 !> Abstract type for polytope actualisation class for elements building blocks
module polytope_actualisation
  use num_types, only : i4
  use polytope, only : polytope_t
  use polytope_topology, only : polytope_topology_t
  use alignment, only : alignment_t
  implicit none
  private

  public :: polytope_actualisation_t, element_object_t

  !> Base type for a polytope actualisation.
  !! @details This is an abstract type building on topology and alignment data
  !! with nonconforming information and operator. Note that vertices can be
  !! hanging but have no transformation operation. This type corresponds to
  !! the realisation of an abstract objects as parts of an existing higher
  !! dimension elements.
  type, abstract :: polytope_actualisation_t
     !> Topology polytope
     class(polytope_topology_t), pointer :: polytope
     !> Is the object aligned
     logical :: ifaligned = .false.
     !> Is the object hanging
     logical :: ifhanging = .false.
     !> Alignment operator
     class(alignment_t), allocatable :: algn_op
   contains
     !> Free aligned polytope and interpolation data
     procedure, pass(this) :: free => polytope_free
     !> Return hanging flag
     procedure, pass(this) :: ifhng => polytope_ifhng_get
     !> Return hanging node information
     procedure, pass(this) :: hng => polytope_hng_get
  end type polytope_actualisation_t

  !> Single element object allocatable space
  type :: element_object_t
     class(polytope_actualisation_t), allocatable :: obj
  end type element_object_t

  !> Abstract interface to initialise a polytope with alignment information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   hng    hanging information
  abstract interface
     subroutine polytope_actualisation_init(this, pltp, algn, hng)
       import i4
       import polytope_t
       import polytope_actualisation_t
       class(polytope_actualisation_t), intent(inout) :: this
       class(polytope_t), target, intent(in) :: pltp
       integer(i4), intent(in) :: algn
     end subroutine polytope_actualisation_init
  end interface

contains

  !> Free aligned polytope and interpolation data
  subroutine polytope_free(this)
    class(polytope_actualisation_t), intent(inout) :: this
    this%polytope => null()
    this%ifaligned = .false.
    if (allocated(this%algn_op)) deallocate(this%polytope)
    this%ifhanging = .false.
  end subroutine polytope_free

  !> @brief Get hanging flag
  !! @return   ifhng
  pure function polytope_ifhng_get(this) result(ifhng)
    class(polytope_actualisation_t), intent(in) :: this
    logical :: ifhng
       ifhng = this%ifhanging
  end function polytope_ifhng_get

  !> @brief Get hanging information
  !! @return   hng
  pure function polytope_hng_get(this) result(hng)
    class(polytope_actualisation_t), intent(in) :: this
    integer(i4) :: hng
    if (this%ifhanging) then
       hng = -1
    else
       hng = -1
    end if
  end function polytope_hng_get

end module polytope_actualisation
