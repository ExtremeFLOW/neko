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
  !! with nonconforming and interpolation information. Note that vertices can be
  !! hanging but have no interpolation operation. This type corresponds to
  !! the realisation of an abstract objects as parts of an existing higher
  !! dimension elements.
  type, abstract :: polytope_actualisation_t
     !> Topology polytope
     class(polytope_topology_t), pointer :: polytope
     !> Is the object aligned
     logical, private :: ifaligned_ = .false.
     !> Alignment operator
     class(alignment_t), allocatable :: algn_op
     !> Is there any interpolation operator active (excluding identity)
     logical, private :: ifinterpolation_ = .false.
     !> Hanging information
     integer(i4), private :: hanging_ = 0
     !> position in the higher dimension polytope
     integer(i4), private :: position_ = -1
   contains
     !> Free aligned polytope and interpolation data
     procedure, pass(this) :: free => polytope_free
     !> Return alignment flag
     procedure, pass(this) :: ifalgn => polytope_ifalgn_get
     !> Return interpolation flag
     procedure, pass(this) :: ifintp => polytope_ifintp_get
     !> Return hanging node information
     procedure, pass(this) :: hng => polytope_hng_get
     !> Return position in the higher dimension object
     procedure, pass(this) :: pos => polytope_pos_get
     !> Initialise an polytope actualisation
     procedure(polytope_actualisation_init), pass(this), deferred :: init
     !> Return higher dimension object direction for local direction @a r
     procedure(polytope_dir), pass(this), deferred :: dirr
     !> Return higher dimension object direction for local direction @a s
     procedure(polytope_dir), pass(this), deferred :: dirs
     !> Test equality and find alignment
     procedure(polytope_equal_algn), pass(this), deferred :: equal_algn
     !> Test alignment
     procedure(polytope_alignment_test), pass(this), deferred :: test
  end type polytope_actualisation_t

  !> Single element object allocatable space
  type :: element_object_t
     class(polytope_actualisation_t), allocatable :: obj
  end type element_object_t

  !> Abstract interface to initialise a polytope with alignment information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  abstract interface
     subroutine polytope_actualisation_init(this, pltp, algn, ifint, hng, pos)
       import i4
       import polytope_topology_t
       import polytope_actualisation_t
       class(polytope_actualisation_t), intent(inout) :: this
       class(polytope_topology_t), target, intent(in) :: pltp
       integer(i4), intent(in) :: algn, hng, pos
       logical :: ifint
     end subroutine polytope_actualisation_init
  end interface

  !> Return higher dimension object direction
  !! @return dir
  abstract interface
     function polytope_dir(this) result(dir)
       import i4
       import polytope_actualisation_t
       class(polytope_actualisation_t), intent(in) :: this
       integer(i4) :: dir
     end function polytope_dir
  end interface

  !> Abstract interface to check polytope equality and return alignment
  !! @parameter[in]    pltp   polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  abstract interface
     subroutine polytope_equal_algn(this, pltp, equal, algn)
       import i4
       import polytope_t
       import polytope_actualisation_t
       class(polytope_actualisation_t), intent(in) :: this
       class(polytope_t), intent(in) :: pltp
       logical, intent(out) :: equal
       integer(i4), intent(out) :: algn
     end subroutine polytope_equal_algn
  end interface

  !> Test alignment
  !! @parameter[in]   pltp   polytope
  !! @return ifalgn
  abstract interface
     function polytope_alignment_test(this, pltp) result(ifalgn)
       import polytope_t
       import polytope_actualisation_t
       class(polytope_actualisation_t), intent(in) :: this
       class(polytope_t), intent(in) :: pltp
       logical :: ifalgn
     end function polytope_alignment_test
  end interface

contains

  !> Free aligned polytope and interpolation data
  subroutine polytope_free(this)
    class(polytope_actualisation_t), intent(inout) :: this
    this%polytope => null()
    this%ifaligned_ = .false.
    if (allocated(this%algn_op)) deallocate(this%polytope)
    this%ifinterpolation_ = .false.
    this%hanging_ = 0
    this%position_ = -1
  end subroutine polytope_free

  !> @brief Get alignment flag
  !! @return   ifalgn
  pure function polytope_ifalgn_get(this) result(ifalgn)
    class(polytope_actualisation_t), intent(in) :: this
    logical :: ifalgn
    ifalgn = this%ifaligned_
  end function polytope_ifalgn_get

     !> @brief Get interpolation flag
  !! @return   ifintp
  pure function polytope_ifintp_get(this) result(ifintp)
    class(polytope_actualisation_t), intent(in) :: this
    logical :: ifintp
    ifintp = this%ifinterpolation_
  end function polytope_ifintp_get

  !> @brief Get hanging information
  !! @return   hng
  pure function polytope_hng_get(this) result(hng)
    class(polytope_actualisation_t), intent(in) :: this
    integer(i4) :: hng
    if (this%ifinterpolation_) then
       hng = this%hanging_
    else
       hng = -1
    end if
  end function polytope_hng_get

  !> @brief Get position information
  !! @return   pos
  pure function polytope_pos_get(this) result(pos)
    class(polytope_actualisation_t), intent(in) :: this
    integer(i4) :: pos
    pos = this%position_
  end function polytope_pos_get

end module polytope_actualisation
