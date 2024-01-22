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
 !> Abstract type for abstract aligned polytope class
module polytope_aligned
  use num_types, only : i4
  use polytope, only : polytope_t
  use alignment, only : alignment_t
  implicit none
  private

  public :: polytope_aligned_t

  !> Base type for an abstract aligned polytope
  !! @details This is an abstract type combining polytope class and alignment
  !! operator. Note that vertices have no alignment. This type corresponds to
  !! building blocs of the higher dimension abstract objects the topology or
  !! mesh consists of.
  type, abstract :: polytope_aligned_t
     !> Polytope pointer
     class(polytope_t), pointer :: polytope => null()
     !> Is the object aligned
     logical, private :: ifaligned_ = .false.
     !> Alignment operator
     class(alignment_t), allocatable :: algn_op
   contains
     !> Free polytope and aligned data
     procedure, pass(this) :: freea => polytope_free
     !> Initialise alignment data
     procedure, pass(this) :: init_data => polytope_init_data
     !> Return a pointer to the polytope
     procedure, pass(this) :: polyp => polytope_ptr
     !> Return alignment flag
     procedure, pass(this) :: ifalgn => polytope_ifalgn_get
     !> Return alignment value
     procedure, pass(this) :: algn => polytope_algn_get
     !> Test equality
     procedure, pass(this) :: equal => polytope_equal
     !> Test equality and find alignment
     procedure(polytope_aligned_equal_algn), pass(this), deferred :: equal_algn
     !> Test alignment
     procedure(polytope_aligned_test), pass(this), deferred :: test
  end type polytope_aligned_t

  !> Check polytope equality and return alignment
  !! @parameter[in]    other  polytope
  !! @parameter[out]   equal  polytope equality
  !! @parameter[out]   algn   alignment information
  abstract interface
     subroutine polytope_aligned_equal_algn(this, other, equal, algn)
       import i4
       import polytope_t
       import polytope_aligned_t
       class(polytope_aligned_t), intent(in) :: this
       class(polytope_t), intent(in) :: other
       logical, intent(out) :: equal
       integer(i4), intent(out) :: algn
     end subroutine polytope_aligned_equal_algn
  end interface

  !> Test alignment
  !! @parameter[in]   other   polytope
  !! @return ifalgn
  abstract interface
     function polytope_aligned_test(this, other) result(ifalgn)
       import polytope_t
       import polytope_aligned_t
       class(polytope_aligned_t), intent(in) :: this
       class(polytope_t), intent(in) :: other
       logical :: ifalgn
     end function polytope_aligned_test
  end interface

contains

  !> Free polytope and aligned data
  subroutine polytope_free(this)
    class(polytope_aligned_t), intent(inout) :: this
    this%polytope => null()
    this%ifaligned_ = .false.
    if (allocated(this%algn_op)) deallocate(this%algn_op)
  end subroutine polytope_free

  !> Initialise general data
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   ifalgn if non identity alignment
  subroutine polytope_init_data(this, pltp, ifalgn)
    class(polytope_aligned_t), intent(inout) :: this
    class(polytope_t), target, intent(in) :: pltp
    logical, intent(in)  :: ifalgn

    this%polytope => pltp
    this%ifaligned_ = ifalgn
  end subroutine polytope_init_data

  !> @brief Return pointer to the polytope
  !! @parameter[out]  poly   polytope pointer
  subroutine polytope_ptr(this, poly)
    class(polytope_aligned_t), intent(in) :: this
    class(polytope_t), pointer, intent(out) :: poly
    poly => this%polytope
  end subroutine polytope_ptr

  !> @brief Get polytope alignment flag
  !! @return   ifalgn
  pure function polytope_ifalgn_get(this) result(ifalgn)
    class(polytope_aligned_t), intent(in) :: this
    logical :: ifalgn
    ifalgn = this%ifaligned_
  end function polytope_ifalgn_get

  !> @brief Get polytope alignment
  !! @return   algn
  pure function polytope_algn_get(this) result(algn)
    class(polytope_aligned_t), intent(in) :: this
    integer(i4) :: algn
    if (allocated(this%algn_op)) then
       algn = this%algn_op%algn()
    else
       algn = -1
    end if
  end function polytope_algn_get

  !> Test equality
  !! @parameter[in]   other   polytope
  !! @return equal
  function polytope_equal(this, other) result(equal)
    class(polytope_aligned_t), intent(in) :: this
    class(polytope_t), intent(in) :: other
    logical :: equal
    integer(i4) :: itmp
    call this%equal_algn(other, equal, itmp)
  end function polytope_equal

end module polytope_aligned
