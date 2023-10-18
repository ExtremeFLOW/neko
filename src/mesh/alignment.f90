! Copyright (c) 2018-2023, The Neko Authors
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
!> Abstract type for polytope alignment
module alignment
  use num_types, only : i2, i4
  use utils, only : neko_error
  implicit none
  private

  public :: alignment_t

  !> Base type for polytope alignment
  type, abstract :: alignment_t
     !> number of different operations excluding identity
     integer(i2), private :: noperation_ = -1
     !> Relative polytope alignment
     integer(i2), private :: alignment_ = -1
   contains
     !> setter for number of operations
     procedure, pass(this) :: set_nop => noperation_set
     !> getter for number of operations
     procedure, pass(this) :: nop => noperation_get
     !> setter for alignment
     procedure, pass(this) :: set_algn => alignment_set
     !> getter for alignment
     procedure, pass(this) :: algn => alignment_get
     !> initialisation routine setting noperation is required
     procedure(alignment_init), pass(this), deferred :: init
  end type alignment_t

  abstract interface
     subroutine alignment_init(this)
       import :: alignment_t
       class(alignment_t), intent(inout) :: this
     end subroutine alignment_init
  end interface

contains
  !> @brief Set number of different operations excluding identity (marked by 0)
  !! @parameter[in]   noperation    number of operations
  subroutine  noperation_set(this, noperation)
    class(alignment_t), intent(inout) :: this
    integer(i4), intent(in) :: noperation
    this%noperation_ = noperation
    return
  end subroutine noperation_set

  !> @brief Get number of operations
  !! @return   noperation
  pure function noperation_get(this) result(noperation)
    class(alignment_t), intent(in) :: this
    integer(i4) :: noperation
    noperation = this%noperation_
  end function noperation_get

  !> @brief Set relative polytope alignment
  !! @parameter[in]   alignment       relative polytope alignment
  subroutine  alignment_set(this, alignment)
    class(alignment_t), intent(inout) :: this
    integer(i4), intent(in) :: alignment
    if ((alignment < 0).or.(alignment > this%noperation_)) &
         & call neko_error('Not proper alignment.')
    this%alignment_ = alignment
    return
  end subroutine alignment_set

  !> @brief Get relative polytope alignment
  !! @return   alignment
  pure function alignment_get(this) result(alignment)
    class(alignment_t), intent(in) :: this
    integer(i4) :: alignment
    alignment = this%alignment_
  end function alignment_get

end module alignment
