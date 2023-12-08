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
!> Abstract type for nonconforming interpolation operators
module ncnf_interpolation
use num_types, only : i4, rp
  implicit none
  private

  public :: ncnf_interpolation_t

  !> Base type for an nonconforming interpolation operator for a given abstract
  !! polytope
  type, abstract :: ncnf_interpolation_t
     !> Hanging node information
     integer(i4), private :: hanging_ = -1
   contains
     !> Return hanging node information
     procedure, pass(this) :: hng => interpolation_hng_get
     !> Set relative polytope alignment
     procedure, pass(this) :: set_hng => interpolation_hng_set
     !> Patent-child interpolation
     procedure(transform_rp), nopass, deferred :: intp
     !> Transposed interpolation
     procedure(transform_rp), nopass, deferred :: intpT
     !> Initialise interpolation data
     procedure(operator_init), pass(this), deferred :: set_jmat
     !> Free interpolation data
     procedure(operator_free), pass(this), deferred :: free_jmat
  end type ncnf_interpolation_t

  !> Abstract interface for various transformations; real type
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  abstract interface
     pure subroutine transform_rp(vec, n1, n2)
       import i4
       import rp
       integer(i4), intent(in) :: n1, n2
       real(rp), dimension(n1, n2), intent(inout) :: vec
     end subroutine transform_rp
  end interface

  !> Abstract interface to initialisation the interpolation data
  !! @notice It is a common interface for 1D and 2D operations, so there
  !! are two dimensions @a lr and @a ls. This routine can be called
  !! after space_t gets initialised.
  !! @parameter[in]      lr, ls   array sizes for r and s dimensions
  abstract interface
     subroutine operator_init(this, lr, ls)
       import i4
       import rp
       import :: ncnf_interpolation_t
       class(ncnf_interpolation_t), intent(inout) :: this
       integer(i4), intent(in) :: lr, ls
     end subroutine operator_init
  end interface

  !> Abstract interface to free the interpolation data
  abstract interface
     subroutine operator_free(this)
       import :: ncnf_interpolation_t
       class(ncnf_interpolation_t), intent(inout) :: this
     end subroutine operator_free
  end interface

contains

  !> @brief Get hanging information
  !! @return   hng
  pure function interpolation_hng_get(this) result(hng)
    class(ncnf_interpolation_t), intent(in) :: this
    integer(i4) :: hng
    hng = this%hanging_
  end function interpolation_hng_get

  !> @brief Set hanging information
  !! @parameter[in]   hng     hanging information
  pure subroutine interpolation_hng_set(this, hng)
    class(ncnf_interpolation_t), intent(inout) :: this
    integer(i4), intent(in) :: hng
    this%hanging_ = hng
  end subroutine interpolation_hng_set

end module ncnf_interpolation
