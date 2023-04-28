! Copyright (c) 2020-2022, The Neko Authors
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
!> Defines output formatting types for floating-point arrays
!! @details This module defines two data types, vector_output_t and
!! matrix_output_t, which help in writing floating-point arrays
!! (1 or 2 dimensional) in text (.dat) files.
module dat_outputs
  use num_types
  use math
  use logger
  implicit none

  type, public :: vector_output_t
     integer :: n                                    !< Size (number of entries)
     real(kind=rp), pointer, dimension(:) :: data => null() !< Array containing
                                                         ! values of the vector
   contains
     procedure, pass(this) :: init => vector_output_init
     procedure, pass(this) :: free => vector_output_free
  end type vector_output_t

  type, public :: matrix_output_t
     integer :: ni !< Number of rows
     integer :: nj !< Number of columns
     real(kind=rp), pointer, dimension(:,:) :: data => null() !< 2D array containing the values
   contains
     procedure, pass(this) :: init => matrix_output_init
     procedure, pass(this) :: free => matrix_output_free
  end type matrix_output_t

contains

  !> Initializes a vector_output_t with another vector
  subroutine vector_output_init(this, vec, n)
    class(vector_output_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), target, dimension(n), intent(in) :: vec
    
    this%data => vec
    this%n = n

  end subroutine vector_output_init

  !> Deallocate vector_output_t
  subroutine vector_output_free(this)
    class(vector_output_t), intent(inout) :: this
    
    nullify(this%data)
    this%n = -1

  end subroutine vector_output_free

  !> Initializes a matrix_output_t with another matrix
  subroutine matrix_output_init(this, vec, ni, nj)
    class(matrix_output_t), intent(inout) :: this
    integer, intent(in) :: ni, nj
    real(kind=rp), target, dimension(ni,nj), intent(in) :: vec
    
    this%ni = ni
    this%nj = nj
    this%data => vec

  end subroutine matrix_output_init
  
  !> Deallocate matrix_output_t
  subroutine matrix_output_free(this)
    class(matrix_output_t), intent(inout) :: this
    
    nullify(this%data)
    this%ni = -1
    this%nj = -1

  end subroutine matrix_output_free

end module dat_outputs
