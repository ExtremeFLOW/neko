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
  use utils, only : neko_error
  use polytope_aligned, only : polytope_aligned_t
  use ncnf_interpolation, only : ncnf_interpolation_t
  implicit none
  private

  public :: polytope_actualisation_t, element_object_t

  !> Base type for a polytope actualisation.
  !! @details This is an abstract type extending aligned polytope class
  !! with nonconforming information and operator. Note that vertices cna be
  !! hanging but have no transformation operation. This type corresponds to
  !! the realisation of an abstract objects as parts of an existing higher
  !! dimension elements.
  type, extends(polytope_aligned_t), abstract :: polytope_actualisation_t
     !> Is the object hanging
     logical :: ifhanging = .false.
     !> Alignment operator
     class(ncnf_interpolation_t), allocatable :: intp_op
   contains
  end type polytope_actualisation_t

  !> Single element object allocatable space
  type :: element_object_t
     class(polytope_actualisation_t), allocatable :: obj
  end type element_object_t

contains

end module polytope_actualisation
