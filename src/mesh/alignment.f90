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
!> Abstract type for alignment operators
!! @details To understand alignment, consider, two neighbouring quads in a 2D
!! mesh, sharing a single edge. The edge is globally defined by two points, p1
!! and p2. However, in the internal numbering of the two quads, the same point
!! can be numbered as either the first point or the second, depending on which
!! one gets the higher internal number. Similarly, for the other point.
!! Alignment basically refers to the relative orientation of the edge
!! in the internal numbering systems of the two quads it belongs two.
!! For the case of the edge, only two alignments are possible: either p1 is
!! is the "first" point in both quads or in one of the quads it is the "second".
!! For higher dimensional polytopes the situation is more involved, for example,
!! for the quad, there are 8 possible numberings.
!! In the code, the various polytope alignments are assigned an integer number.
!! Each alignment type for each polytope has a designated alignment operator
!! implemented in the code, all descending from `alignment_t` defined in this
!! module. The alignment operators are applied to the data living on the
!! polytope and transform it according to the alignment. So, for example, for
!! the edge there are two operators: the "identity" one, which doesn't do
!! anything, and another one, which simply reverses the order of the data. Each
!! operator type implements both the corresponding transformation and its
!! inverse. Moreover, there are implementations of the same transformation for
!! different types of the underlying data: integers (i4 and i8) and the rp real.
!! It is possible to apply the transformation to only part of the data by
!! specifying the starting array index from which the transformation is applied.
!! Finally, subroutines for checking the relative alignment of two polytopes
!! are provided. Currently, alignment is implemented for edges and quads only.
module alignment
  use num_types, only : i4, i8, rp
  implicit none
  private

  public :: alignment_t

  !> Base type for an alignment operator for a given abstract polytope
  !! @note This class for now works with the square arrays only assuming
  !! @a lx = @a ly = @a lz. However, this can be easily updated to the general
  !! case.
  type, abstract :: alignment_t
     !> Relative polytope alignment
     integer(i4), private :: alignment_ = -1
   contains
     !> Return relative polytope alignment
     procedure, pass(this) :: algn => alignment_algn_get
     !> Set relative polytope alignment
     procedure, pass(this) :: set_algn => alignment_algn_set
     !> Is transformation identity
     procedure(alignment_ifidentity), pass(this), deferred :: ifid
     !> Direct transformation; different types
     procedure(transform_i4), nopass, deferred :: trns_i4
     procedure(transform_i8), nopass, deferred :: trns_i8
     procedure(transform_rp), nopass, deferred :: trns_rp
     !> Inverse transformation; different types
     procedure(transform_i4), nopass, deferred :: trns_inv_i4
     procedure(transform_i8), nopass, deferred :: trns_inv_i8
     procedure(transform_rp), nopass, deferred :: trns_inv_rp
  end type alignment_t

  !> Abstract interface for function returning identity flag
  !! @return   ifid
  abstract interface
     pure function alignment_ifidentity(this) result(ifid)
       import :: alignment_t
       class(alignment_t), intent(in) :: this
       logical :: ifid
     end function alignment_ifidentity
  end interface

  !> Abstract interface for various transformations; single integer type
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations. There is as well a starting
  !! position for operations @a ist to be able to avoid vertices.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  abstract interface
     pure subroutine transform_i4(vec, n1, n2, ist)
       import i4
       integer(i4), intent(in) :: n1, n2, ist
       integer(i4), dimension(n1, n2), intent(inout) :: vec
     end subroutine transform_i4
  end interface

  !> Abstract interface for various transformations; double integer type
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations. There is as well a starting
  !! position for operations @a ist to be able to avoid vertices.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  abstract interface
     pure subroutine transform_i8(vec, n1, n2, ist)
       import i4
       import i8
       integer(i4), intent(in) :: n1, n2, ist
       integer(i8), dimension(n1, n2), intent(inout) :: vec
     end subroutine transform_i8
  end interface

  !> Abstract interface for various transformations; real type
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations. There is as well a starting
  !! position for operations @a ist to be able to avoid vertices.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  abstract interface
     pure subroutine transform_rp(vec, n1, n2, ist)
       import i4
       import rp
       integer(i4), intent(in) :: n1, n2, ist
       real(rp), dimension(n1, n2), intent(inout) :: vec
     end subroutine transform_rp
  end interface

contains

  !> @brief Get polytope alignment
  !! @return   algn
  pure function alignment_algn_get(this) result(algn)
    class(alignment_t), intent(in) :: this
    integer(i4) :: algn
    algn = this%alignment_
  end function alignment_algn_get

  !> @brief Set polytope alignment
  !! @parameter[in]   algn     polytope alignment
  pure subroutine alignment_algn_set(this, algn)
    class(alignment_t), intent(inout) :: this
    integer(i4), intent(in) :: algn
    this%alignment_ = algn
  end subroutine alignment_algn_set

end module alignment
