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
!> Edge alignment operators
module alignment_edge
  use num_types, only : i4, i8, rp
  use utils, only : neko_error
  use alignment, only : alignment_t
  implicit none
  private

  public :: alignment_edge_init, alignment_edge_find, alignment_edge_I_t, &
       & alignment_edge_P_t

  !> number of operations different from identity
  integer(i4), public, parameter :: NEKO_EDGE_NOPERATION = 1

  !> Edge identity (I) transformation type
  type, extends(alignment_t) :: alignment_edge_I_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_edge_I
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_edge_I_i4
     procedure, nopass :: trns_i8 => alignment_edge_I_i8
     procedure, nopass :: trns_rp => alignment_edge_I_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_edge_I_i4
     procedure, nopass :: trns_inv_i8 => alignment_edge_I_i8
     procedure, nopass :: trns_inv_rp => alignment_edge_I_rp
  end type alignment_edge_I_t

  !> Edge row permutation (P) transformation type
  type, extends(alignment_t) :: alignment_edge_P_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_edge_P
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_edge_P_i4
     procedure, nopass :: trns_i8 => alignment_edge_P_i8
     procedure, nopass :: trns_rp => alignment_edge_P_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_edge_P_i4
     procedure, nopass :: trns_inv_i8 => alignment_edge_P_i8
     procedure, nopass :: trns_inv_rp => alignment_edge_P_rp
  end type alignment_edge_P_t

contains

  !> @brief Allocate a single alignment operator
  !! @parameter[in]      algn  relative edge alignment
  !! @parameter[inout]   trns  alignment operator
  subroutine alignment_edge_init(algn, trns)
    integer(i4), intent(in) :: algn
    class(alignment_t), allocatable, intent(inout) :: trns

    if (allocated(trns)) then
       deallocate(trns)
    end if

    select case (algn)
    case(0) ! identity
       allocate(alignment_edge_I_t :: trns)
       call trns%set_algn(algn)
    case(1) ! permutation
       allocate(alignment_edge_P_t :: trns)
       call trns%set_algn(algn)
    case default
       call neko_error('Wrong edge alignment')
    end select

  end subroutine alignment_edge_init

  !> @brief Find relative alignment for edges
  !! @parameter[out]     equal       array equality flag
  !! @parameter[out]     algn        relative edge alignment
  !! @parameter[inout]   edg1, edg2  edges
  !! @parameter[in]      n1, n2      dimensions
  subroutine alignment_edge_find(equal, algn, edg1, edg2, n1, n2)
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    integer(i4), intent(in) :: n1, n2
    integer(i4), dimension(n1, n2), intent(inout) :: edg1, edg2
    type(alignment_edge_P_t) :: op_P

    algn = -1
    equal = all(edg1 .eq. edg2)
    if (equal) then
       algn = 0
       return
    end if

    call op_P%trns_inv_i4(edg1, n1, n2, 1)
    equal = all(edg1 .eq. edg2)
    if (equal) then
       algn = 1
       return
    end if

  end subroutine alignment_edge_find

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_edge_I(this) result(ifid)
    class(alignment_edge_I_t), intent(in) :: this
    logical :: ifid
    ifid = .true.
  end function ifidentity_edge_I

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_edge_P(this) result(ifid)
    class(alignment_edge_P_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_edge_P

  !> @brief Identity transformation, single integer array
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_edge_I_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec
  end subroutine alignment_edge_I_i4

  !> @brief Permutation transformation, single integer
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_edge_P_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: istart, il, itmp1, itmp2
    integer(i4) :: iedg

    itmp1 = n1 + 1
    do il = ist, n1 / 2
       itmp2 = itmp1 - il
       iedg = vec(il, n2)
       vec(il, n2) = vec(itmp2, n2)
       vec(itmp2, n2) = iedg
    end do
  end subroutine alignment_edge_P_i4

  !> @brief Identity transformation, double integer array
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_edge_I_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec
  end subroutine alignment_edge_I_i8

  !> @brief Permutation transformation, double integer, full array
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_edge_P_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, itmp1, itmp2
    integer(i8) :: iedg

    itmp1 = n1 + 1
    do il = ist, n1 / 2
       itmp2 = itmp1 - il
       iedg = vec(il, n2)
       vec(il, n2) = vec(itmp2, n2)
       vec(itmp2, n2) = iedg
    end do
  end subroutine alignment_edge_P_i8

  !> @brief Identity transformation, double precision array
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_edge_I_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec
  end subroutine alignment_edge_I_rp

  !> @brief Permutation transformation, double precision, full array
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_edge_P_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, itmp1, itmp2
    real(rp) :: redg

    itmp1 = n1 + 1
    do il = ist, n1 / 2
       itmp2 = itmp1 - il
       redg = vec(il, n2)
       vec(il, n2) = vec(itmp2, n2)
       vec(itmp2, n2) = redg
    end do
  end subroutine alignment_edge_P_rp

end module alignment_edge
