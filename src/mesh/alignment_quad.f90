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
!> Quad alignment operators
module alignment_quad
  use num_types, only : i4, i8, rp
  use utils, only : neko_error
  use math, only : arreq
  use alignment, only : alignment_t
  implicit none
  private

  public :: alignment_quad_init, alignment_quad_find, alignment_quad_I_t, &
       & alignment_quad_T_t, alignment_quad_PX_t, alignment_quad_PXT_t, &
       & alignment_quad_PYT_t, alignment_quad_PY_t, alignment_quad_PXPYT_t, &
       & alignment_quad_PXPY_t

  !> number of operations different from identity
  integer(i4), public, parameter :: NEKO_QUAD_NOPERATION = 7

  !> Quad identity (I) transformation type
  type, extends(alignment_t) :: alignment_quad_I_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_I
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_I_i4
     procedure, nopass :: trns_i8 => alignment_quad_I_i8
     procedure, nopass :: trns_rp => alignment_quad_I_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_I_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_I_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_I_rp
  end type alignment_quad_I_t

  !> Quad transposition (T) transformation type
  type, extends(alignment_t) :: alignment_quad_T_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_T
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_T_i4
     procedure, nopass :: trns_i8 => alignment_quad_T_i8
     procedure, nopass :: trns_rp => alignment_quad_T_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_T_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_T_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_T_rp
  end type alignment_quad_T_t

  !> Quad row permutation (PX) transformation type
  type, extends(alignment_t) :: alignment_quad_PX_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_PX
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_PX_i4
     procedure, nopass :: trns_i8 => alignment_quad_PX_i8
     procedure, nopass :: trns_rp => alignment_quad_PX_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_PX_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_PX_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_PX_rp
  end type alignment_quad_PX_t

  !> Quad row permutation and transposition (PXT) transformation type
  type, extends(alignment_t) :: alignment_quad_PXT_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_PXT
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_PXT_i4
     procedure, nopass :: trns_i8 => alignment_quad_PXT_i8
     procedure, nopass :: trns_rp => alignment_quad_PXT_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_PYT_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_PYT_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_PYT_rp
  end type alignment_quad_PXT_t

  !> Quad column permutation and transposition (PYT) transformation type
  type, extends(alignment_t) :: alignment_quad_PYT_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_PYT
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_PYT_i4
     procedure, nopass :: trns_i8 => alignment_quad_PYT_i8
     procedure, nopass :: trns_rp => alignment_quad_PYT_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_PXT_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_PXT_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_PXT_rp
  end type alignment_quad_PYT_t

  !> Quad column permutation (PY) transformation type
  type, extends(alignment_t) :: alignment_quad_PY_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_PY
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_PY_i4
     procedure, nopass :: trns_i8 => alignment_quad_PY_i8
     procedure, nopass :: trns_rp => alignment_quad_PY_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_PY_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_PY_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_PY_rp
  end type alignment_quad_PY_t

  !> Quad row, column permutation and transposition (PXPYT) transformation type
  type, extends(alignment_t) :: alignment_quad_PXPYT_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_PXPYT
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_PXPYT_i4
     procedure, nopass :: trns_i8 => alignment_quad_PXPYT_i8
     procedure, nopass :: trns_rp => alignment_quad_PXPYT_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_PXPYT_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_PXPYT_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_PXPYT_rp
  end type alignment_quad_PXPYT_t

  !> Quad row, column permutation (PXPY) transformation type
  type, extends(alignment_t) :: alignment_quad_PXPY_t
   contains
     !> Is transformation identity
     procedure, pass(this) :: ifid => ifidentity_quad_PXPY
     !> Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => alignment_quad_PXPY_i4
     procedure, nopass :: trns_i8 => alignment_quad_PXPY_i8
     procedure, nopass :: trns_rp => alignment_quad_PXPY_rp
     !> Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => alignment_quad_PXPY_i4
     procedure, nopass :: trns_inv_i8 => alignment_quad_PXPY_i8
     procedure, nopass :: trns_inv_rp => alignment_quad_PXPY_rp
  end type alignment_quad_PXPY_t

contains

  !> @brief Allocate a single alignment operator
  !! @parameter[in]      algn  relative edge alignment
  !! @parameter[inout]   trns  alignment operator
  subroutine alignment_quad_init(algn, trns)
    integer(i4), intent(in) :: algn
    class(alignment_t), allocatable, intent(inout) :: trns

    if (allocated(trns) ) then
       deallocate(trns)
    end if

    select case (algn)
    case(0) ! identity
       allocate(alignment_quad_I_t :: trns)
       call trns%set_algn(algn)
    case(1) ! transposition
       allocate(alignment_quad_T_t :: trns)
       call trns%set_algn(algn)
    case(2) ! row permutation
       allocate(alignment_quad_PX_t :: trns)
       call trns%set_algn(algn)
    case(3) ! row permutation and transposition
       allocate(alignment_quad_PXT_t :: trns)
       call trns%set_algn(algn)
    case(4) ! column permutation and transposition
       allocate(alignment_quad_PYT_t :: trns)
       call trns%set_algn(algn)
    case(5) ! column permutation
       allocate(alignment_quad_PY_t :: trns)
       call trns%set_algn(algn)
    case(6) ! row, column permutation and transposition
       allocate(alignment_quad_PXPYT_t :: trns)
       call trns%set_algn(algn)
    case(7) ! row, column permutation
       allocate(alignment_quad_PXPY_t :: trns)
       call trns%set_algn(algn)
    case default
       call neko_error('Wrong quad alignment')
    end select

  end subroutine alignment_quad_init

  !> @brief Find relative alignment for quads
  !! @parameter[out]     equal       array equality flag
  !! @parameter[out]     algn        relative quad alignment
  !! @parameter[inout]   qad1, qad2  quads
  !! @parameter[in]      n1, n2      dimensions
  subroutine alignment_quad_find(equal, algn, qad1, qad2, n1, n2)
    logical, intent(out) :: equal
    integer(i4), intent(out) :: algn
    integer(i4), intent(in) :: n1, n2
    integer(i4), dimension(n1, n2), intent(inout) :: qad1, qad2
    type(alignment_quad_T_t) :: op_T
    type(alignment_quad_PX_t) :: op_PX
    type(alignment_quad_PXT_t) :: op_PXT
    type(alignment_quad_PYT_t) :: op_PYT
    type(alignment_quad_PY_t) :: op_PY
    type(alignment_quad_PXPYT_t) :: op_PXPYT
    type(alignment_quad_PXPY_t) :: op_PXPY

    algn = -1
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 0
       return
    end if

    call op_T%trns_inv_i4(qad1, n1, n2, 1)
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 1
       return
    end if
    call op_T%trns_i4(qad1, n1, n2, 1)

    call op_PX%trns_inv_i4(qad1, n1, n2, 1)
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 2
       return
    end if
    call op_PX%trns_i4(qad1, n1, n2, 1)

    call op_PXT%trns_inv_i4(qad1, n1, n2, 1)
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 3
       return
    end if
    call op_PXT%trns_i4(qad1, n1, n2, 1)

    call op_PYT%trns_inv_i4(qad1, n1, n2, 1)
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 4
       return
    end if
    call op_PYT%trns_i4(qad1, n1, n2, 1)

    call op_PY%trns_inv_i4(qad1, n1, n2, 1)
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 5
       return
    end if
    call op_PY%trns_i4(qad1, n1, n2, 1)

    call op_PXPYT%trns_inv_i4(qad1, n1, n2, 1)
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 6
       return
    end if
    call op_PXPYT%trns_i4(qad1, n1, n2, 1)

    call op_PXPY%trns_inv_i4(qad1, n1, n2, 1)
    equal = arreq(qad1, qad2, n1, n2)
    if (equal) then
       algn = 7
       return
    end if

  end subroutine alignment_quad_find

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_I(this) result(ifid)
    class(alignment_quad_I_t), intent(in) :: this
    logical :: ifid
    ifid = .true.
  end function ifidentity_quad_I

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_T(this) result(ifid)
    class(alignment_quad_T_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_quad_T

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_PX(this) result(ifid)
    class(alignment_quad_PX_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_quad_PX

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_PXT(this) result(ifid)
    class(alignment_quad_PXT_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_quad_PXT

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_PYT(this) result(ifid)
    class(alignment_quad_PYT_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_quad_PYT

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_PY(this) result(ifid)
    class(alignment_quad_PY_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_quad_PY

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_PXPYT(this) result(ifid)
    class(alignment_quad_PXPYT_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_quad_PXPYT

  !> Function returning identity flag
  !! @return   ifid
  pure function ifidentity_quad_PXPY(this) result(ifid)
    class(alignment_quad_PXPY_t), intent(in) :: this
    logical :: ifid
    ifid = .false.
  end function ifidentity_quad_PXPY

  !> @brief Identity transformation, single integer array
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_I_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec
  end subroutine alignment_quad_I_i4

  !> @brief Transpose transformation, single integer
  !! @note This routine works for lx=ly(=lz) only.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_T_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = ist, n2 - (ist - 1)
       do il = ist, jl -1
          iface = vec(il, jl)
          vec(il, jl) = vec(jl, il)
          vec(jl, il) = iface
       end do
    end do
  end subroutine alignment_quad_T_i4

  !> @brief Column permutation transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PX_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = ist, n2 - (ist - 1)
       do il = ist, n1 / 2
          iface = vec(il, jl)
          vec(il, jl) = vec(n1 + 1 - il, jl)
          vec(n1 + 1 - il, jl) = iface
       end do
    end do
  end subroutine alignment_quad_PX_i4

  !> @brief Row permutation transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PY_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = ist, n2 / 2
       do il = ist, n1 - (ist - 1)
          iface = vec(il, jl)
          vec(il, jl) = vec(il, n2 + 1 - jl)
          vec(il, n2 + 1 - jl) = iface
       end do
    end do
  end subroutine alignment_quad_PY_i4

  !> @brief PXT = TPY transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXT_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_i4(vec, n1, n2, ist)
    call alignment_quad_T_i4(vec, n1, n2, ist)
  end subroutine alignment_quad_PXT_i4

  !> @brief PYT = TPX transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PYT_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PY_i4(vec, n1, n2, ist)
    call alignment_quad_T_i4(vec, n1, n2, ist)
  end subroutine alignment_quad_PYT_i4

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXPYT_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_i4(vec, n1, n2, ist)
    call alignment_quad_PY_i4(vec, n1, n2, ist)
    call alignment_quad_T_i4(vec, n1, n2, ist)
  end subroutine alignment_quad_PXPYT_i4

  !> @brief PXPY = PYPX transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXPY_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_i4(vec, n1, n2, ist)
    call alignment_quad_PY_i4(vec, n1, n2, ist)
  end subroutine alignment_quad_PXPY_i4

  !> @brief Identity transformation, double integer array
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_I_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec
  end subroutine alignment_quad_I_i8

  !> @brief Transpose transformation, double integer
  !! @note This routine works for lx=ly(=lz) only.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_T_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = ist, n2 - (ist - 1)
       do il = ist, jl -1
          iface = vec(il, jl)
          vec(il, jl) = vec(jl, il)
          vec(jl, il) = iface
       end do
    end do
  end subroutine alignment_quad_T_i8

  !> @brief Column permutation transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PX_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = ist, n2 - (ist - 1)
       do il = ist, n1 / 2
          iface = vec(il, jl)
          vec(il, jl) = vec(n1 + 1 - il, jl)
          vec(n1 + 1 - il, jl) = iface
       end do
    end do
  end subroutine alignment_quad_PX_i8

  !> @brief Row permutation transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PY_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = ist, n2 / 2
       do il = ist, n1 - (ist - 1)
          iface = vec(il, jl)
          vec(il, jl) = vec(il, n2 + 1 - jl)
          vec(il, n2 + 1 - jl) = iface
       end do
    end do
  end subroutine alignment_quad_PY_i8

  !> @brief PXT = TPY transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXT_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_i8(vec, n1, n2, ist)
    call alignment_quad_T_i8(vec, n1, n2, ist)
  end subroutine alignment_quad_PXT_i8

  !> @brief PYT = TPX transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PYT_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PY_i8(vec, n1, n2, ist)
    call alignment_quad_T_i8(vec, n1, n2, ist)
  end subroutine alignment_quad_PYT_i8

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXPYT_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_i8(vec, n1, n2, ist)
    call alignment_quad_PY_i8(vec, n1, n2, ist)
    call alignment_quad_T_i8(vec, n1, n2, ist)
  end subroutine alignment_quad_PXPYT_i8

  !> @brief PXPY = PYPX transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXPY_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_i8(vec, n1, n2, ist)
    call alignment_quad_PY_i8(vec, n1, n2, ist)
  end subroutine alignment_quad_PXPY_i8

  !> @brief Identity transformation, double precision array
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_I_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec
  end subroutine alignment_quad_I_rp

  !> @brief Transpose transformation, double precision
  !! @note This routine works for lx=ly(=lz) only.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_T_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    real(rp) :: rface

    do jl = ist, n2 - (ist - 1)
       do il = ist, jl -1
          rface = vec(il, jl)
          vec(il, jl) = vec(jl, il)
          vec(jl, il) = rface
       end do
    end do
  end subroutine alignment_quad_T_rp

  !> @brief Column permutation transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PX_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    real(rp) :: rface

    do jl = ist, n2 - (ist - 1)
       do il = ist, n1 / 2
          rface = vec(il, jl)
          vec(il, jl) = vec(n1 + 1 - il, jl)
          vec(n1 + 1 - il, jl) = rface
       end do
    end do
  end subroutine alignment_quad_PX_rp

  !> @brief Row permutation transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PY_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    real(rp) :: rface

    do jl = ist, n2 / 2
       do il = ist, n1 - (ist - 1)
          rface = vec(il, jl)
          vec(il, jl) = vec(il, n2 + 1 - jl)
          vec(il, n2 + 1 - jl) = rface
       end do
    end do
  end subroutine alignment_quad_PY_rp

  !> @brief PXT = TPY transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXT_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_rp(vec, n1, n2, ist)
    call alignment_quad_T_rp(vec, n1, n2, ist)
  end subroutine alignment_quad_PXT_rp

  !> @brief PYT = TPX transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PYT_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PY_rp(vec, n1, n2, ist)
    call alignment_quad_T_rp(vec, n1, n2, ist)
  end subroutine alignment_quad_PYT_rp

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXPYT_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_rp(vec, n1, n2, ist)
    call alignment_quad_PY_rp(vec, n1, n2, ist)
    call alignment_quad_T_rp(vec, n1, n2, ist)
  end subroutine alignment_quad_PXPYT_rp

  !> @brief PXPY = PYPX transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine alignment_quad_PXPY_rp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(rp), dimension(n1, n2), intent(inout) :: vec

    call alignment_quad_PX_rp(vec, n1, n2, ist)
    call alignment_quad_PY_rp(vec, n1, n2, ist)
  end subroutine alignment_quad_PXPY_rp

end module alignment_quad
