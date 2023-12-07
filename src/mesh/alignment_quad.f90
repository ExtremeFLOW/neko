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
!> Quad alignment operators
module alignment_quad
  use num_types, only : i4, i8, dp
  use utils, only : neko_error
  use alignment, only : alignment_t, alignment_set_t
  implicit none
  private

  public :: alignment_quad_init, alignment_quad_I_t, alignment_quad_T_t, &
       & alignment_quad_PX_t, alignment_quad_PXT_t, alignment_quad_PYT_t, &
       & alignment_quad_PY_t, alignment_quad_PXPYT_t, alignment_quad_PXPY_t, &
       & alignment_quad_set_t

  !> number of operations different from identity
  integer(i4), public, parameter :: NEKO_QUAD_NOPERATION = 7

  !> Quad identity (I) transformation type
  type, extends(alignment_t) :: alignment_quad_I_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_I_i4
     procedure, nopass :: trns_i8 => transform_quad_I_i8
     procedure, nopass :: trns_dp => transform_quad_I_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_I_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_I_i8
     procedure, nopass :: trns_inv_dp => transform_quad_I_dp
  end type alignment_quad_I_t

  !> Quad transposition (T) transformation type
  type, extends(alignment_t) :: alignment_quad_T_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_T_i4
     procedure, nopass :: trns_i8 => transform_quad_T_i8
     procedure, nopass :: trns_dp => transform_quad_T_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_T_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_T_i8
     procedure, nopass :: trns_inv_dp => transform_quad_T_dp
  end type alignment_quad_T_t

  !> Quad row permutation (PX) transformation type
  type, extends(alignment_t) :: alignment_quad_PX_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_PX_i4
     procedure, nopass :: trns_i8 => transform_quad_PX_i8
     procedure, nopass :: trns_dp => transform_quad_PX_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_PX_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_PX_i8
     procedure, nopass :: trns_inv_dp => transform_quad_PX_dp
  end type alignment_quad_PX_t

  !> Quad row permutation and transposition (PXT) transformation type
  type, extends(alignment_t) :: alignment_quad_PXT_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_PXT_i4
     procedure, nopass :: trns_i8 => transform_quad_PXT_i8
     procedure, nopass :: trns_dp => transform_quad_PXT_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_PYT_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_PYT_i8
     procedure, nopass :: trns_inv_dp => transform_quad_PYT_dp
  end type alignment_quad_PXT_t

  !> Quad column permutation and transposition (PYT) transformation type
  type, extends(alignment_t) :: alignment_quad_PYT_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_PYT_i4
     procedure, nopass :: trns_i8 => transform_quad_PYT_i8
     procedure, nopass :: trns_dp => transform_quad_PYT_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_PXT_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_PXT_i8
     procedure, nopass :: trns_inv_dp => transform_quad_PXT_dp
  end type alignment_quad_PYT_t

  !> Quad column permutation (PY) transformation type
  type, extends(alignment_t) :: alignment_quad_PY_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_PY_i4
     procedure, nopass :: trns_i8 => transform_quad_PY_i8
     procedure, nopass :: trns_dp => transform_quad_PY_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_PY_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_PY_i8
     procedure, nopass :: trns_inv_dp => transform_quad_PY_dp
  end type alignment_quad_PY_t

  !> Quad row, column permutation and transposition (PXPYT) transformation type
  type, extends(alignment_t) :: alignment_quad_PXPYT_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_PXPYT_i4
     procedure, nopass :: trns_i8 => transform_quad_PXPYT_i8
     procedure, nopass :: trns_dp => transform_quad_PXPYT_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_PXPYT_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_PXPYT_i8
     procedure, nopass :: trns_inv_dp => transform_quad_PXPYT_dp
  end type alignment_quad_PXPYT_t

  !> Quad row, column permutation (PXPY) transformation type
  type, extends(alignment_t) :: alignment_quad_PXPY_t
   contains
     ! Direct transformation of full array, different types
     procedure, nopass :: trns_i4 => transform_quad_PXPY_i4
     procedure, nopass :: trns_i8 => transform_quad_PXPY_i8
     procedure, nopass :: trns_dp => transform_quad_PXPY_dp
     ! Inverse transformation of full array, different types
     procedure, nopass :: trns_inv_i4 => transform_quad_PXPY_i4
     procedure, nopass :: trns_inv_i8 => transform_quad_PXPY_i8
     procedure, nopass :: trns_inv_dp => transform_quad_PXPY_dp
  end type alignment_quad_PXPY_t

  !> Type containing set of quad alignment operators
  !! @details There are four main operations : identity (I), column
  !! permutation (PX), row permutation (PY) and transposition (T).
  !! They are combined into 8 allowed quad transformations: I, T, PX, PXT,
  !! PYT, PY, PXPYT, PXPY
  !! @note The identity operator is not really needed, but i keep it for
  !! completeness.
  type, extends(alignment_set_t) :: alignment_quad_set_t
   contains
     !> Initialise type
     procedure, pass(this) :: init => quad_set_init
  end type alignment_quad_set_t

contains

  !> @brief Allocate a single alignment operator
  !! @parameter[in]      algn  relative edge alignment
  !! @parameter[inout]   trns  alignment operator
  subroutine alignment_quad_init(algn, trns)
    integer(i4), intent(in) :: algn
    class(alignment_t), allocatable, intent(inout) :: trns

    if (allocated(trns)) then
       deallocate(trns)
    end if

    select case(algn)
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

  !> @brief Initialise alignment operator set
  subroutine quad_set_init(this)
    class(alignment_quad_set_t), intent(inout) :: this
    integer(i4) :: il

    call this%free()

    call this%set_nop(NEKO_QUAD_NOPERATION)

    allocate(this%trns(0 : NEKO_QUAD_NOPERATION))

    do il = 0, NEKO_QUAD_NOPERATION
       call alignment_quad_init(il, this%trns(il)%op)
    end do

  end subroutine quad_set_init

  !> @brief Identity transformation, single integer array
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_I_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec
  end subroutine transform_quad_I_i4

  !> @brief Transpose transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_T_i4(vec, n1, n2, ist)
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
  end subroutine transform_quad_T_i4

  !> @brief Column permutation transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PX_i4(vec, n1, n2, ist)
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
  end subroutine transform_quad_PX_i4

  !> @brief Row permutation transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PY_i4(vec, n1, n2, ist)
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
  end subroutine transform_quad_PY_i4

  !> @brief PXT = TPY transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXT_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_i4(vec, n1, n2, ist)
    call transform_quad_T_i4(vec, n1, n2, ist)
  end subroutine transform_quad_PXT_i4

  !> @brief PYT = TPX transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PYT_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PY_i4(vec, n1, n2, ist)
    call transform_quad_T_i4(vec, n1, n2, ist)
  end subroutine transform_quad_PYT_i4

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXPYT_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_i4(vec, n1, n2, ist)
    call transform_quad_PY_i4(vec, n1, n2, ist)
    call transform_quad_T_i4(vec, n1, n2, ist)
  end subroutine transform_quad_PXPYT_i4

  !> @brief PXPY = PYPX transformation, single integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXPY_i4(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i4), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_i4(vec, n1, n2, ist)
    call transform_quad_PY_i4(vec, n1, n2, ist)
  end subroutine transform_quad_PXPY_i4

  !> @brief Identity transformation, double integer array
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_I_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec
  end subroutine transform_quad_I_i8

  !> @brief Transpose transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_T_i8(vec, n1, n2, ist)
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
  end subroutine transform_quad_T_i8

  !> @brief Column permutation transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PX_i8(vec, n1, n2, ist)
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
  end subroutine transform_quad_PX_i8

  !> @brief Row permutation transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PY_i8(vec, n1, n2, ist)
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
  end subroutine transform_quad_PY_i8

  !> @brief PXT = TPY transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXT_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_i8(vec, n1, n2, ist)
    call transform_quad_T_i8(vec, n1, n2, ist)
  end subroutine transform_quad_PXT_i8

  !> @brief PYT = TPX transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PYT_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PY_i8(vec, n1, n2, ist)
    call transform_quad_T_i8(vec, n1, n2, ist)
  end subroutine transform_quad_PYT_i8

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXPYT_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_i8(vec, n1, n2, ist)
    call transform_quad_PY_i8(vec, n1, n2, ist)
    call transform_quad_T_i8(vec, n1, n2, ist)
  end subroutine transform_quad_PXPYT_i8

  !> @brief PXPY = PYPX transformation, double integer
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXPY_i8(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    integer(i8), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_i8(vec, n1, n2, ist)
    call transform_quad_PY_i8(vec, n1, n2, ist)
  end subroutine transform_quad_PXPY_i8

  !> @brief Identity transformation, double precision array
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_I_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec
  end subroutine transform_quad_I_dp

  !> @brief Transpose transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_T_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = ist, n2 - (ist - 1)
       do il = ist, jl -1
          rface = vec(il, jl)
          vec(il, jl) = vec(jl, il)
          vec(jl, il) = rface
       end do
    end do
  end subroutine transform_quad_T_dp

  !> @brief Column permutation transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PX_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = ist, n2 - (ist - 1)
       do il = ist, n1 / 2
          rface = vec(il, jl)
          vec(il, jl) = vec(n1 + 1 - il, jl)
          vec(n1 + 1 - il, jl) = rface
       end do
    end do
  end subroutine transform_quad_PX_dp

  !> @brief Row permutation transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PY_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = ist, n2 / 2
       do il = ist, n1 - (ist - 1)
          rface = vec(il, jl)
          vec(il, jl) = vec(il, n2 + 1 - jl)
          vec(il, n2 + 1 - jl) = rface
       end do
    end do
  end subroutine transform_quad_PY_dp

  !> @brief PXT = TPY transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXT_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_dp(vec, n1, n2, ist)
    call transform_quad_T_dp(vec, n1, n2, ist)
  end subroutine transform_quad_PXT_dp

  !> @brief PYT = TPX transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PYT_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PY_dp(vec, n1, n2, ist)
    call transform_quad_T_dp(vec, n1, n2, ist)
  end subroutine transform_quad_PYT_dp

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXPYT_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_dp(vec, n1, n2, ist)
    call transform_quad_PY_dp(vec, n1, n2, ist)
    call transform_quad_T_dp(vec, n1, n2, ist)
  end subroutine transform_quad_PXPYT_dp

  !> @brief PXPY = PYPX transformation, double precision
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  !! @parameter[in]      ist      starting position
  pure subroutine transform_quad_PXPY_dp(vec, n1, n2, ist)
    integer(i4), intent(in) :: n1, n2, ist
    real(dp), dimension(n1, n2), intent(inout) :: vec

    call transform_quad_PX_dp(vec, n1, n2, ist)
    call transform_quad_PY_dp(vec, n1, n2, ist)
  end subroutine transform_quad_PXPY_dp

end module alignment_quad
