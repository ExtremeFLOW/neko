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
  implicit none
  private

  public :: alignment_quad_t, alignment_quad_op_set_t

  !> number of operations different from identity
  integer(i4), public, parameter :: NEKO_QUAD_NOPERATION = 7

  !> procedure pointer type; i4
  type :: algn_quad_proc_i4_ptr
     procedure(transform_i4), pointer, nopass :: obj
  end type algn_quad_proc_i4_ptr
  !> procedure pointer type; i8
  type :: algn_quad_proc_i8_ptr
     procedure(transform_i8), pointer, nopass :: obj
  end type algn_quad_proc_i8_ptr
  !> procedure pointer type; dp
  type :: algn_quad_proc_dp_ptr
     procedure(transform_dp), pointer, nopass :: obj
  end type algn_quad_proc_dp_ptr

  !> Type containing set of quad alignment operators
  !! @details There are four main operations : identity (I), column
  !! permutation (PX), row permutation (PY) and transposition (T).
  !! They are combined into 8 allowed quad transformations: I, T, PX, PXT,
  !! PYT, PY, PXPYT, PXPY
  !! @note The identity operator is not really needed, but i keep it for
  !! completeness.
  type :: alignment_quad_t
     !> number of different operations excluding identity
     integer(i4), private :: noperation_ = NEKO_QUAD_NOPERATION
     !> Direct array transformations for whole array
     type(algn_quad_proc_i4_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_f_i4
     type(algn_quad_proc_i8_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_f_i8
     type(algn_quad_proc_dp_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_f_dp
     !> Direct array transformations for array interior
     type(algn_quad_proc_i4_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_i_i4
     type(algn_quad_proc_i8_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_i_i8
     type(algn_quad_proc_dp_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_i_dp
     !> Inverse array transformations for whole array
     type(algn_quad_proc_i4_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_inv_f_i4
     type(algn_quad_proc_i8_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_inv_f_i8
     type(algn_quad_proc_dp_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_inv_f_dp
     !> Inverse array transformations for array interior
     type(algn_quad_proc_i4_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_inv_i_i4
     type(algn_quad_proc_i8_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_inv_i_i8
     type(algn_quad_proc_dp_ptr),&
          & dimension(0 : NEKO_QUAD_NOPERATION) :: trns_inv_i_dp
   contains
     !> Initialise procedure pointers
     procedure, pass(this) :: init => quad_init
     !> Free procedure pointers
     procedure, pass(this) :: free => quad_free
     !> Return number of operations
     procedure, pass(this) :: nop => quad_noperation_get
  end type alignment_quad_t

  !> Type combining a single set of alignment operators
  type :: alignment_quad_op_set_t
     !> Relative quad alignment
     integer(i4) :: alignment = -1
     !> Direct array transformations for whole array
     type(algn_quad_proc_i4_ptr) :: trns_f_i4
     type(algn_quad_proc_i8_ptr) :: trns_f_i8
     type(algn_quad_proc_dp_ptr) :: trns_f_dp
     !> Direct array transformations for array interior
     type(algn_quad_proc_i4_ptr) :: trns_i_i4
     type(algn_quad_proc_i8_ptr) :: trns_i_i8
     type(algn_quad_proc_dp_ptr) :: trns_i_dp
     !> Inverse array transformations for whole array
     type(algn_quad_proc_i4_ptr) :: trns_inv_f_i4
     type(algn_quad_proc_i8_ptr) :: trns_inv_f_i8
     type(algn_quad_proc_dp_ptr) :: trns_inv_f_dp
     !> Inverse array transformations for array interior
     type(algn_quad_proc_i4_ptr) :: trns_inv_i_i4
     type(algn_quad_proc_i8_ptr) :: trns_inv_i_i8
     type(algn_quad_proc_dp_ptr) :: trns_inv_i_dp
   contains
      !> Initialise alignment and procedure pointers
     procedure, pass(this) :: init => quad_op_set_init
     !> Free alignment and pointers
     procedure, pass(this) :: free => quad_op_set_free
  end type alignment_quad_op_set_t

  ! Abstract types for different transformations; various types
  abstract interface
     pure subroutine transform_i4(sz, fcs, work)
       import i4
       integer(i4), intent(in) :: sz
       integer(i4), dimension(sz, sz), intent(inout) :: fcs
       integer(i4), dimension(sz), intent(inout) :: work
     end subroutine transform_i4

     pure subroutine transform_i8(sz, fcs, work)
       import i4
       import i8
       integer(i4), intent(in) :: sz
       integer(i8), dimension(sz, sz), intent(inout) :: fcs
       integer(i8), dimension(sz), intent(inout) :: work
     end subroutine transform_i8

     pure subroutine transform_dp(sz, fcs, work)
       import i4
       import dp
       integer(i4), intent(in) :: sz
       real(dp), dimension(sz, sz), intent(inout) :: fcs
       real(dp), dimension(sz), intent(inout) :: work
     end subroutine transform_dp
  end interface

contains
  !> @brief Initialise procedure pointers
  subroutine quad_init(this)
    class(alignment_quad_t), intent(inout) :: this

    call this%free()

    ! Identity transformation is added for completeness; in general not needed
    ! Direct transformation of full array, different types
    this%trns_f_i4(0)%obj => transform_quad_I_i4
    this%trns_f_i4(1)%obj => transform_quad_T_full_i4
    this%trns_f_i4(2)%obj => transform_quad_PX_full_i4
    this%trns_f_i4(3)%obj => transform_quad_PXT_full_i4
    this%trns_f_i4(4)%obj => transform_quad_PYT_full_i4
    this%trns_f_i4(5)%obj => transform_quad_PY_full_i4
    this%trns_f_i4(6)%obj => transform_quad_PXPYT_full_i4
    this%trns_f_i4(7)%obj => transform_quad_PXPY_full_i4
    this%trns_f_i8(0)%obj => transform_quad_I_i8
    this%trns_f_i8(1)%obj => transform_quad_T_full_i8
    this%trns_f_i8(2)%obj => transform_quad_PX_full_i8
    this%trns_f_i8(3)%obj => transform_quad_PXT_full_i8
    this%trns_f_i8(4)%obj => transform_quad_PYT_full_i8
    this%trns_f_i8(5)%obj => transform_quad_PY_full_i8
    this%trns_f_i8(6)%obj => transform_quad_PXPYT_full_i8
    this%trns_f_i8(7)%obj => transform_quad_PXPY_full_i8
    this%trns_f_dp(0)%obj => transform_quad_I_dp
    this%trns_f_dp(1)%obj => transform_quad_T_full_dp
    this%trns_f_dp(2)%obj => transform_quad_PX_full_dp
    this%trns_f_dp(3)%obj => transform_quad_PXT_full_dp
    this%trns_f_dp(4)%obj => transform_quad_PYT_full_dp
    this%trns_f_dp(5)%obj => transform_quad_PY_full_dp
    this%trns_f_dp(6)%obj => transform_quad_PXPYT_full_dp
    this%trns_f_dp(7)%obj => transform_quad_PXPY_full_dp
    ! Direct transformation of array interior, different types
    this%trns_i_i4(0)%obj => transform_quad_I_i4
    this%trns_i_i4(1)%obj => transform_quad_T_int_i4
    this%trns_i_i4(2)%obj => transform_quad_PX_int_i4
    this%trns_i_i4(3)%obj => transform_quad_PXT_int_i4
    this%trns_i_i4(4)%obj => transform_quad_PYT_int_i4
    this%trns_i_i4(5)%obj => transform_quad_PY_int_i4
    this%trns_i_i4(6)%obj => transform_quad_PXPYT_int_i4
    this%trns_i_i4(7)%obj => transform_quad_PXPY_int_i4
    this%trns_i_i8(0)%obj => transform_quad_I_i8
    this%trns_i_i8(1)%obj => transform_quad_T_int_i8
    this%trns_i_i8(2)%obj => transform_quad_PX_int_i8
    this%trns_i_i8(3)%obj => transform_quad_PXT_int_i8
    this%trns_i_i8(4)%obj => transform_quad_PYT_int_i8
    this%trns_i_i8(5)%obj => transform_quad_PY_int_i8
    this%trns_i_i8(6)%obj => transform_quad_PXPYT_int_i8
    this%trns_i_i8(7)%obj => transform_quad_PXPY_int_i8
    this%trns_i_dp(0)%obj => transform_quad_I_dp
    this%trns_i_dp(1)%obj => transform_quad_T_int_dp
    this%trns_i_dp(2)%obj => transform_quad_PX_int_dp
    this%trns_i_dp(3)%obj => transform_quad_PXT_int_dp
    this%trns_i_dp(4)%obj => transform_quad_PYT_int_dp
    this%trns_i_dp(5)%obj => transform_quad_PY_int_dp
    this%trns_i_dp(6)%obj => transform_quad_PXPYT_int_dp
    this%trns_i_dp(7)%obj => transform_quad_PXPY_int_dp
    ! Inverse transformation of full array, different types
    this%trns_inv_f_i4(0)%obj => transform_quad_I_i4
    this%trns_inv_f_i4(1)%obj => transform_quad_T_full_i4
    this%trns_inv_f_i4(2)%obj => transform_quad_PX_full_i4
    this%trns_inv_f_i4(3)%obj => transform_quad_PYT_full_i4
    this%trns_inv_f_i4(4)%obj => transform_quad_PXT_full_i4
    this%trns_inv_f_i4(5)%obj => transform_quad_PY_full_i4
    this%trns_inv_f_i4(6)%obj => transform_quad_PXPYT_full_i4
    this%trns_inv_f_i4(7)%obj => transform_quad_PXPY_full_i4
    this%trns_inv_f_i8(0)%obj => transform_quad_I_i8
    this%trns_inv_f_i8(1)%obj => transform_quad_T_full_i8
    this%trns_inv_f_i8(2)%obj => transform_quad_PX_full_i8
    this%trns_inv_f_i8(3)%obj => transform_quad_PYT_full_i8
    this%trns_inv_f_i8(4)%obj => transform_quad_PXT_full_i8
    this%trns_inv_f_i8(5)%obj => transform_quad_PY_full_i8
    this%trns_inv_f_i8(6)%obj => transform_quad_PXPYT_full_i8
    this%trns_inv_f_i8(7)%obj => transform_quad_PXPY_full_i8
    this%trns_inv_f_dp(0)%obj => transform_quad_I_dp
    this%trns_inv_f_dp(1)%obj => transform_quad_T_full_dp
    this%trns_inv_f_dp(2)%obj => transform_quad_PX_full_dp
    this%trns_inv_f_dp(3)%obj => transform_quad_PYT_full_dp
    this%trns_inv_f_dp(4)%obj => transform_quad_PXT_full_dp
    this%trns_inv_f_dp(5)%obj => transform_quad_PY_full_dp
    this%trns_inv_f_dp(6)%obj => transform_quad_PXPYT_full_dp
    this%trns_inv_f_dp(7)%obj => transform_quad_PXPY_full_dp
    ! Inverse transformation of array interior, different types
    this%trns_inv_i_i4(0)%obj => transform_quad_I_i4
    this%trns_inv_i_i4(1)%obj => transform_quad_T_int_i4
    this%trns_inv_i_i4(2)%obj => transform_quad_PX_int_i4
    this%trns_inv_i_i4(3)%obj => transform_quad_PYT_int_i4
    this%trns_inv_i_i4(4)%obj => transform_quad_PXT_int_i4
    this%trns_inv_i_i4(5)%obj => transform_quad_PY_int_i4
    this%trns_inv_i_i4(6)%obj => transform_quad_PXPYT_int_i4
    this%trns_inv_i_i4(7)%obj => transform_quad_PXPY_int_i4
    this%trns_inv_i_i8(0)%obj => transform_quad_I_i8
    this%trns_inv_i_i8(1)%obj => transform_quad_T_int_i8
    this%trns_inv_i_i8(2)%obj => transform_quad_PX_int_i8
    this%trns_inv_i_i8(3)%obj => transform_quad_PYT_int_i8
    this%trns_inv_i_i8(4)%obj => transform_quad_PXT_int_i8
    this%trns_inv_i_i8(5)%obj => transform_quad_PY_int_i8
    this%trns_inv_i_i8(6)%obj => transform_quad_PXPYT_int_i8
    this%trns_inv_i_i8(7)%obj => transform_quad_PXPY_int_i8
    this%trns_inv_i_dp(0)%obj => transform_quad_I_dp
    this%trns_inv_i_dp(1)%obj => transform_quad_T_int_dp
    this%trns_inv_i_dp(2)%obj => transform_quad_PX_int_dp
    this%trns_inv_i_dp(3)%obj => transform_quad_PYT_int_dp
    this%trns_inv_i_dp(4)%obj => transform_quad_PXT_int_dp
    this%trns_inv_i_dp(5)%obj => transform_quad_PY_int_dp
    this%trns_inv_i_dp(6)%obj => transform_quad_PXPYT_int_dp
    this%trns_inv_i_dp(7)%obj => transform_quad_PXPY_int_dp

    return
  end subroutine quad_init

  !> @brief Free procedure pointers
  subroutine quad_free(this)
    class(alignment_quad_t), intent(inout) :: this
    integer(i4) :: il

    ! free pointers
    do il = 0, NEKO_QUAD_NOPERATION
       this%trns_f_i4(il)%obj => null()
       this%trns_f_i8(il)%obj => null()
       this%trns_f_dp(il)%obj => null()
       this%trns_i_i4(il)%obj => null()
       this%trns_i_i8(il)%obj => null()
       this%trns_i_dp(il)%obj => null()
       this%trns_inv_f_i4(il)%obj => null()
       this%trns_inv_f_i8(il)%obj => null()
       this%trns_inv_f_dp(il)%obj => null()
       this%trns_inv_i_i4(il)%obj => null()
       this%trns_inv_i_i8(il)%obj => null()
       this%trns_inv_i_dp(il)%obj => null()
    end do

    return
  end subroutine quad_free

  !> @brief Get number of operations
  !! @return   noperation
  pure function quad_noperation_get(this) result(noperation)
    class(alignment_quad_t), intent(in) :: this
    integer(i4) :: noperation
    noperation = this%noperation_
  end function quad_noperation_get

  !> @brief Initialise alignment and procedure pointers
  !! @parameter[in]   algn   alignment
  subroutine quad_op_set_init(this, algn)
    class(alignment_quad_op_set_t), intent(inout) :: this
    integer(i4), intent(in) :: algn

    call this%free()

    ! set relative alignment transformation
    if ((algn >= 0).and.(algn <= NEKO_QUAD_NOPERATION)) then
       this%alignment = algn
    else
       call neko_error('Not proper alignment.')
    end if
    select case(algn)
    case(0)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_I_i4
       this%trns_f_i8%obj => transform_quad_I_i8
       this%trns_f_dp%obj => transform_quad_I_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_I_i4
       this%trns_i_i8%obj => transform_quad_I_i8
       this%trns_i_dp%obj => transform_quad_I_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_I_i4
       this%trns_inv_f_i8%obj => transform_quad_I_i8
       this%trns_inv_f_dp%obj => transform_quad_I_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_I_i4
       this%trns_inv_i_i8%obj => transform_quad_I_i8
       this%trns_inv_i_dp%obj => transform_quad_I_dp
    case(1)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_T_full_i4
       this%trns_f_i8%obj => transform_quad_T_full_i8
       this%trns_f_dp%obj => transform_quad_T_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_T_int_i4
       this%trns_i_i8%obj => transform_quad_T_int_i8
       this%trns_i_dp%obj => transform_quad_T_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_T_full_i4
       this%trns_inv_f_i8%obj => transform_quad_T_full_i8
       this%trns_inv_f_dp%obj => transform_quad_T_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_T_int_i4
       this%trns_inv_i_i8%obj => transform_quad_T_int_i8
       this%trns_inv_i_dp%obj => transform_quad_T_int_dp
    case(2)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_PX_full_i4
       this%trns_f_i8%obj => transform_quad_PX_full_i8
       this%trns_f_dp%obj => transform_quad_PX_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_PX_int_i4
       this%trns_i_i8%obj => transform_quad_PX_int_i8
       this%trns_i_dp%obj => transform_quad_PX_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_PX_full_i4
       this%trns_inv_f_i8%obj => transform_quad_PX_full_i8
       this%trns_inv_f_dp%obj => transform_quad_PX_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_PX_int_i4
       this%trns_inv_i_i8%obj => transform_quad_PX_int_i8
       this%trns_inv_i_dp%obj => transform_quad_PX_int_dp
    case(3)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_PXT_full_i4
       this%trns_f_i8%obj => transform_quad_PXT_full_i8
       this%trns_f_dp%obj => transform_quad_PXT_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_PXT_int_i4
       this%trns_i_i8%obj => transform_quad_PXT_int_i8
       this%trns_i_dp%obj => transform_quad_PXT_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_PYT_full_i4
       this%trns_inv_f_i8%obj => transform_quad_PYT_full_i8
       this%trns_inv_f_dp%obj => transform_quad_PYT_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_PYT_int_i4
       this%trns_inv_i_i8%obj => transform_quad_PYT_int_i8
       this%trns_inv_i_dp%obj => transform_quad_PYT_int_dp
    case(4)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_PYT_full_i4
       this%trns_f_i8%obj => transform_quad_PYT_full_i8
       this%trns_f_dp%obj => transform_quad_PYT_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_PYT_int_i4
       this%trns_i_i8%obj => transform_quad_PYT_int_i8
       this%trns_i_dp%obj => transform_quad_PYT_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_PXT_full_i4
       this%trns_inv_f_i8%obj => transform_quad_PXT_full_i8
       this%trns_inv_f_dp%obj => transform_quad_PXT_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_PXT_int_i4
       this%trns_inv_i_i8%obj => transform_quad_PXT_int_i8
       this%trns_inv_i_dp%obj => transform_quad_PXT_int_dp
    case(5)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_PY_full_i4
       this%trns_f_i8%obj => transform_quad_PY_full_i8
       this%trns_f_dp%obj => transform_quad_PY_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_PY_int_i4
       this%trns_i_i8%obj => transform_quad_PY_int_i8
       this%trns_i_dp%obj => transform_quad_PY_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_PY_full_i4
       this%trns_inv_f_i8%obj => transform_quad_PY_full_i8
       this%trns_inv_f_dp%obj => transform_quad_PY_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_PY_int_i4
       this%trns_inv_i_i8%obj => transform_quad_PY_int_i8
       this%trns_inv_i_dp%obj => transform_quad_PY_int_dp
    case(6)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_PXPYT_full_i4
       this%trns_f_i8%obj => transform_quad_PXPYT_full_i8
       this%trns_f_dp%obj => transform_quad_PXPYT_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_PXPYT_int_i4
       this%trns_i_i8%obj => transform_quad_PXPYT_int_i8
       this%trns_i_dp%obj => transform_quad_PXPYT_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_PXPYT_full_i4
       this%trns_inv_f_i8%obj => transform_quad_PXPYT_full_i8
       this%trns_inv_f_dp%obj => transform_quad_PXPYT_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_PXPYT_int_i4
       this%trns_inv_i_i8%obj => transform_quad_PXPYT_int_i8
       this%trns_inv_i_dp%obj => transform_quad_PXPYT_int_dp
    case(7)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_quad_PXPY_full_i4
       this%trns_f_i8%obj => transform_quad_PXPY_full_i8
       this%trns_f_dp%obj => transform_quad_PXPY_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_quad_PXPY_int_i4
       this%trns_i_i8%obj => transform_quad_PXPY_int_i8
       this%trns_i_dp%obj => transform_quad_PXPY_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_quad_PXPY_full_i4
       this%trns_inv_f_i8%obj => transform_quad_PXPY_full_i8
       this%trns_inv_f_dp%obj => transform_quad_PXPY_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_quad_PXPY_int_i4
       this%trns_inv_i_i8%obj => transform_quad_PXPY_int_i8
       this%trns_inv_i_dp%obj => transform_quad_PXPY_int_dp
    end select

    return
  end subroutine quad_op_set_init

  !> @brief Free alignment and procedure pointers
  subroutine quad_op_set_free(this)
    class(alignment_quad_op_set_t), intent(inout) :: this

    this%alignment = -1
    ! free pointers
    this%trns_f_i4%obj => null()
    this%trns_f_i8%obj => null()
    this%trns_f_dp%obj => null()
    this%trns_i_i4%obj => null()
    this%trns_i_i8%obj => null()
    this%trns_i_dp%obj => null()
    this%trns_inv_f_i4%obj => null()
    this%trns_inv_f_i8%obj => null()
    this%trns_inv_f_dp%obj => null()
    this%trns_inv_i_i4%obj => null()
    this%trns_inv_i_i8%obj => null()
    this%trns_inv_i_dp%obj => null()

    return
  end subroutine quad_op_set_free

  !> @brief Identity transformation, single integer array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_I_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    return
  end subroutine transform_quad_I_i4

  !> @brief Transpose transformation, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_T_full_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = 1, sz
       do il = 1, jl -1
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(jl, il)
          fcs(jl, il) = iface
       end do
    end do

    return
  end subroutine transform_quad_T_full_i4

  !> @brief Transpose transformation, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_T_int_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = 2, sz - 1
       do il = 2, jl -1
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(jl, il)
          fcs(jl, il) = iface
       end do
    end do

    return
  end subroutine transform_quad_T_int_i4

  !> @brief Column permutation transformation, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PX_full_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = 1, sz
       do il = 1, sz/2
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(sz + 1 - il, jl)
          fcs(sz + 1 - il, jl) = iface
       end do
    end do

    return
  end subroutine transform_quad_PX_full_i4

  !> @brief Column permutation transformation, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PX_int_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = 2, sz - 1
       do il = 2, sz/2
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(sz + 1 - il, jl)
          fcs(sz + 1 - il, jl) = iface
       end do
    end do

    return
  end subroutine transform_quad_PX_int_i4

  !> @brief Row permutation transformation, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PY_full_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = 1, sz/2
       work(1: sz) = fcs(1: sz, jl)
       fcs(1: sz, jl) = fcs(1: sz, sz + 1 - jl)
       fcs(1: sz, sz + 1 - jl) = work(1: sz)
    end do
    ! or
    !do jl = 1, sz/2
    !   do il = 1, sz
    !      iface = fcs(il, jl)
    !      fcs(il, jl) = fcs(il, sz + 1 - jl)
    !      fcs(il, sz + 1 - jl) = iface
    !   end do
    !end do

    return
  end subroutine transform_quad_PY_full_i4

  !> @brief Row permutation transformation, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PY_int_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i4) :: iface

    do jl = 2, sz/2
       work(2: sz - 1) = fcs(2: sz - 1, jl)
       fcs(2: sz - 1, jl) = fcs(2: sz - 1, sz + 1 - jl)
       fcs(2: sz - 1, sz + 1 - jl) = work(2: sz - 1)
    end do
    ! or
    !do jl = 2, sz/2
    !   do il = 2, sz - 1
    !      iface = fcs(il, jl)
    !      fcs(il, jl) = fcs(il, sz + 1 - jl)
    !      fcs(il, sz + 1 - jl) = iface
    !   end do
    !end do

    return
  end subroutine transform_quad_PY_int_i4

  !> @brief PXT = TPY transformation, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXT_full_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_i4(sz, fcs, work)
    call transform_quad_T_full_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PXT_full_i4

  !> @brief PXT = TPY transformation, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXT_int_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_i4(sz, fcs, work)
    call transform_quad_T_int_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PXT_int_i4

  !> @brief PYT = TPX transformation, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PYT_full_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PY_full_i4(sz, fcs, work)
    call transform_quad_T_full_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PYT_full_i4

  !> @brief PYT = TPX transformation, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PYT_int_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PY_int_i4(sz, fcs, work)
    call transform_quad_T_int_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PYT_int_i4

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPYT_full_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_i4(sz, fcs, work)
    call transform_quad_PY_full_i4(sz, fcs, work)
    call transform_quad_T_full_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PXPYT_full_i4

  !> @brief  PXPYT=PYPXT=TPYPX=TPXPY, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPYT_int_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_i4(sz, fcs, work)
    call transform_quad_PY_int_i4(sz, fcs, work)
    call transform_quad_T_int_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PXPYT_int_i4

  !> @brief PXPY = PYPX transformation, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPY_full_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_i4(sz, fcs, work)
    call transform_quad_PY_full_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PXPY_full_i4

  !> @brief  PXPY = PYPX transformation, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPY_int_i4(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz, sz), intent(inout) :: fcs
    integer(i4), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_i4(sz, fcs, work)
    call transform_quad_PY_int_i4(sz, fcs, work)

    return
  end subroutine transform_quad_PXPY_int_i4

  !> @brief Identity transformation, double integer array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_I_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    return
  end subroutine transform_quad_I_i8

  !> @brief Transpose transformation, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_T_full_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = 1, sz
       do il = 1, jl -1
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(jl, il)
          fcs(jl, il) = iface
       end do
    end do

    return
  end subroutine transform_quad_T_full_i8

  !> @brief Transpose transformation, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_T_int_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = 2, sz - 1
       do il = 2, jl -1
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(jl, il)
          fcs(jl, il) = iface
       end do
    end do

    return
  end subroutine transform_quad_T_int_i8

  !> @brief Column permutation transformation, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PX_full_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = 1, sz
       do il = 1, sz/2
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(sz + 1 - il, jl)
          fcs(sz + 1 - il, jl) = iface
       end do
    end do

    return
  end subroutine transform_quad_PX_full_i8

  !> @brief Column permutation transformation, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PX_int_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = 2, sz - 1
       do il = 2, sz/2
          iface = fcs(il, jl)
          fcs(il, jl) = fcs(sz + 1 - il, jl)
          fcs(sz + 1 - il, jl) = iface
       end do
    end do

    return
  end subroutine transform_quad_PX_int_i8

  !> @brief Row permutation transformation, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PY_full_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = 1, sz/2
       work(1: sz) = fcs(1: sz, jl)
       fcs(1: sz, jl) = fcs(1: sz, sz + 1 - jl)
       fcs(1: sz, sz + 1 - jl) = work(1: sz)
    end do
    ! or
    !do jl = 1, sz/2
    !   do il = 1, sz
    !      iface = fcs(il, jl)
    !      fcs(il, jl) = fcs(il, sz + 1 - jl)
    !      fcs(il, sz + 1 - jl) = iface
    !   end do
    !end do

    return
  end subroutine transform_quad_PY_full_i8

  !> @brief Row permutation transformation, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PY_int_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    integer(i8) :: iface

    do jl = 2, sz/2
       work(2: sz - 1) = fcs(2: sz - 1, jl)
       fcs(2: sz - 1, jl) = fcs(2: sz - 1, sz + 1 - jl)
       fcs(2: sz - 1, sz + 1 - jl) = work(2: sz - 1)
    end do
    ! or
    !do jl = 2, sz/2
    !   do il = 2, sz - 1
    !      iface = fcs(il, jl)
    !      fcs(il, jl) = fcs(il, sz + 1 - jl)
    !      fcs(il, sz + 1 - jl) = iface
    !   end do
    !end do

    return
  end subroutine transform_quad_PY_int_i8

  !> @brief PXT = TPY transformation, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXT_full_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_i8(sz, fcs, work)
    call transform_quad_T_full_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PXT_full_i8

  !> @brief PXT = TPY transformation, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXT_int_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_i8(sz, fcs, work)
    call transform_quad_T_int_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PXT_int_i8

  !> @brief PYT = TPX transformation, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PYT_full_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PY_full_i8(sz, fcs, work)
    call transform_quad_T_full_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PYT_full_i8

  !> @brief PYT = TPX transformation, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PYT_int_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PY_int_i8(sz, fcs, work)
    call transform_quad_T_int_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PYT_int_i8

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPYT_full_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_i8(sz, fcs, work)
    call transform_quad_PY_full_i8(sz, fcs, work)
    call transform_quad_T_full_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PXPYT_full_i8

  !> @brief  PXPYT=PYPXT=TPYPX=TPXPY, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPYT_int_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_i8(sz, fcs, work)
    call transform_quad_PY_int_i8(sz, fcs, work)
    call transform_quad_T_int_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PXPYT_int_i8

  !> @brief PXPY = PYPX transformation, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPY_full_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_i8(sz, fcs, work)
    call transform_quad_PY_full_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PXPY_full_i8

  !> @brief  PXPY = PYPX transformation, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPY_int_i8(sz, fcs, work)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz, sz), intent(inout) :: fcs
    integer(i8), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_i8(sz, fcs, work)
    call transform_quad_PY_int_i8(sz, fcs, work)

    return
  end subroutine transform_quad_PXPY_int_i8

  !> @brief Identity transformation, double precision array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_I_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    return
  end subroutine transform_quad_I_dp

  !> @brief Transpose transformation, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_T_full_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = 1, sz
       do il = 1, jl -1
          rface = fcs(il, jl)
          fcs(il, jl) = fcs(jl, il)
          fcs(jl, il) = rface
       end do
    end do

    return
  end subroutine transform_quad_T_full_dp

  !> @brief Transpose transformation, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_T_int_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = 2, sz - 1
       do il = 2, jl -1
          rface = fcs(il, jl)
          fcs(il, jl) = fcs(jl, il)
          fcs(jl, il) = rface
       end do
    end do

    return
  end subroutine transform_quad_T_int_dp

  !> @brief Column permutation transformation, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PX_full_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = 1, sz
       do il = 1, sz/2
          rface = fcs(il, jl)
          fcs(il, jl) = fcs(sz + 1 - il, jl)
          fcs(sz + 1 - il, jl) = rface
       end do
    end do

    return
  end subroutine transform_quad_PX_full_dp

  !> @brief Column permutation transformation, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PX_int_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = 2, sz - 1
       do il = 2, sz/2
          rface = fcs(il, jl)
          fcs(il, jl) = fcs(sz + 1 - il, jl)
          fcs(sz + 1 - il, jl) = rface
       end do
    end do

    return
  end subroutine transform_quad_PX_int_dp

  !> @brief Row permutation transformation, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PY_full_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = 1, sz/2
       work(1: sz) = fcs(1: sz, jl)
       fcs(1: sz, jl) = fcs(1: sz, sz + 1 - jl)
       fcs(1: sz, sz + 1 - jl) = work(1: sz)
    end do
    ! or
    !do jl = 1, sz/2
    !   do il = 1, sz
    !      rface = fcs(il, jl)
    !      fcs(il, jl) = fcs(il, sz + 1 - jl)
    !      fcs(il, sz + 1 - jl) = rface
    !   end do
    !end do

    return
  end subroutine transform_quad_PY_full_dp

  !> @brief Row permutation transformation, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PY_int_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work
    ! local variables
    integer(i4) :: il, jl
    real(dp) :: rface

    do jl = 2, sz/2
       work(2: sz - 1) = fcs(2: sz - 1, jl)
       fcs(2: sz - 1, jl) = fcs(2: sz - 1, sz + 1 - jl)
       fcs(2: sz - 1, sz + 1 - jl) = work(2: sz - 1)
    end do
    ! or
    !do jl = 2, sz/2
    !   do il = 2, sz - 1
    !      rface = fcs(il, jl)
    !      fcs(il, jl) = fcs(il, sz + 1 - jl)
    !      fcs(il, sz + 1 - jl) = rface
    !   end do
    !end do

    return
  end subroutine transform_quad_PY_int_dp

  !> @brief PXT = TPY transformation, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXT_full_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_dp(sz, fcs, work)
    call transform_quad_T_full_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PXT_full_dp

  !> @brief PXT = TPY transformation, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXT_int_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_dp(sz, fcs, work)
    call transform_quad_T_int_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PXT_int_dp

  !> @brief PYT = TPX transformation, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PYT_full_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PY_full_dp(sz, fcs, work)
    call transform_quad_T_full_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PYT_full_dp

  !> @brief PYT = TPX transformation, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PYT_int_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PY_int_dp(sz, fcs, work)
    call transform_quad_T_int_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PYT_int_dp

  !> @brief PXPYT=PYPXT=TPYPX=TPXPY, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPYT_full_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_dp(sz, fcs, work)
    call transform_quad_PY_full_dp(sz, fcs, work)
    call transform_quad_T_full_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PXPYT_full_dp

  !> @brief  PXPYT=PYPXT=TPYPX=TPXPY, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPYT_int_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_dp(sz, fcs, work)
    call transform_quad_PY_int_dp(sz, fcs, work)
    call transform_quad_T_int_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PXPYT_int_dp

  !> @brief PXPY = PYPX transformation, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPY_full_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PX_full_dp(sz, fcs, work)
    call transform_quad_PY_full_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PXPY_full_dp

  !> @brief  PXPY = PYPX transformation, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  !! @parameter[inout]  work     work space
  pure subroutine transform_quad_PXPY_int_dp(sz, fcs, work)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs
    real(dp), dimension(sz), intent(inout) :: work

    call transform_quad_PX_int_dp(sz, fcs, work)
    call transform_quad_PY_int_dp(sz, fcs, work)

    return
  end subroutine transform_quad_PXPY_int_dp

end module alignment_quad
