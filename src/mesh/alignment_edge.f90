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
!> Edge alignment operators
module alignment_edge
  use num_types, only : i4, i8, dp
  use utils, only : neko_error
  implicit none
  private

  public :: alignment_edge_t, alignment_edge_op_set_t

  !> number of operations different from identity
  integer(i4), public, parameter :: NEKO_EDGE_NOPERATION = 1

  !> procedure pointer type; i4
  type :: algn_edge_proc_i4_ptr
     procedure(transform_i4), pointer, nopass :: obj
  end type algn_edge_proc_i4_ptr
  !> procedure pointer type; i8
  type :: algn_edge_proc_i8_ptr
     procedure(transform_i8), pointer, nopass :: obj
  end type algn_edge_proc_i8_ptr
  !> procedure pointer type; dp
  type :: algn_edge_proc_dp_ptr
     procedure(transform_dp), pointer, nopass :: obj
  end type algn_edge_proc_dp_ptr

  !> Type containing set of edge alignment operators
  !! @details There are two main operations : identity (I) and row
  !! permutation (P).
  !! @note In this case inverse operators are not necessary, but I keep
  !! code consistent with quad structure. The same about identity operation.
  !! It is not really needed, but i keep it just to have it complete.
  type :: alignment_edge_t
     !> number of different operations excluding identity
     integer(i4), private :: noperation_ = NEKO_EDGE_NOPERATION
     !> Direct array transformations for whole array
     type(algn_edge_proc_i4_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_f_i4
     type(algn_edge_proc_i8_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_f_i8
     type(algn_edge_proc_dp_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_f_dp
     !> Direct array transformations for array interior
     type(algn_edge_proc_i4_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_i_i4
     type(algn_edge_proc_i8_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_i_i8
     type(algn_edge_proc_dp_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_i_dp
     !> Inverse array transformations for whole array
     type(algn_edge_proc_i4_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_inv_f_i4
     type(algn_edge_proc_i8_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_inv_f_i8
     type(algn_edge_proc_dp_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_inv_f_dp
     !> Inverse array transformations for array interior
     type(algn_edge_proc_i4_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_inv_i_i4
     type(algn_edge_proc_i8_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_inv_i_i8
     type(algn_edge_proc_dp_ptr),&
          & dimension(0 : NEKO_EDGE_NOPERATION) :: trns_inv_i_dp
   contains
     !> Initialise procedure pointers
     procedure, pass(this) :: init => edge_init
     !> Free procedure pointers
     procedure, pass(this) :: free => edge_free
     !> Return number of operations
     procedure, pass(this) :: nop => edge_noperation_get
  end type alignment_edge_t

  !> Type combining a single set of alignment operators
  type :: alignment_edge_op_set_t
     !> Relative edge alignment
     integer(i4) :: alignment = -1
     !> Direct array transformations for whole array
     type(algn_edge_proc_i4_ptr) :: trns_f_i4
     type(algn_edge_proc_i8_ptr) :: trns_f_i8
     type(algn_edge_proc_dp_ptr) :: trns_f_dp
     !> Direct array transformations for array interior
     type(algn_edge_proc_i4_ptr) :: trns_i_i4
     type(algn_edge_proc_i8_ptr) :: trns_i_i8
     type(algn_edge_proc_dp_ptr) :: trns_i_dp
     !> Inverse array transformations for whole array
     type(algn_edge_proc_i4_ptr) :: trns_inv_f_i4
     type(algn_edge_proc_i8_ptr) :: trns_inv_f_i8
     type(algn_edge_proc_dp_ptr) :: trns_inv_f_dp
     !> Inverse array transformations for array interior
     type(algn_edge_proc_i4_ptr) :: trns_inv_i_i4
     type(algn_edge_proc_i8_ptr) :: trns_inv_i_i8
     type(algn_edge_proc_dp_ptr) :: trns_inv_i_dp
   contains
      !> Initialise alignment and procedure pointers
     procedure, pass(this) :: init => edge_op_set_init
     !> Free alignment and pointers
     procedure, pass(this) :: free => edge_op_set_free
  end type alignment_edge_op_set_t

  ! Abstract types for different transformations; various types
  abstract interface
     pure subroutine transform_i4(sz, edg)
       import i4
       integer(i4), intent(in) :: sz
       integer(i4), dimension(sz), intent(inout) :: edg
     end subroutine transform_i4

     pure subroutine transform_i8(sz, edg)
       import i4
       import i8
       integer(i4), intent(in) :: sz
       integer(i8), dimension(sz), intent(inout) :: edg
     end subroutine transform_i8

     pure subroutine transform_dp(sz, edg)
       import i4
       import dp
       integer(i4), intent(in) :: sz
       real(dp), dimension(sz), intent(inout) :: edg
     end subroutine transform_dp
  end interface

contains
  !> @brief Initialise procedure pointers
  subroutine edge_init(this)
    class(alignment_edge_t), intent(inout) :: this

    call this%free()

    ! Identity transformation is added for completeness; in general not needed
    ! Direct transformation of full array, different types
    this%trns_f_i4(0)%obj => transform_edge_I_i4
    this%trns_f_i4(1)%obj => transform_edge_P_full_i4
    this%trns_f_i8(0)%obj => transform_edge_I_i8
    this%trns_f_i8(1)%obj => transform_edge_P_full_i8
    this%trns_f_dp(0)%obj => transform_edge_I_dp
    this%trns_f_dp(1)%obj => transform_edge_P_full_dp
    ! Direct transformation of array interior, different types
    this%trns_i_i4(0)%obj => transform_edge_I_i4
    this%trns_i_i4(1)%obj => transform_edge_P_int_i4
    this%trns_i_i8(0)%obj => transform_edge_I_i8
    this%trns_i_i8(1)%obj => transform_edge_P_int_i8
    this%trns_i_dp(0)%obj => transform_edge_I_dp
    this%trns_i_dp(1)%obj => transform_edge_P_int_dp
    ! Inverse transformation of full array, different types
    this%trns_inv_f_i4(0)%obj => transform_edge_I_i4
    this%trns_inv_f_i4(1)%obj => transform_edge_P_full_i4
    this%trns_inv_f_i8(0)%obj => transform_edge_I_i8
    this%trns_inv_f_i8(1)%obj => transform_edge_P_full_i8
    this%trns_inv_f_dp(0)%obj => transform_edge_I_dp
    this%trns_inv_f_dp(1)%obj => transform_edge_P_full_dp
    ! Inverse transformation of array interior, different types
    this%trns_inv_i_i4(0)%obj => transform_edge_I_i4
    this%trns_inv_i_i4(1)%obj => transform_edge_P_int_i4
    this%trns_inv_i_i8(0)%obj => transform_edge_I_i8
    this%trns_inv_i_i8(1)%obj => transform_edge_P_int_i8
    this%trns_inv_i_dp(0)%obj => transform_edge_I_dp
    this%trns_inv_i_dp(1)%obj => transform_edge_P_int_dp

    return
  end subroutine edge_init

  !> @brief Free procedure pointers
  subroutine edge_free(this)
    class(alignment_edge_t), intent(inout) :: this
    integer(i4) :: il

    ! free pointers
    do il = 0, NEKO_EDGE_NOPERATION
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
  end subroutine edge_free

  !> @brief Get number of operations
  !! @return   noperation
  pure function edge_noperation_get(this) result(noperation)
    class(alignment_edge_t), intent(in) :: this
    integer(i4) :: noperation
    noperation = this%noperation_
  end function edge_noperation_get

  !> @brief Initialise alignment and procedure pointers
  !! @parameter[in]   algn   alignment
  subroutine edge_op_set_init(this, algn)
    class(alignment_edge_op_set_t), intent(inout) :: this
    integer(i4), intent(in) :: algn

    call this%free()

    ! set relative alignment transformation
    if ((algn >= 0).and.(algn <= NEKO_EDGE_NOPERATION)) then
       this%alignment = algn
    else
       call neko_error('Not proper alignment.')
    end if
    select case(algn)
    case(0)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_edge_I_i4
       this%trns_f_i8%obj => transform_edge_I_i8
       this%trns_f_dp%obj => transform_edge_I_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_edge_I_i4
       this%trns_i_i8%obj => transform_edge_I_i8
       this%trns_i_dp%obj => transform_edge_I_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_edge_I_i4
       this%trns_inv_f_i8%obj => transform_edge_I_i8
       this%trns_inv_f_dp%obj => transform_edge_I_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_edge_I_i4
       this%trns_inv_i_i8%obj => transform_edge_I_i8
       this%trns_inv_i_dp%obj => transform_edge_I_dp
    case(1)
       ! Direct transformation of full array, different types
       this%trns_f_i4%obj => transform_edge_P_full_i4
       this%trns_f_i8%obj => transform_edge_P_full_i8
       this%trns_f_dp%obj => transform_edge_P_full_dp
       ! Direct transformation of array interior, different types
       this%trns_i_i4%obj => transform_edge_P_int_i4
       this%trns_i_i8%obj => transform_edge_P_int_i8
       this%trns_i_dp%obj => transform_edge_P_int_dp
       ! Inverse transformation of full array, different types
       this%trns_inv_f_i4%obj => transform_edge_P_full_i4
       this%trns_inv_f_i8%obj => transform_edge_P_full_i8
       this%trns_inv_f_dp%obj => transform_edge_P_full_dp
       ! Inverse transformation of array interior, different types
       this%trns_inv_i_i4%obj => transform_edge_P_int_i4
       this%trns_inv_i_i8%obj => transform_edge_P_int_i8
       this%trns_inv_i_dp%obj => transform_edge_P_int_dp
    end select

    return
  end subroutine edge_op_set_init

  !> @brief Free alignment and procedure pointers
  subroutine edge_op_set_free(this)
    class(alignment_edge_op_set_t), intent(inout) :: this

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
  end subroutine edge_op_set_free

  !> @brief Identity transformation, single integer array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_I_i4(sz, edg)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz), intent(inout) :: edg

    return
  end subroutine transform_edge_I_i4

  !> @brief Permutation transformation, single integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_P_full_i4(sz, edg)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz), intent(inout) :: edg
    ! local variables
    integer(i4) :: istart, il, itmp1, itmp2
    integer(i4) :: iedg

    itmp1 = sz + 1
    do il = 1, sz/2
       itmp2 = itmp1 - il
       iedg = edg(il)
       edg(il) = edg(itmp2)
       edg(itmp2) = iedg
    end do

    return
  end subroutine transform_edge_P_full_i4

  !> @brief Permutation transformation, single integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_P_int_i4(sz, edg)
    integer(i4), intent(in) :: sz
    integer(i4), dimension(sz), intent(inout) :: edg
    ! local variables
    integer(i4) :: istart, il, itmp1, itmp2
    integer(i4) :: iedg

    itmp1 = sz + 1
    do il = 2, sz/2
       itmp2 = itmp1 - il
       iedg = edg(il)
       edg(il) = edg(itmp2)
       edg(itmp2) = iedg
    end do

    return
  end subroutine transform_edge_P_int_i4

  !> @brief Identity transformation, double integer array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_I_i8(sz, edg)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz), intent(inout) :: edg

    return
  end subroutine transform_edge_I_i8

  !> @brief Permutation transformation, double integer, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_P_full_i8(sz, edg)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz), intent(inout) :: edg
    ! local variables
    integer(i4) :: il, itmp1, itmp2
    integer(i8) :: iedg

    itmp1 = sz + 1
    do il = 1, sz/2
       itmp2 = itmp1 - il
       iedg = edg(il)
       edg(il) = edg(itmp2)
       edg(itmp2) = iedg
    end do

    return
  end subroutine transform_edge_P_full_i8

  !> @brief Permutation transformation, double integer, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_P_int_i8(sz, edg)
    integer(i4), intent(in) :: sz
    integer(i8), dimension(sz), intent(inout) :: edg
    ! local variables
    integer(i4) :: il, itmp1, itmp2
    integer(i8) :: iedg

    itmp1 = sz + 1
    do il = 2, sz/2
       itmp2 = itmp1 - il
       iedg = edg(il)
       edg(il) = edg(itmp2)
       edg(itmp2) = iedg
    end do

    return
  end subroutine transform_edge_P_int_i8

  !> @brief Identity transformation, double precision array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_I_dp(sz, edg)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz), intent(inout) :: edg

    return
  end subroutine transform_edge_I_dp

  !> @brief Permutation transformation, double precision, full array
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_P_full_dp(sz, edg)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz), intent(inout) :: edg
    ! local variables
    integer(i4) :: il, itmp1, itmp2
    real(dp) :: redg

    itmp1 = sz + 1
    do il = 1, sz/2
       itmp2 = itmp1 - il
       redg = edg(il)
       edg(il) = edg(itmp2)
       edg(itmp2) = redg
    end do

    return
  end subroutine transform_edge_P_full_dp

  !> @brief Permutation transformation, double precision, array interior
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  edg      edge data
  pure subroutine transform_edge_P_int_dp(sz, edg)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz), intent(inout) :: edg
    ! local variables
    integer(i4) :: il, itmp1, itmp2
    real(dp) :: redg

    itmp1 = sz + 1
    do il = 2, sz/2
       itmp2 = itmp1 - il
       redg = edg(il)
       edg(il) = edg(itmp2)
       edg(itmp2) = redg
    end do

    return
  end subroutine transform_edge_P_int_dp

end module alignment_edge
