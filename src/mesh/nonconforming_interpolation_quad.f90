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
!> Interpolation operators for nonconforming quads
module ncnf_interpolation_quad
  use num_types, only : i4, dp
  implicit none
  private

  public :: ncnf_intp_quad_op_set_t

  !> number of various interpolation operators
  integer(i4), public, parameter :: NEKO_INTP_QUAD_NOPERATION = 4

  !> procedure pointer type; dp
  type :: intp_quad_proc_dp_ptr
     procedure(transform_dp), pointer, nopass :: ptr
  end type intp_quad_proc_dp_ptr

  !> Type combining a single set of interpolation operators
  type :: ncnf_intp_quad_op_set_t
     !> Hanging information
     integer(i4) :: hanging = -1
     !> Direct interpolation
     type(intp_quad_proc_dp_ptr) :: intp
     !> Transposed interpolation
     type(intp_quad_proc_dp_ptr) :: intpT
   contains
     !> Initialise hanging info and procedure pointers
     procedure, pass(this) :: init => quad_op_set_init
     !> Free hanging info and pointers
     procedure, pass(this) :: free => quad_op_set_free
  end type ncnf_intp_quad_op_set_t

  ! Abstract types for various transformations
  abstract interface
     pure subroutine transform_dp(sz, fcs)
       import i4
       import dp
       integer(i4), intent(in) :: sz
       real(dp), dimension(sz, sz), intent(inout) :: fcs
     end subroutine transform_dp
  end interface

contains

  !> @brief Initialise hanging info and procedure pointers
  !! @parameter[in]   hng   hanging information
  subroutine quad_op_set_init(this, hng)
    class(ncnf_intp_quad_op_set_t), intent(inout) :: this
    integer(i4), intent(in) :: hng

    call this%free()
    this%hanging = hng
    select case(hng)
    case(1) ! first parent face corner
       ! for now just identity
       this%intp%ptr => transform_quad_I
       this%intpT%ptr => transform_quad_I
    case(2) ! second parent face corner
       ! for now just identity
       this%intp%ptr => transform_quad_I
       this%intpT%ptr => transform_quad_I
    case(3) ! third parent face corner
       ! for now just identity
       this%intp%ptr => transform_quad_I
       this%intpT%ptr => transform_quad_I
    case(4) ! fourth parent face corner
       ! for now just identity
       this%intp%ptr => transform_quad_I
       this%intpT%ptr => transform_quad_I
    case default ! any other option is just identity operation
       this%intp%ptr => transform_quad_I
       this%intpT%ptr => transform_quad_I
    end select
  end subroutine quad_op_set_init

  !> @brief Free hanging info and procedure pointers
  subroutine quad_op_set_free(this)
    class(ncnf_intp_quad_op_set_t), intent(inout) :: this

    this%hanging = -1
    ! free pointers
    this%intp%ptr => null()
    this%intpT%ptr => null()

  end subroutine quad_op_set_free

  !> @brief Identity transformation
  !! @parameter[in]     sz       array size
  !! @parameter[inout]  fcs      face data
  pure subroutine transform_quad_I(sz, fcs)
    integer(i4), intent(in) :: sz
    real(dp), dimension(sz, sz), intent(inout) :: fcs

  end subroutine transform_quad_I

end module ncnf_interpolation_quad
