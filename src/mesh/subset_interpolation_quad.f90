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
module subset_interpolation_quad
  use num_types, only : i4, rp
  use subset_interpolation, only : subset_interpolation_t
  implicit none
  private

  public :: subset_intp_quad_init, subset_interpolation_quad_dmy_t, &
       & subset_interpolation_quad_hng_t

  !> number of various interpolation operators
  integer(i4), public, parameter :: NEKO_INTP_QUAD_NOPERATION = 4

  !> Dummy interpolation operators
  type, extends(subset_interpolation_t) :: subset_interpolation_quad_dmy_t
     !> Array size in local "r" direction
     integer(i4), private :: lr_ = -1
     !> Array size in local "s" direction
     integer(i4), private :: ls_ = -1
   contains
     !> Is transformation a dummy one
     procedure, pass(this) :: ifdmy => ifdummy_quad_dmy
     !> Patent-child interpolation
     procedure, nopass :: intp  => transform_quad_dmy
     !> Transposed interpolation
     procedure, nopass :: intpT => transform_quad_dmy
     !> Initialise interpolation data
     procedure, pass(this) :: set_jmat  => quad_set_jmat_dmy
     !> Free interpolation data
     procedure, pass(this) :: free_jmat => quad_free_jmat_dmy
  end type subset_interpolation_quad_dmy_t

  !> Quad interpolation operators for hanging values 1 to 4
  type, extends(subset_interpolation_t) :: subset_interpolation_quad_hng_t
     !> Array size in local "r" direction
     integer(i4), private :: lr_ = -1
     !> Array size in local "s" direction
     integer(i4), private :: ls_ = -1
     !> Interpolation data "r"
     real(rp), private, dimension(:, :), allocatable :: jmatr_
     !> Interpolation data "s"
     real(rp), private, dimension(:, :), allocatable :: jmats_
   contains
     !> Is transformation a dummy one
     procedure, pass(this) :: ifdmy => ifdummy_quad_hng
     !> Patent-child interpolation
     procedure, nopass :: intp  => transform_quad_hng
     !> Transposed interpolation
     procedure, nopass :: intpT => transform_trn_quad_hng
     !> Initialise interpolation data
     procedure, pass(this) :: set_jmat  => quad_set_jmat_hng
     !> Free interpolation data
     procedure, pass(this) :: free_jmat => quad_free_jmat_hng
  end type subset_interpolation_quad_hng_t

contains

  !> @brief Allocate a single interpolation operator
  !! @parameter[in]      hng   hanging information
  !! @parameter[inout]   trns  interpolation operator
  subroutine subset_intp_quad_init(hng, trns)
    integer(i4), intent(in) :: hng
    class(subset_interpolation_t), allocatable, intent(inout) :: trns

    if (allocated(trns)) then
       call trns%free_jmat()
       deallocate(trns)
    end if

    select case(hng)
    case(1:4) ! the quarter corresponding to the parent vertex
       allocate(subset_interpolation_quad_hng_t :: trns)
       call trns%set_hng(hng)
    case default ! any other option is just a dummy operation
       allocate(subset_interpolation_quad_dmy_t :: trns)
       call trns%set_hng(hng)
    end select

  end subroutine subset_intp_quad_init

  !> Function returning dummy operation flag
  !! @return   ifdmy
  pure function ifdummy_quad_dmy(this) result(ifdmy)
    class(subset_interpolation_quad_dmy_t), intent(in) :: this
    logical :: ifdmy
    ifdmy = .true.
  end function ifdummy_quad_dmy

  !> Function returning dummy operation flag
  !! @return   ifdmy
  pure function ifdummy_quad_hng(this) result(ifdmy)
    class(subset_interpolation_quad_hng_t), intent(in) :: this
    logical :: ifdmy
    ifdmy = .false.
  end function ifdummy_quad_hng

  !> Dummy interpolation operator, do nothing
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  pure subroutine transform_quad_dmy(vec, n1, n2)
    integer(i4), intent(in) :: n1, n2
    real(rp), dimension(n1, n2), intent(inout) :: vec
  end subroutine transform_quad_dmy

  !> Dummy initialisation of the interpolation data
  !! @notice  This routine can be called after space_t gets initialised.
  !! @parameter[in]      lr, ls   array sizes for r and s dimensions
  subroutine quad_set_jmat_dmy(this, lr, ls)
    class(subset_interpolation_quad_dmy_t), intent(inout) :: this
    integer(i4), intent(in) :: lr, ls

    call this%free_jmat()

    this%lr_ = lr
    this%ls_ = ls
  end subroutine quad_set_jmat_dmy

  !> Free the dummy interpolation data
  subroutine quad_free_jmat_dmy(this)
    class(subset_interpolation_quad_dmy_t), intent(inout) :: this
    this%lr_ = -1
    this%ls_ = -1
  end subroutine quad_free_jmat_dmy

  !> Parent-child interpolation
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  pure subroutine transform_quad_hng(vec, n1, n2)
    integer(i4), intent(in) :: n1, n2
    real(rp), dimension(n1, n2), intent(inout) :: vec
    ! will be added later
  end subroutine transform_quad_hng

  !> Transposed interpolation
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  pure subroutine transform_trn_quad_hng(vec, n1, n2)
    integer(i4), intent(in) :: n1, n2
    real(rp), dimension(n1, n2), intent(inout) :: vec
    ! will be added later
  end subroutine transform_trn_quad_hng

  !> Dummy initialisation of the interpolation data
  !! @notice This routine can be called after space_t gets initialised.
  !! @parameter[in]      lr, ls   array sizes for r and s dimensions
  subroutine quad_set_jmat_hng(this, lr, ls)
    class(subset_interpolation_quad_hng_t), intent(inout) :: this
    integer(i4), intent(in) :: lr, ls
    integer(i4) :: hng

    call this%free_jmat()

    this%lr_ = lr
    this%ls_ = ls
    allocate(this%jmatr_(lr, lr), this%jmats_(ls, ls))
    hng = this%hng()
    select case(hng)
    case(1) ! the first parent face corner
       ! will be added later
    case(2) ! the second parent face corner
       ! will be added later
    case(3) ! the third parent face corner
       ! will be added later
    case(4) ! the fourth parent face corner
       ! will be added later
    end select
  end subroutine quad_set_jmat_hng

  !> Free the interpolation data
  subroutine quad_free_jmat_hng(this)
    class(subset_interpolation_quad_hng_t), intent(inout) :: this
    if (allocated(this%jmatr_)) deallocate(this%jmatr_)
    if (allocated(this%jmats_)) deallocate(this%jmats_)
    this%lr_ = -1
    this%ls_ = -1
  end subroutine quad_free_jmat_hng

end module subset_interpolation_quad
