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
!> Interpolation operators for nonconforming edges
module ncnf_interpolation_edge
  use num_types, only : i4, dp
  use ncnf_interpolation, only : ncnf_interpolation_t
  implicit none
  private

  public :: ncnf_intp_edge_init, ncnf_interpolation_edge_dmy_t, &
       & ncnf_interpolation_edge_hng_t

  !> number of various interpolation operators
  integer(i4), public, parameter :: NEKO_INTP_EDGE_NOPERATION = 2

  !> Dummy interpolation operators
  type, extends(ncnf_interpolation_t) :: ncnf_interpolation_edge_dmy_t
     !> Array size
     integer(i4), private :: lr_ = -1
   contains
     !> Patent-child interpolation
     procedure, nopass :: intp  => transform_edge_dmy
     !> Transposed interpolation
     procedure, nopass :: intpT => transform_edge_dmy
     !> Initialise interpolation data
     procedure, pass(this) :: set_jmat  => edge_set_jmat_dmy
     !> Free interpolation data
     procedure, pass(this) :: free_jmat => edge_free_jmat_dmy
  end type ncnf_interpolation_edge_dmy_t

  !> Edge interpolation operators for hanging values 1 and 2
  type, extends(ncnf_interpolation_t) :: ncnf_interpolation_edge_hng_t
     !> Array size
     integer(i4), private :: lr_ = -1
     !> Interpolation data
     real(dp), private, dimension(:, :), allocatable :: jmatr_
   contains
     !> Patent-child interpolation
     procedure, nopass :: intp  => transform_edge_hng
     !> Transposed interpolation
     procedure, nopass :: intpT => transform_trn_edge_hng
     !> Initialise interpolation data
     procedure, pass(this) :: set_jmat  => edge_set_jmat_hng
     !> Free interpolation data
     procedure, pass(this) :: free_jmat => edge_free_jmat_hng
  end type ncnf_interpolation_edge_hng_t

contains

  !> @brief Allocate a single interpolation operator
  !! @parameter[in]      hng   hanging information
  !! @parameter[inout]   trns  interpolation operator
  subroutine ncnf_intp_edge_init(hng, trns)
    integer(i4), intent(in) :: hng
    class(ncnf_interpolation_t), allocatable, intent(inout) :: trns

    if (allocated(trns)) then
       call trns%free_jmat()
       deallocate(trns)
    end if

    select case(hng)
    case(1:2) ! the first and the second half of the full edge
       allocate(ncnf_interpolation_edge_hng_t :: trns)
       call trns%set_hng(hng)
    case default ! any other option is just a dummy operation
       allocate(ncnf_interpolation_edge_dmy_t :: trns)
       call trns%set_hng(hng)
    end select

  end subroutine ncnf_intp_edge_init

  !> Dummy interpolation operator, do nothing
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  pure subroutine transform_edge_dmy(vec, n1, n2)
    integer(i4), intent(in) :: n1, n2
    real(dp), dimension(n1, n2), intent(inout) :: vec
  end subroutine transform_edge_dmy

  !> Dummy initialisation of the interpolation data
  !! @notice It is a common interface for 1D and 2D operations, so there
  !! are two dimensions @a lr and @a ls. This routine can be called
  !! after space_t gets initialised.
  !! @parameter[in]      lr, ls   array sizes for r and s dimensions
  subroutine edge_set_jmat_dmy(this, lr, ls)
    class(ncnf_interpolation_edge_dmy_t), intent(inout) :: this
    integer(i4), intent(in) :: lr, ls

    call this%free_jmat()

    this%lr_ = lr
  end subroutine edge_set_jmat_dmy

  !> Free the dummy interpolation data
  subroutine edge_free_jmat_dmy(this)
    class(ncnf_interpolation_edge_dmy_t), intent(inout) :: this
    this%lr_ = -1
  end subroutine edge_free_jmat_dmy

  !> Parent-child interpolation
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  pure subroutine transform_edge_hng(vec, n1, n2)
    integer(i4), intent(in) :: n1, n2
    real(dp), dimension(n1, n2), intent(inout) :: vec
    ! will be added later
  end subroutine transform_edge_hng

  !> Transposed interpolation
  !! @notice It is a common interface for 1D and 2D operations, so the data
  !! array @a vec is rank 2 even for 1D operations.
  !! @parameter[inout]   vec      data vector
  !! @parameter[in]      n1, n2   dimensions
  pure subroutine transform_trn_edge_hng(vec, n1, n2)
    integer(i4), intent(in) :: n1, n2
    real(dp), dimension(n1, n2), intent(inout) :: vec
    ! will be added later
  end subroutine transform_trn_edge_hng

  !> Dummy initialisation of the interpolation data
  !! @notice It is a common interface for 1D and 2D operations, so there
  !! are two dimensions @a lr and @a ls. This routine can be called
  !! after space_t gets initialised.
  !! @parameter[in]      lr, ls   array sizes for r and s dimensions
  subroutine edge_set_jmat_hng(this, lr, ls)
    class(ncnf_interpolation_edge_hng_t), intent(inout) :: this
    integer(i4), intent(in) :: lr, ls
    integer(i4) :: hng

    call this%free_jmat()

    this%lr_ = lr
    allocate(this%jmatr_(lr, lr))
    hng = this%hng()
    select case(hng)
    case(1) ! the first half of the full edge
       ! will be added later
    case(2) ! the second half of the full edge
       ! will be added later
    end select
  end subroutine edge_set_jmat_hng

  !> Free the interpolation data
  subroutine edge_free_jmat_hng(this)
    class(ncnf_interpolation_edge_hng_t), intent(inout) :: this
    if (allocated(this%jmatr_)) deallocate(this%jmatr_)
    this%lr_ = -1
  end subroutine edge_free_jmat_hng

end module ncnf_interpolation_edge
