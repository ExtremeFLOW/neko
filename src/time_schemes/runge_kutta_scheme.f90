! Copyright (c) 2025, The Neko Authors
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
module runge_kutta_time_scheme
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use time_scheme, only: time_scheme_t
  use math, only : rzero
  use utils, only : neko_error
  use device, only : device_free, device_map, device_memcpy, HOST_TO_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr, c_associated, C_NULL_PTR
  implicit none
  private

  type, public :: runge_kutta_time_scheme_t
     real(kind=rp), allocatable :: coeffs_A(:, :)
     real(kind=rp), allocatable :: coeffs_b(:)
     real(kind=rp), allocatable :: coeffs_c(:)
     type(c_ptr) :: coeffs_A_d = C_NULL_PTR
     type(c_ptr) :: coeffs_b_d = C_NULL_PTR
     type(c_ptr) :: coeffs_c_d = C_NULL_PTR
     integer :: order
   contains
     !> Constructor
     procedure, pass(this) :: init => runge_kutta_scheme_coeffs_init
     !> Destructor
     procedure, pass(this) :: free => runge_kutta_scheme_coeffs_free
  end type runge_kutta_time_scheme_t

contains
  !> Constructor for Runge-Kutta time scheme
  !> @param this The Runge-Kutta scheme object
  !> @param order Order of accuracy (1-4), determines coefficients:
  !>              1: Forward Euler
  !>              2: Heun's method
  !>              3: SSPRK3
  !>              4: Classic RK4
  subroutine runge_kutta_scheme_coeffs_init(this, order)
    class(runge_kutta_time_scheme_t), intent(inout) :: this
    integer, intent(in) :: order
    integer :: s

    ! Assume number of stages is equal to the order
    s = order
    this%order = order
    allocate(this%coeffs_A(order, order))
    allocate(this%coeffs_b(order))
    allocate(this%coeffs_c(order))

    associate( &
         coeffs_A => this%coeffs_A, &
         coeffs_b => this%coeffs_b, &
         coeffs_c => this%coeffs_c, &
         coeffs_A_d => this%coeffs_A_d, &
         coeffs_b_d => this%coeffs_b_d, &
         coeffs_c_d => this%coeffs_c_d)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_map(coeffs_A, coeffs_A_d, s)
         call device_map(coeffs_b, coeffs_b_d, s)
         call device_map(coeffs_c, coeffs_c_d, s)
      end if

      select case (order)
      case (1)
         !> Forward Euler
         coeffs_A(1, 1) = 0.0_rp
         coeffs_b(1) = 1.0_rp
         coeffs_c(1) = 0.0_rp
      case (2)
         !> Heun's method
         coeffs_A = 0.0_rp
         coeffs_A(2, 1) = 1.0_rp

         coeffs_b(1) = 0.5_rp
         coeffs_b(2) = 0.5_rp

         coeffs_c(1) = 0.0_rp
         coeffs_c(2) = 1.0_rp
      case (3)
         !> SSPRK3
         coeffs_A = 0.0_rp
         coeffs_A(2, 1) = 1.0_rp
         coeffs_A(3, 1) = 0.25_rp
         coeffs_A(3, 2) = 0.25_rp

         coeffs_b(1) = 1.0_rp / 6.0_rp
         coeffs_b(2) = 1.0_rp / 6.0_rp
         coeffs_b(3) = 2.0_rp / 3.0_rp

         coeffs_c(1) = 0.0_rp
         coeffs_c(2) = 1.0_rp
         coeffs_c(3) = 0.5_rp
      case (4)
         !> RK4
         coeffs_A = 0.0_rp
         coeffs_A(2, 1) = 0.5_rp
         coeffs_A(3, 2) = 0.5_rp
         coeffs_A(4, 3) = 1.0_rp

         coeffs_b(1) = 1.0_rp / 6.0_rp
         coeffs_b(2) = 1.0_rp / 3.0_rp
         coeffs_b(3) = 1.0_rp / 3.0_rp
         coeffs_b(4) = 1.0_rp / 6.0_rp

         coeffs_c(1) = 0.0_rp
         coeffs_c(2) = 0.5_rp
         coeffs_c(3) = 0.5_rp
         coeffs_c(4) = 1.0_rp
      case default
         call neko_error("The time order must be 1 to 4.")
      end select

      if (c_associated(coeffs_A_d)) then
         call device_memcpy(coeffs_A, coeffs_A_d, s, &
              HOST_TO_DEVICE, sync = .false.)
      end if

      if (c_associated(coeffs_b_d)) then
         call device_memcpy(coeffs_b, coeffs_b_d, s, &
              HOST_TO_DEVICE, sync = .false.)
      end if

      if (c_associated(coeffs_c_d)) then
         call device_memcpy(coeffs_c, coeffs_c_d, s, &
              HOST_TO_DEVICE, sync = .false.)
      end if
    end associate

  end subroutine runge_kutta_scheme_coeffs_init

  !> Destructor - deallocates device memory
  !> @param this The Runge-Kutta scheme object to be destroyed
  subroutine runge_kutta_scheme_coeffs_free(this)
    implicit none
    class(runge_kutta_time_scheme_t) :: this

    if (c_associated(this%coeffs_A_d)) then
       call device_free(this%coeffs_A_d)
    end if
    if (c_associated(this%coeffs_b_d)) then
       call device_free(this%coeffs_b_d)
    end if
    if (c_associated(this%coeffs_c_d)) then
       call device_free(this%coeffs_c_d)
    end if

    if (allocated(this%coeffs_A)) then
       deallocate(this%coeffs_A)
    end if
    if (allocated(this%coeffs_b)) then
       deallocate(this%coeffs_b)
    end if
    if (allocated(this%coeffs_c)) then
       deallocate(this%coeffs_c)
    end if
  end subroutine runge_kutta_scheme_coeffs_free

end module runge_kutta_time_scheme
