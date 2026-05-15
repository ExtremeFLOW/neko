! Copyright (c) 2026, The Neko Authors
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
! Implemented as shown in:
! Brandon E. Merrill, Yulia T. Peet, Paul F. Fischer, & James W. Lottes (2016).
! A spectrally accurate method for overlapping grid solution of incompressible
! Navier–Stokes equations. Journal of Computational Physics, 307, 60-93.
module iextm_time_scheme
  use num_types, only : rp
  use time_scheme, only : time_scheme_t
  use math, only : rzero
  use utils, only : neko_error
  implicit none
  private

  !> Explicit interface extrapolation scheme for overset grids.
  type, public, extends(time_scheme_t) :: iextm_time_scheme_t
   contains
     !> Compute the scheme coefficients
     procedure, nopass :: compute_coeffs => iextm_time_scheme_compute_coeffs
  end type iextm_time_scheme_t

contains

  !> Compute the scheme coefficients
  !! @param dt Timestep values, first element is the current timestep.
  !! @param order Order the scheme.
  !! @note dt(1) is the current timestep, dt(2) is the previous timestep, etc.
  subroutine iextm_time_scheme_compute_coeffs(coeffs, dt, order)
    implicit none
    real(kind=rp), intent(out) :: coeffs(4)
    real(kind=rp), intent(in) :: dt(10)
    integer, intent(in) :: order
    real(kind=rp) :: x(0:3), xtest
    integer :: i, j
    real(kind=rp) :: basis

    call rzero(coeffs, 4)

    ! This is the direct way, working for constant tstep
    !select case (order)
    !case (1)
    !   coeffs(1) = 1.0_rp
    !case (2)
    !   coeffs(2) = -1.0_rp
    !   coeffs(1) = 2.0_rp
    !case (3)
    !   coeffs(3) = 1.0_rp
    !   coeffs(2) = -3.0_rp
    !   coeffs(1) = 3.0_rp
    !case default
    !   call neko_error("The order of the IEXTm time scheme must be 1 to 3.")
    !end select

    if (order .lt. 1 .or. order .gt. 3) then
       call neko_error("The order of the IEXTm time scheme must be 1 to 3.")
    end if

    ! Create the nodes to build the interpolation basis
    ! This simply places old points as nodes in a 1d grid
    x(0) = dt(1) ! Current timestep
    do i = 1, order
       x(i) = x(i-1) - dt(i)
    end do
    ! Define the test point
    xtest = x(0)

    !> Compute the Lagrange interpolation coefficients for xtest
    do i = 1, order
       basis = 1.0_rp
       do j = 1, order
          if (j .ne. i) then
             basis = basis * (xtest - x(j)) / (x(i) - x(j))
          end if
       end do
       coeffs(i) = basis
    end do

  end subroutine iextm_time_scheme_compute_coeffs

end module iextm_time_scheme
