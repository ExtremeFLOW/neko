
! Copyright (c) 2023, The Neko Authors
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
!> Explicit extrapolation scheme for time integration.
module ext_time_scheme
  use neko_config
  use num_types, only : rp
  use time_scheme, only: time_scheme_t
  use math, only : rzero
  use utils, only : neko_warning
  use device, only : HOST_TO_DEVICE, device_memcpy, device_free
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  !> Explicit extrapolation scheme for time integration.
  !! @details
  !! Compute the value at the current time-step by evaluating a polynomial
  !! built using the values on previous time-steps.
  !!
  !! For a constant time-step corresponds to the following schemes of order
  !! 1 to 3:
  !! - Order 1: \f$  u^{n+1} = u^n \f$
  !! - Order 2: \f$  u^{n+1} = 2u^n - u^{n-1} \f$, linear extrapolation
  !! - Order 3: \f$  u^{n+1} = 3u^n - 3u^{n-1} + u^{n-2} \f$
  type, public, extends(time_scheme_t) :: ext_time_scheme_t 
     contains
       !> Compute the scheme coefficients
       procedure, pass(this) :: set_coeffs => ext_time_scheme_set_coeffs
  end type ext_time_scheme_t

  contains 

  !> Compute the scheme coefficients
  !! @param t Timestep values, first element is the current timestep.
  subroutine ext_time_scheme_set_coeffs(this, dt)
    class(ext_time_scheme_t), intent(inout)  :: this
    real(kind=rp), intent(inout), dimension(10) :: dt
    real(kind=rp), dimension(4) :: coeffs_old
    associate(n => this%n, coeffs => this%coeffs, coeffs_d => this%coeffs_d)
      
      ! To check whether the coefficients changed
      coeffs_old = coeffs
      
      ! Increment the order of the scheme if below time_order
      n = n + 1
      n = min(n, this%time_order)
      
      call rzero(coeffs, 4)
      
      if (n .eq. 1) then
         coeffs(1) = 1.0_rp
      else if (n .eq. 2) then
         coeffs(2) = -dt(1) / dt(2)
         coeffs(1) =  1.0_rp - coeffs(2)
      else if (n .eq. 3) then
         coeffs(3) =  dt(1) / (dt(2) + dt(3)) * (dt(1) + dt(2)) / dt(3)
         coeffs(2) = - dt(1) / dt(2) * (1.0_rp + dt(2) / dt(3) + dt(1) / dt(3))
         coeffs(1) =  1.0_rp - coeffs(2) - coeffs(3)
      endif

      if (c_associated(coeffs_d)) then
         if (maxval(abs(coeffs - coeffs_old)) .gt. 1e-10_rp) then
            call device_memcpy(coeffs, coeffs_d, 4, HOST_TO_DEVICE)
         end if
      end if
    end associate
    
  end subroutine ext_time_scheme_set_coeffs
end module ext_time_scheme