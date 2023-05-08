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
!> Backward-differencing scheme for time integration.
module bdf_time_scheme
  use neko_config
  use num_types, only : rp
  use time_scheme, only: time_scheme_t
  use math, only : rzero
  use utils, only : neko_warning
  use device, only : HOST_TO_DEVICE, device_memcpy, device_free
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  !> Implicit backward-differencing scheme for time integration.
  !! @details
  !! The explicit forumlas for the coefficients are taken from the following
  !! techincal note:
  !! "Derivation of BDF2/BDF3 for Variable Step Size" by Hiroaki Nishikawa,
  !! which can be found on ResearchGate.
  !!
  !! For a contant time-step this corresponds to the following schemes for
  !! order 1 to 3:
  !! - Order 1: \f$ \frac{1}{\Delta t} u^{n+1} - \frac{1}{\Delta t} u^{n} \f$ 
  !! - Order 2: \f$ \frac{3}{2\Delta t} u^{n+1} - \frac{4}{2\Delta t} u^{n}  + 
  !!                \frac{1}{2\Delta t} u^{n-1}\f$ 
  !! - Order 3: \f$ \frac{11}{6\Delta t} u^{n+1} - \frac{18}{6\Delta t} u^{n}  + 
  !!                \frac{9}{6\Delta t} u^{n-1} - \frac{2}{6\Delta t} u^{n-2}\f$ 
  !!
  !! It is assumed that all the coefficients but the first one premultiply terms
  !! that go to the right-hand side of the equation. 
  !! Accordingly, the signs of these coefficients are reversed in the `coeffs`
  !! array. This is taken into account, for example, in the implemeation of the
  !! `rhs_maker` class.
  !!
  !! Another important convention is that the coefficients are meant to be
  !! later divided by the current value of the timestep.
  !!
  !! In line with the above assumptions, the first order scheme always returns
  !! the array \f$[1, 1]\f$, and **not** \f$[1/\Delta t, -1/\Delta t]\f$, as one
  !! might expect. Similar for the second and third order.
  !!
  !! @remark 
  !! The current implementation can be easily extended to schemes of arbitrary
  !! order, by using the `fd_weights_full` subroutine to compute the
  !! coefficients. A demonstration of this is implemented in a test in
  !! `tests/ext_bdf_scheme/test_bdf.pf`
  type, public, extends(time_scheme_t) :: bdf_time_scheme_t 
     contains
       !> Constructor
    !    procedure, pass(this) :: init => bdf_time_scheme_init
       !> Destructor
    !    procedure, pass(this) :: free => bdf_time_scheme_free
       !> Compute the scheme coefficients
       procedure, pass(this) :: set_coeffs => bdf_time_scheme_set_coeffs
  end type bdf_time_scheme_t
  
  contains

!   !> Constructor
!   !! @param torder Desired order of the scheme: 1, 2 or 3.
!   subroutine bdf_time_scheme_init(this, torder)
!     class(bdf_time_scheme_t), intent(inout) :: this
!     integer, intent(in) :: torder
!     call this%time_scheme_t%init(torder)
     
!   end subroutine bdf_time_scheme_init

  !> Compute the scheme coefficients
  !! @param t Timestep values, first element is the current timestep.
  subroutine bdf_time_scheme_set_coeffs(this, dt)
    implicit none
    class(bdf_time_scheme_t), intent(inout) :: this
    real(kind=rp), intent(inout), dimension(10) :: dt
    real(kind=rp), dimension(4) :: coeffs_old

    associate(n => this%n, coeffs => this%coeffs, coeffs_d => this%coeffs_d)
      ! will be used to check whether the coefficients changed
      coeffs_old = coeffs
      
      ! Increment the order of the scheme if below time_order
      n = n + 1
      n = min(n, this%time_order)

      call rzero(coeffs, 4)
      
      ! Note, these are true coeffs, multiplied by dt(1)
      select case (n)
      case (1)
         coeffs(1) = 1.0_rp
         coeffs(2) = 1.0_rp
      case (2)
         coeffs(1) = (1 + dt(1)/(dt(1) + dt(2)))
         coeffs(3) = -dt(1)**2/dt(2)/(dt(1) + dt(2))
         coeffs(2) =  coeffs(1) - coeffs(3)
      case (3)
         coeffs(2) = (dt(1) + dt(2)) * (dt(1) + dt(2) + dt(3)) / &
                     (dt(1) * dt(2) * (dt(2) + dt(3)))
         coeffs(3) = -dt(1) * (dt(1) + dt(2) + dt(3)) / &
                     (dt(2) * dt(3) * (dt(1) + dt(2)))
         coeffs(4) = dt(1) * (dt(1) + dt(2)) / &
                     (dt(3) * (dt(2) + dt(3)) * (dt(1) + dt(2) + dt(3)))
         coeffs(1) = coeffs(2) + coeffs(3) + coeffs(4)
         coeffs = coeffs * dt(1)
      end select
      
      if (c_associated(coeffs_d)) then
         if (maxval(abs(coeffs - coeffs_old)) .gt. 1e-10_rp) then
            call device_memcpy(coeffs, coeffs_d, 4, HOST_TO_DEVICE)
         end if
      end if
    end associate
    
  end subroutine bdf_time_scheme_set_coeffs

end module bdf_time_scheme