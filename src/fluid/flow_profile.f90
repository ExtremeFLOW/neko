! Copyright (c) 2021, The Neko Authors
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
!> Defines a flow profile 
module flow_profile
  use num_types
  implicit none

  !> Abstract interface for computing a Blasius flow profile
  abstract interface
     function blasius_profile(y, delta, u)
       import rp
       real(kind=rp), intent(in) :: y, delta, u
       real(kind=rp) :: blasius_profile
     end function blasius_profile
  end interface

contains

  !> Linear approximate Blasius profile
  !! \f$ \frac{u}{U} = \frac{y}{\delta} \f$
  function blasius_linear(y, delta, u)
    real(kind=rp), intent(in) :: y, delta, u
    real(kind=rp) :: blasius_linear
    real(kind=rp) :: arg

    arg = y /delta
    if (arg .gt. 1.0_rp) then
       blasius_linear = u
    else
       blasius_linear = u * (y / delta)
    end if
    
  end function blasius_linear

  !> Quadratic approximate Blasius Profile
  !! \f$ \frac{u}{U} = 2 \frac{y}{\delta} - \frac{y}{\delta}^2 \f$
  function blasius_quadratic(y, delta, u)
    real(kind=rp), intent(in) :: y, delta, u
    real(kind=rp) :: blasius_quadratic
    real(kind=rp) :: arg

    arg = ( 2.0_rp * (y / delta) - (y / delta)**2 )

    if (arg .gt. 1.0_rp) then
       blasius_quadratic = u
    else
       blasius_quadratic = u * arg
    end if
    
  end function blasius_quadratic

  !> Cubic approximate Blasius Profile
  !! \f$ \frac{u}{U} = 3/2 \frac{y}{\delta} - 1/2\frac{y}{\delta}^3 \f$
  function blasius_cubic(y, delta, u)
    real(kind=rp), intent(in) :: y, delta, u
    real(kind=rp) :: blasius_cubic
    real(kind=rp) :: arg

    arg = ( 3.0_rp / 2.0_rp * (y / delta) - 0.5_rp * (y / delta)**3 )

    if (arg .gt. 1.0_rp) then
       blasius_cubic = 1.0_rp
    else
       blasius_cubic = u * arg
    end if
    
  end function blasius_cubic

  !> Quartic approximate Blasius Profile
  !! \f$ \frac{u}{U} = 2 \frac{y}{\delta} - 2\frac{y}{\delta}^3 +
  !! \frac{y}{\delta}^4 \f$
  function blasius_quartic(y, delta, u)
    real(kind=rp), intent(in) :: y, delta, u
    real(kind=rp) :: blasius_quartic
    real(kind=rp) :: arg

    arg = 2.0_rp * (y / delta) - 2.0_rp * (y / delta)**3 + (y / delta)**4

    if (arg .gt. 1.0_rp) then
       blasius_quartic = u
    else
       blasius_quartic = u * arg
    end if
    
  end function blasius_quartic

  !> Sinusoidal approximate Blasius Profile
  !! \f$ \frac{u}{U} = \sin(\frac{\pi}{2}\frac{y}{\delta}) \f$
  function blasius_sin(y, delta, u)
    real(kind=rp), intent(in) :: y, delta, u
    real(kind=rp) :: blasius_sin
    real(kind=rp), parameter :: PI = 4.0_rp * atan(1.0_rp)
    real(kind=rp) :: arg

    arg = (PI / 2.0_rp) * (y/delta)

    if (arg .gt. 0.5_rp * PI) then
       blasius_sin = u
    else
       blasius_sin = u * sin(arg)
    end if
    
  end function blasius_sin
    
end module flow_profile
