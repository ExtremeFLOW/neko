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
  !! \f$ \fract{u}{U} = 2 \frac{y}{\delta} - \frac{y}{delta}^2 \f$
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
  !! \f$ \fract{u}{U} = 3/2 \frac{y}{\delta} - 1/2\frac{y}{delta}^3 \f$
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
  !! \f$ \fract{u}{U} = 2 \frac{y}{\delta} - 2\frac{y}{delta}^3 +
  !! frac{y}{delta}^4 \f$
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
  !! \f$ \frac{u}{U} = \sin(\frac{\pi}{2}\frac{y}{\delta})
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
