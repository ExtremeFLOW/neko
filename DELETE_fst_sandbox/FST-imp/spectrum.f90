module spectrum
    use num_types, only: rp
    implicit none
    
contains 

  !> Computes the Von Karman spectrum
  !! @param k Wave number.
  !! @param L Integral length scale.
  !! @param q Turbulence intensity.
  function ek(k,L,q) result(E)
    real(kind=rp), intent(in) :: k, L, q
    real(kind=rp) :: E

    E = 2._rp/3._rp*q*1.606_rp * (k*L)**4._rp * L / &
      (1.350_rp+(k*L)**2._rp)**(17._rp/6._rp)
  end function ek

end module spectrum
