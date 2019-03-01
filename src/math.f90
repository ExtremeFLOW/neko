module math
  use num_types
  implicit none

  real(kind=dp), parameter :: NEKO_EPS = epsilon(1d0)

contains

  !> Return absolute comparison \f$ | x - y | < \epsilon \f$
  pure function abscmp(x, y)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: y
    logical, intent(out) ::  abscmp 

    abscmp = abs(x - y) .lt. NEKO_EPS

  end function abscmp
  
end module math
