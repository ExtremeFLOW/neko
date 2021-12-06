module mathops
  use num_types
  implicit none

contains

  !> \f$ a = -a \f$
  subroutine opchsign (a1,a2,a3,gdim,n)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    integer :: i

    if (gdim .eq. 3) then
       do i=1,n
          a1(i)=-a1(i)
          a2(i)=-a2(i)
          a3(i)=-a3(i)
       end do
    else
       do i=1,n
          a1(i)=-a1(i)
          a2(i)=-a2(i)
       end do
    end if

  end subroutine opchsign
  
  !> \f$ a = a * c \f$
  subroutine opcolv (a1,a2,a3,c,gdim,n)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: c
    integer :: i

    if (gdim .eq. 3) then
       do i=1,n
          a1(i)=a1(i)*c(i)
          a2(i)=a2(i)*c(i)
          a3(i)=a3(i)*c(i)
       end do
    else
       do i=1,n
          a1(i)=a1(i)*c(i)
          a2(i)=a2(i)*c(i)
       end do
    end if

  end subroutine opcolv

  !> \f$ a(i) = b(i) * c(i) * d \f$ 
  subroutine opcolv3c(a1, a2, a3, b1, b2, b3, c, d, n, gdim)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=rp), intent(in) :: c(n), d
    integer :: i

    if (gdim .eq. 3) then
       do i=1,n
          a1(i) = b1(i)*c(i)*d
          a2(i) = b2(i)*c(i)*d
          a3(i) = b3(i)*c(i)*d
       end do
    else
       do i=1,n
          a1(i) =  b1(i)*c(i)*d
          a2(i) =  b2(i)*c(i)*d
       end do
    endif

  end subroutine opcolv3c

  !> \f$ a(i) = a + b(i) * c \f$ 
  subroutine opadd2cm(a1, a2, a3, b1, b2, b3, c, n, gdim)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=rp), intent(in) :: c
    integer :: i

    if (gdim .eq. 3) then
       do i = 1,n
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
          a3(i) = a3(i) + b3(i)*c
       end do
    else
       do i = 1, n
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
       end do
    endif

  end subroutine opadd2cm

  !> \f$ a(i) = a + b(i) * c(i) \f$
  subroutine opadd2col(a1, a2, a3, b1, b2, b3, c, n, gdim)
    integer, intent(in) :: n, gdim
    real(kind=rp), dimension(n), intent(inout) :: a1, a2, a3
    real(kind=rp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=rp), intent(in) :: c(n)
    integer :: i
    
    if (gdim .eq. 3) then
       do i=1,n
          a1(i) = a1(i) + b1(i)*c(i)
          a2(i) = a2(i) + b2(i)*c(i)
          a3(i) = a3(i) + b3(i)*c(i)
       end do
    else
       do i=1,n
          a1(i) = a1(i) + b1(i)*c(i)
          a2(i) = a2(i) + b2(i)*c(i)
       end do
    endif
    
  end subroutine opadd2col
  
end module mathops
