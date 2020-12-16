module mathops
  use num_types
  use comm
  implicit none
 
contains

  subroutine opchsign (a1,a2,a3,gdim,n)
    integer, intent(in) :: n, gdim
    real(kind=dp), dimension(n), intent(inout) :: a1, a2, a3
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
  subroutine opcolv (a1,a2,a3,c,gdim,n)
    integer, intent(in) :: n, gdim
    real(kind=dp), dimension(n), intent(inout) :: a1, a2, a3, c
    integer :: i
    if (gdim .eq. 3) then
       do i=1,n
         a1(I)=a1(I)*c(I)
         a2(I)=a2(I)*c(I)
         a3(I)=a3(I)*c(I)
       end do 
    else
       do i=1,n
         a1(I)=a1(I)*c(I)
         a2(I)=a2(I)*c(I)
       end do
    end if
  end subroutine opcolv
  subroutine opcolv3c (a1,a2,a3,b1,b2,b3,c,d,n,gdim)
    integer, intent(in) :: n, gdim
    real(kind=dp), dimension(n), intent(inout) :: a1,a2,a3
    real(kind=dp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=dp), intent(in) :: c(n), d
    integer :: i
    if (gdim .eq. 3) then
       do i=1,n
          a1(i) = b1(i)*c(i)*d
          a2(i) = b2(i)*c(i)*d
          a3(i) = b3(i)*c(i)*d
       enddo
    else
       do i=1,n
          a1(i) =  b1(i)*c(i)*d
          a2(i) =  b2(i)*c(i)*d
       enddo
    endif
  end subroutine opcolv3c
 subroutine opadd2cm (a1,a2,a3,b1,b2,b3,c,n,gdim)
    integer, intent(in) :: n, gdim
    real(kind=dp), dimension(n), intent(inout) :: a1,a2,a3
    real(kind=dp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=dp), intent(in) :: c
    integer :: i
    if (gdim .eq. 3) then
       do i=1,n
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
          a3(i) = a3(i) + b3(i)*c
       enddo
    else
       do i=1,n
          a1(i) = a1(i) + b1(i)*c
          a2(i) = a2(i) + b2(i)*c
       enddo
    endif
  end subroutine opadd2cm
 subroutine opadd2col (a1,a2,a3,b1,b2,b3,c,n,gdim)
    integer, intent(in) :: n, gdim
    real(kind=dp), dimension(n), intent(inout) :: a1,a2,a3
    real(kind=dp), dimension(n), intent(in) :: b1, b2, b3
    real(kind=dp), intent(in) :: c(n)
    integer :: i
    if (gdim .eq. 3) then
       do i=1,n
          a1(i) = a1(i) + b1(i)*c(i)
          a2(i) = a2(i) + b2(i)*c(i)
          a3(i) = a3(i) + b3(i)*c(i)
       enddo
    else
       do i=1,n
          a1(i) = a1(i) + b1(i)*c(i)
          a2(i) = a2(i) + b2(i)*c(i)
       enddo
    endif
  end subroutine opadd2col
end module mathops
