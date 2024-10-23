! Setup fake data
subroutine set_data(u, v, n)
  use gather_scatter
  use num_types
  implicit none
  
  real(kind=dp), intent(inout), dimension(n) :: u
  real(kind=dp), intent(inout), dimension(n) :: v
  integer,  intent(inout) :: n
  real(kind=dp) :: arg
  integer :: i

  do i = 1, n
     arg = (i * i)
     arg =  cos(arg)
     u(i) = sin(arg)
     v(i) = sin(arg)
  end do

end subroutine set_data

! Setup geom terms
subroutine setup_g(g, w, lx, ly, lz, n)
  use num_types
  implicit none
  
  real(kind=dp), intent(inout), dimension(6, lx, ly, lz, n) :: g
  real(kind=dp), intent(inout), dimension(lx) :: w
  integer, intent(in) :: lx, ly, lz, n
  integer :: i, j, k, l

  g = 0d0
  
  do i = 1, n
     do l = 1, lz
        do k = 1, ly
           do j = 1, lx
              g(1, j, k, l, i) = w(j) * w(k) * w(l)
              g(4, j, k, l, i) = w(j) * w(k) * w(l)
              g(6, j, k, l, i) = w(j) * w(k) * w(l)
           end do
        end do
     end do
  end do

end subroutine setup_g
  
