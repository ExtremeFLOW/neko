subroutine ax(w, u, gxyz, ur, us, ut, wk, n)
  use num_types
  use field
  implicit none

  type(field_t), intent(inout) :: w
  real(kind=dp), intent(inout) :: u(w%Xh%lx**3, w%msh%nelv)
  real(kind=dp), intent(inout) :: ur(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: us(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: ut(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: wk(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: gxyz(6, w%Xh%lx**3, w%msh%nelv)
  integer, intent(inout) :: n
  integer :: e

  do e = 1, w%msh%nelv
     call ax_e(w%x(1,1,1,e), u(1,e), gxyz(1,1,e), &
          ur, us, ut, wk, w%Xh%lx, w%Xh%dx, w%Xh%dxt)
  end do
  
end subroutine ax

subroutine ax_e(w, u, g, ur, us, ut, wk, lx, D, Dt)
  use num_types
  implicit none
  
  real(kind=dp), intent(inout) :: w(lx**3)
  real(kind=dp), intent(inout) :: u(lx**3)
  real(kind=dp), intent(inout) :: ur(lx**3)
  real(kind=dp), intent(inout) :: us(lx**3)
  real(kind=dp), intent(inout) :: ut(lx**3)
  real(kind=dp), intent(inout) :: wk(lx**3)
  real(kind=dp), intent(inout) :: g(6, lx**3)
  real(kind=dp), intent(inout) :: D(lx, lx)
  real(kind=dp), intent(inout) :: Dt(lx, lx)
  integer, intent(inout) :: lx
  real(kind=dp) :: wr, ws, wt
  integer :: i, n

  n = lx - 1
  call local_grad3(ur, us, ut, u, n, D, Dt)

  do i=1, lx**3
     wr = g(1,i)*ur(i) + g(2,i)*us(i) + g(3,i)*ut(i)
     ws = g(2,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
     wt = g(3,i)*ur(i) + g(5,i)*us(i) + g(6,i)*ut(i)
     ur(i) = wr
     us(i) = ws
     ut(i) = wt
  enddo

  call local_grad3_t(w, ur, us, ut, n, D, Dt, wk)

end subroutine ax_e

subroutine local_grad3(ur, us, ut, u, n, D, Dt)
  use num_types
  use mxm_wrapper
  implicit none
  
  real(kind=dp), intent(inout) :: ur(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: us(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: ut(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: u(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: D(0:n, 0:n)
  real(kind=dp), intent(inout) :: Dt(0:n, 0:n)
  integer, intent(inout) :: n
  integer :: m1, m2, k

  m1 = n + 1
  m2 = m1*m1

  call mxm(D ,m1,u,m1,ur,m2)
  do k=0,n
     call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
  enddo
  call mxm(u,m2,Dt,m1,ut,m1)
  
end subroutine local_grad3

subroutine local_grad3_t(u,ur,us,ut,N,D,Dt,w)
  use num_types
  use mxm_wrapper
  use math
  implicit none

  real(kind=dp), intent(inout) :: ur(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: us(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: ut(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: u(0:n, 0:n, 0:n)
  real(kind=dp), intent(inout) :: D(0:n, 0:n)
  real(kind=dp), intent(inout) :: Dt(0:n, 0:n)
  real(kind=dp), intent(inout) :: w(0:n, 0:n, 0:n)
  integer, intent(inout) :: n
  integer :: m1, m2, m3, k

  m1 = n + 1
  m2 = m1**2
  m3 = m1**3
  
  call mxm(Dt,m1,ur,m1,u,m2)
  
  do k=0,N
     call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
  enddo
  call add2(u,w,m3)
  
  call mxm(ut,m2,D ,m1,w,m1)
  call add2(u,w,m3)
  
end subroutine local_grad3_t
