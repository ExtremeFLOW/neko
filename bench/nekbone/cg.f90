subroutine cg(x, f, g, c, r, w, p, z, n, msk, niter, gs_h)
  use num_types
  use gather_scatter
  use field
  use comm  
  implicit none

  type(field_t), intent(inout) :: x
  type(field_t), intent(inout) :: w
  type(field_t), intent(inout) :: msk
  real(kind=dp), intent(inout) :: g(6, x%Xh%lx, x%Xh%ly, x%Xh%lz, x%msh%nelv)
  real(kind=dp), intent(inout), dimension(n) :: f
  real(kind=dp), intent(inout), dimension(n) :: c
  real(kind=dp), intent(inout), dimension(n) :: r
  real(kind=dp), intent(inout), dimension(n) :: p
  real(kind=dp), intent(inout), dimension(n) :: z
  integer, intent(inout) :: n
  integer, intent(inout) :: niter
  type(gs_t), intent(inout) :: gs_h

  real(kind=dp) :: rnorm, rtr, rtr0, rtz2, rtz1, beta, pap, alpha, alphm, eps
  integer :: iter

  real(kind=dp) :: ur(x%Xh%lx**3), us(x%Xh%lx**3)
  real(kind=dp) :: ut(x%Xh%lx**3), wk(x%Xh%lx**3)
  
  eps = 1d-20
  pap = 0d0
  rtz1 = 1d0
  
  call rzero(x%x, n)
  call copy(r, f, n)
  call col2(r, msk%x, n)
  
  rnorm = sqrt(glsc3(r, c, r, n))
  iter = 0
  if (pe_rank .eq. 0) write(6,6) iter, rnorm
  
  do iter = 1, niter
     call solveM(z, r, n)

     rtz2 = rtz1
     rtz1 = glsc3(r, c, z, n)

     beta = rtz1 / rtz2
     if (iter .eq. 1) beta = 0d0
     call add2s1(p, z, beta, n)
     
     call ax(w, p, g, ur, us, ut, wk, n, msk, gs_h)
     pap = glsc3(w%x, c, p, n)

     alpha = rtz1/pap
     alphm = -alpha
     call add2s2(x%x, p, alpha, n)
     call add2s2(r, w%x, alphm, n)

     rtr = glsc3(r, c, r, n)
     if (iter .eq. 1) rtr0  = rtr
     rnorm = sqrt(rtr)
6    format('cg:',i4,1p4e12.4)
  end do
  
  if (pe_rank .eq. 0) write(6,6) iter,rnorm,alpha,beta,pap
  
end subroutine cg

subroutine solveM(z, r, n)
  use num_types
  use math
  implicit none
  real(kind=dp), intent(inout), dimension(n) :: z
  real(kind=dp), intent(inout), dimension(n) :: r
  integer, intent(inout) :: n

  call copy(z, r, n)
  
end subroutine solveM

subroutine ax(w, u, gxyz, ur, us, ut, wk, n, msk, gs_h)
  use num_types
  use gather_scatter
  use field
  implicit none

  type(field_t), intent(inout) :: w
  real(kind=dp), intent(inout) :: u(w%Xh%lx**3, w%msh%nelv)
  real(kind=dp), intent(inout) :: ur(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: us(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: ut(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: wk(w%Xh%lx**3)
  real(kind=dp), intent(inout) :: gxyz(6, w%Xh%lx**3, w%msh%nelv)
  real(kind=dp), dimension(w%Xh%lx, w%Xh%lx) :: D, Dt
  integer, intent(inout) :: n
  type(gs_t), intent(inout) :: gs_h
  type(field_t), intent(inout) :: msk
  integer :: e

  D = dble(w%Xh%dx)
  Dt = dble(w%Xh%dxt)
  do e = 1, w%msh%nelv
     call ax_e(w%x(1,1,1,e), u(1, e), gxyz(1, 1, e), &
          ur, us, ut, wk, w%Xh%lx, D, Dt)
  end do

  call gs_op(gs_h, w, GS_OP_ADD)
  call add2s2(w%x, u, 1d-1, n)
  call col2(w%x, msk%x, n)
  
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
