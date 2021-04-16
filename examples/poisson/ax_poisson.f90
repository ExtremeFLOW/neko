module ax_poisson
  use ax_product
  implicit none

  type, public, extends(ax_t) :: ax_poisson_t
  contains
     procedure, nopass :: compute => ax_poisson_compute
  end type ax_poisson_t

contains 
  
  subroutine ax_poisson_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
  
    real(kind=rp) :: ur(Xh%lx**3)
    real(kind=rp) :: us(Xh%lx**3)
    real(kind=rp) :: ut(Xh%lx**3)
    real(kind=rp) :: wk(Xh%lx**3)
    integer :: e
  
    do e = 1, msh%nelv
       call ax_e(w(1,1,1,e), u(1,1,1,e), &
            coef%g1(1,1,1,e), coef%g2(1,1,1,e), coef%g3(1,1,1,e), &
            coef%g4(1,1,1,e), coef%g5(1,1,1,e), coef%g6(1,1,1,e), &
            ur, us, ut, wk, Xh%lx, Xh%dx, Xh%dxt)
    end do

  end subroutine ax_poisson_compute
  
  subroutine ax_e(w, u, g1, g2, g3, g4, g5, g6, ur, us, ut, wk, lx, D, Dt)
    real(kind=rp), intent(inout) :: w(lx**3)
    real(kind=rp), intent(inout) :: u(lx**3)
    real(kind=rp), intent(inout) :: g1(lx**3)
    real(kind=rp), intent(inout) :: g2(lx**3)
    real(kind=rp), intent(inout) :: g3(lx**3)
    real(kind=rp), intent(inout) :: g4(lx**3)
    real(kind=rp), intent(inout) :: g5(lx**3)
    real(kind=rp), intent(inout) :: g6(lx**3)
    real(kind=rp), intent(inout) :: ur(lx**3)
    real(kind=rp), intent(inout) :: us(lx**3)
    real(kind=rp), intent(inout) :: ut(lx**3)
    real(kind=rp), intent(inout) :: wk(lx**3)
    real(kind=rp), intent(inout) :: D(lx, lx)
    real(kind=rp), intent(inout) :: Dt(lx, lx)
    integer, intent(inout) :: lx
    real(kind=rp) :: wr, ws, wt
    integer :: i, n
  
    n = lx - 1
    call local_grad3(ur, us, ut, u, n, D, Dt)
  
    do i=1, lx**3
       wr = g1(i)*ur(i) 
       ws = g2(i)*us(i) 
       wt = g3(i)*ut(i) 
       ur(i) = wr
       us(i) = ws
       ut(i) = wt
    enddo
  
    call local_grad3_t(w, ur, us, ut, n, D, Dt, wk)
  
  end subroutine ax_e
  
  subroutine local_grad3(ur, us, ut, u, n, D, Dt)
    
    real(kind=rp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: D(0:n, 0:n)
    real(kind=rp), intent(inout) :: Dt(0:n, 0:n)
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
  
    real(kind=rp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=rp), intent(inout) :: D(0:n, 0:n)
    real(kind=rp), intent(inout) :: Dt(0:n, 0:n)
    real(kind=rp), intent(inout) :: w(0:n, 0:n, 0:n)
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
   
end module ax_poisson
