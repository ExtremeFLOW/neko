module ax_helm
  use ax_product
  implicit none

  type, public, extends(ax_t) :: ax_helm_t
   contains
     procedure, nopass :: compute => ax_helm_compute
  end type ax_helm_t

contains 

  subroutine ax_helm_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    
    
    select case(Xh%lx)
    case (12)
       call ax_helm_lx12(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case (11)
       call ax_helm_lx11(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(10)
       call ax_helm_lx10(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(9)
       call ax_helm_lx9(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(8)
       call ax_helm_lx8(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(7)
       call ax_helm_lx7(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(6)
       call ax_helm_lx6(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(5)
       call ax_helm_lx5(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(4)
       call ax_helm_lx4(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(3)
       call ax_helm_lx3(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    case(2)
       call ax_helm_lx2(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    end select
    
    if (coef%ifh2) call addcol4 (w,coef%h2,coef%B,u,coef%dof%n_dofs)
    
 
  end subroutine ax_helm_compute
  
  subroutine ax_helm_lx12(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e) &
                        + Dx(i,11) * u(11,j,1,e) &
                        + Dx(i,12) * u(12,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e) &
                           + Dy(j,11) * u(i,11,k,e) &
                           + Dy(j,12) * u(i,12,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e) &
                        + Dz(k,11) * u(i,1,11,e) &
                        + Dz(k,12) * u(i,1,12,e)
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) &
                        + Dxt(i,6) * ur(6,j,1) &
                        + Dxt(i,7) * ur(7,j,1) &
                        + Dxt(i,8) * ur(8,j,1) &
                        + Dxt(i,9) * ur(9,j,1) &
                        + Dxt(i,10) * ur(10,j,1) &
                        + Dxt(i,11) * ur(11,j,1) &
                        + Dxt(i,12) * ur(12,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) &
                           + Dyt(j,6) * us(i,6,k) &
                           + Dyt(j,7) * us(i,7,k) &
                           + Dyt(j,8) * us(i,8,k) &
                           + Dyt(j,9) * us(i,9,k) &
                           + Dyt(j,10) * us(i,10,k) &
                           + Dyt(j,11) * us(i,11,k) &
                           + Dyt(j,12) * us(i,12,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) &
                        + Dzt(k,6) * ut(i,1,6) &
                        + Dzt(k,7) * ut(i,1,7) &
                        + Dzt(k,8) * ut(i,1,8) &
                        + Dzt(k,9) * ut(i,1,9) &
                        + Dzt(k,10) * ut(i,1,10) &
                        + Dzt(k,11) * ut(i,1,11) &
                        + Dzt(k,12) * ut(i,1,12) 
          end do
       end do

    end do
  end subroutine ax_helm_lx12

  subroutine ax_helm_lx11(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e) &
                        + Dx(i,11) * u(11,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e) &
                           + Dy(j,11) * u(i,11,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e) &
                        + Dz(k,11) * u(i,1,11,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) &
                        + Dxt(i,6) * ur(6,j,1) &
                        + Dxt(i,7) * ur(7,j,1) &
                        + Dxt(i,8) * ur(8,j,1) &
                        + Dxt(i,9) * ur(9,j,1) &
                        + Dxt(i,10) * ur(10,j,1) &
                        + Dxt(i,11) * ur(11,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) &
                           + Dyt(j,6) * us(i,6,k) &
                           + Dyt(j,7) * us(i,7,k) &
                           + Dyt(j,8) * us(i,8,k) &
                           + Dyt(j,9) * us(i,9,k) &
                           + Dyt(j,10) * us(i,10,k) &
                           + Dyt(j,11) * us(i,11,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) &
                        + Dzt(k,6) * ut(i,1,6) &
                        + Dzt(k,7) * ut(i,1,7) &
                        + Dzt(k,8) * ut(i,1,8) &
                        + Dzt(k,9) * ut(i,1,9) &
                        + Dzt(k,10) * ut(i,1,10) &
                        + Dzt(k,11) * ut(i,1,11) 
          end do
       end do

    end do
  end subroutine ax_helm_lx11

  subroutine ax_helm_lx10(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) &
                        + Dxt(i,6) * ur(6,j,1) &
                        + Dxt(i,7) * ur(7,j,1) &
                        + Dxt(i,8) * ur(8,j,1) &
                        + Dxt(i,9) * ur(9,j,1) &
                        + Dxt(i,10) * ur(10,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) &
                           + Dyt(j,6) * us(i,6,k) &
                           + Dyt(j,7) * us(i,7,k) &
                           + Dyt(j,8) * us(i,8,k) &
                           + Dyt(j,9) * us(i,9,k) &
                           + Dyt(j,10) * us(i,10,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) &
                        + Dzt(k,6) * ut(i,1,6) &
                        + Dzt(k,7) * ut(i,1,7) &
                        + Dzt(k,8) * ut(i,1,8) &
                        + Dzt(k,9) * ut(i,1,9) &
                        + Dzt(k,10) * ut(i,1,10) 
          end do
       end do

    end do
  end subroutine ax_helm_lx10

  subroutine ax_helm_lx9(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) &
                        + Dxt(i,6) * ur(6,j,1) &
                        + Dxt(i,7) * ur(7,j,1) &
                        + Dxt(i,8) * ur(8,j,1) &
                        + Dxt(i,9) * ur(9,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) &
                           + Dyt(j,6) * us(i,6,k) &
                           + Dyt(j,7) * us(i,7,k) &
                           + Dyt(j,8) * us(i,8,k) &
                           + Dyt(j,9) * us(i,9,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) &
                        + Dzt(k,6) * ut(i,1,6) &
                        + Dzt(k,7) * ut(i,1,7) &
                        + Dzt(k,8) * ut(i,1,8) &
                        + Dzt(k,9) * ut(i,1,9) 
          end do
       end do

    end do
  end subroutine ax_helm_lx9

  subroutine ax_helm_lx8(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) &
                        + Dxt(i,6) * ur(6,j,1) &
                        + Dxt(i,7) * ur(7,j,1) &
                        + Dxt(i,8) * ur(8,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) &
                           + Dyt(j,6) * us(i,6,k) &
                           + Dyt(j,7) * us(i,7,k) &
                           + Dyt(j,8) * us(i,8,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) &
                        + Dzt(k,6) * ut(i,1,6) &
                        + Dzt(k,7) * ut(i,1,7) &
                        + Dzt(k,8) * ut(i,1,8) 
          end do
       end do

    end do
  end subroutine ax_helm_lx8

  subroutine ax_helm_lx7(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) &
                        + Dxt(i,6) * ur(6,j,1) &
                        + Dxt(i,7) * ur(7,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) &
                           + Dyt(j,6) * us(i,6,k) &
                           + Dyt(j,7) * us(i,7,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) &
                        + Dzt(k,6) * ut(i,1,6) &
                        + Dzt(k,7) * ut(i,1,7) 
          end do
       end do

    end do
  end subroutine ax_helm_lx7

  subroutine ax_helm_lx6(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) &
                        + Dxt(i,6) * ur(6,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) &
                           + Dyt(j,6) * us(i,6,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) &
                        + Dzt(k,6) * ut(i,1,6) 
          end do
       end do

    end do
  end subroutine ax_helm_lx6

  subroutine ax_helm_lx5(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) &
                        + Dxt(i,5) * ur(5,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) &
                           + Dyt(j,5) * us(i,5,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) &
                        + Dzt(k,5) * ut(i,1,5) 
          end do
       end do

    end do
  end subroutine ax_helm_lx5

  subroutine ax_helm_lx4(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) &
                        + Dxt(i,4) * ur(4,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) &
                           + Dyt(j,4) * us(i,4,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) &
                        + Dzt(k,4) * ut(i,1,4) 
          end do
       end do

    end do
  end subroutine ax_helm_lx4

  subroutine ax_helm_lx3(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) &
                        + Dxt(i,3) * ur(3,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) &
                           + Dyt(j,3) * us(i,3,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) &
                        + Dzt(k,3) * ut(i,1,3) 
          end do
       end do

    end do
  end subroutine ax_helm_lx3

  subroutine ax_helm_lx2(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) 
          end do
       end do

       do i = 1, lx*lx*lx          
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             w(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                        + Dxt(i,2) * ur(2,j,1) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                w(i,j,k,e) = w(i,j,k,e) &
                           + Dyt(j,1) * us(i,1,k) &
                           + Dyt(j,2) * us(i,2,k) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             w(i,1,k,e) = w(i,1,k,e) &
                        + Dzt(k,1) * ut(i,1,1) &
                        + Dzt(k,2) * ut(i,1,2) 
          end do
       end do

    end do
  end subroutine ax_helm_lx2

end module ax_helm
