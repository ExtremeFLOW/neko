module ax_helm_sx
  use ax_product
  use num_types
  implicit none
  private

  type, public, extends(ax_t) :: ax_helm_sx_t
   contains
     procedure, nopass :: compute => ax_helm_sx_compute
  end type ax_helm_sx_t

contains 

  subroutine ax_helm_sx_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    if (Xh%lx .eq. 12) then
       call sx_ax_helm_lx12(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    else if (Xh%lx .eq. 10) then
       call sx_ax_helm_lx10(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    else if (Xh%lx .eq. 8) then
       call sx_ax_helm_lx8(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    else if (Xh%lx .eq. 4) then
       call sx_ax_helm_lx4(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    else if (Xh%lx .eq. 2) then
       call sx_ax_helm_lx2(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, msh%nelv)
    end if

    if (coef%ifh2) call addcol4 (w,coef%h2,coef%B,u,coef%dof%n_dofs)

  end subroutine ax_helm_sx_compute

  subroutine sx_ax_helm_lx12(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
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

    integer :: e, i, j, k, jj, kk    
    real(kind=rp) :: ur(lx,lx,lx,n)
    real(kind=rp) :: us(lx,lx,lx,n)
    real(kind=rp) :: ut(lx,lx,lx,n)
    real(kind=rp) :: uur(lx,lx,lx,n)
    real(kind=rp) :: uus(lx,lx,lx,n)
    real(kind=rp) :: uut(lx,lx,lx,n)
    real(kind=rp) :: wr, ws, wt

    do i=1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk=1,lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k=1,lx
       do i=1,lx                
          do j=1,lx       
             do e = 1,n                                   
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, n * lx * lx * lx
       uur(i,1,1,1) = h1(i,1,1,1) * &
            ( G11(i,1,1,1) * ur(i,1,1,1) &
            + G12(i,1,1,1) * us(i,1,1,1) &
            + G13(i,1,1,1) * ut(i,1,1,1))

       uus(i,1,1,1) = h1(i,1,1,1) * &
            ( G22(i,1,1,1) * us(i,1,1,1) &
            + G12(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * ut(i,1,1,1) )

       uut(i,1,1,1) = h1(i,1,1,1) * &
            ( G33(i,1,1,1) * ut(i,1,1,1) &
            + G13(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * us(i,1,1,1))
    end do

    do i=1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk=1,lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,n
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt 
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_lx12

  subroutine sx_ax_helm_lx10(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
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

    integer :: e, i, j, k, jj, kk    
    real(kind=rp) :: ur(lx,lx,lx,n)
    real(kind=rp) :: us(lx,lx,lx,n)
    real(kind=rp) :: ut(lx,lx,lx,n)
    real(kind=rp) :: uur(lx,lx,lx,n)
    real(kind=rp) :: uus(lx,lx,lx,n)
    real(kind=rp) :: uut(lx,lx,lx,n)
    real(kind=rp) :: wr, ws, wt

    do i=1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk=1,lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k=1,lx
       do i=1,lx                
          do j=1,lx       
             do e = 1,n                                   
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, n * lx * lx * lx
       uur(i,1,1,1) = h1(i,1,1,1) * &
            ( G11(i,1,1,1) * ur(i,1,1,1) &
            + G12(i,1,1,1) * us(i,1,1,1) &
            + G13(i,1,1,1) * ut(i,1,1,1))

       uus(i,1,1,1) = h1(i,1,1,1) * &
            ( G22(i,1,1,1) * us(i,1,1,1) &
            + G12(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * ut(i,1,1,1) )

       uut(i,1,1,1) = h1(i,1,1,1) * &
            ( G33(i,1,1,1) * ut(i,1,1,1) &
            + G13(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * us(i,1,1,1))
    end do

    do i=1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk=1,lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,n
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt 
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_lx10
  
  subroutine sx_ax_helm_lx8(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
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

    integer :: e, i, j, k, jj, kk    
    real(kind=rp) :: ur(lx,lx,lx,n)
    real(kind=rp) :: us(lx,lx,lx,n)
    real(kind=rp) :: ut(lx,lx,lx,n)
    real(kind=rp) :: uur(lx,lx,lx,n)
    real(kind=rp) :: uus(lx,lx,lx,n)
    real(kind=rp) :: uut(lx,lx,lx,n)
    real(kind=rp) :: wr, ws, wt

    do i=1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk=1,lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k=1,lx
       do i=1,lx                
          do j=1,lx       
             do e = 1,n                                   
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, n * lx * lx * lx
       uur(i,1,1,1) = h1(i,1,1,1) * &
            ( G11(i,1,1,1) * ur(i,1,1,1) &
            + G12(i,1,1,1) * us(i,1,1,1) &
            + G13(i,1,1,1) * ut(i,1,1,1))

       uus(i,1,1,1) = h1(i,1,1,1) * &
            ( G22(i,1,1,1) * us(i,1,1,1) &
            + G12(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * ut(i,1,1,1) )

       uut(i,1,1,1) = h1(i,1,1,1) * &
            ( G33(i,1,1,1) * ut(i,1,1,1) &
            + G13(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * us(i,1,1,1))
    end do

    do i=1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk=1,lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,n
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt 
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_lx8

  subroutine sx_ax_helm_lx4(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
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

    integer :: e, i, j, k, jj, kk
    real(kind=rp) :: ur(lx,lx,lx,n)
    real(kind=rp) :: us(lx,lx,lx,n)
    real(kind=rp) :: ut(lx,lx,lx,n)
    real(kind=rp) :: uur(lx,lx,lx,n)
    real(kind=rp) :: uus(lx,lx,lx,n)
    real(kind=rp) :: uut(lx,lx,lx,n)

    do i=1,lx
       do jj = 1, lx * lx * n
          ur(i,jj,1,1) = Dx(i,1)*u(1,jj,1,1) &
                       + Dx(i,2)*u(2,jj,1,1) &
                       + Dx(i,3)*u(3,jj,1,1) &
                       + Dx(i,4)*u(4,jj,1,1) 
       end do
    end do


    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1, n
                us(i,j,k,e) = Dy(j,1) * u(i,1,k,e) &
                            + Dy(j,2) * u(i,2,k,e) &
                            + Dy(j,3) * u(i,3,k,e) &
                            + Dy(j,4) * u(i,4,k,e) 
             end do
          end do
       end do
    end do



    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1, n
                ut(i,j,k,e) = Dz(k,1) * u(i,j,1,e) &
                            + Dz(k,2) * u(i,j,2,e) &
                            + Dz(k,3) * u(i,j,3,e) &
                            + Dz(k,4) * u(i,j,4,e) 
             end do
          end do
       end do
    end do

    do i = 1, n * lx * lx * lx
       uur(i,1,1,1) = h1(i,1,1,1) * &
            ( G11(i,1,1,1) * ur(i,1,1,1) &
            + G12(i,1,1,1) * us(i,1,1,1) &
            + G13(i,1,1,1) * ut(i,1,1,1))

       uus(i,1,1,1) = h1(i,1,1,1) * &
            ( G22(i,1,1,1) * us(i,1,1,1) &
            + G12(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * ut(i,1,1,1) )

       uut(i,1,1,1) = h1(i,1,1,1) * &
            ( G33(i,1,1,1) * ut(i,1,1,1) &
            + G13(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * us(i,1,1,1))
    end do

    do i=1,lx
       do jj = 1, lx * lx * n
          w(i,jj,1,1) = Dxt(i,1) * uur(1,jj,1,1) &
                      + Dxt(i,2) * uur(2,jj,1,1) &
                      + Dxt(i,3) * uur(3,jj,1,1) &
                      + Dxt(i,4) * uur(4,jj,1,1) 
       end do
    end do

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,n
                w(i,j,k,e) = w(i,j,k,e) + Dyt(j,1) * uus(i,1,k,e) &
                                        + Dyt(j,2) * uus(i,2,k,e) &
                                        + Dyt(j,3) * uus(i,3,k,e) &
                                        + Dyt(j,4) * uus(i,4,k,e) 
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,n
                w(i,j,k,e) = w(i,j,k,e) + Dzt(k,1) * uut(i,j,1,e) &
                                        + Dzt(k,2) * uut(i,j,2,e) &
                                        + Dzt(k,3) * uut(i,j,3,e) &
                                        + Dzt(k,4) * uut(i,j,4,e) 
             end do
          end do
       end do
    end do

  end subroutine sx_ax_helm_lx4

  subroutine sx_ax_helm_lx2(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
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

    integer :: e, i, j, k, jj, kk
    real(kind=rp) :: ur(lx,lx,lx,n)
    real(kind=rp) :: us(lx,lx,lx,n)
    real(kind=rp) :: ut(lx,lx,lx,n)
    real(kind=rp) :: uur(lx,lx,lx,n)
    real(kind=rp) :: uus(lx,lx,lx,n)
    real(kind=rp) :: uut(lx,lx,lx,n)

    do i=1,lx
       do jj = 1, lx * lx * n
          ur(i,jj,1,1) = Dx(i,1) * u(1,jj,1,1) &
                       + Dx(i,2) * u(2,jj,1,1)
       end do
    end do

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,n
                us(i,j,k,e) = Dy(j,1) * u(i,1,k,e) &
                            + Dy(j,2) * u(i,2,k,e)
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,n
                ut(i,j,k,e) = Dz(k,1) * u(i,j,1,e) &
                            + Dz(k,2) * u(i,j,2,e)
             end do
          end do
       end do
    end do

    do i = 1, n * lx * lx * lx
       uur(i,1,1,1) = h1(i,1,1,1) * &
            ( G11(i,1,1,1) * ur(i,1,1,1) &
            + G12(i,1,1,1) * us(i,1,1,1) &
            + G13(i,1,1,1) * ut(i,1,1,1))

       uus(i,1,1,1) = h1(i,1,1,1) * &
            ( G22(i,1,1,1) * us(i,1,1,1) &
            + G12(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * ut(i,1,1,1) )

       uut(i,1,1,1) = h1(i,1,1,1) * &
            ( G33(i,1,1,1) * ut(i,1,1,1) &
            + G13(i,1,1,1) * ur(i,1,1,1) &
            + G23(i,1,1,1) * us(i,1,1,1))
    end do

    do i=1,lx
       do jj = 1, lx * lx * n
          w(i,jj,1,1) = Dxt(i,1) * uur(1,jj,1,1) &
                      + Dxt(i,2) * uur(2,jj,1,1)
       end do
    end do

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e=1,n
                w(i,j,k,e) = w(i,j,k,e) + Dyt(j,1) * uus(i,1,k,e) &
                                        + Dyt(j,2) * uus(i,2,k,e)
             end do
          end do
       end do
    end do

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e=1,n
                w(i,j,k,e) = w(i,j,k,e) + Dzt(k, 1) * uut(i,j,1,e) &
                                        + Dzt(k, 2) * uut(i,j,2,e)
             end do
          end do
       end do
    end do

  end subroutine sx_ax_helm_lx2

end module ax_helm_sx
