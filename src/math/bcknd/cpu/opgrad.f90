!> Gradient kernels
module opgrad
  use num_types
  implicit none

contains

  subroutine opgrad_lx12(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e) &
                        + dx(i,8) * u(8,j,1,e) &
                        + dx(i,9) * u(9,j,1,e) &
                        + dx(i,10) * u(10,j,1,e) &
                        + dx(i,11) * u(11,j,1,e) &
                        + dx(i,12) * u(12,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) &
                          + dy(j,7) * u(i,7,k,e) &
                          + dy(j,8) * u(i,8,k,e) &
                          + dy(j,9) * u(i,9,k,e) &
                          + dy(j,10) * u(i,10,k,e) &
                          + dy(j,11) * u(i,11,k,e) &
                          + dy(j,12) * u(i,12,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) &
                       + dz(k,7) * u(i,1,7,e) &
                       + dz(k,8) * u(i,1,8,e) &
                       + dz(k,9) * u(i,1,9,e) &
                       + dz(k,10) * u(i,1,10,e) &
                       + dz(k,11) * u(i,1,11,e) &
                       + dz(k,12) * u(i,1,12,e)
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx12

  subroutine opgrad_lx11(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e) &
                        + dx(i,8) * u(8,j,1,e) &
                        + dx(i,9) * u(9,j,1,e) &
                        + dx(i,10) * u(10,j,1,e) &
                        + dx(i,11) * u(11,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) &
                          + dy(j,7) * u(i,7,k,e) &
                          + dy(j,8) * u(i,8,k,e) &
                          + dy(j,9) * u(i,9,k,e) &
                          + dy(j,10) * u(i,10,k,e) &
                          + dy(j,11) * u(i,11,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) &
                       + dz(k,7) * u(i,1,7,e) &
                       + dz(k,8) * u(i,1,8,e) &
                       + dz(k,9) * u(i,1,9,e) &
                       + dz(k,10) * u(i,1,10,e) &
                       + dz(k,11) * u(i,1,11,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx11

  subroutine opgrad_lx10(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e) &
                        + dx(i,8) * u(8,j,1,e) &
                        + dx(i,9) * u(9,j,1,e) &
                        + dx(i,10) * u(10,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) &
                          + dy(j,7) * u(i,7,k,e) &
                          + dy(j,8) * u(i,8,k,e) &
                          + dy(j,9) * u(i,9,k,e) &
                          + dy(j,10) * u(i,10,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) &
                       + dz(k,7) * u(i,1,7,e) &
                       + dz(k,8) * u(i,1,8,e) &
                       + dz(k,9) * u(i,1,9,e) &
                       + dz(k,10) * u(i,1,10,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                         + dsdx(i,1,1,e) * us(i,1,1) &
                         + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx10

  subroutine opgrad_lx9(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e) &
                        + dx(i,8) * u(8,j,1,e) &
                        + dx(i,9) * u(9,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) &
                          + dy(j,7) * u(i,7,k,e) &
                          + dy(j,8) * u(i,8,k,e) &
                          + dy(j,9) * u(i,9,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) &
                       + dz(k,7) * u(i,1,7,e) &
                       + dz(k,8) * u(i,1,8,e) &
                       + dz(k,9) * u(i,1,9,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx9

  subroutine opgrad_lx8(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e) &
                        + dx(i,8) * u(8,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) &
                          + dy(j,7) * u(i,7,k,e) &
                          + dy(j,8) * u(i,8,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) &
                       + dz(k,7) * u(i,1,7,e) &
                       + dz(k,8) * u(i,1,8,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx8

  subroutine opgrad_lx7(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) &
                          + dy(j,7) * u(i,7,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) &
                       + dz(k,7) * u(i,1,7,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx7

  subroutine opgrad_lx6(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx6

  subroutine opgrad_lx5(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx5

  subroutine opgrad_lx4(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx4

  subroutine opgrad_lx3(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx3

  subroutine opgrad_lx2(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx)
    real(kind=rp) :: us(lx,lx,lx)
    real(kind=rp) :: ut(lx,lx,lx)
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) 
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) 
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) 
          end do
       end do
    
       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )          
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine opgrad_lx2

end module opgrad
