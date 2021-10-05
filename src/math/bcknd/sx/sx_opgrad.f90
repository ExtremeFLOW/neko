!> Gradient kernels for SX-Aurora
module sx_opgrad
  use num_types
  implicit none

contains

  subroutine sx_opgrad_lx12(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx12

  subroutine sx_opgrad_lx11(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx11

  subroutine sx_opgrad_lx10(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx10

  subroutine sx_opgrad_lx9(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx9
  
  subroutine sx_opgrad_lx8(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx8

    subroutine sx_opgrad_lx7(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx7
  
  subroutine sx_opgrad_lx6(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx6

  subroutine sx_opgrad_lx5(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx5

  subroutine sx_opgrad_lx4(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx4
  
  subroutine sx_opgrad_lx3(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx3

  subroutine sx_opgrad_lx2(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz    
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    real(kind=rp) :: ur(lx,lx,lx, n)
    real(kind=rp) :: us(lx,lx,lx, n)
    real(kind=rp) :: ut(lx,lx,lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, ii, jj, kk

    do i = 1,lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do
 
    do k = 1,lx
       do i = 1,lx
          do j = 1,lx
             do e = 1,n  
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1,lx
       do i = 1,lx
          do k = 1,lx
             do e = 1,n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do
    
    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )
          
          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx2
  
end module sx_opgrad
