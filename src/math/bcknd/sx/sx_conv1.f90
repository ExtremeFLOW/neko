!> conv1 SX-Aurora kernels
module sx_conv1
  use num_types
  implicit none

contains

  subroutine sx_conv1_lx(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim, lx)
    integer, intent(in) :: nelv, gdim, lx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx
  
  subroutine sx_conv1_lx14(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 14
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx14

  subroutine sx_conv1_lx13(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 13
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx13
  
  subroutine sx_conv1_lx12(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 12
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx12

  subroutine sx_conv1_lx11(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 11
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx11
  
  subroutine sx_conv1_lx10(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 10
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx10

  subroutine sx_conv1_lx9(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 9
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx9
  
  subroutine sx_conv1_lx8(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 8
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx8

  subroutine sx_conv1_lx7(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 7
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx7
  
  subroutine sx_conv1_lx6(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 6
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx6
  
  subroutine sx_conv1_lx5(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 5
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx5

  subroutine sx_conv1_lx4(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 4
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx4

  subroutine sx_conv1_lx3(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 3
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j=1, lx
       do i=1, lx
          do k=1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx3

  subroutine sx_conv1_lx2(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, gdim)
    integer, parameter :: lx = 2
    integer, intent(in) :: nelv, gdim
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) ::  u, vx, vy, vz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz        
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudr
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  duds
    real(kind=rp), dimension(lx,lx,lx,nelv) ::  dudt
    real(kind=rp) :: wr, ws, wt, www
    integer :: e, i, j, k, ii, jj, kk    

    do i = 1, lx
       do jj = 1, lx * lx * nelv     
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          dudr(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv     
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do
    
  end subroutine sx_conv1_lx2
  
end module sx_conv1
