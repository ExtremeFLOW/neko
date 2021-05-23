!> Derivative kernels for SX-Aurora
module sx_dudxyz
  use num_types
  use math
  implicit none

contains

  subroutine sx_dudxyz_lx12(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx12

  subroutine sx_dudxyz_lx11(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx11

  subroutine sx_dudxyz_lx10(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx10

  subroutine sx_dudxyz_lx9(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx9

  subroutine sx_dudxyz_lx8(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx8

  subroutine sx_dudxyz_lx7(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx7

  subroutine sx_dudxyz_lx6(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx6

  subroutine sx_dudxyz_lx5(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx5

  subroutine sx_dudxyz_lx4(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx4

  subroutine sx_dudxyz_lx3(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx3

  subroutine sx_dudxyz_lx2(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=rp) :: wr, ws, wt, www

    do i=1,lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk=1,lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx2

end module sx_dudxyz
