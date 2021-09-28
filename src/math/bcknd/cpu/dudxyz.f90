!> Derivative kernels
module dudxyz
  use num_types
  use math
  implicit none

contains
  
  subroutine dudxyz_lx12(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k, l

    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx12

   subroutine dudxyz_lx11(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 11
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx11

   subroutine dudxyz_lx10(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 10
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx10

   subroutine dudxyz_lx9(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 9
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx9

   subroutine dudxyz_lx8(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 8
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                          + dx(i,2) * u(2,j,1,e) &
                          + dx(i,3) * u(3,j,1,e) &
                          + dx(i,4) * u(4,j,1,e) &
                          + dx(i,5) * u(5,j,1,e) &
                          + dx(i,6) * u(6,j,1,e) &
                          + dx(i,7) * u(7,j,1,e) &
                          + dx(i,8) * u(8,j,1,e) 
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx8

   subroutine dudxyz_lx7(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 7
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                          + dx(i,2) * u(2,j,1,e) &
                          + dx(i,3) * u(3,j,1,e) &
                          + dx(i,4) * u(4,j,1,e) &
                          + dx(i,5) * u(5,j,1,e) &
                          + dx(i,6) * u(6,j,1,e) &
                          + dx(i,7) * u(7,j,1,e) 
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                             + dy(j,2) * u(i,2,k,e) &
                             + dy(j,3) * u(i,3,k,e) &
                             + dy(j,4) * u(i,4,k,e) &
                             + dy(j,5) * u(i,5,k,e) &
                             + dy(j,6) * u(i,6,k,e) &
                             + dy(j,7) * u(i,7,k,e)
              end do
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                          + dz(k,2) * u(i,1,2,e) &
                          + dz(k,3) * u(i,1,3,e) &
                          + dz(k,4) * u(i,1,4,e) &
                          + dz(k,5) * u(i,1,5,e) &
                          + dz(k,6) * u(i,1,6,e) &
                          + dz(k,7) * u(i,1,7,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx7

   subroutine dudxyz_lx6(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 6
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                          + dx(i,2) * u(2,j,1,e) &
                          + dx(i,3) * u(3,j,1,e) &
                          + dx(i,4) * u(4,j,1,e) &
                          + dx(i,5) * u(5,j,1,e) &
                          + dx(i,6) * u(6,j,1,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                             + dy(j,2) * u(i,2,k,e) &
                             + dy(j,3) * u(i,3,k,e) &
                             + dy(j,4) * u(i,4,k,e) &
                             + dy(j,5) * u(i,5,k,e) &
                             + dy(j,6) * u(i,6,k,e)
              end do
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                          + dz(k,2) * u(i,1,2,e) &
                          + dz(k,3) * u(i,1,3,e) &
                          + dz(k,4) * u(i,1,4,e) &
                          + dz(k,5) * u(i,1,5,e) &
                          + dz(k,6) * u(i,1,6,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx6

   subroutine dudxyz_lx5(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 5
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                          + dx(i,2) * u(2,j,1,e) &
                          + dx(i,3) * u(3,j,1,e) &
                          + dx(i,4) * u(4,j,1,e) &
                          + dx(i,5) * u(5,j,1,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                             + dy(j,2) * u(i,2,k,e) &
                             + dy(j,3) * u(i,3,k,e) &
                             + dy(j,4) * u(i,4,k,e) &
                             + dy(j,5) * u(i,5,k,e)
              end do
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                          + dz(k,2) * u(i,1,2,e) &
                          + dz(k,3) * u(i,1,3,e) &
                          + dz(k,4) * u(i,1,4,e) &
                          + dz(k,5) * u(i,1,5,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx5

   subroutine dudxyz_lx4(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 4
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                             + dy(j,2) * u(i,2,k,e) &
                             + dy(j,3) * u(i,3,k,e) &
                             + dy(j,4) * u(i,4,k,e)
              end do
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                          + dz(k,2) * u(i,1,2,e) &
                          + dz(k,3) * u(i,1,3,e) &
                          + dz(k,4) * u(i,1,4,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx4

   subroutine dudxyz_lx3(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 3
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                             + dy(j,2) * u(i,2,k,e) &
                             + dy(j,3) * u(i,3,k,e)
              end do
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                          + dz(k,2) * u(i,1,2,e) &
                          + dz(k,3) * u(i,1,3,e)
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx3

   subroutine dudxyz_lx2(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
     integer, parameter :: lx = 2
     integer, intent(in) :: nel
     real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
     real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
     real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
     real(kind=rp), dimension(lx,lx,lx) :: drst
     integer :: e, i, j, k, l

     do e = 1, nel
        do j = 1, lx * lx
           do i = 1, lx
              du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                          + dx(i,2) * u(2,j,1,e) 
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
        end do

        do k = 1, lx
           do j = 1, lx
              do i = 1, lx
                 drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                             + dy(j,2) * u(i,2,k,e)
              end do
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
        end do

        do k = 1, lx
           do i = 1, lx*lx
              drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                          + dz(k,2) * u(i,1,2,e) 
           end do
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
        end do

        do i = 1, lx * lx * lx
           du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
        end do

     end do

   end subroutine dudxyz_lx2

end module dudxyz
