!> conv1 kernels
module cpu_conv1
  use num_types
  implicit none

contains

  subroutine cpu_conv1_lx12(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx12

   subroutine cpu_conv1_lx11(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx11
   
   subroutine cpu_conv1_lx10(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx10
   
   subroutine cpu_conv1_lx9(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx9

   subroutine cpu_conv1_lx8(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx8

   subroutine cpu_conv1_lx7(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e)
           end do
        end do
       
        do i = 1, lx * lx * lx
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx7

   subroutine cpu_conv1_lx6(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k
    
    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) 
           end do
        end do
       
        do i = 1, lx * lx * lx
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx6

   subroutine cpu_conv1_lx5(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) 
          end do
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) 
             end do
          end do
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) 
           end do
        end do
       
        do i = 1, lx * lx * lx
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx5

   subroutine cpu_conv1_lx4(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) 
          end do
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) 
             end do
          end do
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) 
           end do
        end do
       
        do i = 1, lx * lx * lx
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx4

   subroutine cpu_conv1_lx3(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) 
          end do
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) 
             end do
          end do
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) 
           end do
        end do
       
        do i = 1, lx * lx * lx
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx3

   subroutine cpu_conv1_lx2(du, u, vx, vy, vz, dx, dy, dz, &
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
    real(kind=rp), dimension(lx,lx,lx) ::  dudr
    real(kind=rp), dimension(lx,lx,lx) ::  duds
    real(kind=rp), dimension(lx,lx,lx) ::  dudt
    integer :: e, i, j, k

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             dudr(i,j,1) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) 
          end do
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                duds(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) 
             end do
          end do
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dudt(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) 
           end do
        end do
       
        do i = 1, lx * lx * lx
           du(i,1,1,e) = jacinv(i,1,1,e) &
                       * ( vx(i,1,1,e) &
                         * ( drdx(i,1,1,e) * dudr(i,1,1) &
                           + dsdx(i,1,1,e) * duds(i,1,1) &
                           + dtdx(i,1,1,e) * dudt(i,1,1) ) &
                         + vy(i,1,1,e) &
                         * ( drdy(i,1,1,e) * dudr(i,1,1) &
                           + dsdy(i,1,1,e) * duds(i,1,1) &
                           + dtdy(i,1,1,e) * dudt(i,1,1) ) &
                         + vz(i,1,1,e) &
                         * ( drdz(i,1,1,e) * dudr(i,1,1) &
                           + dsdz(i,1,1,e) * duds(i,1,1) &
                           + dtdz(i,1,1,e) * dudt(i,1,1) ) )
        end do
     end do
    
   end subroutine cpu_conv1_lx2

end module cpu_conv1
