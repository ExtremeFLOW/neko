!> Derivative kernels for SX-Aurora
module sx_dudxyz
  use num_types
  use space
  use coefs
  use field
  use math
  implicit none

contains

  subroutine sx_dudxyz_lx12(du, u, dr, ds, dt, coef)
    integer, parameter :: lx = 12
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  u
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dr
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  ds
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dt
    real(kind=dp) :: drst(lx,lx,lx,coef%msh%nelv)
    real(kind=dp) :: dx(lx,lx), dy(lx, lx), dz(lx, lx)
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=dp) :: wr, ws, wt, www

    dx = coef%Xh%dx
    dy = coef%Xh%dy
    dz = coef%Xh%dz
    nelv = coef%msh%nelv

    do i=1,lx
       do jj = 1, lx*lx*coef%msh%nelv
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, coef%dof%n_dofs)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, ds, coef%dof%n_dofs)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, dt, coef%dof%n_dofs)
    call col2 (du,coef%jacinv,coef%dof%n_dofs)
  end subroutine sx_dudxyz_lx12

  subroutine sx_dudxyz_lx11(du, u, dr, ds, dt, coef)
    integer, parameter :: lx = 11
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  u
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dr
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  ds
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dt
    real(kind=dp) :: drst(lx,lx,lx,coef%msh%nelv)
    real(kind=dp) :: dx(lx,lx), dy(lx, lx), dz(lx, lx)
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=dp) :: wr, ws, wt, www

    dx = coef%Xh%dx
    dy = coef%Xh%dy
    dz = coef%Xh%dz
    nelv = coef%msh%nelv

    do i=1,lx
       do jj = 1, lx*lx*coef%msh%nelv
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, coef%dof%n_dofs)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, ds, coef%dof%n_dofs)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, dt, coef%dof%n_dofs)
    call col2 (du,coef%jacinv,coef%dof%n_dofs)
  end subroutine sx_dudxyz_lx11

  subroutine sx_dudxyz_lx10(du, u, dr, ds, dt, coef)
    integer, parameter :: lx = 10
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  u
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dr
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  ds
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dt
    real(kind=dp) :: drst(lx,lx,lx,coef%msh%nelv)
    real(kind=dp) :: dx(lx,lx), dy(lx, lx), dz(lx, lx)
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=dp) :: wr, ws, wt, www

    dx = coef%Xh%dx
    dy = coef%Xh%dy
    dz = coef%Xh%dz
    nelv = coef%msh%nelv

    do i=1,lx
       do jj = 1, lx*lx*coef%msh%nelv
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, coef%dof%n_dofs)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, ds, coef%dof%n_dofs)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, dt, coef%dof%n_dofs)
    call col2 (du,coef%jacinv,coef%dof%n_dofs)
  end subroutine sx_dudxyz_lx10

  subroutine sx_dudxyz_lx9(du, u, dr, ds, dt, coef)
    integer, parameter :: lx = 9
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  u
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dr
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  ds
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dt
    real(kind=dp) :: drst(lx,lx,lx,coef%msh%nelv)
    real(kind=dp) :: dx(lx,lx), dy(lx, lx), dz(lx, lx)
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=dp) :: wr, ws, wt, www

    dx = coef%Xh%dx
    dy = coef%Xh%dy
    dz = coef%Xh%dz
    nelv = coef%msh%nelv

    do i=1,lx
       do jj = 1, lx*lx*coef%msh%nelv
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, coef%dof%n_dofs)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, ds, coef%dof%n_dofs)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, dt, coef%dof%n_dofs)
    call col2 (du,coef%jacinv,coef%dof%n_dofs)
  end subroutine sx_dudxyz_lx9

  subroutine sx_dudxyz_lx8(du, u, dr, ds, dt, coef)
    integer, parameter :: lx = 8
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  u
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dr
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  ds
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dt
    real(kind=dp) :: drst(lx,lx,lx,coef%msh%nelv)
    real(kind=dp) :: dx(lx,lx), dy(lx, lx), dz(lx, lx)
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=dp) :: wr, ws, wt, www

    dx = coef%Xh%dx
    dy = coef%Xh%dy
    dz = coef%Xh%dz
    nelv = coef%msh%nelv

    do i=1,lx
       do jj = 1, lx*lx*coef%msh%nelv
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, coef%dof%n_dofs)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, ds, coef%dof%n_dofs)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, dt, coef%dof%n_dofs)
    call col2 (du,coef%jacinv,coef%dof%n_dofs)
  end subroutine sx_dudxyz_lx8

  subroutine sx_dudxyz_lx6(du, u, dr, ds, dt, coef)
    integer, parameter :: lx = 6
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  u
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dr
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  ds
    real(kind=dp), dimension(lx,lx,lx,coef%msh%nelv), intent(inout) ::  dt
    real(kind=dp) :: drst(lx,lx,lx,coef%msh%nelv)
    real(kind=dp) :: dx(lx,lx), dy(lx, lx), dz(lx, lx)
    integer :: e, k, lxy, lyz, lxyz
    integer :: i, j, ii, jj, kk, nelv 
    real(kind=dp) :: wr, ws, wt, www

    dx = coef%Xh%dx
    dy = coef%Xh%dy
    dz = coef%Xh%dz
    nelv = coef%msh%nelv

    do i=1,lx
       do jj = 1, lx*lx*coef%msh%nelv
          wr = 0d0
          do kk=1,lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, coef%dof%n_dofs)

    do k=1,lx
       do i=1,lx
          do j=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, ds, coef%dof%n_dofs)

    do j=1,lx
       do i=1,lx
          do k=1,lx
             do e = 1,nelv     
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

    call addcol3(du, drst, dt, coef%dof%n_dofs)
    call col2 (du,coef%jacinv,coef%dof%n_dofs)
  end subroutine sx_dudxyz_lx6

end module sx_dudxyz
