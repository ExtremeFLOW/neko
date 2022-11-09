!> Operators SX-Aurora backend
module opr_sx
  use sx_dudxyz
  use sx_opgrad
  use sx_conv1
  use sx_cdtp
  use sx_cfl
  use gather_scatter
  use num_types
  use space
  use coefs
  use math
  use mesh
  use field
  use mathops
  implicit none
  private

  public :: opr_sx_dudxyz, opr_sx_opgrad, opr_sx_cdtp, opr_sx_conv1, &
       opr_sx_curl, opr_sx_cfl

contains

  subroutine opr_sx_dudxyz(du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(in) ::  u, dr, ds, dt

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)
      select case(coef%Xh%lx)
      case(14)
         call sx_dudxyz_lx14(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(13)
         call sx_dudxyz_lx13(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(12)
         call sx_dudxyz_lx12(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(11)
         call sx_dudxyz_lx11(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(10)
         call sx_dudxyz_lx10(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(9)
         call sx_dudxyz_lx9(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(8)
         call sx_dudxyz_lx8(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(7)
         call sx_dudxyz_lx7(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(6)
         call sx_dudxyz_lx6(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(5)
         call sx_dudxyz_lx5(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(4)
         call sx_dudxyz_lx4(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(3)
         call sx_dudxyz_lx3(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case(2)
         call sx_dudxyz_lx2(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size())
      case default
         call sx_dudxyz_lx(du, u, dr, ds, dt, Xh%dx, Xh%dy, Xh%dz, &
              coef%jacinv, msh%nelv, dof%size(), Xh%lx)
      end select
    end associate

  end subroutine opr_sx_dudxyz

  subroutine opr_sx_opgrad(ux,uy,uz,u,coef) 
    type(coef_t), intent(in) :: coef  
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: u

    associate(Xh => coef%Xh, msh => coef%msh)
      select case(Xh%lx)
      case(18)
         call sx_opgrad_lx18(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(17)
         call sx_opgrad_lx17(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(16)
         call sx_opgrad_lx16(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(15)
         call sx_opgrad_lx15(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(14)
         call sx_opgrad_lx14(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(13)
         call sx_opgrad_lx13(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(12)
         call sx_opgrad_lx12(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(11)
         call sx_opgrad_lx11(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(10)
         call sx_opgrad_lx10(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(9)
         call sx_opgrad_lx9(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(8)
         call sx_opgrad_lx8(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(7)
         call sx_opgrad_lx7(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(6)
         call sx_opgrad_lx6(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(5)
         call sx_opgrad_lx5(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(4)
         call sx_opgrad_lx4(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(3)
         call sx_opgrad_lx3(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case(2)
         call sx_opgrad_lx2(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
      case default
         call sx_opgrad_lx(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv, Xh%lx)
      end select
    end associate

  end subroutine opr_sx_opgrad

  subroutine opr_sx_cdtp(dtx,x,dr,ds,dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dt

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)
      select case(Xh%lx)
      case(14)
         call sx_cdtp_lx14(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(13)
         call sx_cdtp_lx13(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(12)
         call sx_cdtp_lx12(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(11)
         call sx_cdtp_lx11(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(10)
         call sx_cdtp_lx10(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(9)
         call sx_cdtp_lx9(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(8)
         call sx_cdtp_lx8(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(7)
         call sx_cdtp_lx7(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(6)
         call sx_cdtp_lx6(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(5)
         call sx_cdtp_lx5(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(4)
         call sx_cdtp_lx4(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(3)
         call sx_cdtp_lx3(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case(2)
         call sx_cdtp_lx2(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size())
      case default
         call sx_cdtp_lx(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, &
              coef%B, coef%jac, msh%nelv, dof%size(), Xh%lx)
      end select
    end associate

  end subroutine opr_sx_cdtp

  subroutine opr_sx_conv1(du,u, vx, vy, vz, Xh, coef, nelv, gdim)  
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=rp), intent(inout) ::  du(Xh%lxyz,nelv)
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vx
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vy
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vz

    select case(Xh%lx)
    case(14)
       call sx_conv1_lx14(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(13)
       call sx_conv1_lx13(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)   
    case(12)
       call sx_conv1_lx12(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(11)
       call sx_conv1_lx11(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(10)
       call sx_conv1_lx10(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(9)
       call sx_conv1_lx9(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(8)
       call sx_conv1_lx8(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(7)
       call sx_conv1_lx7(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(6)
       call sx_conv1_lx6(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(5)
       call sx_conv1_lx5(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(4)
       call sx_conv1_lx4(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(3)
       call sx_conv1_lx3(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case(2)
       call sx_conv1_lx2(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim)
    case default
       call sx_conv1_lx(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, gdim, Xh%lx)
    end select

  end subroutine opr_sx_conv1

  subroutine opr_sx_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in)  :: c_Xh
    integer :: gdim, n

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim

    !     this%work1=dw/dy ; this%work2=dv/dz
    call opr_sx_dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call opr_sx_dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    end if
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call opr_sx_dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call opr_sx_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    else
       call rzero (work1%x, n)
       call opr_sx_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    end if
    !     this%work1=dv/dx ; this%work2=du/dy
    call opr_sx_dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call opr_sx_dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call sub3(w3%x, work1%x, work2%x, n)
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x,w2%x,w3%x,c_Xh%B, gdim, n)
    call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w3, GS_OP_ADD) 
    call opcolv  (w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)

  end subroutine opr_sx_curl

  function opr_sx_cfl(dt, u, v, w, Xh, coef, nelv, gdim) result(cfl)
    type(space_t) :: Xh
    type(coef_t) :: coef
    integer :: nelv, gdim
    real(kind=rp) :: dt
    real(kind=rp), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u, v, w
    real(kind=rp) :: cfl

   select case(Xh%lx)
   case (14)
      cfl = sx_cfl_lx14(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (13)
      cfl = sx_cfl_lx13(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (12)
      cfl = sx_cfl_lx12(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (11)
      cfl = sx_cfl_lx11(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (10)
      cfl = sx_cfl_lx10(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (9)
      cfl = sx_cfl_lx9(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (8)
      cfl = sx_cfl_lx8(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (7)
      cfl = sx_cfl_lx7(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (6)
      cfl = sx_cfl_lx6(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (5)
      cfl = sx_cfl_lx5(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (4)
      cfl = sx_cfl_lx4(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (3)
      cfl = sx_cfl_lx3(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case (2)
      cfl = sx_cfl_lx2(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim)
   case default
      cfl = sx_cfl_lx(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, gdim, Xh%lx)
   end select

  end function opr_sx_cfl

end module opr_sx
