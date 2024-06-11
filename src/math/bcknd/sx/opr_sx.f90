!> Operators SX-Aurora backend
module opr_sx
   use sx_dudxyz
   use sx_opgrad
   use sx_conv1
   use sx_cdtp
   use sx_cfl
   use sx_lambda2
   use sx_conv_fst_3d
   use sx_set_convect_new
   use gather_scatter
   use interpolation
   use num_types, only : rp
   use space, only : space_t
   use coefs, only : coef_t
   use math
   use field, only : field_t
   use mathops
  implicit none
  private

  public :: opr_sx_dudxyz, opr_sx_opgrad, opr_sx_cdtp, opr_sx_conv1, &
       opr_sx_curl, opr_sx_cfl, opr_sx_lambda2, opr_sx_conv_fst_3d, opr_sx_set_convect_new

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

  subroutine opr_sx_conv1(du,u, vx, vy, vz, Xh, coef, nelv)
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: nelv
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
            coef%jacinv, nelv)
    case(13)
       call sx_conv1_lx13(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(12)
       call sx_conv1_lx12(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(11)
       call sx_conv1_lx11(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(10)
       call sx_conv1_lx10(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(9)
       call sx_conv1_lx9(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(8)
       call sx_conv1_lx8(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(7)
       call sx_conv1_lx7(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(6)
       call sx_conv1_lx6(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(5)
       call sx_conv1_lx5(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(4)
       call sx_conv1_lx4(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(3)
       call sx_conv1_lx3(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case(2)
       call sx_conv1_lx2(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
    case default
       call sx_conv1_lx(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, Xh%lx)
    end select

  end subroutine opr_sx_conv1

  subroutine opr_sx_conv_fst_3d(du, u, c, Xh_GLL, Xh_GL, coef_GLL, coef_GL, GLL_to_GL)
    type(space_t), intent(in) :: Xh_GL
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(coef_t), intent(in) :: coef_GL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%ly, Xh_GLL%lz, coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: u(Xh_GL%lx, Xh_GL%lx, Xh_GL%lx, coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: c(Xh_GL%lxyz,coef_GL%msh%nelv,3)
    associate(dx => Xh_GL%dx, dy => Xh_GL%dy, dz => Xh_GL%dz, &
         lx => Xh_GL%lx, nelv => coef_GL%msh%nelv)
      
      select case(lx)
      case(18)
         call sx_conv_fst_3d_lx18(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(17)
         call sx_conv_fst_3d_lx17(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(16)
         call sx_conv_fst_3d_lx16(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(15)
         call sx_conv_fst_3d_lx15(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(14)
         call sx_conv_fst_3d_lx14(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(13)
         call sx_conv_fst_3d_lx13(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(12)
         call sx_conv_fst_3d_lx12(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(11)
         call sx_conv_fst_3d_lx11(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(10)
         call sx_conv_fst_3d_lx10(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(9)
         call sx_conv_fst_3d_lx9(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(8)
         call sx_conv_fst_3d_lx8(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(7)
         call sx_conv_fst_3d_lx7(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(6)
         call sx_conv_fst_3d_lx6(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(5)
         call sx_conv_fst_3d_lx5(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(4)
         call sx_conv_fst_3d_lx4(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(3)
         call sx_conv_fst_3d_lx3(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case(2)
         call sx_conv_fst_3d_lx2(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case default
         call sx_conv_fst_3d_lx(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv, lx)
      end select
    end associate

  end subroutine opr_sx_conv_fst_3d

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
    call c_Xh%gs_h%op(w1, GS_OP_ADD)
    call c_Xh%gs_h%op(w2, GS_OP_ADD)
    call c_Xh%gs_h%op(w3, GS_OP_ADD)
    call opcolv(w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)

  end subroutine opr_sx_curl

  function opr_sx_cfl(dt, u, v, w, Xh, coef, nelv) result(cfl)
    type(space_t) :: Xh
    type(coef_t) :: coef
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u, v, w
    real(kind=rp) :: cfl

    select case(Xh%lx)
    case (14)
       cfl = sx_cfl_lx14(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (13)
       cfl = sx_cfl_lx13(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (12)
       cfl = sx_cfl_lx12(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (11)
       cfl = sx_cfl_lx11(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (10)
       cfl = sx_cfl_lx10(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (9)
       cfl = sx_cfl_lx9(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (8)
       cfl = sx_cfl_lx8(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (7)
       cfl = sx_cfl_lx7(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (6)
       cfl = sx_cfl_lx6(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (5)
       cfl = sx_cfl_lx5(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (4)
       cfl = sx_cfl_lx4(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (3)
       cfl = sx_cfl_lx3(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (2)
       cfl = sx_cfl_lx2(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case default
       cfl = sx_cfl_lx(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, Xh%lx)
    end select

  end function opr_sx_cfl

  subroutine opr_sx_lambda2(lambda2, u, v, w, coef)
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: lambda2
    type(field_t), intent(in) :: u, v, w

    associate(Xh => coef%Xh, msh => coef%msh)
      select case(Xh%lx)
      case (18)
         call sx_lambda2_lx18(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (17)
         call sx_lambda2_lx17(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (16)
         call sx_lambda2_lx16(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (15)
         call sx_lambda2_lx15(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (14)
         call sx_lambda2_lx14(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (13)
         call sx_lambda2_lx13(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (12)
         call sx_lambda2_lx12(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (11)
         call sx_lambda2_lx11(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (10)
         call sx_lambda2_lx10(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (9)
         call sx_lambda2_lx9(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (8)
         call sx_lambda2_lx8(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (7)
         call sx_lambda2_lx7(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (6)
         call sx_lambda2_lx6(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (5)
         call sx_lambda2_lx5(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (4)
         call sx_lambda2_lx4(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (3)
         call sx_lambda2_lx3(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (2)
         call sx_lambda2_lx2(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case default
         call sx_lambda2_lx(lambda2%x, u%x, v%x, w%x, &
                            Xh%dx, Xh%dy, Xh%dz, &
                            coef%drdx, coef%dsdx, coef%dtdx, &
                            coef%drdy, coef%dsdy, coef%dtdy, &
                            coef%drdz, coef%dsdz, coef%dtdz, &
                            Xh%w3, coef%B, msh%nelv, Xh%lx)
      end select
    end associate
    
  end subroutine opr_sx_lambda2  

  subroutine opr_sx_set_convect_new(cr, cs, ct, cx, cy, cz, Xh, coef)   
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), intent(in) :: cx, cy, cz
    associate(drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
      dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
      dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz, &  
      nelv => coef%msh%nelv, lx=>Xh%lx, w3 => Xh%w3)
     
      select case(lx)
      case(18)
         call sx_set_convect_new_lx18(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(17)
         call sx_set_convect_new_lx17(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(16)
         call sx_set_convect_new_lx16(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(15)
         call sx_set_convect_new_lx15(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(14)
         call sx_set_convect_new_lx14(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(13)
         call sx_set_convect_new_lx13(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(12)
         call sx_set_convect_new_lx12(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(11)
         call sx_set_convect_new_lx11(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(10)
         call sx_set_convect_new_lx10(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(9)
         call sx_set_convect_new_lx9(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(8)
         call sx_set_convect_new_lx8(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(7)
         call sx_set_convect_new_lx7(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(6)
         call sx_set_convect_new_lx6(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(5)
         call sx_set_convect_new_lx5(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(4)
         call sx_set_convect_new_lx4(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(3)
         call sx_set_convect_new_lx3(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case(2)
         call sx_set_convect_new_lx2(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case default
         call sx_set_convect_new_lx(cr, cs, ct, cx, cy, cz, &
              drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv, lx)
      end select
    end associate
 end subroutine opr_sx_set_convect_new

end module opr_sx
