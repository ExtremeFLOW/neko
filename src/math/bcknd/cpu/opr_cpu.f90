!> Operators CPU backend
module opr_cpu
  use cpu_dudxyz
  use cpu_opgrad
  use cpu_cdtp
  use cpu_conv1
  use num_types
  use space
  use coefs
  use math
  use mesh
  use field
  use gather_scatter
  use mathops
  implicit none
  private

  public :: opr_cpu_dudxyz, opr_cpu_opgrad, opr_cpu_cdtp, &
       opr_cpu_conv1, opr_cpu_curl

  
contains

  subroutine opr_cpu_dudxyz(du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(in) ::  u, dr, ds, dt

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)
      select case(coef%Xh%lx)
      case(12)
         call cpu_dudxyz_lx12(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(11)
         call cpu_dudxyz_lx11(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(10)
         call cpu_dudxyz_lx10(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(9)
         call cpu_dudxyz_lx9(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(8)
         call cpu_dudxyz_lx8(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(7)
         call cpu_dudxyz_lx7(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(6)
         call cpu_dudxyz_lx6(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(5)
         call cpu_dudxyz_lx5(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(4)
         call cpu_dudxyz_lx4(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(3)
         call cpu_dudxyz_lx3(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      case(2)
         call cpu_dudxyz_lx2(du, u, dr, ds, dt, & 
              Xh%dx, Xh%dy, Xh%dz, coef%jacinv, msh%nelv)
      end select

    end associate

  end subroutine opr_cpu_dudxyz

  subroutine opr_cpu_opgrad(ux, uy, uz, u, coef, e_start, e_end) 
    type(coef_t), intent(in) :: coef  
    integer, intent(in) :: e_start, e_end
    real(kind=rp), dimension(coef%Xh%lxyz,e_end-e_start+1), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz,e_end-e_start+1), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz,e_end-e_start+1), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz,e_end-e_start+1), intent(in) :: u
    integer :: e_len
    e_len = e_end-e_start+1
    associate(Xh => coef%Xh, msh => coef%msh, &
         drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
         dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
         dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz)
      
      select case(Xh%lx)
      case(12)
         call cpu_opgrad_lx12(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case(11)
         call cpu_opgrad_lx11(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case(10)
         call cpu_opgrad_lx10(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)

      case(9)
         call cpu_opgrad_lx9(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)  
      case(8)
         call cpu_opgrad_lx8(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)  
      case(7)
         call cpu_opgrad_lx7(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)  
      case(6)
         call cpu_opgrad_lx6(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)  
      case(5)
         call cpu_opgrad_lx5(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)  
      case(4)
         call cpu_opgrad_lx4(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)  
      case(3)
         call cpu_opgrad_lx3(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len) 
      case(2)
         call cpu_opgrad_lx2(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len) 
      end select
    end associate

  end subroutine opr_cpu_opgrad

  subroutine opr_cpu_cdtp(dtx, x, dr, ds, dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dt

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)
      select case(Xh%lx)
      case(12)
         call cpu_cdtp_lx12(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(11)
         call cpu_cdtp_lx11(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(10)
         call cpu_cdtp_lx10(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(9)
         call cpu_cdtp_lx9(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(8)
         call cpu_cdtp_lx8(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(7)
         call cpu_cdtp_lx7(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(6)
         call cpu_cdtp_lx6(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(5)
         call cpu_cdtp_lx5(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(4)
         call cpu_cdtp_lx4(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(3)
         call cpu_cdtp_lx3(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      case(2)
         call cpu_cdtp_lx2(dtx, x, dr, ds, dt, &
              Xh%dxt, Xh%dyt, Xh%dzt, coef%B, coef%jac, msh%nelv)
      end select
    end associate

  end subroutine opr_cpu_cdtp

  subroutine opr_cpu_conv1(du, u, vx, vy, vz, Xh, coef, e_start, e_end)
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e_start, e_end
    real(kind=rp), intent(inout) ::  du(Xh%lxyz,e_end-e_start+1)
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,e_end-e_start+1) ::  u
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,e_end-e_start+1) ::  vx
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,e_end-e_start+1) ::  vy
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,e_end-e_start+1) ::  vz
    integer :: e_len

    e_len = e_end-e_start+1
    associate(drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
         dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
         dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz, &
         jacinv => coef%jacinv)
      select case(Xh%lx)
      case(12)
         call cpu_conv1_lx12(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)
      case(11)
         call cpu_conv1_lx11(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)
      case(10)
         call cpu_conv1_lx10(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)         
      case(9)
         call cpu_conv1_lx9(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)         
      case(8)
         call cpu_conv1_lx8(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)         
      case(7)
         call cpu_conv1_lx7(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)         
      case(6)
         call cpu_conv1_lx6(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)
      case(5)
         call cpu_conv1_lx5(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)
      case(4)
         call cpu_conv1_lx4(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)
      case(3)
         call cpu_conv1_lx3(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)
      case(2)
         call cpu_conv1_lx2(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              jacinv(1,1,1,e_start), e_len)
      end select
    end associate
    
  end subroutine opr_cpu_conv1

  subroutine opr_cpu_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
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
    call opr_cpu_dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call opr_cpu_dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    end if
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call opr_cpu_dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
       call opr_cpu_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    else
       call rzero (work1%x, n)
       call opr_cpu_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    end if
    !     this%work1=dv/dx ; this%work2=du/dy
    call opr_cpu_dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call opr_cpu_dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call sub3(w3%x, work1%x, work2%x, n)
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x,w2%x,w3%x,c_Xh%B, gdim, n)
    call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w3, GS_OP_ADD) 
    call opcolv  (w1%x,w2%x,w3%x,c_Xh%Binv, gdim, n)

  end subroutine opr_cpu_curl



end module opr_cpu
