!> Operators SX-Aurora backend
module opr_sx
  use gather_scatter, only : gs_t, GS_OP_ADD
  use interpolation, only : interpolator_t
  use num_types, only : rp
  use space, only : space_t
  use coefs, only : coef_t
  use math, only : sub3, copy, rzero
  use field, only : field_t
  use mathops, only : opcolv
  implicit none
  private

  public :: opr_sx_dudxyz, opr_sx_opgrad, opr_sx_cdtp, opr_sx_conv1, &
       opr_sx_curl, opr_sx_cfl, opr_sx_lambda2, opr_sx_convect_scalar, &
       opr_sx_set_convect_rst


  interface
     module subroutine opr_sx_dudxyz(du, u, dr, ds, dt, coef)
       type(coef_t), intent(in), target :: coef
       real(kind=rp), intent(inout), &
            dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv) ::  du
       real(kind=rp), intent(in), &
            dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv) :: &
            u, dr, ds, dt
     end subroutine opr_sx_dudxyz

     module subroutine opr_sx_opgrad(ux, uy, uz, u, coef)
       type(coef_t), intent(in) :: coef
       real(kind=rp), intent(inout) :: ux(coef%Xh%lxyz, coef%msh%nelv)
       real(kind=rp), intent(inout) :: uy(coef%Xh%lxyz, coef%msh%nelv)
       real(kind=rp), intent(inout) :: uz(coef%Xh%lxyz, coef%msh%nelv)
       real(kind=rp), intent(in) :: u(coef%Xh%lxyz, coef%msh%nelv)
     end subroutine opr_sx_opgrad

     module subroutine opr_sx_cdtp(dtx, x, dr, ds, dt, coef)
       type(coef_t), intent(in) :: coef
       real(kind=rp), intent(inout) :: dtx(coef%Xh%lxyz, coef%msh%nelv)
       real(kind=rp), intent(inout) :: x(coef%Xh%lxyz, coef%msh%nelv)
       real(kind=rp), intent(in) :: dr(coef%Xh%lxyz, coef%msh%nelv)
       real(kind=rp), intent(in) :: ds(coef%Xh%lxyz, coef%msh%nelv)
       real(kind=rp), intent(in) :: dt(coef%Xh%lxyz, coef%msh%nelv)
     end subroutine opr_sx_cdtp

     module subroutine opr_sx_conv1(du, u, vx, vy, vz, Xh, coef, nelv)
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       integer, intent(in) :: nelv
       real(kind=rp), intent(inout) ::  du(Xh%lxyz, nelv)
       real(kind=rp), intent(inout) ::  u(Xh%lx, Xh%ly, Xh%lz, nelv)
       real(kind=rp), intent(inout) ::  vx(Xh%lx, Xh%ly, Xh%lz, nelv)
       real(kind=rp), intent(inout) ::  vy(Xh%lx, Xh%ly, Xh%lz, nelv)
       real(kind=rp), intent(inout) ::  vz(Xh%lx, Xh%ly, Xh%lz, nelv)
     end subroutine opr_sx_conv1

     module subroutine opr_sx_convect_scalar(du, u, c, Xh_GLL, Xh_GL, &
                                             coef_GLL, coef_GL, GLL_to_GL)
       type(space_t), intent(in) :: Xh_GL
       type(space_t), intent(in) :: Xh_GLL
       type(coef_t), intent(in) :: coef_GLL
       type(coef_t), intent(in) :: coef_GL
       type(interpolator_t), intent(inout) :: GLL_to_GL
       real(kind=rp), intent(inout) :: &
                      du(Xh_GLL%lx, Xh_GLL%ly, Xh_GLL%lz, coef_GL%msh%nelv)
       real(kind=rp), intent(inout) :: &
                      u(Xh_GL%lx, Xh_GL%lx, Xh_GL%lx, coef_GL%msh%nelv)
       real(kind=rp), intent(inout) :: c(Xh_GL%lxyz, coef_GL%msh%nelv, 3)

     end subroutine opr_sx_convect_scalar

     module function opr_sx_cfl(dt, u, v, w, Xh, coef, nelv) result(cfl)
       type(space_t), intent(in) :: Xh
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: nelv
       real(kind=rp), intent(in) :: dt
       real(kind=rp), dimension(Xh%lx, Xh%ly, Xh%lz, nelv) ::  u, v, w
       real(kind=rp) :: cfl
     end function opr_sx_cfl

     module subroutine opr_sx_lambda2(lambda2, u, v, w, coef)
       type(coef_t), intent(in) :: coef
       type(field_t), intent(inout) :: lambda2
       type(field_t), intent(in) :: u, v, w
     end subroutine opr_sx_lambda2

     module subroutine opr_sx_set_convect_rst(cr, cs, ct, cx, cy, cz, Xh, coef)
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                      intent(inout) :: cr, cs, ct
       real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                      intent(in) :: cx, cy, cz
     end subroutine opr_sx_set_convect_rst
  end interface

contains

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

    call opcolv(w1%x, w2%x, w3%x, c_Xh%B, gdim, n)
    call c_Xh%gs_h%op(w1, GS_OP_ADD)
    call c_Xh%gs_h%op(w2, GS_OP_ADD)
    call c_Xh%gs_h%op(w3, GS_OP_ADD)
    call opcolv(w1%x, w2%x, w3%x, c_Xh%Binv, gdim, n)

  end subroutine opr_sx_curl





end module opr_sx
