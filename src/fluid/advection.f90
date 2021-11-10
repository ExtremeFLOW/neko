!> Adams-Bashforth coefs for Backward Differentiation schemes
module advection
  use num_types
  use math
  use utils
  use space
  use field
  use coefs
  use neko_config
  use operators
  use interpolation
  implicit none

  type, public, abstract :: advection_t
  contains
     procedure(apply_adv), pass(this), deferred :: apply
  end type advection_t

  type, public, extends(advection_t) :: adv_no_dealias_t
  contains
    procedure, pass(this) :: apply => advab
  end type adv_no_dealias_t

  type, public, extends(advection_t) :: adv_dealias_t
    type(coef_t) :: coef_GL
    type(coef_t), pointer :: coef_GLL
    type(interpolator_t) :: GLL_to_GL
    type(space_t) :: Xh_GL
    type(space_t), pointer :: Xh_GLL
  contains
    procedure, pass(this) :: apply => apply_adv_dealias
    procedure, pass(this) :: init => init_dealias
  end type adv_dealias_t


  abstract interface
     subroutine apply_adv(this, ta1, ta2, ta3, vx, vy, vz, bfx, bfy, bfz, Xh, coef, nelv, n, gdim)
        import :: advection_t
        import :: coef_t
        import :: space_t
        import :: field_t
        import :: rp
        class(advection_t), intent(inout) :: this
        type(space_t), intent(inout) :: Xh
        type(coef_t), intent(inout) :: coef
        type(field_t), intent(inout) :: ta1, ta2, ta3, vx, vy, vz
        integer, intent(inout) :: nelv, n, gdim
        real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz
     end subroutine apply_adv
  end interface

contains
  subroutine advection_factory(this, coef, dealias, lxd)
    class(advection_t), allocatable, intent(inout) :: this
    type(coef_t) :: coef
    logical, intent(in) :: dealias
    integer, intent(in), optional :: lxd
     
    if (allocated(this)) then
       deallocate(this)
    end if

    if (dealias) then
       allocate(adv_dealias_t::this)
    else
       allocate(adv_no_dealias_t::this)
    end if
    select type(adv => this)
       type is(adv_dealias_t)
          if (present(lxd)) then
             call init_dealias(adv,lxd,coef) 
          else
             call init_dealias(adv,coef%Xh%lx*3/2, coef)
          end if
       end select

  end subroutine advection_factory

  subroutine init_dealias(this, lxd, coef)
    class(adv_dealias_t), intent(inout) :: this
    integer, intent(in) :: lxd
    type(coef_t), intent(inout), target :: coef
    integer :: nel
     
    call space_init(this%Xh_GL,GL, lxd, lxd, lxd)
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef
    call this%GLL_to_GL%init(this%Xh_GL, this%Xh_GLL)

    call coef_empty_init(this%coef_GL,this%Xh_GL,coef%msh)

    nel = coef%msh%nelv
    call this%GLL_to_GL%map(this%coef_GL%drdx, coef%drdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdx, coef%dsdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdx, coef%dtdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdy, coef%drdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdy, coef%dsdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdy, coef%dtdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdz, coef%drdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdz, coef%dsdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdz, coef%dtdz, nel, this%Xh_GL)

  end subroutine init_dealias
  !> Eulerian scheme, add convection term to forcing function
  !! at current time step.
  subroutine apply_adv_dealias(this, ta1, ta2, ta3, vx, vy, vz, bfx, bfy, bfz, Xh, coef, nelv, n, gdim)
    class(adv_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: ta1, ta2, ta3, vx, vy, vz
    integer, intent(inout) :: nelv, n, gdim
    real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tx, ty, tz, tr
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: ts, tt, tbfx, tbfy 
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tbfz, vr, vs, vt
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz
    integer :: e, i, idx, NEKO_NEL_SIZE = 1
    associate(c_GL => this%coef_GL)
    do e = 1, nelv, NEKO_NEL_SIZE
       call this%GLL_to_GL%map(tx, vx%x(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(ty, vy%x(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(tz, vz%x(1,1,1,e), 1, this%Xh_GL)
       do i = 1, this%Xh_GL%lxyz
          tr(i) = c_GL%drdx(i,1,1,e) * tx(i) + c_GL%drdy(i,1,1,e)*ty(i) + c_GL%drdz(i,1,1,e)*tz(i)
          ts(i) = c_GL%dsdx(i,1,1,e) * tx(i) + c_GL%dsdy(i,1,1,e)*ty(i) + c_GL%dsdz(i,1,1,e)*tz(i)
          tt(i) = c_GL%dtdx(i,1,1,e) * tx(i) + c_GL%dtdy(i,1,1,e)*ty(i) + c_GL%dtdz(i,1,1,e)*tz(i)
       end do
       call opgrad(vr, vs, vt, tx, c_GL, e, e)
       do i = 1, this%Xh_GL%lxyz
          tbfx(i) = tr(i)*vr(i) + ts(i)*vs(i) + tt(i)*vt(i)
       end do
       call opgrad(vr, vs, vt, ty, c_GL, e, e)
       do i = 1, this%Xh_GL%lxyz
          tbfy(i) = tr(i)*vr(i) + ts(i)*vs(i) + tt(i)*vt(i)
       end do
       call opgrad(vr, vs, vt, tz, c_GL, e, e)
       do i = 1, this%Xh_GL%lxyz
          tbfz(i) = tr(i)*vr(i) + ts(i)*vs(i) + tt(i)*vt(i)
       end do
       idx = (e-1)*this%Xh_GLL%lxyz+1
       call this%GLL_to_GL%map(tempx, tbfx, 1, this%Xh_GLL)
       call this%GLL_to_GL%map(tempy, tbfy, 1, this%Xh_GLL)
       call this%GLL_to_GL%map(tempz, tbfz, 1, this%Xh_GLL)
       call subcol3(bfx(idx),this%coef_GLL%B(1,1,1,e), tempx,this%Xh_GLL%lxyz)
       call subcol3(bfy(idx),this%coef_GLL%B(1,1,1,e), tempy,this%Xh_GLL%lxyz)
       call subcol3(bfz(idx),this%coef_GLL%B(1,1,1,e), tempz,this%Xh_GLL%lxyz)

    end do
    end associate

  end subroutine apply_adv_dealias



  !> Eulerian scheme, add convection term to forcing function
  !! at current time step.
  subroutine advab(this, ta1, ta2, ta3, vx, vy, vz, bfx, bfy, bfz, Xh, coef, nelv, n, gdim)
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: ta1, ta2, ta3, vx, vy, vz
    integer, intent(inout) :: nelv, n, gdim
    real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz

    call conv1(ta1%x, vx%x, vx%x, vy%x, vz%x, Xh, coef, nelv, gdim)
    call conv1(ta2%x, vy%x, vx%x, vy%x, vz%x, Xh, coef, nelv, gdim)
    call subcol3 (bfx, coef%B, ta1%x, n)
    call subcol3 (bfy, coef%B, ta2%x, n)
    if (gdim .eq. 2) then
       call rzero (ta3%x, n)
    else
       call conv1(ta3%x, vz%x, vx%x, vy%x, vz%x, Xh, coef, nelv, gdim)
       call subcol3(bfz, coef%B, ta3%x, n)
    end if
  end subroutine advab
  
end module advection
