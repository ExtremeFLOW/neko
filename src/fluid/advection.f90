!> Subroutines to apply advection to RHS
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
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public, abstract :: advection_t
   contains
     procedure(apply_adv), pass(this), deferred :: apply
  end type advection_t

  type, public, extends(advection_t) :: adv_no_dealias_t
     real(kind=rp), allocatable :: temp(:)
     type(c_ptr) :: temp_d = C_NULL_PTR
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
     subroutine apply_adv(this, vx, vy, vz, bfx, bfy, bfz, Xh, coef, n)
       import :: advection_t
       import :: coef_t
       import :: space_t
       import :: field_t
       import :: rp
       class(advection_t), intent(inout) :: this
       type(space_t), intent(inout) :: Xh
       type(coef_t), intent(inout) :: coef
       type(field_t), intent(inout) :: vx, vy, vz
       integer, intent(inout) :: n
       real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz
     end subroutine apply_adv
  end interface

  public :: advection_factory

contains
  
  subroutine advection_factory(this, coef, dealias, lxd)
    class(advection_t), allocatable, intent(inout) :: this
    type(coef_t) :: coef
    logical, intent(in) :: dealias
    integer, intent(in) :: lxd

    if (allocated(this)) then
       select type(adv => this)
       type is(adv_no_dealias_t)
          if (allocated(adv%temp)) then
             deallocate(adv%temp)
          end if
          if (c_associated(adv%temp_d)) then
             call device_free(adv%temp_d)
          end if
       end select
       deallocate(this)
    end if

    if (dealias) then
       allocate(adv_dealias_t::this)
    else
       allocate(adv_no_dealias_t::this)
    end if
    
    select type(adv => this)
    type is(adv_dealias_t)
       if (lxd .gt. 0) then
          call init_dealias(adv,lxd,coef) 
       else
          call init_dealias(adv,coef%Xh%lx*3/2, coef)
       end if
    type is(adv_no_dealias_t)
       call init_no_dealias(adv,coef)
    end select

  end subroutine advection_factory

  subroutine init_no_dealias(this, coef)
    class(adv_no_dealias_t) :: this
    type(coef_t) :: coef

    allocate(this%temp(coef%dof%n_dofs))

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_map(this%temp, this%temp_d, coef%dof%n_dofs)
    end if

  end subroutine init_no_dealias

  subroutine init_dealias(this, lxd, coef)
    class(adv_dealias_t), intent(inout) :: this
    integer, intent(in) :: lxd
    type(coef_t), intent(inout), target :: coef
    integer :: nel

    call space_init(this%Xh_GL,GL, lxd, lxd, lxd)
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef
    call this%GLL_to_GL%init(this%Xh_GL, this%Xh_GLL)

    call coef_init(this%coef_GL, this%Xh_GL,coef%msh)

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
  subroutine apply_adv_dealias(this, vx, vy, vz, bfx, bfy, bfz, Xh, coef, n)
    class(adv_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    integer, intent(inout) :: n
    real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tx, ty, tz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tbfx, tbfy, tbfz 
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vr, vs, vt
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz
    integer :: e, i, idx
    associate(c_GL => this%coef_GL)

      do e = 1, coef%msh%nelv
         call this%GLL_to_GL%map(tx, vx%x(1,1,1,e), 1, this%Xh_GL)
         call this%GLL_to_GL%map(ty, vy%x(1,1,1,e), 1, this%Xh_GL)
         call this%GLL_to_GL%map(tz, vz%x(1,1,1,e), 1, this%Xh_GL)

         call opgrad(vr, vs, vt, tx, c_GL, e, e)
         do i = 1, this%Xh_GL%lxyz
            tbfx(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
         end do

         call opgrad(vr, vs, vt, ty, c_GL, e, e)
         do i = 1, this%Xh_GL%lxyz
            tbfy(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
         end do

         call opgrad(vr, vs, vt, tz, c_GL, e, e)
         do i = 1, this%Xh_GL%lxyz
            tbfz(i) = tx(i)*vr(i) + ty(i)*vs(i) + tz(i)*vt(i)
         end do

         call this%GLL_to_GL%map(tempx, tbfx, 1, this%Xh_GLL)
         call this%GLL_to_GL%map(tempy, tbfy, 1, this%Xh_GLL)
         call this%GLL_to_GL%map(tempz, tbfz, 1, this%Xh_GLL)

         idx = (e-1)*this%Xh_GLL%lxyz+1
         call sub2(bfx(idx), tempx,this%Xh_GLL%lxyz)
         call sub2(bfy(idx), tempy,this%Xh_GLL%lxyz)
         call sub2(bfz(idx), tempz,this%Xh_GLL%lxyz)
      end do
    end associate

  end subroutine apply_adv_dealias



  !> Eulerian scheme, add convection term to forcing function
  !! at current time step.
  subroutine advab(this, vx, vy, vz, bfx, bfy, bfz, Xh, coef, n)
    class(adv_no_dealias_t), intent(inout) :: this
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    integer, intent(inout) :: n
    real(kind=rp), intent(inout), dimension(n) :: bfx, bfy, bfz

    call conv1(this%temp, vx%x, vx%x, vy%x, vz%x, Xh, coef)
    call subcol3 (bfx, coef%B, this%temp, n)
    call conv1(this%temp, vy%x, vx%x, vy%x, vz%x, Xh, coef)
    call subcol3 (bfy, coef%B, this%temp, n)
    if (coef%Xh%lz .eq. 1) then
       call rzero (this%temp, n)
    else
       call conv1(this%temp, vz%x, vx%x, vy%x, vz%x, Xh, coef)
       call subcol3(bfz, coef%B, this%temp, n)
    end if
  end subroutine advab

end module advection
