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
  implicit none

  type, public, abstract :: advection_t
  contains
     procedure(apply_adv), pass(this), deferred :: apply
  end type advection_t

  type, public, extends(advection_t) :: adv_no_dealias_t
  contains
    procedure, pass(this) :: apply => advab
  end type adv_no_dealias_t

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
  subroutine advection_factory(this)
     class(advection_t), allocatable, intent(inout) :: this
     
    if (allocated(this)) then
       deallocate(this)
    end if

    allocate(adv_no_dealias_t::this)

  end subroutine advection_factory
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
