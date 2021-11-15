!> Initial flow condition
module flow_ic
  use flow_profile
  use parameters
  use field
  use utils
  implicit none
  private

  public :: set_flow_ic
contains

  !> Set initial flow condition
  subroutine set_flow_ic(u, v, w, p, type, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    character(len=*) :: type
    type(param_t), intent(in) :: params

    if (trim(type) .eq. 'uniform') then       
       call set_flow_ic_uniform(u, v, w, params%uinf)
    else if (trim(type) .eq. 'blasius') then
       call set_flow_ic_blasius(u, v, w, &
            params%delta, params%uinf, params%blasius_approx)
    else
       call neko_error('Invalid initial condition')
    end if
    
  end subroutine set_flow_ic

  !> Uniform initial condition
  subroutine set_flow_ic_uniform(u, v, w, uinf)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    real(kind=rp), intent(in) :: uinf(3)
    u = uinf(1)
    v = uinf(2)
    w = uinf(3)
  end subroutine set_flow_ic_uniform

  !> Set a Blasius profile as initial condition
  !! @note currently limited to axis aligned flow
  subroutine set_flow_ic_blasius(u, v, w, delta, uinf, type)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    real(kind=dp), intent(in) :: delta
    real(kind=rp), intent(in) :: uinf(3)
    character(len=*), intent(in) :: type
    procedure(blasius_profile), pointer :: bla => null()
    integer :: i

    select case(trim(type))
    case('linear')
       bla => blasius_linear
    case('quadratic')
       bla => blasius_quadratic
    case('cubic')
       bla => blasius_cubic
    case('quartic')
       bla => blasius_quartic
    case('sin')
       bla => blasius_sin
    case default
       call neko_error('Invalid Blasius approximation')
    end select
    
    if ((uinf(1) .gt. 0.0_rp) .and. (uinf(2) .eq. 0.0_rp) &
         .and. (uinf(3) .eq. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = bla(u%dof%z(i,1,1,1), delta, uinf(1))
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .eq. 0.0_rp) .and. (uinf(2) .gt. 0.0_rp) &
         .and. (uinf(3) .eq. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = 0.0_rp
          v%x(i,1,1,1) = bla(u%dof%x(i,1,1,1), delta, uinf(2))
          w%x(i,1,1,1) = 0.0_rp
       end do
    else if ((uinf(1) .eq. 0.0_rp) .and. (uinf(2) .eq. 0.0_rp) &
         .and. (uinf(3) .gt. 0.0_rp)) then
       do i = 1, u%dof%size()
          u%x(i,1,1,1) = 0.0_rp
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = bla(u%dof%y(i,1,1,1), delta, uinf(3))
       end do       
    end if
    
  end subroutine set_flow_ic_blasius
  
end module flow_ic
