!> Initial flow condition
module flow_ic
  use gather_scatter
  use neko_config
  use flow_profile
  use device_math
  use parameters
  use field
  use utils
  use coefs
  use math
  use user_intf, only : useric
  implicit none
  private

  interface set_flow_ic
     module procedure set_flow_ic_int, set_flow_ic_usr
  end interface

  public :: set_flow_ic
  
contains

  !> Set initial flow condition (builtin)
  subroutine set_flow_ic_int(u, v, w, p, coef, gs, type, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
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
    
    call set_flow_ic_common(u, v, w, p, coef, gs)
    
  end subroutine set_flow_ic_int

  !> Set intial flow condition (user defined)
  subroutine set_flow_ic_usr(u, v, w, p, coef, gs, usr_ic, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    procedure(useric) :: usr_ic
    type(param_t), intent(inout) :: params

    call usr_ic(u, v, w, p, params)
    
    call set_flow_ic_common(u, v, w, p, coef, gs)
    
  end subroutine set_flow_ic_usr

  subroutine set_flow_ic_common(u, v, w, p, coef, gs)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(u%x, u%x_d, u%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(v%x, v%x_d, v%dof%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(w%x, w%x_d, w%dof%n_dofs, HOST_TO_DEVICE)
    end if
    
    ! Ensure continuity across elements for initial conditions
    call gs_op_vector(gs, u%x, u%dof%n_dofs, GS_OP_ADD) 
    call gs_op_vector(gs, v%x, v%dof%n_dofs, GS_OP_ADD) 
    call gs_op_vector(gs, w%x, w%dof%n_dofs, GS_OP_ADD) 

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_col2(u%x_d, coef%mult_d, u%dof%n_dofs)
       call device_col2(v%x_d, coef%mult_d, v%dof%n_dofs)
       call device_col2(w%x_d, coef%mult_d, w%dof%n_dofs)
    else
       call col2(u%x, coef%mult, u%dof%n_dofs)
       call col2(v%x, coef%mult, v%dof%n_dofs)
       call col2(w%x, coef%mult, w%dof%n_dofs)
    end if
    
  end subroutine set_flow_ic_common

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
    real(kind=rp), intent(in) :: delta
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
