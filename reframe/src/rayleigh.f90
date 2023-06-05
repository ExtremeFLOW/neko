module user
  use neko
  implicit none

  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Pr = 0
  real(kind=rp) :: ta2 = 0

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => set_Pr
    u%fluid_user_ic => set_ic
    u%fluid_user_f_vector => forcing
  end subroutine user_setup

  !> Dummy user initial condition
  subroutine set_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i
    real(kind=rp) :: rand
    s => neko_field_registry%get_field('s')

    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())
    
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1) + 0.1*sin(4*pi/4.5*s%dof%x(i,1,1,1)) &
                 * sin(4*pi/4.5*s%dof%y(i,1,1,1))*sqrt((0.5**2-(0.5_rp-s%dof%z(i,1,1,1))**2))
    end do
    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x,s%x_d,s%dof%size(),HOST_TO_DEVICE)
    end if


  end subroutine set_ic

  subroutine set_Pr(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    ! Reset the relevant nondimensional parameters
    ! Pr = input Pr
    ! Ra = input Re
    ! Re = 1/Pr
    Pr = params%Pr
    Ra = params%Re
    params%Re = 1._rp / Pr
  end subroutine set_Pr



  !> Forcing
  subroutine forcing(f, t)
    class(source_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')
    rapr = Ra*Pr
    ta2pr = ta2*Pr

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_cmult2(f%u_d,v%x_d,Ta2Pr,f%dm%size())
       call device_cmult2(f%v_d,u%x_d,Ta2Pr,f%dm%size())
       call device_cmult2(f%w_d,s%x_d,rapr,f%dm%size())
    else
       call cmult2(f%u,v%x,Ta2Pr,f%dm%size())
       call cmult2(f%v,u%x,Ta2Pr,f%dm%size())
       call cmult2(f%w,s%x,rapr,f%dm%size())
    end if
  end subroutine forcing
end module user
