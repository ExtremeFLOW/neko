module user
  use neko
  implicit none

  real(kind=rp) :: Ra
  real(kind=rp) :: Pr
  real(kind=rp) :: ta2 = 0

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%scalar_user_ic => set_ic
    user%source_term => source_term
    user%scalar_user_bc => scalar_bc
    user%startup => startup
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    real(kind=rp) :: rho, mu, cp, lambda, Re

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.scalar.Pr", Pr)

    Re = 1.0_rp / Pr
    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp

    call params%add("case.fluid.mu", mu)
    call params%add("case.fluid.rho", rho)
    call params%add("case.scalar.lambda", lambda)
    call params%add("case.scalar.cp", cp)
  end subroutine startup

  subroutine scalar_bc(scalar_name, s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    character(len=*), intent(in) :: scalar_name
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    ! If we set scalar_bcs(*) = 'user' instead
    ! this will be used instead on that zone
    s = 1.0_rp-z
  end subroutine scalar_bc

  !> User initial condition
  subroutine set_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: rand, z

    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1)
    end do
    ! perturb not on element boundaries
    ! Maybe not necessary, but lets be safe
    do e = 1, s%msh%nelv
       do k = 2,s%Xh%lx-1
          do j = 2,s%Xh%lx-1
             do i = 2,s%Xh%lx-1

                !call random_number(rand)
                !Somewhat random
                rand = cos(real(e+s%msh%offset_el,rp)*real(i*j*k,rp))
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1-z + 0.0001* rand*&
                     sin(4*pi/4.5*s%dof%x(i,j,k,e)) &
                     * sin(4*pi/4.5*s%dof%y(i,j,k,e))

             end do
          end do
       end do
    end do

    if ((NEKO_BCKND_DEVICE .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
            HOST_TO_DEVICE, sync=.false.)
    end if


  end subroutine set_ic

  !> Forcing
  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    integer :: i
    type(field_t), pointer :: u, v, w, s, rhs_u, rhs_v, rhs_w
    real(kind=rp) :: rapr, ta2pr

    if (scheme_name .eq. 'fluid') then
       u => neko_field_registry%get_field('u')
       v => neko_field_registry%get_field('v')
       w => neko_field_registry%get_field('w')
       s => neko_field_registry%get_field('temperature')

       rhs_u => rhs%get_by_index(1)
       rhs_v => rhs%get_by_index(2)
       rhs_w => rhs%get_by_index(3)

       rapr = Ra*Pr
       ta2pr = ta2*Pr

       call field_cmult2(rhs_u, v, Ta2Pr)
       call field_cmult2(rhs_v, u, Ta2Pr)
       call field_cmult2(rhs_w, s, rapr)
    end if
  end subroutine source_term
end module user
