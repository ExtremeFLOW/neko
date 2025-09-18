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
    user%initial_conditions => initial_conditions
    user%source_term => source_term
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

  !> User initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    integer :: i, e, k, j
    real(kind=rp) :: rand, z
    type(field_t), pointer :: s

    ! See scalar.name in the case file, makes sure that we only
    ! run this for the scalar field.
    if (scheme_name .ne. 'temperature') return

    s => fields%items(1)%ptr

    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1 - s%dof%z(i,1,1,1)
    end do

    ! perturb not on element boundaries
    ! Maybe not necessary, but lets be safe
    do e = 1, s%msh%nelv
       do k = 2, s%Xh%lx-1
          do j = 2, s%Xh%lx-1
             do i = 2, s%Xh%lx-1

                !call random_number(rand)
                !Somewhat random
                rand = cos(real(e + s%msh%offset_el, rp) * real(i*j*k, rp))
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 1 - z + 0.0001*rand* &
                     sin(4*pi/4.5 * s%dof%x(i,j,k,e)) &
                     * sin(4*pi/4.5 * s%dof%y(i,j,k,e))

             end do
          end do
       end do
    end do

  end subroutine initial_conditions

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
