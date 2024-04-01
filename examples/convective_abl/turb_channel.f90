! Martin Karp 13/3-2023
module user
  use neko
  implicit none

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_ic => fluid_ic
    u%scalar_user_ic => scalar_ic
    u%fluid_user_f_vector => set_bousinesq_forcing_term
  end subroutine user_setup

  ! User defined initial condition
  subroutine fluid_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i
    real(kind=rp) :: uvw(3)

    do i = 1, u%dof%size()
       u%x(i,1,1,1) = 10_rp
       v%x(i,1,1,1) = 0_rp
       w%x(i,1,1,1) = 0_rp
    end do
  end subroutine fluid_ic

  subroutine scalar_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y, z

    do i = 1, s%dof%size()
       y = s%dof%y(i,1,1,1)

       if (y < 937) then
          s%x(i,1,1,1) = 300
       else if (y > 1063) then
          s%x(i,1,1,1) = 308 + 3e-3_rp*(y - 1063)
       else
          s%x(i,1,1,1) = 300 + 8.0_rp/126.0_rp*(y - 937)
       end if

    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if
  end subroutine scalar_ic

  subroutine set_bousinesq_forcing_term(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: g, tref

    g = 9.81
    tref = 300

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_rzero(f%w_d,f%dm%size())
    else
       call rzero(f%u, f%dm%size())
       call copy(f%v, g*(s%x - tref)/tref, f%dm%size())
       call rzero(f%w, f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term
end module user
