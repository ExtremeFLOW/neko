module user
  use neko
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initialize => user_initialize
  end subroutine user_setup

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t)
    type(time_state_t), intent(in) :: t

    type(field_t), pointer :: u, fringe
    real(kind=rp) :: x, y, xmin1, delta_rise1, xmin2, delta_rise2
    integer :: i, imask

    fringe => null()
    u => null()

    !
    ! 1. Add the "sponge_field" to the field registry.
    !    NOTE: The name of the fringe field in the registry
    !    can be changed with the parameter `fringe_registry_name`.
    !
    !
    u => neko_field_registry%get_field("u")
    call neko_field_registry%add_field(u%dof,"sponge_fringe")
    fringe => neko_field_registry%get_field("sponge_fringe")

    !
    ! 2. Set the function f(x,y,z) from 0 to 1.
    !

    ! Bottom boundary
    xmin1 = 3.0_rp
    delta_rise1 = 3.0_rp

    ! Top boundary
    xmin2 = 20.0_rp
    delta_rise2 = 7.0_rp

    fringe = 0.0_rp
    do i = 1, fringe%size()
       x = fringe%dof%x(i,1,1,1)
       y = fringe%dof%y(i,1,1,1)

       ! Bottom boundary
       if ( (y .lt. 0.0_rp) .and. (x .gt. xmin1)) then
          fringe%x(i,1,1,1) = S( (x - xmin1)/delta_rise1 )

          ! Top boundary
       else if ( (y .gt. 0.0_rp) .and. (x .gt. xmin2)) then
          fringe%x(i,1,1,1) = S( (x - xmin2)/delta_rise2 )

       end if
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(fringe%x, fringe%x_d, fringe%size(), &
            HOST_TO_DEVICE, .true.)
    end if

    ! NOTE: You can dump the fringe field to file using the `dump_fields`
    ! parameter. The fringe field will be stored under `pressure`.

    nullify(fringe)
    nullify(u)

  end subroutine user_initialize

  ! Smooth step function, 0 if x <= 0, 1 if x >= 1, 1/exp(1/(x-1) + 1/x) between 0 and 1
  function S(x) result(y)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: y

    if ( x .le. 0.0_rp ) then
       y = 0.0_rp
    else if ( x .ge. 1.0_rp ) then
       y = 1.0_rp
    else
       y = 1.0_rp / (1.0_rp + exp( 1.0_rp / (x - 1.0_rp) + 1.0_rp / x))
    end if
  end function S

end module user
