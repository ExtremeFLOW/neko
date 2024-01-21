module user
  use neko
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%scalar_user_ic => set_s_ic
    u%fluid_user_ic => set_velocity
  end subroutine user_setup

  !> User initial condition for the scalar
  subroutine set_s_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: cone_radius, mux, muy, x, y, r, theta

    ! Center of the cone
    mux = 1
    muy = 0

    cone_radius = 0.5

    do i = 1, s%dof%size()
       x = s%dof%x(i,1,1,1) - mux
       y = s%dof%y(i,1,1,1) - muy

       r = sqrt(x**2 + y**2)
       theta = atan2(y, x)

       ! Check if the point is inside the cone's base
       if (r > cone_radius) then
         s%x(i,1,1,1) = 0.0
       else
         s%x(i,1,1,1) = 1.0 - r / cone_radius
       endif
    end do

    if ((NEKO_BCKND_DEVICE .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if
  end subroutine set_s_ic

  !> Set the advecting velocity field.
  subroutine set_velocity(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y

    do i = 1, u%dof%size()
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)

       ! Angular velocity is pi, giving a full rotation in 2 sec
       u%x(i,1,1,1) = -y*pi
       v%x(i,1,1,1) = x*pi
       w%x(i,1,1,1) = 0
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if
  end subroutine set_velocity

end module user
