module user
  use neko
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%scalar_user_ic => set_scalars_ic
    u%fluid_user_ic => set_velocity
    u%scalar_user_f_vector => set_source_vector
    u%scalar_user_f => set_source
  end subroutine user_setup

  !> User initial condition for the scalars
  subroutine set_scalars_ic(field_name, s, params)
    CHARACTER(len=*), INTENT(IN) :: field_name
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

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if
  end subroutine set_scalars_ic

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

  !> Set source term vector
  subroutine set_source_vector(field_name, f, t)
    character(len=*), intent(in) :: field_name
    class(scalar_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: x, y
    integer :: i

    do i = 1, f%dm%size()
       x = f%dm%x(i,1,1,1)
       y = f%dm%y(i,1,1,1)

       ! 0.01 is the viscosity
       f%s(i,1,1,1) = cos(x) - 0.01 * sin(x) - 1.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f%s, f%s_d, f%dm%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_source_vector

  !> Set source term vector
  subroutine set_source(field_name, s, j, k, l, e, t)
    CHARACTER(len=*), INTENT(IN) :: field_name
    real(kind=rp), intent(inout) :: s
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, intent(in) :: e
    real(kind=rp), intent(in) :: t

    if (field_name == "s1") then
      s = 0.01_rp
    else
      s = 0.0_rp
    end if

  end subroutine set_source

end module user
