module user
  use neko
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%scalar_user_ic => set_s_ic
    u%fluid_user_ic => set_velocity
    u%scalar_user_f_vector => set_source
    u%user_mesh_setup => user_mesh_scale
  end subroutine user_setup

  ! Stretch bounds to 2pi
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = msh%points(i)%x(1)* pi *2
       msh%points(i)%x(2) = msh%points(i)%x(2)* pi *2
    end do
  end subroutine user_mesh_scale

  !> Set source term
  subroutine set_source(f, t)
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

  end subroutine set_source

  !> User initial condition for the scalar
  subroutine set_s_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y

    do i = 1, s%dof%size()
       x = s%dof%x(i,1,1,1)
       y = s%dof%y(i,1,1,1)

       s%x(i,1,1,1) = sin(x)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
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

      ! Simple advection in x
       u%x(i,1,1,1) = 1
       v%x(i,1,1,1) = 0
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
