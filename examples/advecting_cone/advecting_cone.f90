module user
  use neko
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  !> User initial condition for the scalar
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    integer :: i, e, k, j
    real(kind=rp) :: cone_radius, mux, muy, x, y, r, theta

    type(dofmap_t), pointer :: dof
    type (field_t), pointer :: u, v, w, s

    dof => fields%dof(1)
    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       do i = 1, u%dof%size()
          x = u%dof%x(i,1,1,1)
          y = u%dof%y(i,1,1,1)

          ! Angular velocity is pi, giving a full rotation in 2 sec
          u%x(i,1,1,1) = -y*pi
          v%x(i,1,1,1) = x*pi
          w%x(i,1,1,1) = 0
       end do
    else
       s => fields%get("s")
       ! Center of the cone
       mux = 1
       muy = 0

       cone_radius = 0.5

       do i = 1, s%dof%size()
          x = dof%x(i,1,1,1) - mux
          y = dof%y(i,1,1,1) - muy

          r = sqrt(x**2 + y**2)
          theta = atan2(y, x)

          ! Check if the point is inside the cone's base
          if (r > cone_radius) then
             s%x(i,1,1,1) = 0.0
          else
             s%x(i,1,1,1) = 1.0 - r / cone_radius
          end if
       end do
    end if

  end subroutine initial_conditions

  !> Set the advecting velocity field.
  subroutine set_velocity(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y


  end subroutine set_velocity

end module user
