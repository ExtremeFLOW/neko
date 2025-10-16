module user
  use neko
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%source_term => source_term
    user%mesh_setup => user_mesh_scale
  end subroutine user_setup

  ! Stretch bounds to 2pi
  subroutine user_mesh_scale(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time
    integer :: i, p, nvert

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = msh%points(i)%x(1)* pi *2
       msh%points(i)%x(2) = msh%points(i)%x(2)* pi *2
    end do
  end subroutine user_mesh_scale

  !> Set source term
  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    real(kind=rp) :: x
    type(field_t), pointer :: f
    integer :: i

    if (scheme_name .eq. 'fluid') return

    f => rhs%items(1)%ptr
    do i = 1, f%size()
       x = f%dof%x(i,1,1,1)

       ! 0.01 is the viscosity
       f%x(i,1,1,1) = cos(x) - 0.01 * sin(x) - 1.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f%x, f%x_d, f%size(), &
            HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine source_term

  !> User initial condition for the scalar
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    integer :: i, e, k, j
    real(kind=rp) :: x, y
    type (field_t), pointer :: u, v, w, s
    type(dofmap_t), pointer :: dof

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       call field_rone(u)
       call field_rzero(v)
       call field_rzero(w)
    else !scalar
       s => fields%get("s")
       do i = 1, s%dof%size()
          x = s%dof%x(i,1,1,1)
          y = s%dof%y(i,1,1,1)
          s%x(i,1,1,1) = sin(x)
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

    do i = 1, u%dof%size()

       ! Simple advection in x
       u%x(i,1,1,1) = 1
       v%x(i,1,1,1) = 0
       w%x(i,1,1,1) = 0
    end do

  end subroutine set_velocity

end module user
