module user
  use neko
  implicit none

  logical :: skip_pressure = .false.
  logical :: skip_temperature = .false.
  logical :: skip_velocity = .false.
  logical :: is_list = .false.
  integer :: n_items = 0

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => compute
    user%startup => startup
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, "skip_pressure", skip_pressure)
    call json_get(params, "skip_temperature", skip_temperature)
    call json_get(params, "skip_velocity", skip_velocity)
    call json_get(params, "is_list", is_list)
    call json_get(params, "n_items", n_items)
  end subroutine startup

  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type (field_t), pointer :: u, v, w, s
    integer :: i
    type(field_list_t) :: list
    type(file_t) :: file

    u => neko_registry%get_field("u")

    call list%init(n_items)
    do i = 1 , n_items
       call list%assign(i, u)
    end do

    call file%init("./tests/test_fld_file/test.fld")

    select type (ft => file%file_type)
    type is (fld_file_t)
       ft%skip_pressure = skip_pressure
       ft%skip_temperature = skip_temperature
       ft%skip_velocity = skip_velocity
    end select

    if (is_list) then
       call file%write(list)
    else
       call file%write(u)
    end if
    call file%free()
  end subroutine compute
end module user
