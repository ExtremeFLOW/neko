module user
  use neko
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => compute
  end subroutine user_setup

  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type (field_t), pointer :: u, v, w, s
    integer :: i
    type(field_list_t) :: list
    type(file_t) :: file

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    call list%init(3)
    call list%assign(1, u)
    call list%assign(2, v)
    call list%assign(3, w)

    call file%init("test.fld")

    call file%write(list)
    call file%free()

  end subroutine compute
end module user
