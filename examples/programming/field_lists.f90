! This tutorial demonstrates how to work with field lists.

module user
  use neko
  implicit none


  ! We declare some fields that we will pack into a list.
  type(field_t) :: my_field1, my_field2

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup
    user%initialize => initialize
    user%finalize => finalize
  end subroutine user_setup

  ! We will use the user_startup routine to manipulate the end time.
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call params%add("case.end_time", 0.0_rp)
  end subroutine startup

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(field_list_t) :: field_list

    ! Initialize the field so that we have something to work with.
    call my_field1%init(neko_user_access%case%fluid%dm_Xh, "my_field1")
    call my_field2%init(neko_user_access%case%fluid%dm_Xh, "my_field2")

    ! To initialize the field list, we provide some initial size for it.
    call field_list%init(2)


  end subroutine initialize


  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    call my_field1%free()
    call my_field2%free()
  end subroutine finalize

end module user
