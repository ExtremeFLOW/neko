module user
  use neko
  implicit none

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%initial_conditions => initial_conditions
    user%initialize => initialize
    user%compute => report_solution
  end subroutine user_setup

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, p, s

    if (scheme_name .eq. "fluid") then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")
       p => fields%get("p")
       call field_rzero(u)
       call field_rzero(v)
       call field_rzero(w)
       call field_rzero(p)
    else
       s => fields%get("s")
       call field_cfill(s, 1.0_rp)
    end if
  end subroutine initial_conditions

  !> The source is a s + b with a = -90 and b = 2.
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: a, b

    a => neko_registry%get_field("scalar_source_a")
    b => neko_registry%get_field("scalar_source_b")

    call field_cfill(a, -90.0_rp)
    call field_cfill(b, 2.0_rp)
  end subroutine initialize

  subroutine report_solution(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: s

    s => neko_registry%get_field("s")
    if (pe_rank .eq. 0) then
       write(*, '(A,1X,ES24.16,1X,ES24.16)') "implicit_source_scalar", &
            field_glmin(s), field_glmax(s)
    end if
  end subroutine report_solution

end module user
