!> Defines an output for a mean flow field
module mean_flow_output
  use mean_flow
  use num_types
  use output
  implicit none

  type, extends(output_t) :: mean_flow_output_t
     type(mean_flow_t), pointer :: mf
     real(kind=rp) :: T_begin
   contains
     procedure, pass(this) :: sample => mean_flow_output_sample
  end type mean_flow_output_t

  interface mean_flow_output_t
     module procedure mean_flow_output_init
  end interface mean_flow_output_t

contains
  
  function mean_flow_output_init(mf, T_begin, name, path) result(this)
    type(mean_flow_t), intent(in), target ::mf
    real(kind=rp), intent(in) :: T_begin
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(mean_flow_output_t) :: this
    character(len=80) :: fname

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '.fld'
    else if (present(name)) then
       fname = trim(name) // '.fld'
    else if (present(path)) then
       fname = trim(path) // 'mean_field.fld'
    else
       fname = 'mean_field.fld'
    end if

    call output_init(this, fname)
    this%mf => mf
    this%T_begin = T_begin
  end function mean_flow_output_init

  !> Sample a mean flow field at time @a t
  subroutine mean_flow_output_sample(this, t)
    class(mean_flow_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    if (t .ge. this%T_begin) then
       call this%file_%write(this%mf, t)
    end if

  end subroutine mean_flow_output_sample
  
end module mean_flow_output


