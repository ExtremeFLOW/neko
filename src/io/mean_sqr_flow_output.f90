!> Defines an output for a mean squared flow field
module mean_sqr_flow_output
  use mean_sqr_flow
  use num_types
  use output
  implicit none

  type, extends(output_t) :: mean_sqr_flow_output_t
     type(mean_sqr_flow_t), pointer :: msqrf
     real(kind=rp) :: T_begin
   contains
     procedure, pass(this) :: sample => mean_sqr_flow_output_sample
  end type mean_sqr_flow_output_t

  interface mean_sqr_flow_output_t
     module procedure mean_sqr_flow_output_init
  end interface mean_sqr_flow_output_t

contains
  
  function mean_sqr_flow_output_init(msqrf, T_begin, name) result(this)
    type(mean_sqr_flow_t), intent(in), target ::msqrf
    real(kind=rp), intent(in) :: T_begin
    character(len=*), intent(in), optional :: name
    type(mean_sqr_flow_output_t) :: this
    character(len=80) :: fname

    if (present(name)) then
       fname = trim(name) // '.fld'
    else
       fname = 'mean_sqr_field.fld'
    end if

    call output_init(this, fname)
    this%msqrf => msqrf
    this%T_begin = T_begin
  end function mean_sqr_flow_output_init

  !> Sample a mean squared flow field at time @a t
  subroutine mean_sqr_flow_output_sample(this, t)
    class(mean_sqr_flow_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    if (t .ge. this%T_begin) then
       call this%file_%write(this%msqrf, t)
    end if

  end subroutine mean_sqr_flow_output_sample
  
end module mean_sqr_flow_output


