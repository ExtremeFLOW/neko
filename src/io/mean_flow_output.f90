!> Defines an output for a mean flow field
module mean_flow_output
  use mean_flow
  use output
  implicit none

  type, extends(output_t) :: mean_flow_output_t
     type(mean_flow_t), pointer :: mf
   contains
     procedure, pass(this) :: sample => mean_flow_output_sample
  end type mean_flow_output_t

  interface mean_flow_output_t
     module procedure mean_flow_output_init
  end interface mean_flow_output_t

contains
  
  function mean_flow_output_init(mf, name) result(this)
    type(mean_flow_t), intent(in), target ::mf
    character(len=*), intent(in), optional :: name
    type(mean_flow_output_t) :: this
    character(len=80) :: fname

    if (present(name)) then
       fname = trim(name) // '.fld'
    else
       fname = 'mean_field.chkp'
    end if

    call output_init(this, fname)
    this%mf => mf
  end function mean_flow_output_init

  !> Sample a mean flow field at time @a t
  subroutine mean_flow_output_sample(this, t)
    class(mean_flow_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%file_%write(this%mf, t)

  end subroutine mean_flow_output_sample
  
end module mean_flow_output


