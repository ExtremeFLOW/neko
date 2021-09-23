!> Defines an output for a checkpoint
module chkp_output
  use checkpoint
  use output
  implicit none

  type, extends(output_t) :: chkp_output_t
     type(chkp_t), pointer :: chkp
   contains
     procedure, pass(this) :: sample => chkp_output_sample
  end type chkp_output_t

  interface chkp_output_t
     module procedure chkp_output_init
  end interface chkp_output_t
  
contains

  function chkp_output_init(chkp, name, path) result(this)
    type(chkp_t), intent(in), target :: chkp
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(chkp_output_t) :: this
    character(len=80) :: fname

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '.chkp'
    else if (present(name)) then
       fname = trim(name) // '.chkp'
    else if (present(path)) then
       fname = trim(path) // 'fluid.chkp'
    else
       fname = 'fluid.chkp'
    end if

    call output_init(this, fname)
    this%chkp => chkp
  end function chkp_output_init

  !> Sample a checkpoint at time @a t
  subroutine chkp_output_sample(this, t)
    class(chkp_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%file_%write(this%chkp, t)

  end subroutine chkp_output_sample
  
end module chkp_output
