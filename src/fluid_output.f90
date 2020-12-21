  !> Defines an output for a fluid
module fluid_output
  use fluid_method
  use output
  implicit none

  !> Fluid output
  type, extends(output_t) :: fluid_output_t
     class(fluid_scheme_t), pointer :: fluid
   contains
     procedure, pass(this) :: sample => fluid_output_sample
  end type fluid_output_t

  interface fluid_output_t
     module procedure fluid_output_init
  end interface fluid_output_t

contains

  function fluid_output_init(fluid, name) result(this)
    class(fluid_scheme_t), intent(in), target :: fluid
    character(len=*), intent(in), optional :: name
    type(fluid_output_t) :: this
    character(len=80) :: fname

    if (present(name)) then
       fname = trim(name) // '.fld'
    else       
       fname = 'fluid.fld'
    end if
    
    call output_init(this, fname)    
    this%fluid => fluid
  end function fluid_output_init

  subroutine fluid_output_sample(this)
    class(fluid_output_t), intent(inout) :: this
    call this%file_%write(this%fluid)
  end subroutine fluid_output_sample
  
end module fluid_output
