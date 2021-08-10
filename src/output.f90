!> Defines an output
module output
  use num_types
  use file, only : file_t
  implicit none

  !> Abstract type defining an output type
  type, abstract :: output_t
     type(file_t) :: file_
   contains
     procedure, pass(this) :: init => output_init
     procedure, pass(this) :: set_counter => output_set_counter
     procedure(output_sample), pass(this), deferred :: sample
  end type output_t

  !> Abstract interface for sampling an output type at time @a t
  abstract interface
     subroutine output_sample(this, t)
       import :: output_t
       import rp
       class(output_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
     end subroutine output_sample
  end interface

contains

  !> Output constructor
  subroutine output_init(this, fname)
    class(output_t), intent(inout) :: this
    character(len=*), intent(inout) :: fname

    this%file_ = file_t(fname)
    
  end subroutine output_init

  !> Update the output's file counter
  subroutine output_set_counter(this, n)
    class(output_t), intent(inout) :: this
    integer, intent(in) :: n
    call this%file_%set_counter(n)
  end subroutine output_set_counter
  
end module output
