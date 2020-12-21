!> Defines an output
module output
  use file, only : file_t
  implicit none

  !> Abstract type defining an output type
  type, abstract :: output_t
     type(file_t) :: file_
   contains
     procedure, pass(this) :: init => output_init
     procedure(output_sample), pass(this), deferred :: sample
  end type output_t

  abstract interface
     subroutine output_sample(this)
       import :: output_t
       class(output_t), intent(inout) :: this
     end subroutine output_sample
  end interface

contains

  !> Output constructor
  subroutine output_init(this, fname)
    class(output_t), intent(inout) :: this
    character(len=*), intent(inout) :: fname

    this%file_ = file_t(fname)
    
  end subroutine output_init
  
end module output
