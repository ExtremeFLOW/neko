module generic_file
  use num_types
  implicit none
  
  !> A generic file handler
  type, abstract :: generic_file_t
     character(len=80) :: fname
     integer :: counter
   contains
     procedure :: init => generic_file_init           !< Constructor
     procedure(generic_file_write), deferred :: write !< Write method
     procedure(generic_file_read), deferred :: read   !< Read method
     procedure :: set_counter => generic_file_set_counter !< Update counter
  end type generic_file_t

  abstract interface
     subroutine generic_file_write(this, data, t)
       import :: generic_file_t
       import :: rp
       class(generic_file_t), intent(inout) :: this
       class(*), target, intent(in) :: data
       real(kind=rp), intent(in), optional :: t
     end subroutine generic_file_write
  end interface
  
  abstract interface
     subroutine generic_file_read(this, data)
       import :: generic_file_t
       class(generic_file_t) :: this
       class(*), target, intent(inout) :: data
     end subroutine generic_file_read
  end interface

contains
  
  !> Generic file constructor
  !! @param fname Filename
  subroutine generic_file_init(this, fname)
    class(generic_file_t) :: this
    character(len=*) :: fname
    
    this%fname = fname
    this%counter = 0
    
  end subroutine generic_file_init

  !> Set the file counter to @a n
  subroutine generic_file_set_counter(this, n)
    class(generic_file_t), intent(inout) :: this
    integer, intent(in) :: n
    this%counter = n    
  end subroutine generic_file_set_counter

end module generic_file
