module error_redirection
  use utils, only: throw_error, throw_warning, throw_intf
  use funit, only: SourceLocation, throw

  implicit none
  private

  public :: redirect_errors

contains

  subroutine redirect_errors()
    throw_error => fail_with_pfunit
    throw_warning => fail_with_pfunit
  end subroutine redirect_errors

  subroutine fail_with_pfunit(filename, line, message)
    character(*), intent(in) :: filename
    integer, intent(in) :: line
    character(*), optional, intent(in) :: message

    character(len=:), allocatable :: msg

    if (present(message)) then
       msg = message
    else
       msg = '<no message>'
    end if

    call throw(msg, SourceLocation(filename, line))
  end subroutine fail_with_pfunit

end module error_redirection
