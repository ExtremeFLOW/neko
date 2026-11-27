!> Redirects Neko errors and warnings to pFUnit exceptions.
!! @details Point throw_error and throw_warning in utils at a routine that
!! raises a pFUnit exception instead of stopping execution, so that tests can
!! check for error emission with @assertExceptionRaised().
!!
!! Errors and warnings are switched independently. Call redirect_errors from
!! a suite's PFUNIT_EXTRA_INITIALIZE hook (see tests/unit/errors/Makefile.in)
!! or from setUp, and restore_errors from tearDown if the suite needs the
!! default behaviour again. Warnings are common on valid code paths, so only
!! redirect them inside the test that expects one: call redirect_warnings at
!! the start of the test and restore_warnings before it returns.
module error_redirection
  use utils, only : throw_error, throw_warning, throw_intf, &
       default_throw_error, default_throw_warning
  use funit, only : SourceLocation, throw

  implicit none
  private

  public :: redirect_errors, restore_errors, &
       redirect_warnings, restore_warnings

contains

  !> Raise a pFUnit exception whenever neko_error is called.
  subroutine redirect_errors()
    throw_error => fail_with_pfunit
  end subroutine redirect_errors

  !> Restore the default behaviour of neko_error, which stops execution.
  subroutine restore_errors()
    throw_error => default_throw_error
  end subroutine restore_errors

  !> Raise a pFUnit exception whenever neko_warning is called.
  subroutine redirect_warnings()
    throw_warning => fail_with_pfunit
  end subroutine redirect_warnings

  !> Restore the default behaviour of neko_warning, which only prints.
  subroutine restore_warnings()
    throw_warning => default_throw_warning
  end subroutine restore_warnings

  subroutine fail_with_pfunit(filename, line, message)
    character(*), intent(in) :: filename
    integer, intent(in) :: line
    character(*), optional, intent(in) :: message

    character(len=:), allocatable :: msg

    if (present(message)) then
       msg = message
    else
       msg = '(no message)'
    end if

    call throw(msg, SourceLocation(filename, line))
  end subroutine fail_with_pfunit

end module error_redirection
