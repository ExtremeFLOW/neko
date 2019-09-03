!> Utilities
!! @details Various utility functions
module utils
  implicit none

  integer, parameter :: NEKO_FNAME_LEN = 1024

  interface neko_error
     module procedure neko_error_plain, neko_error_msg
  end interface neko_error
  
contains
  
  subroutine neko_error_plain(error_code)
    integer, optional :: error_code
    integer :: code

    if (present(error_code)) then
       write(*,*) '*** ERROR ***', error_code
       stop
    else
       write(*,*) '*** ERROR ***'
       stop
    end if

  end subroutine neko_error_plain

  subroutine neko_error_msg(error_msg)
    character(len=*) :: error_msg
    write(*,*) '*** ERROR: ', error_msg,' ***'
    stop 
  end subroutine neko_error_msg

  subroutine neko_warning(warning_msg)
    character(len=*) :: warning_msg
    write(*,*) '*** WARNING: ', warning_msg,' ***'
  end subroutine neko_warning
    
end module utils
