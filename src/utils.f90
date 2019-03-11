!> Utilities
!! @details Various utility functions
!
module utils
  implicit none

  integer, parameter :: NEKO_FNAME_LEN = 1024

  interface neko_error
     module procedure neko_error_plain, neko_error_msg
  end interface neko_error
  
contains
  
  subroutine neko_error_plain(error_code)
    integer, optional :: error_code

    write(*,*) '*** ERROR ***'

    if (present(error_code)) then
       stop error_code
    else       
       stop 
    end if

  end subroutine neko_error_plain

  subroutine neko_error_msg(error_msg)
    character(len=*) :: error_msg
    write(*,*) '*** ERROR: ', error_msg,' ***'
    stop 
  end subroutine neko_error_msg
    
end module utils
