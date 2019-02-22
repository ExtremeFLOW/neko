!> Utilities
!! @details Various utility functions
!
module utils
  implicit none

  integer, parameter :: NEKO_FNAME_LEN = 1024

contains

  subroutine neko_error(error_msg)
    character(len=*) :: error_msg
    write(*,*) '*** ERROR: ', error_msg,' ***'
    stop 
  end subroutine neko_error
    
end module utils
