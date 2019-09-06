!> Utilities
!! @details Various utility functions
module utils
  implicit none

  integer, parameter :: NEKO_FNAME_LEN = 1024

  interface neko_error
     module procedure neko_error_plain, neko_error_msg
  end interface neko_error
  
contains
  
  !> Find position (in the string) of a filename's suffix
  pure function filename_suffix_pos(fname) result(suffix_pos)
    character(len=*), intent(in) :: fname
    integer :: suffix_pos
    suffix_pos = scan(trim(fname), '.', back=.true.)    
  end function filename_suffix_pos

  !> Extract a filename's suffix
  subroutine filename_suffix(fname, suffix)
    character(len=*) :: fname
    character(len=*) :: suffix
    suffix = trim(fname(filename_suffix_pos(fname) + 1:len_trim(fname)))
  end subroutine filename_suffix

  !> Change a filename's suffix
  subroutine filename_chsuffix(fname, new_fname, new_suffix)
    character(len=*) :: fname
    character(len=*) :: new_fname
    character(len=*) :: new_suffix
    integer :: suffix_pos

    suffix_pos = filename_suffix_pos(fname)
    new_fname = trim(fname(1:suffix_pos))//new_suffix

  end subroutine filename_chsuffix

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
