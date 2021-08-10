!> Interface to signal handler
module signal
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

  interface
     integer (c_int8_t) function sighdl_timeout() &
          bind(c, name='sighdl_timeout')
       use, intrinsic :: iso_c_binding
     end function sighdl_timeout
  end interface

  interface
     integer (c_int8_t) function sighdl_usr() &
          bind(c, name='sighdl_usr')
       use, intrinsic :: iso_c_binding
     end function sighdl_usr
  end interface
  
  interface
     integer (c_int) function sighdl_set_timeout(sec) &
          bind(c, name='sighdl_set_timeout')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: sec
     end function sighdl_set_timeout
  end interface
  
  interface
     integer (c_int) function sighdl_trap_cpulimit() &
          bind(c, name='sighdl_trap_cpulimit')
       use, intrinsic :: iso_c_binding
     end function sighdl_trap_cpulimit
  end interface

  interface
     integer (c_int) function sighdl_trap_usr() &
          bind(c, name='sighdl_trap_usr')
       use, intrinsic :: iso_c_binding
     end function sighdl_trap_usr
  end interface

contains
  
  !> Check if any timeout has occurred (either SIGXCPU or SIGALRM)
  function signal_timeout() result(timeout)
    logical :: timeout

    if (sighdl_timeout() .gt. 0)  then
       timeout = .true.
    else
       timeout = .false.
    end if
    
  end function signal_timeout

  !> Check if a user signal has been raised
  function signal_usr(usr) result(raised)
    integer, intent(in) :: usr
    logical :: raised
    integer(kind=c_int8_t) :: usr12

    if (usr .gt. 2) then
       call neko_error('Invalid usr signal')
    end if

    usr12 = sighdl_timeout()
    
    if (bge(usr12, usr)) then
       raised = .true.
    else
       raised = .false.
    end if
    
  end function signal_usr
  
  !> Set a timeout after @a seconds
  subroutine signal_set_timeout(sec)
    integer(kind=c_int) :: sec

    if (sighdl_set_timeout(sec) .lt. 0) then
       call neko_error('sighdl failed to set SIGALRM')
    end if
    
  end subroutine signal_set_timeout
  
  !> Initialize signal handler to trap SIGXCPU
  subroutine signal_trap_cpulimit()
    logical, save :: initialized = .false.

    if (.not. initialized) then
       if (sighdl_trap_cpulimit() .lt. 0) then
          call neko_error('sighdl failed to trap SIGXCPU')
       end if
       initialized = .true.
    end if
    
  end subroutine signal_trap_cpulimit

  !> Initialize signal handler to trap SIGUSR1 and SIGUSR2
  subroutine signal_trap_usr()
    logical, save :: initialized = .false.

    if (.not. initialized) then
       if (sighdl_trap_usr() .lt. 0) then
          call neko_error('sighdl failed to trap SIGUSR*')
       end if
       initialized = .true.
    end if
    
  end subroutine signal_trap_usr
  
end module signal
