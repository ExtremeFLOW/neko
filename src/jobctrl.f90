!> Job control
module jobctrl
  use num_types
  use signal
  use utils
  use comm
  use log
  implicit none

contains

  !> Initialize jobctrl
  subroutine jobctrl_init()
    real(kind=dp) :: jobtime

    ! Start the job clock
    jobtime = jobctrl_jobtime()

    ! Catch SIGXCPU
    call signal_trap_cpulimit()

    ! Catch SIGUSR1 and SIGUSR2
    call signal_trap_usr()

  end subroutine jobctrl_init

  !> Set a job's time limit (in seconds)
  subroutine jobctrl_set_time_limit(sec)
    integer :: sec
    integer :: jstop_sec

    jstop_sec = sec - jobctrl_jobtime()
    call signal_set_timeout(jstop_sec)
    
  end subroutine jobctrl_set_time_limit
  
  !> Check if the job's time limit has been reached
  function jobctrl_time_limit() result(jstop)
    logical :: jstop
    integer :: ierr
    character(len=LOG_SIZE) :: log_buf

    jstop = signal_timeout()

    if (jstop) then
       write(log_buf, '(A)') '! stop at job limit >>>'
       call neko_log%message(log_buf)
    end if

    ! Let rank zero decide if we should stop
    call MPI_Bcast(jstop, 1, MPI_LOGICAL, 0, NEKO_COMM, ierr)
    
  end function jobctrl_time_limit

  !> Returns a job's time in seconds relative to the first call
  function jobctrl_jobtime() result(jobtime)
    real(kind=rp), save :: stime
    real(kind=rp) :: jobtime
    logical, save :: init = .false.
    
    if (.not. init) then
       stime = MPI_WTIME()
       init = .true.
    end if
    
    jobtime = MPI_WTIME() - stime
  end function jobctrl_jobtime
  
end module jobctrl
