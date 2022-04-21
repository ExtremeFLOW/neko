! Copyright (c) 2021, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Job control
module jobctrl
  use num_types
  use signal
  use utils
  use comm
  use logger
  implicit none

  interface jobctrl_set_time_limit
     module procedure jobctrl_set_time_limit_sec, jobctrl_set_time_limit_str
  end interface jobctrl_set_time_limit

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

  !> Set a job's time limit (in walltime 'HH:MM:SS')
  subroutine jobctrl_set_time_limit_str(limit_str)
    character(len=*) limit_str
    integer :: str_len, sep_h, h, m, s

    str_len = len_trim(limit_str)
    
    if (str_len .lt. 8) then
       call neko_error('Invalid job limit')
    end if

    ! hour
    sep_h = scan(trim(limit_str), ':')
    read(limit_str(1:sep_h-1), *) h

    !min
    read(limit_str(sep_h+1:sep_h+2), *) m

    !sec
    read(limit_str(sep_h+4:str_len), *) s

    call jobctrl_set_time_limit_sec(h*3600 + m * 60 + s)
    
  end subroutine jobctrl_set_time_limit_str
  
  !> Set a job's time limit (in seconds)
  subroutine jobctrl_set_time_limit_sec(sec)
    integer :: sec
    integer :: jstop_sec

    if (sec .gt. 0) then
       jstop_sec = sec - jobctrl_jobtime()
       call signal_set_timeout(jstop_sec)
    end if
    
  end subroutine jobctrl_set_time_limit_sec
  
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
