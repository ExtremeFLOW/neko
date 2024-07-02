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
!> Interface to signal handler
module signal
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: signal_timeout, signal_set_timeout, signal_trap_cpulimit,&
            signal_trap_usr, signal_usr

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

    usr12 = sighdl_usr()

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
