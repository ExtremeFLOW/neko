! Copyright (c) 2025, The Neko Authors
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

!> Module with things related to the simulation time
module time_state
  use num_types, only : rp
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_QUIET
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  implicit none
  private

  !> A struct that contains all info about the time, expand as needed
  type, public :: time_state_t
     real(kind=rp), dimension(10) :: tlag = 0.0_rp !< Old times
     real(kind=rp), dimension(10) :: dtlag = 0.0_rp !< Old dts
     ! @note Might be good to make these of kind=xp in the future.
     ! Would ensure that the time is not subject to stagnation.
     real(kind=rp) :: t = 0.0_rp !< Current time
     real(kind=rp) :: dt = 0.0_rp !< Current dt
     real(kind=rp) :: start_time = 0.0_rp !< Start time
     real(kind=rp) :: end_time = 0.0_rp !< End time
     integer :: tstep = 0 !< Current timestep

   contains
     generic :: init => init_from_components, init_from_json
     procedure, pass(this) :: init_from_components => &
          time_state_init_from_components
     procedure, pass(this) :: init_from_json => time_state_init_from_json
     procedure, pass(this) :: reset => time_state_reset
     procedure, pass(this) :: status => time_state_status
  end type time_state_t

contains

  !> Initialize time state
  subroutine time_state_init_from_components(this, start_time, end_time, step)
    class(time_state_t), intent(inout) :: this
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time
    real(kind=rp), intent(in) :: step

    this%start_time = start_time
    this%end_time = end_time
    this%dt = step

    this%t = start_time

  end subroutine time_state_init_from_components

  !> Initialize time state from JSON
  subroutine time_state_init_from_json(this, params)
    class(time_state_t), intent(inout) :: this
    type(json_file), intent(inout) :: params

    real(kind=rp) :: time_step
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time
    logical :: is_variable

    call json_get_or_default(params, 'variable_timestep', is_variable, .false.)
    if (.not. is_variable) then
       call json_get(params, 'timestep', time_step)
    else
       ! randomly set an initial dt to get cfl when dt is variable
       time_step = 1.0_rp
    end if

    call json_get_or_default(params, 'start_time', start_time, 0.0_rp)
    call json_get(params, 'end_time', end_time)

    call this%init_from_components(start_time, end_time, time_step)

  end subroutine time_state_init_from_json

  !> Reset time state
  subroutine time_state_reset(this)
    class(time_state_t), intent(inout) :: this

    this%t = this%start_time
    this%tstep = 0
    this%dtlag = 0.0_rp
    this%tlag = 0.0_rp

  end subroutine time_state_reset

  !> Write status banner
  subroutine time_state_status(this)
    class(time_state_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf
    character(len=38) :: log_fmt
    real(kind=rp) :: t_prog

    t_prog = 100.0_rp * (this%t - this%start_time) / &
         (this%end_time - this%start_time)

    write(log_fmt, '(A,I2,A)') &
         '(A7,1X,I10,1X,A4,E15.7,', LOG_SIZE - 49, 'X,A2,F6.2,A3)'
    write(log_buf, log_fmt) 'Step = ', this%tstep, 't = ', this%t, &
         '[ ', t_prog, '% ]'

    call neko_log%message(repeat('-', LOG_SIZE))
    call neko_log%message(log_buf, NEKO_LOG_QUIET)
    call neko_log%message(repeat('-', LOG_SIZE))

  end subroutine time_state_status

end module time_state
