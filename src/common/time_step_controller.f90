! Copyright (c) 2022, The Neko Authors
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
!> Implements type time_step_controller.
module time_step_controller
  use num_types, only : rp
  use case
  use json_utils, only : json_get_or_default
  implicit none
  private

  !> Provides a tool to set time step dt
  type, public :: time_step_controller_t
     !> Components recording time stepping info
     logical :: if_variable_dt = .false.
     logical :: if_restart = .false.
     real(kind=rp) :: set_cfl = 0.0_rp
     real(kind=rp) :: max_dt = 0.0_rp
     integer :: max_update_frequency = 1
   contains
     !> Initialize object.
     procedure, pass(this) :: init => time_step_controller_init
     !> Set time stepping
     procedure, pass(this) :: set_dt => time_step_controller_set_dt

  end type time_step_controller_t

contains

  !> Constructor
  !! @param order order of the interpolation
  subroutine time_step_controller_init(this, C)
    class(time_step_controller_t), intent(inout) :: this
    type(case_t), intent(inout) :: C
    logical :: found
    character(len=:), allocatable :: restart_file

    call C%params%get('case.timestep', this%max_dt, found)
    call C%params%get('case.constant_cfl', this%set_cfl, this%if_variable_dt)
    call C%params%get('case.restart_file', restart_file, this%if_restart)
    call json_get_or_default(C%params, 'case.cfl_max_update_frequency',&
                                    this%max_update_frequency, 1)

  end subroutine time_step_controller_init

  !> Set new dt based on cfl if requested
  !! @param C case type.
  !! @param cfl courant number of current iteration.
  !! @param cfl_avrg average Courant number.
  !! @param dt_last_change time step since last dt change.
  !! @param tstep the current time step.
  !! @Algorithm: 
  !! 1. If the simulation is not a restarted one,
  !! set the first time step such that cfl is the set one;
  !! 2. If the simulation is a restarted one,
  !! set the first time step from the restart file.
  !! 3. During time-stepping, adjust dt when cfl_avrg is offset by 20%.
  subroutine time_step_controller_set_dt(this, C, cfl, cfl_avrg, dt_last_change, tstep)
    implicit none
    class(time_step_controller_t), intent(inout) :: this
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(in) :: cfl
    real(kind=rp), intent(inout) :: cfl_avrg
    integer, intent(inout) :: dt_last_change
    real(kind=rp) :: dt_old, scaling_factor
    character(len=LOG_SIZE) :: log_buf    
    real(kind=rp) :: alpha = 0.5_rp !coefficient of running average
    integer, intent(in):: tstep

    if (this%if_variable_dt .eqv. .true.) then
       if ((tstep .eq. 1) .and. .not. this%if_restart) then
          ! set the first dt for desired cfl
          C%dt = min(this%set_cfl/cfl*C%dt, this%max_dt)
       else if ((tstep .eq. 1) .and. this%if_restart) then
          C%dt = C%dtlag(1)
       else
          ! Calculate the average of cfl over the desired interval
          cfl_avrg = alpha * cfl + (1-alpha) * cfl_avrg

          if (abs(cfl_avrg - this%set_cfl) .ge. 0.2*this%set_cfl .and. &
             dt_last_change .ge. this%max_update_frequency) then

             if (this%set_cfl/cfl .ge. 1) then 
                scaling_factor = min(1.2_rp, this%set_cfl/cfl) 
             else
                scaling_factor = max(0.8_rp, this%set_cfl/cfl) 
             end if

             dt_old = C%dt
             C%dt = scaling_factor * dt_old
             C%dt = min(C%dt, this%max_dt)

             write(log_buf, '(A,E15.7,1x,A,E15.7)') 'Avrg CFL:', cfl_avrg, &
                         'set_cfl:', this%set_cfl
             call neko_log%message(log_buf)

             write(log_buf, '(A,E15.7,1x,A,E15.7)') 'old dt:', dt_old, &
                         'new dt:', C%dt
             call neko_log%message(log_buf)

             dt_last_change = 0

          else
             dt_last_change = dt_last_change + 1
          end if
       end if

    end if

  end subroutine time_step_controller_set_dt




end module time_step_controller
