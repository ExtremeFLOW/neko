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
  use logger, only : neko_log, LOG_SIZE
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  implicit none
  private

  !> Provides a tool to set time step dt
  type, public :: time_step_controller_t
     !> Components recording time stepping info
     logical :: if_variable_dt
     real(kind=rp) :: set_cfl
     real(kind=rp) :: max_dt
     integer :: max_update_frequency
     integer :: dt_last_change
     real(kind=rp) :: alpha !coefficient of running average
     real(kind=rp) :: max_dt_increase_factor, min_dt_decrease_factor
   contains
     !> Initialize object.
     procedure, pass(this) :: init => time_step_controller_init
     !> Set time stepping
     procedure, pass(this) :: set_dt => time_step_controller_set_dt

  end type time_step_controller_t

contains

  !> Constructor
  !! @param order order of the interpolation
  subroutine time_step_controller_init(this, params)
    class(time_step_controller_t), intent(inout) :: this
    type(json_file), intent(inout) :: params

    this%dt_last_change = 0
    call json_get_or_default(params, 'case.variable_timestep',&
                                    this%if_variable_dt, .false.)
    call json_get_or_default(params, 'case.target_cfl',&
                                    this%set_cfl, 0.4_rp)
    call json_get_or_default(params, 'case.max_timestep',&
                                    this%max_dt, huge(0.0_rp))
    call json_get_or_default(params, 'case.cfl_max_update_frequency',&
                                    this%max_update_frequency, 0)
    call json_get_or_default(params, 'case.cfl_running_avg_coeff',&
                                    this%alpha, 0.5_rp)
    call json_get_or_default(params, 'case.max_dt_increase_factor',&
                                    this%max_dt_increase_factor, 1.2_rp)
    call json_get_or_default(params, 'case.min_dt_decrease_factor',&
                                    this%min_dt_decrease_factor, 0.5_rp)

  end subroutine time_step_controller_init

  !> Set new dt based on cfl if requested
  !! @param dt time step in case_t.
  !! @param cfl courant number of current iteration.
  !! @param cfl_avrg average Courant number.
  !! @param tstep the current time step.
  !! @Algorithm:
  !! 1. Set the first time step such that cfl is the set one;
  !! 2. During time-stepping, adjust dt when cfl_avrg is offset by 20%.
  subroutine time_step_controller_set_dt(this, dt, cfl, cfl_avrg, tstep)
    implicit none
    class(time_step_controller_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: dt
    real(kind=rp), intent(in) :: cfl
    real(kind=rp), intent(inout) :: cfl_avrg
    real(kind=rp) :: dt_old, scaling_factor
    character(len=LOG_SIZE) :: log_buf
    integer, intent(in):: tstep

    if (this%if_variable_dt .eqv. .true.) then
       if (tstep .eq. 1) then
          ! set the first dt for desired cfl
          dt = min(this%set_cfl/cfl*dt, this%max_dt)
       else
          ! Calculate the average of cfl over the desired interval
          cfl_avrg = this%alpha * cfl + (1-this%alpha) * cfl_avrg

          if (abs(cfl_avrg - this%set_cfl) .ge. 0.2*this%set_cfl .and. &
             this%dt_last_change .ge. this%max_update_frequency) then

             if (this%set_cfl/cfl .ge. 1) then
                ! increase of time step
                scaling_factor = min(this%max_dt_increase_factor, this%set_cfl/cfl)
             else
                ! reduction of time step
                scaling_factor = max(this%min_dt_decrease_factor, this%set_cfl/cfl)
             end if

             dt_old = dt
             dt = scaling_factor * dt_old
             dt = min(dt, this%max_dt)

             write(log_buf, '(A,E15.7,1x,A,E15.7)') 'Avrg CFL:', cfl_avrg, &
                         'set_cfl:', this%set_cfl
             call neko_log%message(log_buf)

             write(log_buf, '(A,E15.7,1x,A,E15.7)') 'old dt:', dt_old, &
                         'new dt:', dt
             call neko_log%message(log_buf)

             this%dt_last_change = 0

          else
             this%dt_last_change = this%dt_last_change + 1
          end if
       end if

    end if

  end subroutine time_step_controller_set_dt




end module time_step_controller
