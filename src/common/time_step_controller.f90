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
  use time_state, only : time_state_t
  implicit none
  private

  !> Provides a tool to set time step dt
  type, public :: time_step_controller_t
     logical :: is_variable_dt
     real(kind=rp) :: cfl_trg
     real(kind=rp) :: cfl_avg
     real(kind=rp) :: max_dt, min_dt
     integer :: max_update_frequency, min_update_frequency
     integer :: dt_last_change
     real(kind=rp) :: alpha !< coefficient of running average
     real(kind=rp) :: max_dt_increase_factor, min_dt_decrease_factor
     real(kind=rp) :: dev_tol
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

    this%dt_last_change = -1
    call json_get_or_default(params, 'variable_timestep', &
         this%is_variable_dt, .false.)
    call json_get_or_default(params, 'target_cfl', &
         this%cfl_trg, 0.4_rp)
    call json_get_or_default(params, 'max_timestep', &
         this%max_dt, huge(0.0_rp))
    call json_get_or_default(params, 'min_timestep', &
         this%min_dt, 0.0_rp)
    call json_get_or_default(params, 'max_update_frequency',&
         this%max_update_frequency, 0)
    call json_get_or_default(params, 'min_update_frequency',&
         this%min_update_frequency, huge(0))
    call json_get_or_default(params, 'cfl_running_avg_coeff', &
         this%alpha, 0.5_rp)
    call json_get_or_default(params, 'max_dt_increase_factor', &
         this%max_dt_increase_factor, 1.2_rp)
    call json_get_or_default(params, 'min_dt_decrease_factor', &
         this%min_dt_decrease_factor, 0.5_rp)
    call json_get_or_default(params, 'cfl_deviation_tolerance', &
         this%dev_tol, 0.2_rp)

  end subroutine time_step_controller_init

  !> Set new dt based on cfl if requested
  !! @param dt time step in case_t.
  !! @param cfl courant number of current iteration.
  !! @param tstep the current time step.
  !! @Algorithm:
  !! 1. Set the first time step such that cfl is the set one;
  !! 2. During time-stepping, adjust dt when cfl_avg is offset by 20%.
  subroutine time_step_controller_set_dt(this, time, cfl)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=rp), intent(in) :: cfl
    real(kind=rp) :: dt_old, scaling_factor
    character(len=LOG_SIZE) :: log_buf

    ! Check if variable dt is requested
    if (.not. this%is_variable_dt) return

    ! Reset the average cfl if it is the first time step since the last change
    if (this%dt_last_change .eq. 0) then
       this%cfl_avg = cfl
    end if

    if (this%dt_last_change .eq. -1) then

       ! set the first dt for desired cfl
       time%dt = max(min(this%cfl_trg / cfl * time%dt, &
            this%max_dt), this%min_dt)
       this%dt_last_change = 0
       this%cfl_avg = cfl

    else
       ! Calculate the average of cfl over the desired interval
       this%cfl_avg = this%alpha * cfl + (1 - this%alpha) * this%cfl_avg

       if (abs(this%cfl_avg - this%cfl_trg) .ge. this%dev_tol * this%cfl_trg &
            .and. this%dt_last_change .ge. this%max_update_frequency &
            .or. this%dt_last_change .ge. this%min_update_frequency) then

          if (this%cfl_trg/cfl .ge. 1) then
             ! increase of time step
             scaling_factor = min(this%max_dt_increase_factor, this%cfl_trg/cfl)
          else
             ! reduction of time step
             scaling_factor = max(this%min_dt_decrease_factor, this%cfl_trg/cfl)
          end if

          dt_old = time%dt
          time%dt = scaling_factor * dt_old
          time%dt = max(min(time%dt, this%max_dt), this%min_dt)

          write(log_buf, '(A,E15.7,1x,A,E15.7)') &
               'Average CFL:', this%cfl_avg, &
               'Target  CFL:', this%cfl_trg
          call neko_log%message(log_buf)

          write(log_buf, '(A,E15.7,1x,A,E15.7)') 'Old dt:', dt_old, &
               'New dt:', time%dt
          call neko_log%message(log_buf)

          this%dt_last_change = 0

       else
          this%dt_last_change = this%dt_last_change + 1
       end if
    end if

  end subroutine time_step_controller_set_dt




end module time_step_controller
