! Copyright (c) 2023, The Neko Authors
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
!> Contains the `time_based_controller_t` type.
module time_based_controller
  use num_types, only : rp
  use utils, only : neko_error
  implicit none
  private

  !> A utility type for determening whether an action should be executed based
  !! on the current time value. Used to e.g. control whether we should write a
  !! file or execute a simcomp.
  !! Note that the nexecutions variable should be incremented externally by
  !! calling the `register_execution` procedure.
  !! This is to allow running the the `check` multiple times at the same time
  !! step.
  type, public :: time_based_controller_t
     !> Frequency of execution.
     real(kind=rp) :: frequency = 0
     !> Time interval between executions.
     real(kind=rp) :: time_interval = 0
     !> Number of timesteps in between executions.
     integer :: nsteps = 0
     !> Simulation end time.
     real(kind=rp) :: end_time = 0
     !> Number of times already executed.
     integer :: nexecutions = 0
     !> Whether to never output.
     logical :: never = .false.

   contains
    !> Constructor.
     procedure, pass(this) :: init => time_based_controller_init
    !> Check if the execution should be performed.
     procedure, pass(this) :: check => time_based_controller_check
    !> Increment `nexectutions`.
     procedure, pass(this) :: register_execution => &
       time_based_controller_register_execution

  end type time_based_controller_t

  interface assignment(=)
     module procedure time_based_controller_assignment
  end interface assignment(=)

contains
  
  !> Constructor.
  !! @param end_time The final simulation time.
  !! @param control_mode The way to interpret the `control_value` parameter.
  !! @param control_value The value definining the execution frequency.
  subroutine time_based_controller_init(this, end_time, control_mode, &
     control_value)
    class(time_based_controller_t), intent(inout) :: this
    real(kind=rp), intent(in) :: end_time
    character(len=*), intent(in) :: control_mode
    real(kind=rp), intent(in) :: control_value

    this%end_time = end_time

    if (trim(control_mode) .eq. 'simulationtime') then
        this%time_interval = control_value
        this%frequency = 1/this%time_interval
        this%nsteps = 0 
    else if (trim(control_mode) .eq. 'nsamples') then
        this%frequency = control_value / end_time
        this%time_interval = 1.0_rp / this%frequency
        this%nsteps = 0
    else if (trim(control_mode) .eq. 'tsteps') then 
        this%nsteps = control_value
        ! if the timestep will be variable, we cannot compute these.
        this%frequency = 0
        this%time_interval = 0
    else if (trim(control_mode) .eq. 'never') then 
        this%never = .true.
    else
        call neko_error("The control parameter must be simulationtime, nsamples&
          & or tsteps, but received "//trim(control_mode))
    end if
  end subroutine time_based_controller_init

  !> Check if the execution should be performed.
  !! @param t Time value.
  !! @param tstep Current timestep.
  !! @param force Whether to force returning true. Optional.
  !! @note In the logic, `nsteps` being zero corresponds to us not knowing the
  !! number of time-steps between executions and thus having to rely on
  !! `nexecutions`. This is done in anticipation of having a variable timestep.
  function time_based_controller_check(this, t, tstep, force) result(check)
    class(time_based_controller_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    logical, intent(in), optional :: force
    logical :: check
    logical :: ifforce

    if (present(force)) then
       ifforce = force
    else
       ifforce = .false.
    end if

    check = .false.
    if (ifforce) then
        check = .true.
    else if (this%never) then
        check = .false.
    else if ( (this%nsteps .eq. 0) .and. &
              (t .ge. this%nexecutions * this%time_interval) ) then
        check = .true.
    else if (this%nsteps .gt. 0) then
       if (mod(tstep, this%nsteps) .eq. 0) then
        check = .true.
       end if
    end if
  end function

  !> Assignment operator. Simply copies attribute values.
  !! @param ctrl1 Left-hand side.
  !! @param ctrl2 Right-hand side.
  subroutine time_based_controller_assignment(ctrl1, ctrl2)
    type(time_based_controller_t), intent(inout) :: ctrl1
    type(time_based_controller_t), intent(in) :: ctrl2

    ctrl1%end_time = ctrl2%end_time
    ctrl1%frequency = ctrl2%frequency
    ctrl1%nsteps = ctrl2%nsteps
    ctrl1%time_interval = ctrl2%time_interval
    ctrl1%nexecutions = ctrl2%nexecutions

  end subroutine time_based_controller_assignment

  !> Increment `nexectutions`.
  subroutine time_based_controller_register_execution(this)
    class(time_based_controller_t), intent(inout) :: this

    this%nexecutions = this%nexecutions + 1

  end subroutine time_based_controller_register_execution



end module time_based_controller