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
  use num_types, only : dp
  use utils, only : neko_error
  use time_state, only : time_state_t
  implicit none
  private

  !> A utility type for determining whether an action should be executed based
  !! on the current time value. Used to e.g. control whether we should write a
  !! file or execute a simcomp.
  !! Note that the nexecutions variable should be incremented externally by
  !! calling the `register_execution` procedure.
  !! This is to allow running the the `check` multiple times at the same time
  !! step.
  type, public :: time_based_controller_t
     !> Frequency of execution.
     real(kind=dp) :: frequency = 0.0_dp
     !> Time interval between executions.
     real(kind=dp) :: time_interval = 0.0_dp
     !> Number of time steps in between executions.
     integer :: nsteps = 0
     !> Simulation start time.
     real(kind=dp) :: start_time = 0.0_dp
     !> Simulation end time.
     real(kind=dp) :: end_time = 0.0_dp
     !> Number of times already executed.
     integer :: nexecutions = 0
     !> Whether to never output.
     logical :: never = .false.
     !> Control mode defining the meaning of `control_value`.
     !> Can be `simulationtime`, `tsteps`, `nsamples` or `never`.
     character(len=:), allocatable :: control_mode
     !> Defines the frequency of writes.
     real(kind=dp) :: control_value

   contains
     !> Constructor.
     procedure, pass(this) :: init => time_based_controller_init
     !> Destructor.
     procedure, pass(this) :: free => time_based_controller_free
     !> Check if the execution should be performed.
     procedure, pass(this) :: check => time_based_controller_check
     !> Increment `nexectutions`.
     procedure, pass(this) :: register_execution => &
          time_based_controller_register_execution
     !> Set the counter based on a time (for restarts)
     procedure, pass(this) :: set_counter => &
          time_based_controller_set_counter
     !> The next time at which `check` will return true.
     procedure, pass(this) :: next_time => time_based_controller_next_time

  end type time_based_controller_t

  !> The execution schedule of a controller, stripped down to what is needed
  !! to predict when it will fire next.
  type :: schedule_t
     !> Time interval between executions.
     real(kind=dp) :: time_interval = 0.0_dp
     !> Time after which the controller stops executing.
     real(kind=dp) :: end_time = 0.0_dp
  end type schedule_t

  !> The distinct schedules of all time-based controllers constructed so far.
  !! @note Kept as a module variable rather than as a list of pointers to the
  !! controllers themselves, since the latter live in arrays that are
  !! reallocated as outputs and simcomps are added.
  type(schedule_t), allocatable :: registered_schedules(:)

  !> Number of used entries in `registered_schedules`.
  integer :: n_registered_schedules = 0

  interface assignment(=)
     module procedure time_based_controller_assignment
  end interface assignment(=)

  public :: time_based_controller_next_scheduled_time, &
       time_based_controller_reset_schedules

contains

  !> Constructor.
  !! @param end_time The final simulation time.
  !! @param control_mode The way to interpret the `control_value` parameter.
  !! @param control_value The value defining the execution frequency.
  subroutine time_based_controller_init(this, start_time, end_time, &
       control_mode, control_value)
    class(time_based_controller_t), intent(inout) :: this
    real(kind=dp), intent(in) :: start_time
    real(kind=dp), intent(in) :: end_time
    character(len=*), intent(in) :: control_mode
    real(kind=dp), intent(in) :: control_value

    this%start_time = start_time
    this%end_time = end_time
    this%control_mode = control_mode
    this%control_value = control_value

    if (trim(control_mode) .eq. 'simulationtime') then
       this%time_interval = control_value
       this%frequency = 1/this%time_interval
       this%nsteps = 0
    else if (trim(control_mode) .eq. 'nsamples') then
       if (control_value .le. 0) then
          call neko_error("nsamples must be positive")
       end if

       this%frequency = control_value / (end_time - start_time)
       this%time_interval = 1.0_dp / this%frequency
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
       & tsteps, or never, but received "//trim(control_mode))
    end if

    ! Make the schedule known globally, so that the time-step controller can
    ! shrink dt to land exactly on the execution times, if asked to.
    if (.not. this%never .and. this%time_interval .gt. 0.0_dp) then
       call register_schedule(this%time_interval, this%end_time)
    end if
  end subroutine time_based_controller_init

  !> Destructor.
  subroutine time_based_controller_free(this)
    class(time_based_controller_t), intent(inout) :: this

    if (allocated(this%control_mode)) then
       deallocate(this%control_mode)
    end if

    this%frequency = 0.0_dp
    this%time_interval = 0.0_dp
    this%nsteps = 0
    this%start_time = 0.0_dp
    this%end_time = 0.0_dp
    this%nexecutions = 0
    this%never = .false.
  end subroutine time_based_controller_free

  !> Check if the execution should be performed.
  !! @param t Time value.
  !! @param tstep Current timestep.
  !! @param dt Timestep size.
  !! @param force Whether to force returning true. Optional.
  !! @note In the logic, `nsteps` being zero corresponds to us not knowing the
  !! number of time-steps between executions and thus having to rely on
  !! `nexecutions`. This is done in anticipation of having a variable timestep.
  !! A fraction of the time step (10 percent) is used as a tolerance.
  function time_based_controller_check(this, time, force) result(check)
    class(time_based_controller_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    logical, intent(in), optional :: force
    real(kind=dp) :: t
    integer :: tstep
    real(kind=dp) :: dt
    logical :: check
    logical :: ifforce

    t = time%t - time%start_time
    dt = time%dt
    tstep = time%tstep

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
    else if (time%t - this%start_time .gt. this%end_time - this%start_time) then
       check = .false.
    else if ( (this%nsteps .eq. 0) .and. &
         (t .ge. this%nexecutions * this%time_interval - 0.1_dp * dt) ) then
       check = .true.
    else if (this%nsteps .gt. 0) then
       if (mod(tstep, this%nsteps) .eq. 0) then
          check = .true.
       end if
    end if
  end function time_based_controller_check

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

  !> Set the counter based on a time (for restarts)
  !! @param time Current time.
  subroutine time_based_controller_set_counter(this, time)
    class(time_based_controller_t), intent(inout) :: this
    type(time_state_t) :: time

    if (this%nsteps .eq. 0) then
       this%nexecutions = int(((time%t - time%start_time) + 0.1_dp*time%dt) &
            / this%time_interval) + 1
    end if

  end subroutine time_based_controller_set_counter

  !> The next time at which this controller is scheduled to execute.
  !! @param time Current time.
  !! @return The next execution time, or `huge(0.0_dp)` if the controller has
  !! no time-based schedule, or will never execute again.
  pure function time_based_controller_next_time(this, time) result(next_time)
    class(time_based_controller_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    real(kind=dp) :: next_time

    if (this%never .or. this%nsteps .gt. 0) then
       next_time = huge(0.0_dp)
    else
       next_time = schedule_next_time(this%time_interval, this%end_time, time)
    end if

  end function time_based_controller_next_time

  !> The earliest time at which any of the constructed controllers is
  !! scheduled to execute.
  !! @param time Current time.
  !! @return The next execution time, or `huge(0.0_dp)` if no controller has a
  !! time-based schedule left to run.
  pure function time_based_controller_next_scheduled_time(time) &
       result(next_time)
    type(time_state_t), intent(in) :: time
    real(kind=dp) :: next_time
    integer :: i

    next_time = huge(0.0_dp)
    do i = 1, n_registered_schedules
       next_time = min(next_time, &
            schedule_next_time(registered_schedules(i)%time_interval, &
            registered_schedules(i)%end_time, time))
    end do

  end function time_based_controller_next_scheduled_time

  !> Forget all registered schedules, e.g. when setting up a new case.
  subroutine time_based_controller_reset_schedules()

    if (allocated(registered_schedules)) deallocate(registered_schedules)
    n_registered_schedules = 0

  end subroutine time_based_controller_reset_schedules

  !> Add a schedule to `registered_schedules`, unless an identical one is
  !! already there.
  !! @param time_interval Time interval between executions.
  !! @param end_time Time after which the controller stops executing.
  subroutine register_schedule(time_interval, end_time)
    real(kind=dp), intent(in) :: time_interval
    real(kind=dp), intent(in) :: end_time
    type(schedule_t), allocatable :: tmp(:)
    integer :: i

    if (.not. allocated(registered_schedules)) then
       allocate(registered_schedules(8))
    end if

    do i = 1, n_registered_schedules
       if (registered_schedules(i)%time_interval .eq. time_interval .and. &
            registered_schedules(i)%end_time .eq. end_time) return
    end do

    if (n_registered_schedules .eq. size(registered_schedules)) then
       allocate(tmp(2 * size(registered_schedules)))
       tmp(1:n_registered_schedules) = &
            registered_schedules(1:n_registered_schedules)
       call move_alloc(tmp, registered_schedules)
    end if

    n_registered_schedules = n_registered_schedules + 1
    registered_schedules(n_registered_schedules)%time_interval = time_interval
    registered_schedules(n_registered_schedules)%end_time = end_time

  end subroutine register_schedule

  !> The next time a periodic schedule fires, strictly after the current time.
  !! @details The times at which `check` returns true are the grid
  !! `time%start_time + k * time_interval`, so this is the first grid point
  !! ahead of `time%t`. A tolerance of a few units in the last place is used,
  !! such that the point we have just landed on is not returned again.
  !! @param time_interval Time interval between executions.
  !! @param end_time Time after which the controller stops executing.
  !! @param time Current time.
  !! @return The next execution time, or `huge(0.0_dp)` if there is none.
  pure function schedule_next_time(time_interval, end_time, time) &
       result(next_time)
    real(kind=dp), intent(in) :: time_interval
    real(kind=dp), intent(in) :: end_time
    type(time_state_t), intent(in) :: time
    real(kind=dp) :: next_time
    real(kind=dp) :: nintervals, tol

    next_time = huge(0.0_dp)
    if (time_interval .le. 0.0_dp) return
    if (time%t .gt. end_time) return

    ! Number of whole intervals since the start, rounded towards minus
    ! infinity. Kept in floating point to not overflow for tiny intervals.
    nintervals = aint((time%t - time%start_time) / time_interval)
    if (nintervals .gt. (time%t - time%start_time) / time_interval) then
       nintervals = nintervals - 1.0_dp
    end if

    tol = 10.0_dp * spacing(max(abs(time%t), time_interval))

    next_time = time%start_time + (nintervals + 1.0_dp) * time_interval
    if (next_time .le. time%t + tol) next_time = next_time + time_interval

    if (next_time .gt. end_time) next_time = huge(0.0_dp)

  end function schedule_next_time


end module time_based_controller
