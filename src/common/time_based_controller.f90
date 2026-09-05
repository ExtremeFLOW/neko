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
  use num_types, only : dp, i8
  use utils, only : neko_error
  use time_state, only : time_state_t
  implicit none
  private

  !> Tolerance used when comparing the simulation time to a scheduled
  !! execution time, as a fraction of the current time-step size (or of the
  !! output interval, whichever is smaller).
  !! A scheduled execution is performed at the first time step for which
  !! `t >= t_scheduled - TIME_TOL * dt`, so that a step landing a hair short
  !! of the scheduled time (round-off in the accumulation of `t`) still
  !! triggers it, instead of postponing it by a full step. Capping the
  !! tolerance by the output interval keeps a time step that is much larger
  !! than the interval from dragging the first scheduled execution to before
  !! the start of the simulation.
  real(kind=dp), public, parameter :: TIME_TOL = 0.1_dp

  !> Relative tolerance used when deciding whether a scheduled time still
  !! falls inside the interval covered by the controller. Only meant to
  !! absorb round-off in `k * time_interval`.
  real(kind=dp), public, parameter :: SPAN_TOL = 1.0e-9_dp

  !> Largest number of intervals between the anchor of a schedule and the
  !! start of the simulation for which anchoring is still meaningful. Beyond
  !! it the interval is comparable to the resolution of the time itself, and
  !! the schedule is anchored to the start of the simulation instead.
  real(kind=dp), public, parameter :: MAX_ANCHOR_INDEX = 1.0e15_dp

  !> A utility type for determining whether an action should be executed based
  !! on the current time value. Used to e.g. control whether we should write a
  !! file or execute a simcomp.
  !!
  !! @details
  !! The controller defines a *schedule*: the sequence of times (or time
  !! steps) at which the action should be performed. For the time based
  !! control modes the schedule is
  !! \f$ t_k = t_{anchor} + k \Delta t_{out} \f$, restricted to the
  !! \f$ t_k \f$ that fall inside \f$ [t_{start}, t_{end}] \f$, plus, if
  !! `write_at_start` is set, one execution at \f$ t_{start} \f$ itself.
  !!
  !! A `simulationtime` schedule is anchored at zero, so that an output asked
  !! for every \f$ \Delta t_{out} \f$ lands on whole multiples of it no
  !! matter what time the simulation was started from. An `nsamples` schedule
  !! divides the simulated interval instead, and is therefore anchored at
  !! \f$ t_{start} \f$.
  !!
  !! The schedule is an absolute property of the case: it does not depend on
  !! how many executions have been performed, nor on when the run was
  !! started. This is what makes it possible to restart a simulation without
  !! either losing or duplicating an output: `set_counter` simply moves the
  !! cursor to the first scheduled time that has not been reached yet, using
  !! exactly the same predicate that `check` uses.
  !!
  !! `check` is free of observable side effects and idempotent within a time
  !! step: it may be called any number of times per step, and returns
  !! `.false.` once an execution has been registered for the current step.
  !! The `nexecutions` counter is incremented externally by calling
  !! `register_execution`.
  type, public :: time_based_controller_t
     !> Frequency of execution.
     real(kind=dp) :: frequency = 0.0_dp
     !> Time interval between executions. Always positive.
     real(kind=dp) :: time_interval = 0.0_dp
     !> Number of time steps in between executions.
     integer :: nsteps = 0
     !> Time at which the schedule starts (and which it is anchored to).
     real(kind=dp) :: start_time = 0.0_dp
     !> Time after which nothing is scheduled anymore.
     real(kind=dp) :: end_time = 0.0_dp
     !> Number of times already executed.
     integer :: nexecutions = 0
     !> Whether to never execute.
     logical :: never = .false.
     !> Control mode defining the meaning of `control_value`.
     !> Can be `simulationtime`, `tsteps`, `nsamples` or `never`.
     character(len=:), allocatable :: control_mode
     !> Defines the frequency of writes.
     real(kind=dp) :: control_value = 0.0_dp
     !> Time the schedule is anchored to. Zero for `simulationtime`, so that
     !! the executions land on whole multiples of the interval, and the start
     !! time for `nsamples`, which divides the simulated interval.
     real(kind=dp) :: anchor_time = 0.0_dp
     !> Whether an execution is scheduled at the start time itself, on top of
     !! the ones the anchored schedule prescribes. Should be `.false.` for
     !! outputs for which an execution at the very first time step is
     !! meaningless, such as checkpoints and running statistics.
     logical :: write_at_start = .true.
     !> Index, relative to the anchor, of the first scheduled execution.
     integer(kind=i8) :: first_index = 0
     !> Whether the start time is itself one of the anchored schedule times.
     logical :: start_is_scheduled = .false.
     !> Whether the execution at the start time is still to be performed.
     logical :: start_pending = .false.
     !> Index of the next scheduled execution, counted from `first_index`.
     integer :: next_index = 0
     !> Direction of time, +1 for a forward and -1 for a backward run.
     real(kind=dp) :: direction = 1.0_dp
     !> Value of `tstep` at which the current run started. Only used by the
     !! `tsteps` control mode, for which the schedule cannot be anchored in
     !! absolute time.
     integer :: tstep_offset = 0
     !> Time step at which an execution was last registered, -1 if none.
     !! Guarantees that at most one execution is performed per time step,
     !! also when a forced execution coincides with a scheduled one.
     integer :: last_tstep = -1
     !> Index the cursor is moved to when the pending execution is
     !! registered, together with the time step it was computed for.
     !! @note Private handshake between `check` and `register_execution`.
     integer :: pending_index = -1
     integer :: pending_tstep = -1

   contains
     !> Constructor.
     procedure, pass(this) :: init => time_based_controller_init
     !> Destructor.
     procedure, pass(this) :: free => time_based_controller_free
     !> Check if the execution should be performed.
     procedure, pass(this) :: check => time_based_controller_check
     !> Increment `nexecutions` and advance the schedule.
     procedure, pass(this) :: register_execution => &
          time_based_controller_register_execution
     !> Move the schedule cursor to the current time (for restarts).
     procedure, pass(this) :: set_counter => &
          time_based_controller_set_counter
     !> The time of the next scheduled execution.
     procedure, pass(this) :: next_time => time_based_controller_next_time
     !> The tolerance used when comparing times.
     procedure, pass(this) :: tolerance => time_based_controller_tolerance
  end type time_based_controller_t

contains

  !> Constructor.
  !! @param start_time The time at which the schedule starts. Both the first
  !! possible execution time and the anchor of the schedule.
  !! @param end_time The final simulation time.
  !! @param control_mode The way to interpret the `control_value` parameter.
  !! @param control_value The value defining the execution frequency.
  !! @param write_at_start Whether an execution is scheduled at `start_time`
  !! itself, on top of the ones the anchored schedule prescribes. Optional,
  !! defaults to `.true.`.
  !! @param anchor_time The time the schedule is anchored to. Optional,
  !! defaults to zero for `simulationtime` and to `start_time` for
  !! `nsamples`, which divides the simulated interval rather than tiling it.
  !! @param direction The direction of time of the simulation, +1 forwards
  !! and -1 backwards. Optional, defaults to the direction from `start_time`
  !! to `end_time`, which is the direction of the simulation whenever the
  !! output covers the whole of it. Pass it for an output that starts later,
  !! since the two directions differ when its start time lies beyond the end
  !! of the simulation.
  subroutine time_based_controller_init(this, start_time, end_time, &
       control_mode, control_value, write_at_start, anchor_time, direction)
    class(time_based_controller_t), intent(inout) :: this
    real(kind=dp), intent(in) :: start_time
    real(kind=dp), intent(in) :: end_time
    character(len=*), intent(in) :: control_mode
    real(kind=dp), intent(in) :: control_value
    logical, intent(in), optional :: write_at_start
    real(kind=dp), intent(in), optional :: anchor_time
    real(kind=dp), intent(in), optional :: direction
    real(kind=dp) :: span, offset

    call this%free()

    this%start_time = start_time
    this%end_time = end_time
    this%control_mode = control_mode
    this%control_value = control_value

    if (present(write_at_start)) then
       this%write_at_start = write_at_start
    else
       this%write_at_start = .true.
    end if

    ! The schedule is expressed in terms of the progress made in the
    ! direction the simulation marches in, which is positive also for a
    ! simulation marching backwards in time.
    if (present(direction)) then
       this%direction = sign(1.0_dp, direction)
    else
       this%direction = sign(1.0_dp, end_time - start_time)
    end if
    span = this%direction * (end_time - start_time)

    if (trim(control_mode) .eq. 'simulationtime') then
       if (control_value .le. 0.0_dp) then
          call neko_error("The output interval must be positive")
       end if
       this%time_interval = control_value
       this%frequency = 1.0_dp / this%time_interval
       this%nsteps = 0
    else if (trim(control_mode) .eq. 'nsamples') then
       if (control_value .le. 0.0_dp) then
          call neko_error("nsamples must be positive")
       end if
       if (span .le. 0.0_dp) then
          call neko_error("nsamples requires the output to start before the &
          &end of the simulation")
       end if

       this%frequency = control_value / span
       this%time_interval = 1.0_dp / this%frequency
       this%nsteps = 0
       ! The samples divide the simulated interval, so they are counted from
       ! its start rather than from a global grid.
       this%anchor_time = start_time
    else if (trim(control_mode) .eq. 'tsteps') then
       if (control_value .lt. 1.0_dp) then
          call neko_error("The output interval in time steps must be at &
          &least 1")
       end if
       this%nsteps = int(control_value)
       ! if the timestep will be variable, we cannot compute these.
       this%frequency = 0.0_dp
       this%time_interval = 0.0_dp
    else if (trim(control_mode) .eq. 'never') then
       this%never = .true.
    else
       call neko_error("The control parameter must be simulationtime, nsamples&
       & tsteps, or never, but received "//trim(control_mode))
    end if

    if (present(anchor_time)) this%anchor_time = anchor_time

    ! Point the cursor at the first scheduled execution, and work out whether
    ! the start time is one of the scheduled times or a point of its own.
    this%first_index = 0
    this%start_is_scheduled = .true.
    this%start_pending = .false.

    if (this%time_interval .gt. 0.0_dp) then
       offset = this%direction * (start_time - this%anchor_time) / &
            this%time_interval
       if (abs(offset) .gt. MAX_ANCHOR_INDEX) then
          ! The interval is too fine to resolve against this anchor.
          this%anchor_time = start_time
          offset = 0.0_dp
       end if
       this%first_index = ceiling(offset - SPAN_TOL * max(1.0_dp, &
            abs(offset)), kind = i8)
       this%start_is_scheduled = abs(real(this%first_index, dp) - offset) &
            .le. SPAN_TOL * max(1.0_dp, abs(offset))
       if (this%start_is_scheduled .and. .not. this%write_at_start) then
          this%first_index = this%first_index + 1_i8
       end if
       this%start_pending = this%write_at_start .and. &
            .not. this%start_is_scheduled
    end if

    if (this%nsteps .gt. 0 .and. .not. this%write_at_start) then
       ! A step based schedule is anchored at the first step of the run,
       ! which is the counterpart of the start time here.
       this%next_index = 1
    else
       this%next_index = 0
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
    this%control_value = 0.0_dp
    this%write_at_start = .true.
    this%anchor_time = 0.0_dp
    this%first_index = 0
    this%start_is_scheduled = .false.
    this%start_pending = .false.
    this%next_index = 0
    this%direction = 1.0_dp
    this%tstep_offset = 0
    this%last_tstep = -1
    this%pending_index = -1
    this%pending_tstep = -1
  end subroutine time_based_controller_free

  !> Check if the execution should be performed.
  !! @param time The current time state.
  !! @param force Whether to execute irrespective of the schedule. Optional,
  !! defaults to `.false.`.
  !! @details The result is `.true.` at the first time step at which the next
  !! scheduled time has been reached, up to a tolerance of `TIME_TOL * dt`.
  !! A forced execution is performed unless the simulation has not reached the
  !! start time of the schedule, or an execution has already been registered
  !! for the current time step. The latter is what keeps the final, forced,
  !! execution of a simulation from duplicating a scheduled one that lands on
  !! the very same step. Forcing does override a `never` control, which is
  !! how `output_at_end` writes an output that is otherwise never written.
  function time_based_controller_check(this, time, force) result(check)
    class(time_based_controller_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    logical, intent(in), optional :: force
    logical :: check
    logical :: ifforce
    real(kind=dp) :: progress, tol, t_next, t_start, t_end
    integer :: nstep

    if (present(force)) then
       ifforce = force
    else
       ifforce = .false.
    end if

    check = .false.

    ! At most one execution per time step.
    if (this%last_tstep .eq. time%tstep) return

    ! Nothing is scheduled, but an execution can still be forced.
    if (this%never .and. .not. ifforce) return

    progress = this%direction * (time%t - this%anchor_time)
    t_start = this%direction * (this%start_time - this%anchor_time)
    tol = this%tolerance(time%dt)

    ! Nothing is ever executed before the start of the schedule.
    if (progress .lt. t_start - tol) return

    if (ifforce) then
       check = .true.
    else if (this%nsteps .gt. 0) then
       nstep = time%tstep - this%tstep_offset
       check = nstep .ge. this%next_index * this%nsteps
    else if (this%start_pending) then
       ! The execution at the start time, which the anchored schedule does
       ! not prescribe by itself.
       check = .true.
    else
       t_end = this%direction * (this%end_time - this%anchor_time)
       t_next = real(this%first_index + this%next_index, dp) * &
            this%time_interval
       ! The schedule stops at end_time. Note that this is a condition on the
       ! *scheduled* time and not on the current time: a step that overshoots
       ! end_time still performs the execution scheduled for end_time.
       if (t_next .gt. t_end + SPAN_TOL * max(abs(t_end), &
            this%time_interval)) return
       check = progress .ge. t_next - tol
    end if

    if (check) then
       this%pending_index = next_index_after(this, time)
       this%pending_tstep = time%tstep
    end if

  end function time_based_controller_check

  !> The index of the first scheduled execution that lies strictly ahead of
  !! the current time state.
  !! @note All scheduled times that have already been reached are skipped, so
  !! that a run using a time step larger than the output interval does not
  !! accumulate a backlog of executions.
  pure function next_index_after(this, time) result(index)
    class(time_based_controller_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    integer :: index
    real(kind=dp) :: progress, tol

    if (this%nsteps .gt. 0) then
       index = (time%tstep - this%tstep_offset) / this%nsteps + 1
    else if (this%time_interval .gt. 0.0_dp) then
       progress = this%direction * (time%t - this%anchor_time)
       tol = this%tolerance(time%dt)
       index = int(floor((progress + tol) / this%time_interval, kind = i8) &
            - this%first_index) + 1
    else
       ! Nothing is scheduled, so there is nothing to move past.
       index = this%next_index
    end if

    ! No clamping against `next_index` is done on purpose: a scheduled
    ! execution always yields an index larger than the current one, while a
    ! forced execution that happens before the next scheduled time must
    ! leave the schedule untouched.
    index = max(index, 0)

  end function next_index_after

  !> Increment `nexecutions` and advance the schedule past the current time.
  !! @param time The current time state. Optional, but should be passed
  !! whenever it is available: without it the controller falls back on the
  !! state recorded by the preceding call to `check`.
  subroutine time_based_controller_register_execution(this, time)
    class(time_based_controller_t), intent(inout) :: this
    type(time_state_t), intent(in), optional :: time

    this%nexecutions = this%nexecutions + 1
    this%start_pending = .false.

    if (present(time)) then
       this%next_index = next_index_after(this, time)
       this%last_tstep = time%tstep
    else if (this%pending_tstep .ge. 0) then
       this%next_index = this%pending_index
       this%last_tstep = this%pending_tstep
    else
       this%next_index = this%next_index + 1
    end if

    this%pending_index = -1
    this%pending_tstep = -1

  end subroutine time_based_controller_register_execution

  !> Move the schedule cursor to the current time, and set `nexecutions` to
  !! the number of executions that the schedule has already prescribed.
  !! Called when restarting a simulation.
  !! @param time The current time.
  !! @details The predicate used here is exactly the one used by `check`, so
  !! that the executions performed by the run that wrote the checkpoint are
  !! the ones considered done: nothing is repeated and, more importantly,
  !! nothing is skipped.
  subroutine time_based_controller_set_counter(this, time)
    class(time_based_controller_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=dp) :: progress, tol, dt, t_start
    integer :: n_passed

    if (this%never) return

    if (this%nsteps .gt. 0) then
       ! `tstep` is not stored in the checkpoint, so a step based schedule
       ! cannot be reconstructed. Restart the cadence from the restart step,
       ! without executing at the restart step itself.
       this%tstep_offset = time%tstep
       this%next_index = 1
       this%start_pending = .false.
       this%last_tstep = time%tstep
       return
    end if

    ! The size of the step that produced the checkpoint, which is the one
    ! `check` would have used, and not necessarily the one of the new run.
    dt = time%dt
    if (abs(time%dtlag(1)) .gt. 0.0_dp) dt = time%dtlag(1)

    progress = this%direction * (time%t - this%anchor_time)
    t_start = this%direction * (this%start_time - this%anchor_time)
    tol = this%tolerance(dt)

    if (progress .lt. t_start - tol) then
       ! The schedule has not started yet, so neither has the run.
       this%next_index = 0
       this%nexecutions = 0
    else
       n_passed = int(floor((progress + tol) / this%time_interval, &
            kind = i8) - this%first_index) + 1
       this%next_index = max(n_passed, 0)
       ! The execution at the start time, where the schedule prescribes one
       ! of its own, was performed by the run that reached this time.
       this%nexecutions = this%next_index
       if (this%start_pending) this%nexecutions = this%nexecutions + 1
       this%start_pending = .false.
    end if

    ! Whatever was scheduled for the restart time has been done by the run
    ! that wrote the checkpoint.
    this%last_tstep = time%tstep

  end subroutine time_based_controller_set_counter

  !> The tolerance used when comparing the simulation time to a scheduled
  !! execution time.
  !! @param dt The size of the time step under consideration.
  pure function time_based_controller_tolerance(this, dt) result(tol)
    class(time_based_controller_t), intent(in) :: this
    real(kind=dp), intent(in) :: dt
    real(kind=dp) :: tol

    tol = TIME_TOL * abs(dt)
    if (this%time_interval .gt. 0.0_dp) then
       tol = min(tol, TIME_TOL * this%time_interval)
    end if

  end function time_based_controller_tolerance

  !> The time of the next scheduled execution.
  !! @note Meaningless for the `tsteps` and `never` control modes, for which
  !! the end time of the simulation is returned.
  pure function time_based_controller_next_time(this) result(t)
    class(time_based_controller_t), intent(in) :: this
    real(kind=dp) :: t

    if (this%never .or. this%nsteps .gt. 0) then
       t = this%end_time
    else if (this%start_pending) then
       t = this%start_time
    else
       t = this%anchor_time + this%direction * &
            real(this%first_index + this%next_index, dp) * this%time_interval
    end if

  end function time_based_controller_next_time

end module time_based_controller
