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

  !> Tolerance for comparing the simulation time to a scheduled time, as a
  !! fraction of the time step. A scheduled time `t_k` is considered reached
  !! at the first step with `t >= t_k - TIME_TOL * dt`, so that round-off in
  !! the accumulation of `t` does not postpone the execution by a full step.
  !! The tolerance is limited to the same fraction of the output interval, so
  !! that a time step longer than the interval (in particular the `dt = 1`
  !! placeholder a variable time step run starts from) cannot make a
  !! scheduled time count as reached before the previous one has passed.
  real(kind=dp), public, parameter :: TIME_TOL = 0.1_dp

  !> Relative tolerance for round-off in `k * time_interval` when comparing a
  !! scheduled time to `start_time` and `end_time`.
  real(kind=dp), public, parameter :: SPAN_TOL = 1.0e-9_dp

  !> A utility type for determining whether an action should be executed based
  !! on the current time value. Used to e.g. control whether we should write a
  !! file or execute a simcomp.
  !!
  !! @details
  !! The controller defines a schedule: the times (or time steps) at which
  !! the action is executed. For the time based control modes these are
  !! \f$ t_k = t_{anchor} + k \Delta t_{out} \f$ for the \f$ k \f$ with
  !! \f$ t_k \f$ in \f$ [t_{start}, t_{end}] \f$ and, if `write_at_start`
  !! is set, \f$ t_{start} \f$ itself.
  !!
  !! A `simulationtime` schedule has `anchor_time = 0`, so that the
  !! executions are at whole multiples of the interval regardless of the
  !! start time. An `nsamples` schedule divides the interval from
  !! `start_time` to `end_time`, so its anchor is `start_time`.
  !!
  !! The schedule depends only on the case parameters, not on the number of
  !! executions performed or on the time the run was started from. On a
  !! restart, `set_counter` sets `next_index` to the first scheduled time not
  !! yet reached, using the same comparison as `check`, so that no execution
  !! is repeated or skipped.
  !!
  !! `check` does not modify the controller and returns `.false.` for a time
  !! step at which an execution has already been registered, so it can be
  !! called several times per step. `register_execution` increments
  !! `nexecutions` and advances `next_index`.
  type, public :: time_based_controller_t
     !> Frequency of execution.
     real(kind=dp) :: frequency = 0.0_dp
     !> Time interval between executions. Always positive.
     real(kind=dp) :: time_interval = 0.0_dp
     !> Number of time steps in between executions.
     integer :: nsteps = 0
     !> First time at which an execution can be performed.
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
     !> Time the scheduled times are counted from. Zero for `simulationtime`
     !! and `start_time` for `nsamples`.
     real(kind=dp) :: anchor_time = 0.0_dp
     !> Whether an execution is performed at `start_time`, whether or not it
     !! is one of the scheduled times. Should be `.false.` for outputs for
     !! which an execution at the very first time step is meaningless, such as
     !! checkpoints and running statistics.
     logical :: write_at_start = .true.
     !> Index `k` of the first scheduled time.
     integer(kind=i8) :: first_index = 0
     !> Whether `start_time` is one of the scheduled times.
     logical :: start_is_scheduled = .false.
     !> Whether the execution at `start_time` is still to be performed. Only
     !! set when `write_at_start` is true and `start_time` is not one of the
     !! scheduled times.
     logical :: start_pending = .false.
     !> Number of scheduled times already passed. The next scheduled time is
     !! `t_k` with `k = first_index + next_index`.
     integer :: next_index = 0
     !> Value of `tstep` at which the current run started. Only used by the
     !! `tsteps` control mode, whose schedule is relative to the run.
     integer :: tstep_offset = 0
     !> Time step at which an execution was last registered, -1 if none.
     !! Guarantees that at most one execution is performed per time step,
     !! also when a forced execution coincides with a scheduled one.
     integer :: last_tstep = -1

   contains
     !> Constructor.
     procedure, pass(this) :: init => time_based_controller_init
     !> Destructor.
     procedure, pass(this) :: free => time_based_controller_free
     !> Check if the execution should be performed.
     procedure, pass(this) :: check => time_based_controller_check
     !> Increment `nexecutions` and advance `next_index`.
     procedure, pass(this) :: register_execution => &
          time_based_controller_register_execution
     !> Set `next_index` and `nexecutions` for the current time (restarts).
     procedure, pass(this) :: set_counter => &
          time_based_controller_set_counter
     !> The time of the next scheduled execution.
     procedure, pass(this) :: next_time => time_based_controller_next_time
     !> The tolerance used when comparing times.
     procedure, pass(this) :: tolerance => time_based_controller_tolerance
  end type time_based_controller_t

contains

  !> Constructor.
  !! @param start_time The first time at which an execution can be performed.
  !! @param end_time The final simulation time.
  !! @param control_mode The way to interpret the `control_value` parameter.
  !! @param control_value The value defining the execution frequency.
  !! @param write_at_start Whether an execution is performed at `start_time`,
  !! whether or not it is one of the scheduled times. Optional, defaults to
  !! `.true.`.
  !! @param anchor_time The time the scheduled times are counted from.
  !! Optional, defaults to zero for `simulationtime` and to `start_time` for
  !! `nsamples`.
  subroutine time_based_controller_init(this, start_time, end_time, &
       control_mode, control_value, write_at_start, anchor_time)
    class(time_based_controller_t), intent(inout) :: this
    real(kind=dp), intent(in) :: start_time
    real(kind=dp), intent(in) :: end_time
    character(len=*), intent(in) :: control_mode
    real(kind=dp), intent(in) :: control_value
    logical, intent(in), optional :: write_at_start
    real(kind=dp), intent(in), optional :: anchor_time
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

    span = end_time - start_time

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
       ! The samples divide the interval from start_time to end_time, so
       ! they are counted from start_time.
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

    ! Find the index of the first scheduled time, and whether start_time is
    ! itself one of the scheduled times.
    this%first_index = 0
    this%start_is_scheduled = .true.
    this%start_pending = .false.

    if (this%time_interval .gt. 0.0_dp) then
       offset = (start_time - this%anchor_time) / this%time_interval
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
       ! The first scheduled step of a step based schedule is the first step
       ! of the run. Skip it, as start_time is skipped above.
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
    this%tstep_offset = 0
    this%last_tstep = -1
  end subroutine time_based_controller_free

  !> Check if the execution should be performed.
  !! @param time The current time state.
  !! @param force Whether to execute irrespective of the schedule. Optional,
  !! defaults to `.false.`.
  !! @details The result is `.true.` at the first time step at which the next
  !! scheduled time has been reached, within a tolerance of `TIME_TOL * dt`.
  !! A forced execution is performed unless the time is before `start_time`,
  !! or an execution has already been registered for the current time step.
  !! The latter keeps the forced execution at the end of a simulation from
  !! repeating a scheduled one at the same step. Forcing does override a
  !! `never` control, which is how `output_at_end` writes an output that is
  !! otherwise never written.
  function time_based_controller_check(this, time, force) result(check)
    class(time_based_controller_t), intent(in) :: this
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

    progress = time%t - this%anchor_time
    t_start = this%start_time - this%anchor_time
    tol = this%tolerance(time%dt)

    ! Nothing is executed before start_time.
    if (progress .lt. t_start - tol) return

    if (ifforce) then
       check = .true.
    else if (this%nsteps .gt. 0) then
       nstep = time%tstep - this%tstep_offset
       check = nstep .ge. this%next_index * this%nsteps
    else if (this%start_pending) then
       ! The execution at start_time, which is not one of the scheduled
       ! times.
       check = .true.
    else
       t_end = this%end_time - this%anchor_time
       t_next = real(this%first_index + this%next_index, dp) * &
            this%time_interval
       ! Nothing is scheduled after end_time. The condition is on the
       ! scheduled time and not on the current time, so a step that
       ! overshoots end_time still performs the execution scheduled for it.
       if (t_next .gt. t_end + SPAN_TOL * max(abs(t_end), &
            this%time_interval)) return
       check = progress .ge. t_next - tol
    end if

  end function time_based_controller_check

  !> The number of scheduled times at or before the given time, counted from
  !! `first_index`. Assigning it to `next_index` after an execution makes the
  !! next scheduled time the first one after the given time.
  !! @note Scheduled times passed within a single step are all counted, so a
  !! time step larger than the output interval gives one execution per step
  !! and not one per passed scheduled time.
  pure function next_index_after(this, time) result(index)
    class(time_based_controller_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    integer :: index
    real(kind=dp) :: progress, tol

    if (this%nsteps .gt. 0) then
       index = (time%tstep - this%tstep_offset) / this%nsteps + 1
    else if (this%time_interval .gt. 0.0_dp) then
       progress = time%t - this%anchor_time
       tol = this%tolerance(time%dt)
       index = int(floor((progress + tol) / this%time_interval, kind = i8) &
            - this%first_index) + 1
    else
       ! Nothing is scheduled.
       index = this%next_index
    end if

    ! The result is not limited from below by the current `next_index` on
    ! purpose: a scheduled execution always gives a larger index, and a
    ! forced execution before the next scheduled time must leave
    ! `next_index` unchanged.
    index = max(index, 0)

  end function next_index_after

  !> Increment `nexecutions` and advance `next_index` past the current time.
  !! @param time The current time state.
  subroutine time_based_controller_register_execution(this, time)
    class(time_based_controller_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    this%nexecutions = this%nexecutions + 1
    this%start_pending = .false.
    this%next_index = next_index_after(this, time)
    this%last_tstep = time%tstep

  end subroutine time_based_controller_register_execution

  !> Set `next_index` to the first scheduled time not yet reached, and
  !! `nexecutions` to the number of executions performed up to the current
  !! time. Called when restarting a simulation.
  !! @param time The current time.
  !! @details The comparison is the same as in `check`, so that exactly the
  !! executions performed by the run that wrote the checkpoint are counted as
  !! done, and none is repeated or skipped.
  subroutine time_based_controller_set_counter(this, time)
    class(time_based_controller_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=dp) :: progress, tol, dt, t_start, t_first
    integer :: n_passed

    if (this%never) return

    if (this%nsteps .gt. 0) then
       ! `tstep` is not stored in the checkpoint, so a step based schedule
       ! cannot be continued. Start a new one from the restart step, without
       ! executing at the restart step itself.
       this%tstep_offset = time%tstep
       this%next_index = 1
       this%start_pending = .false.
       this%last_tstep = time%tstep
       return
    end if

    ! The size of the step that produced the checkpoint, which is the step
    ! `check` used at that time, and not the step of the new run.
    dt = time%dt
    if (abs(time%dtlag(1)) .gt. 0.0_dp) dt = time%dtlag(1)

    progress = time%t - this%anchor_time
    t_start = this%start_time - this%anchor_time
    tol = this%tolerance(dt)

    if (progress .lt. t_start - tol) then
       ! The run has not reached start_time.
       this%next_index = 0
       this%nexecutions = 0
    else
       n_passed = int(floor((progress + tol) / this%time_interval, &
            kind = i8) - this%first_index) + 1
       this%next_index = max(n_passed, 0)
       this%nexecutions = this%next_index
       if (this%start_pending) then
          ! The execution at start_time was performed by the previous run.
          ! It is counted separately unless the first scheduled time was
          ! within tolerance of start_time, in which case the same execution
          ! was registered for both (see next_index_after).
          t_first = real(this%first_index, dp) * this%time_interval
          if (t_first .gt. t_start + tol) then
             this%nexecutions = this%nexecutions + 1
          end if
       end if
       this%start_pending = .false.
    end if

    ! Anything scheduled for the restart time was executed by the run that
    ! wrote the checkpoint.
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
       t = this%anchor_time + &
            real(this%first_index + this%next_index, dp) * this%time_interval
    end if

  end function time_based_controller_next_time

end module time_based_controller
