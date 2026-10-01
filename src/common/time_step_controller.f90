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
  use num_types, only : dp, i8
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get_or_lookup_or_default
  use time_state, only : time_state_t
  use comm, only : NEKO_GLOBAL_COMM, is_mpmd
  use mpi_f08, only : MPI_MIN, MPI_MAX, MPI_IN_PLACE, MPI_Allreduce, &
       MPI_DOUBLE_PRECISION, MPI_INTEGER
  implicit none
  private

  !> Relative tolerance on the time step when landing on a scheduled time.
  !! A remaining time that is a whole number of steps within this fraction of
  !! a step is taken as such, so that round-off in the accumulation of the
  !! time does not split it into one step more, and a change of the step by
  !! less than this fraction is not registered as a change.
  real(kind=dp), public, parameter :: LANDING_TOL = 1.0e-6_dp

  !> Largest number of steps up to a scheduled time that is still fitted.
  !! Beyond it the quantisation of the step is below round-off anyway.
  real(kind=dp), parameter :: MAX_LANDING_STEPS = 1.0e15_dp

  !> Number of steps over which the CFL controller may ask in vain for a
  !! change that the landing on the scheduled times denies before the run
  !! stops: the scheduled times are then too dense to follow the controller.
  integer, parameter :: HELD_STEPS_MAX = 100

  !> Provides a tool to set time step dt
  type, public :: time_step_controller_t
     logical :: is_variable_dt = .false.
     real(kind=dp) :: cfl_trg = 0.0_dp
     real(kind=dp) :: cfl_avg = 0.0_dp
     real(kind=dp) :: init_dt = huge(0.0_dp)
     real(kind=dp) :: max_dt = 0.0_dp
     real(kind=dp) :: min_dt = 0.0_dp
     integer :: max_update_frequency = 0
     integer :: min_update_frequency = 0
     !> Number of time steps since the time step last changed, 0 at the step
     !! it changed, -1 until the first time step.
     integer :: dt_last_change = -1
     real(kind=dp) :: alpha = 0.0_dp !< coefficient of running average
     real(kind=dp) :: max_dt_increase_factor = 0.0_dp
     real(kind=dp) :: min_dt_decrease_factor = 0.0_dp
     real(kind=dp) :: dev_tol = 0.0_dp
     !> Whether to adjust the time step so that the scheduled sampling and
     !! output times are reached exactly. Requires a variable time step.
     logical :: exact_output_time = .false.
     !> The time step taken last, which the one about to be taken is bounded
     !! by and compared against.
     real(kind=dp) :: dt_previous = 0.0_dp
     !> The CFL number of the time step taken last.
     real(kind=dp) :: cfl_previous = 0.0_dp
     !> `dt_last_change` on entry to `set_dt`, to undo a change the
     !! controller registered for a step that the landing left unchanged.
     integer :: dt_last_change_prev = -1
     !> Whether the step about to be taken is fitted to a scheduled time.
     logical :: dt_fitted = .false.
     !> The step planned after the one about to be taken, when the fit
     !! changes the step in two stages, and the time it lands on.
     real(kind=dp) :: dt_plateau = 0.0_dp
     real(kind=dp) :: plateau_target = huge(0.0_dp)
     !> Number of steps over which the controller has asked in vain for a
     !! change that only the bounds on the ratio of the steps deny.
     integer :: held_steps = 0
     !> Whether the last change asked for was denied by those bounds.
     logical :: held_by_ratio = .false.
     !> A scheduled time so close to the end that the end is landed on
     !! instead, see `land`.
     real(kind=dp) :: skipped_for_end = huge(0.0_dp)
     !> The time being landed on, or passed, by the step about to be taken.
     real(kind=dp) :: landing_target = huge(0.0_dp)
   contains
     !> Initialize object.
     procedure, pass(this) :: init => time_step_controller_init
     !> Set time stepping
     procedure, pass(this) :: set_dt => time_step_controller_set_dt
     !> The first time step of a variable time step run.
     procedure, pass(this) :: first_dt => time_step_controller_first_dt
     !> Adjust the time step to land on the next scheduled time.
     procedure, pass(this) :: land => time_step_controller_land
     !> The shortest interval, in steps, that is landed on.
     procedure, pass(this) :: min_interval_steps => &
          time_step_controller_min_interval_steps
     !> The shortest interval, in time, that is landed on.
     procedure, pass(this) :: landing_min_interval => &
          time_step_controller_landing_min_interval
     !> The CFL number at the geometric centre of the controller's band.
     procedure, pass(this) :: cfl_centre => time_step_controller_cfl_centre
     !> Divide the time up to a scheduled time into steps within the bounds.
     procedure, pass(this) :: fit_step => time_step_controller_fit_step

  end type time_step_controller_t

contains

  !> Constructor
  !! @param order order of the interpolation
  subroutine time_step_controller_init(this, params)
    class(time_step_controller_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    integer :: flags_any(2), flags_all(2), ierr

    this%dt_last_change = -1
    call json_get_or_default(params, 'variable_timestep', &
         this%is_variable_dt, .false.)
    if (this%is_variable_dt) then
       call json_get_or_lookup_or_default(params, 'target_cfl', &
            this%cfl_trg, 0.4_dp)
       call json_get_or_lookup_or_default(params, 'timestep', &
            this%init_dt, huge(0.0_dp))
       call json_get_or_lookup_or_default(params, 'max_timestep', &
            this%max_dt, huge(0.0_dp))
       call json_get_or_lookup_or_default(params, 'min_timestep', &
            this%min_dt, 0.0_dp)
       call json_get_or_lookup_or_default(params, 'max_update_frequency',&
            this%max_update_frequency, 0)
       call json_get_or_lookup_or_default(params, 'min_update_frequency',&
            this%min_update_frequency, huge(0))
       call json_get_or_lookup_or_default(params, 'cfl_running_avg_coeff', &
            this%alpha, 0.5_dp)
       call json_get_or_lookup_or_default(params, 'max_dt_increase_factor', &
            this%max_dt_increase_factor, 1.2_dp)
       call json_get_or_lookup_or_default(params, 'min_dt_decrease_factor', &
            this%min_dt_decrease_factor, 0.5_dp)
       call json_get_or_lookup_or_default(params, 'cfl_deviation_tolerance', &
            this%dev_tol, 0.2_dp)
    end if

    ! Landing on the scheduled times adjusts the step within the bounds of
    ! the CFL controller. With a fixed step it only checks the schedules.
    call json_get_or_default(params, 'exact_output_time', &
         this%exact_output_time, .false.)
    if (this%exact_output_time .and. .not. this%is_variable_dt) then
       call neko_warning('exact_output_time with a fixed timestep only &
       &checks that the sampling and output times are whole numbers &
       &of steps away, and stops the run otherwise. Use tsteps to &
       &sample every so many steps, or variable_timestep to have the &
       &step fitted.')
    end if

    ! A variable time step takes a collective at every step in an MPMD run,
    ! and the landing on the scheduled times another one, so the coupled
    ! simulations have to agree on both, or their collectives would not
    ! match up.
    if (is_mpmd()) then
       flags_any = merge(1, 0, [this%is_variable_dt, this%exact_output_time])
       flags_all = flags_any
       call MPI_Allreduce(MPI_IN_PLACE, flags_any, 2, MPI_INTEGER, &
            MPI_MAX, NEKO_GLOBAL_COMM, ierr)
       call MPI_Allreduce(MPI_IN_PLACE, flags_all, 2, MPI_INTEGER, &
            MPI_MIN, NEKO_GLOBAL_COMM, ierr)
       if (any(flags_any .ne. flags_all)) then
          call neko_error('variable_timestep and exact_output_time must &
          &each be set in all the coupled cases of an MPMD run, or in &
          &none')
       end if
    end if

  end subroutine time_step_controller_init

  !> The first time step of a variable time step run: the one giving the
  !! target CFL number, limited by `timestep`, `max_timestep` and
  !! `min_timestep` where given.
  !! @param time The time state, whose `dt` is the step the CFL number was
  !! computed with.
  !! @param cfl The CFL number of that step.
  pure function time_step_controller_first_dt(this, time, cfl) result(dt)
    class(time_step_controller_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    real(kind=dp), intent(in) :: cfl
    real(kind=dp) :: dt

    dt = this%init_dt
    if (cfl .gt. 0.0_dp) dt = min(this%cfl_trg / cfl * time%dt, dt)
    dt = max(min(dt, this%max_dt), this%min_dt)

  end function time_step_controller_first_dt

  !> The shortest interval of a schedule, in time steps, that the landing
  !! divides into steps: the number of steps `n` from which on the step of
  !! `n` steps is within `max_dt_increase_factor` squared of the step of
  !! `n + 1` steps, so that the CFL controller can move the step from one
  !! whole number of steps per interval to the next in the two stages the
  !! fit allows. A shorter interval could only be divided into steps that
  !! the controller cannot climb back out of, so its schedule is executed
  !! at the first step past each scheduled time instead, as it is without
  !! the landing. Three steps with the default factor of 1.2.
  pure function time_step_controller_min_interval_steps(this) result(steps)
    class(time_step_controller_t), intent(in) :: this
    real(kind=dp) :: steps

    steps = huge(0.0_dp)
    if (this%max_dt_increase_factor .gt. 1.0_dp) then
       steps = 1.0_dp / (this%max_dt_increase_factor**2 - 1.0_dp) - &
            LANDING_TOL
       steps = real(ceiling(min(steps, MAX_LANDING_STEPS), kind = i8), dp)
    end if
    steps = max(steps, 2.0_dp)

  end function time_step_controller_min_interval_steps

  !> The CFL number at the geometric centre of the controller's band,
  !! `target_cfl * sqrt(1 - cfl_deviation_tolerance**2)`: the CFL number
  !! the controller is the least likely to move away from.
  pure function time_step_controller_cfl_centre(this) result(cfl_mid)
    class(time_step_controller_t), intent(in) :: this
    real(kind=dp) :: cfl_mid

    cfl_mid = this%cfl_trg * sqrt(max(1.0_dp - this%dev_tol**2, 0.0_dp))
    if (.not. cfl_mid .gt. 0.0_dp) cfl_mid = this%cfl_trg

  end function time_step_controller_cfl_centre

  !> The shortest interval of a schedule that the landing divides into
  !! steps: `min_interval_steps` times the step the controller aims at, the
  !! one that would give the CFL number at the centre of the band, as
  !! estimated from the step taken last and its CFL number, within the
  !! minimum and maximum step. That estimate depends on the flow and the
  !! case only, not on the step actually taken, so a schedule does not flip
  !! in and out of the landing as the step is fitted.
  !! @param time The time state, whose `dt` is used before the first step.
  pure function time_step_controller_landing_min_interval(this, time) &
       result(min_interval)
    class(time_step_controller_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    real(kind=dp) :: min_interval, dt_estimate

    dt_estimate = abs(time%dt)
    if (.not. this%is_variable_dt) then
       min_interval = dt_estimate
       return
    end if
    if (this%cfl_previous .gt. 0.0_dp .and. &
         abs(this%dt_previous) .gt. 0.0_dp) then
       dt_estimate = abs(this%dt_previous) * this%cfl_centre() / &
            this%cfl_previous
    end if
    dt_estimate = max(min(dt_estimate, this%max_dt), this%min_dt)
    min_interval = this%min_interval_steps() * dt_estimate

  end function time_step_controller_landing_min_interval

  !> Set new dt based on cfl if requested
  !! @param time The time state, whose `dt` is the step taken last on entry
  !! and the step about to be taken on exit.
  !! @param cfl courant number of current iteration.
  !! @Algorithm:
  !! 1. Set the first time step such that cfl is the set one;
  !! 2. During time-stepping, adjust dt when cfl_avg is offset by 20%.
  !! With `exact_output_time` the step is then fitted to the next scheduled
  !! time by `land`, which is to be called right after.
  subroutine time_step_controller_set_dt(this, time, cfl)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=dp), intent(in) :: cfl
    real(kind=dp) :: dt_old, scaling_factor, global_min_dt
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr

    this%dt_previous = time%dt
    this%cfl_previous = cfl
    this%dt_last_change_prev = this%dt_last_change

    if (this%is_variable_dt) then

       ! Reset the average cfl if it is the first time step since the last
       ! change
       if (this%dt_last_change .eq. 0) then
          this%cfl_avg = cfl
       end if

       if (this%dt_last_change .eq. -1) then

          ! Set the first dt for desired cfl, or use the provided initial dt if
          ! it is smaller. Then clamp between max and min dt if provided.
          time%dt = this%first_dt(time, cfl)
          this%dt_last_change = 0
          this%cfl_avg = cfl
          ! There is no step before the first one, which the landing then
          ! bounds around itself rather than around the placeholder.
          this%dt_previous = time%dt

       else
          ! Calculate the average of cfl over the desired interval
          this%cfl_avg = this%alpha * cfl + (1 - this%alpha) * this%cfl_avg

          if (abs(this%cfl_avg - this%cfl_trg) .ge. this%dev_tol*this%cfl_trg &
               .and. this%dt_last_change .ge. this%max_update_frequency &
               .or. this%dt_last_change .ge. this%min_update_frequency) then

             if (this%cfl_trg/cfl .ge. 1) then
                ! increase of time step
                scaling_factor = min(this%max_dt_increase_factor, &
                     this%cfl_trg/cfl)
             else
                ! reduction of time step
                scaling_factor = max(this%min_dt_decrease_factor, &
                     this%cfl_trg/cfl)
             end if

             dt_old = time%dt
             time%dt = scaling_factor * dt_old
             time%dt = max(min(time%dt, this%max_dt), this%min_dt)

             ! With the landing on the scheduled times, the step is fitted
             ! and reported by `land`.
             if (.not. this%exact_output_time) then
                write(log_buf, '(A,E15.7,1x,A,E15.7)') &
                     'Average CFL:', this%cfl_avg, &
                     'Target  CFL:', this%cfl_trg
                call neko_log%message(log_buf)

                write(log_buf, '(A,E15.7,1x,A,E15.7)') 'Old dt:', dt_old, &
                     'New dt:', time%dt
                call neko_log%message(log_buf)
             end if

             this%dt_last_change = 0

          else
             this%dt_last_change = this%dt_last_change + 1
          end if
       end if

       ! If running in mpmd, the new dt is the minimum across simulations
       if (is_mpmd()) then
          global_min_dt = time%dt
          call MPI_Allreduce(MPI_IN_PLACE, global_min_dt, 1, &
               MPI_DOUBLE_PRECISION, MPI_MIN, NEKO_GLOBAL_COMM, ierr)

          ! If my dt is larger that the global min, mark a change
          if (time%dt .gt. global_min_dt) then
             time%dt = global_min_dt
             this%dt_last_change = 0
          end if
       end if

    end if

  end subroutine time_step_controller_set_dt

  !> Adjust the time step so that the next scheduled sampling or output
  !! time, and the end of the simulation, are reached exactly.
  !! @param time The time state, whose `dt` is the step the CFL controller
  !! chose for the step about to be taken, and is adjusted in place.
  !! @param time_to_next The time until the next scheduled time, as
  !! `time_to_next` of the controllers gives it when called with
  !! `landing_min_interval` of this controller, `huge(0.0_dp)` if there is
  !! none. The end of the simulation is a target of its own and needs not
  !! be included.
  !! @details To be called right after `set_dt`. The time remaining up to
  !! the next scheduled time, or to `end_time` if that is closer, is divided
  !! into a whole number of equal steps, see `fit_step`. The step is fitted
  !! when the CFL controller changes it and at each scheduled time, where
  !! the whole interval up to the next one is divided at once; in between
  !! it is kept, as it divides the remaining time already.
  !!
  !! Every step taken is within `min_dt_decrease_factor` and
  !! `max_dt_increase_factor` of the step before, and within `min_timestep`
  !! and `max_timestep`, so the fit is in effect a quantisation of the
  !! controller's own step. A scheduled time that cannot be reached with a
  !! step within these bounds, which with the default factors takes a
  !! scheduled time closer than half a step to the one before it, is passed
  !! and executed at the first step past it, as it is without the option.
  !!
  !! Scheduled times of several schedules that follow each other only a few
  !! steps apart all the time can hold the step at a whole fraction of
  !! their gaps that leaves the CFL number outside the controller's band,
  !! with the controller asking in vain for a change that only the bounds
  !! on the ratio of the steps deny. After `HELD_STEPS_MAX` such steps the
  !! run stops with an error, as the schedules are too dense to follow the
  !! controller.
  !!
  !! With a fixed time step nothing is adjusted: the next scheduled time is
  !! checked to be a whole number of steps away, and the run stops with an
  !! error otherwise.
  !!
  !! In an MPMD run the fitted step is the smallest one over the
  !! simulations, so that they keep advancing in lockstep.
  subroutine time_step_controller_land(this, time, time_to_next)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=dp), intent(in) :: time_to_next
    real(kind=dp) :: dt_prev, dt_new, direction, remaining, global_min_dt
    real(kind=dp) :: lo, hi, remaining_end, target_time, nearest, plateau, x
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr
    logical :: reachable, fitted, changed, on_end, ratio_bound, step_changed

    if (.not. this%exact_output_time) return

    ! A fixed step is only checked against the next scheduled time
    if (.not. this%is_variable_dt) then
       if (time_to_next .lt. huge(0.0_dp) .and. abs(time%dt) .gt. 0.0_dp) &
            then
          x = time_to_next / abs(time%dt)
          if (x .lt. MAX_LANDING_STEPS) then
             if (abs(x - real(nint(x, kind = i8), dp)) .gt. LANDING_TOL) then
                write(log_buf, '(A,E15.7,A,F0.3,A)') 'The scheduled time ', &
                     time%t + sign(1.0_dp, time%dt) * time_to_next, ' is ', &
                     x, ' steps away'
                call neko_error(trim(log_buf) // ', so the fixed timestep &
                &cannot reach it exactly. Use tsteps for that sampling &
                &or output, or disable exact_output_time, before &
                &restarting.')
             end if
          end if
       end if
       return
    end if

    direction = sign(1.0_dp, time%dt)
    dt_prev = abs(this%dt_previous)
    if (.not. dt_prev .gt. 0.0_dp) dt_prev = abs(time%dt)
    dt_new = abs(time%dt)
    plateau = dt_new
    reachable = .true.
    fitted = .false.
    ratio_bound = .false.
    target_time = huge(0.0_dp)
    ! Whether the controller changed the step (or it is the first step)
    changed = this%dt_last_change .eq. 0

    ! The bounds on the step about to be taken, as the CFL controller has
    ! them: within the factors of the step taken last, and within the
    ! minimum and maximum step.
    lo = min(this%min_dt_decrease_factor, 1.0_dp) * dt_prev
    hi = max(this%max_dt_increase_factor, 1.0_dp) * dt_prev
    lo = max(min(lo, this%max_dt), this%min_dt, LANDING_TOL * dt_prev)
    hi = max(min(hi, this%max_dt), this%min_dt, lo)

    ! The next time to land on, the end of the simulation included. The end
    ! is the target whenever it is within the shortest step allowed of the
    ! next scheduled time, which is then executed at the end rather than
    ! landed on and followed by a step too short to be allowed, that would
    ! leave the end to be passed. A scheduled time once skipped for the end
    ! stays skipped, as the shortest step allowed varies with the step.
    remaining_end = direction * (time%end_time - time%t)
    nearest = time%t + direction * time_to_next
    on_end = remaining_end .le. time_to_next + lo
    if (time_to_next .lt. remaining_end .and. time_to_next .lt. huge(0.0_dp)) &
         then
       if (abs(nearest - this%skipped_for_end) .le. LANDING_TOL * dt_prev) &
            on_end = .true.
       if (on_end) this%skipped_for_end = nearest
    end if
    if (on_end) then
       remaining = remaining_end
    else
       remaining = time_to_next
    end if

    if (dt_prev .gt. 0.0_dp .and. remaining .gt. 0.0_dp .and. &
         remaining .lt. huge(0.0_dp)) then
       target_time = time%t + direction * remaining

       ! Keep the step as long as it divides the remaining time and the
       ! controller did not change it, or move on to the step planned after
       ! a first step of a different length. Dividing the remaining time
       ! anew absorbs the round-off of the accumulated time, so the last
       ! step lands exactly.
       if (.not. changed) then
          call keep_step(remaining, dt_prev, dt_new, fitted)
          if (.not. fitted .and. &
               this%dt_plateau .ge. lo * (1.0_dp - LANDING_TOL) .and. &
               this%dt_plateau .le. hi * (1.0_dp + LANDING_TOL) .and. &
               abs(target_time - this%plateau_target) .le. &
               LANDING_TOL * dt_prev) then
             call keep_step(remaining, this%dt_plateau, dt_new, fitted)
          end if
          plateau = dt_new
       end if

       ! Otherwise fit the step anew, or find the scheduled time out of reach
       if (.not. fitted) then
          call this%fit_step(remaining, dt_prev, abs(time%dt), lo, hi, &
               dt_new, plateau, fitted, ratio_bound)
          reachable = fitted
       end if
    end if

    ! If running in mpmd, the fitted step is the minimum across simulations
    if (is_mpmd()) then
       global_min_dt = dt_new
       call MPI_Allreduce(MPI_IN_PLACE, global_min_dt, 1, &
            MPI_DOUBLE_PRECISION, MPI_MIN, NEKO_GLOBAL_COMM, ierr)
       dt_new = global_min_dt
    end if

    time%dt = direction * dt_new
    this%dt_fitted = fitted
    this%dt_plateau = 0.0_dp
    this%plateau_target = huge(0.0_dp)
    if (fitted .and. abs(plateau - dt_new) .gt. LANDING_TOL * dt_new) then
       this%dt_plateau = plateau
       this%plateau_target = target_time
    end if

    ! The projection spaces watch for a change of the step, whatever its
    ! reason. A change by round-off only is none, and a change the
    ! controller asked for that the fit did not make is none either.
    step_changed = abs(dt_new - dt_prev) .gt. LANDING_TOL * dt_prev
    if (step_changed) then
       this%dt_last_change = 0
       write(log_buf, '(A,E15.7,1x,A,E15.7)') 'Old dt:', &
            direction * dt_prev, 'New dt:', time%dt
       call neko_log%message(log_buf)
    else if (this%dt_last_change .eq. 0 .and. &
         this%dt_last_change_prev .ge. 0) then
       this%dt_last_change = this%dt_last_change_prev + 1
    end if

    ! The controller asking in vain, for too long, for a change that only
    ! the bounds on the ratio of the steps deny
    if (step_changed) then
       this%held_steps = 0
       this%held_by_ratio = .false.
    else
       if (changed) this%held_by_ratio = ratio_bound
       if (abs(this%cfl_avg - this%cfl_trg) .ge. this%dev_tol * this%cfl_trg &
            .and. this%held_by_ratio) then
          this%held_steps = this%held_steps + 1
       else if (abs(this%cfl_avg - this%cfl_trg) .lt. &
            this%dev_tol * this%cfl_trg) then
          this%held_steps = 0
       end if
    end if
    if (this%held_steps .ge. HELD_STEPS_MAX) then
       this%held_steps = -huge(0) / 2
       write(log_buf, '(A,E15.7,A,I0,A)') 'exact_output_time has held dt at ', &
            time%dt, ' for ', HELD_STEPS_MAX, ' steps'
       call neko_error(trim(log_buf) // ' while the CFL controller asked for &
       &a change: the scheduled times are too dense to follow it. &
       &Spread them out, use tsteps for the dense schedule, or disable &
       &exact_output_time, before restarting.')
    end if

    ! Report each scheduled time once
    if (target_time .lt. huge(0.0_dp) .and. &
         abs(target_time - this%landing_target) .gt. LANDING_TOL * dt_prev) &
         then
       this%landing_target = target_time
       if (reachable) then
          write(log_buf, '(A,E15.7,1x,A,E15.7)') &
               'Landing on scheduled time:', target_time, 'dt:', time%dt
       else
          write(log_buf, '(A,E15.7)') &
               'Scheduled time out of the step bounds, passed:', target_time
       end if
       call neko_log%message(log_buf)
    end if

  end subroutine time_step_controller_land

  !> Keep the step if it divides the remaining time up to the scheduled
  !! time, dividing it anew to absorb round-off.
  !! @param remaining The time up to the scheduled time.
  !! @param dt_prev The step taken last.
  !! @param dt_new The step to take, set if the step is kept.
  !! @param kept Whether the step is kept.
  pure subroutine keep_step(remaining, dt_prev, dt_new, kept)
    real(kind=dp), intent(in) :: remaining, dt_prev
    real(kind=dp), intent(inout) :: dt_new
    logical, intent(out) :: kept
    real(kind=dp) :: x
    integer(kind=i8) :: n

    kept = .false.
    x = remaining / dt_prev
    if (x .ge. MAX_LANDING_STEPS) return
    n = nint(x, kind = i8)
    if (n .ge. 1_i8 .and. abs(x - real(n, dp)) .le. LANDING_TOL) then
       dt_new = remaining / real(n, dp)
       kept = .true.
    end if

  end subroutine keep_step

  !> Divide the remaining time up to the scheduled time into the whole
  !! number of equal steps within the bounds whose step is the closest, in
  !! ratio, to the step giving the CFL number at the centre of the band.
  !! Should that leave the CFL number outside the band, the first step may
  !! differ from the rest, each within the bounds of the one before, which
  !! moves the step in two stages where one is not enough.
  !! @param remaining The time up to the scheduled time.
  !! @param dt_prev The step taken last, whose CFL number is `cfl_previous`.
  !! @param dt_ask The step the CFL controller asked for.
  !! @param lo The shortest step allowed.
  !! @param hi The longest step allowed.
  !! @param dt_new The step to take, set if the scheduled time is reachable.
  !! @param plateau The steps after it, `dt_new` unless the step is moved
  !! in two stages.
  !! @param fitted Whether the scheduled time is reachable within the bounds.
  !! @param ratio_bound Whether the step is held back from the one asked for
  !! by the bounds on the ratio of the steps alone: the next number of
  !! steps in that direction is within the minimum and maximum step, but
  !! not within the ratio bounds.
  pure subroutine time_step_controller_fit_step(this, remaining, dt_prev, &
       dt_ask, lo, hi, dt_new, plateau, fitted, ratio_bound)
    class(time_step_controller_t), intent(in) :: this
    real(kind=dp), intent(in) :: remaining, dt_prev, dt_ask, lo, hi
    real(kind=dp), intent(inout) :: dt_new
    real(kind=dp), intent(out) :: plateau
    logical, intent(out) :: fitted, ratio_bound
    real(kind=dp) :: dt_target, dt_aim, x, dt_next, p, pmin, pmax, dist
    real(kind=dp) :: best_dist, fmin, fmax, cfl_lo, cfl_hi, cfl_p
    integer(kind=i8) :: n, n_lo, n_hi, m, m_lo, m_hi
    logical :: in_band

    fitted = .false.
    ratio_bound = .false.
    plateau = dt_new

    ! The step aimed for, and its clamp to the bounds
    dt_aim = dt_prev
    if (this%cfl_previous .gt. 0.0_dp) then
       dt_aim = dt_prev * this%cfl_centre() / this%cfl_previous
    end if
    dt_target = max(min(dt_aim, hi), lo)

    ! The numbers of steps whose step is within the bounds
    if (.not. hi .gt. 0.0_dp) return
    if (remaining / hi .ge. MAX_LANDING_STEPS) return
    n_lo = max(1_i8, ceiling(remaining / hi - LANDING_TOL, kind = i8))
    n_hi = huge(1_i8)
    if (remaining / lo .lt. MAX_LANDING_STEPS) then
       n_hi = floor(remaining / lo + LANDING_TOL, kind = i8)
    end if
    if (n_lo .gt. n_hi) return

    ! The number of steps closest, in ratio, to the step aimed for
    x = remaining / dt_target
    n = max(1_i8, floor(min(x, MAX_LANDING_STEPS), kind = i8))
    if (x * x .gt. real(n, dp) * real(n + 1_i8, dp)) n = n + 1_i8
    n = max(n_lo, min(n_hi, n))
    dt_new = remaining / real(n, dp)
    plateau = dt_new
    fitted = .true.

    ! Whether the ratio bounds alone hold the step back from the one asked
    if (dt_ask .gt. (1.0_dp + LANDING_TOL) * dt_prev .and. n .gt. 1_i8) then
       dt_next = remaining / real(n - 1_i8, dp)
       ratio_bound = dt_next .gt. hi .and. dt_next .le. this%max_dt
    else if (dt_ask .lt. (1.0_dp - LANDING_TOL) * dt_prev) then
       dt_next = remaining / real(n + 1_i8, dp)
       ratio_bound = dt_next .lt. lo .and. dt_next .ge. this%min_dt
    end if

    ! Good enough if the CFL number of that step is within the band
    cfl_lo = this%cfl_trg * (1.0_dp - this%dev_tol)
    cfl_hi = this%cfl_trg * (1.0_dp + this%dev_tol)
    if (.not. this%cfl_previous .gt. 0.0_dp) return
    cfl_p = this%cfl_previous * dt_new / dt_prev
    if (cfl_p .ge. cfl_lo .and. cfl_p .le. cfl_hi) return

    ! Otherwise a first step within the bounds, then m equal steps within
    ! the bounds of the first: the number of steps whose step is the closest
    ! to the one aimed for among those within the band, or else closer than
    ! the equal steps by a band's width
    fmin = min(this%min_dt_decrease_factor, 1.0_dp)
    fmax = max(this%max_dt_increase_factor, 1.0_dp)
    best_dist = abs(log(dt_new / dt_aim))
    x = remaining / dt_aim
    m_lo = max(1_i8, floor(min(x, MAX_LANDING_STEPS), kind = i8) - 8_i8)
    m_hi = min(ceiling(min(x, MAX_LANDING_STEPS), kind = i8) + 8_i8, &
         n_hi + 8_i8)
    do m = m_lo, m_hi
       pmin = max((remaining - hi) / real(m, dp), &
            fmin * remaining / (1.0_dp + real(m, dp) * fmin), this%min_dt)
       pmax = min((remaining - lo) / real(m, dp), &
            fmax * remaining / (1.0_dp + real(m, dp) * fmax), this%max_dt)
       if (pmin .gt. pmax) cycle
       p = max(min(dt_aim, pmax), pmin)
       dist = abs(log(p / dt_aim))
       cfl_p = this%cfl_previous * p / dt_prev
       in_band = cfl_p .ge. cfl_lo .and. cfl_p .le. cfl_hi
       if ((in_band .and. dist .lt. best_dist - LANDING_TOL) .or. &
            dist .lt. best_dist - log(1.0_dp + this%dev_tol)) then
          best_dist = dist
          dt_new = remaining - real(m, dp) * p
          plateau = p
       end if
    end do

  end subroutine time_step_controller_fit_step

end module time_step_controller
