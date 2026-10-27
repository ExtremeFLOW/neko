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
  use utils, only : neko_error
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get_or_lookup_or_default
  use time_state, only : time_state_t
  use math, only : abscmp
  use comm, only : NEKO_GLOBAL_COMM, is_mpmd
  use mpi_f08, only : MPI_MIN, MPI_MAX, MPI_IN_PLACE, MPI_Allreduce, &
       MPI_DOUBLE_PRECISION, MPI_INTEGER
  implicit none
  private

  !> Relative tolerance on the number of steps up to a scheduled time, so
  !! that round-off in the accumulated time does not add a step, and on the
  !! change of the step that counts as one.
  real(kind=dp), public, parameter :: LANDING_TOL = 1.0e-6_dp

  !> Largest number of steps up to a scheduled time that is still fitted.
  real(kind=dp), parameter :: MAX_LANDING_STEPS = 1.0e15_dp

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
     !> Whether to fit the time step so that the scheduled sampling and
     !! output times, and the end time, are reached exactly.
     logical :: exact_output_time = .false.
     !> The step taken last.
     real(kind=dp) :: dt_previous = 0.0_dp
   contains
     !> Initialize object.
     procedure, pass(this) :: init => time_step_controller_init
     !> Set time stepping
     procedure, pass(this) :: set_dt => time_step_controller_set_dt
     !> The first time step of a variable time step run.
     procedure, pass(this) :: first_dt => time_step_controller_first_dt
     !> Fit the time step to land on the next scheduled time.
     procedure, pass(this) :: land => time_step_controller_land

  end type time_step_controller_t

contains

  !> Constructor
  !! @param params The time settings of the case file, `case.time`.
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

    call json_get_or_default(params, 'exact_output_time', &
         this%exact_output_time, .false.)
    if (this%exact_output_time .and. .not. this%is_variable_dt) then
       call neko_error('exact_output_time requires variable_timestep')
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

  !> Set new dt based on cfl if requested
  !! @param time The time state, whose `dt` is the step taken last on entry
  !! and the step about to be taken on exit.
  !! @param cfl courant number of current iteration.
  !! @Algorithm:
  !! 1. Set the first time step such that cfl is the set one;
  !! 2. During time-stepping, adjust dt when cfl_avg is offset by 20%.
  !! 3. An exact zero cfl, e.g. from a zero velocity field, puts no limit on
  !! dt, so the first time step is the initial one and a later one is
  !! increased by the maximum increase factor.
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

    ! Check if variable dt is requested
    if (.not. this%is_variable_dt) return

    ! Reset the average cfl if it is the first time step since the last change
    if (this%dt_last_change .eq. 0) then
       this%cfl_avg = cfl
    end if

    if (this%dt_last_change .eq. -1) then

       ! Set the first dt for desired cfl, or use the provided initial dt if it
       ! is smaller. Then clamp between max and min dt if provided.
       time%dt = this%first_dt(time, cfl)
       this%dt_last_change = 0
       this%cfl_avg = cfl
       ! There is no step before the first one
       this%dt_previous = time%dt

    else
       ! Calculate the average of cfl over the desired interval
       this%cfl_avg = this%alpha * cfl + (1 - this%alpha) * this%cfl_avg

       if (abs(this%cfl_avg - this%cfl_trg) .ge. this%dev_tol * this%cfl_trg &
            .and. this%dt_last_change .ge. this%max_update_frequency &
            .or. this%dt_last_change .ge. this%min_update_frequency) then

          ! Test for an exact zero cfl first, so that it is never divided by.
          if (abscmp(cfl, 0.0_dp)) then
             ! no limit on dt, so increase it as much as allowed
             scaling_factor = this%max_dt_increase_factor
          else if (this%cfl_trg/cfl .ge. 1) then
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


  end subroutine time_step_controller_set_dt

  !> Fit the time step so that the next scheduled time, or the end time,
  !! is reached exactly.
  !! @param time The time state, whose `dt` is the step proposed by the CFL
  !! controller, adjusted in place.
  !! @param time_to_next The time until the next scheduled time, as
  !! `time_to_next` of the controllers gives it, `huge(0.0_dp)` if none.
  !! @details To be called right after `set_dt`. The time up to the next
  !! target is divided into the number of equal steps nearest, in ratio, to
  !! the proposed step, within `min_dt_decrease_factor` and
  !! `max_dt_increase_factor` of the previous step and below `max_timestep`.
  !! As long as the CFL controller keeps the step and the step divides the
  !! remaining time, which it does between scheduled times and at those of
  !! equal intervals, the step is returned unchanged. The end is the target
  !! if the next scheduled time is within `min_dt_decrease_factor` times the
  !! previous step of it, and that scheduled time is executed at the end. A
  !! scheduled time closer than that is passed, and executed at the first
  !! step past it.
  !!
  !! Designed for output intervals of several steps, the same for all
  !! controllers or whole multiples of each other. Schedules that produce
  !! scheduled times close to each other make the step shorter than needed.
  subroutine time_step_controller_land(this, time, time_to_next)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=dp), intent(in) :: time_to_next
    real(kind=dp) :: direction, dt_prev, lo, remaining, remaining_end, x
    integer(kind=i8) :: n
    integer :: ierr

    if (.not. (this%exact_output_time .and. this%is_variable_dt)) return

    direction = sign(1.0_dp, time%dt)
    dt_prev = abs(this%dt_previous)
    lo = this%min_dt_decrease_factor * dt_prev

    ! The end, or the next scheduled time unless it is too close
    remaining_end = direction * (time%end_time - time%t)
    if (remaining_end .le. time_to_next + lo) then
       remaining = remaining_end
    else if (time_to_next .ge. lo) then
       remaining = time_to_next
    else
       remaining = 0.0_dp
    end if

    x = remaining / abs(time%dt)
    if (remaining .gt. 0.0_dp .and. x .lt. MAX_LANDING_STEPS) then
       ! Candidate count: x proposed steps fit into the remaining time,
       ! round down, but take at least one step
       n = max(1_i8, floor(x, kind = i8))
       ! Round up instead if x is closer to n + 1 in ratio, that is, if
       ! (n + 1) / x < x / n. The boundary is the geometric mean
       ! sqrt(n (n + 1)), not n + 0.5, so that a step twice too long and one
       ! half as long count as equally far off
       if (x * x .gt. real(n, dp) * real(n + 1_i8, dp)) n = n + 1_i8
       ! Lower bound on the step, min_dt_decrease_factor times the previous
       ! one: at most remaining / lo steps, and at least one. The tolerance
       ! keeps round-off from dropping a count that fits exactly
       n = min(n, max(1_i8, floor(remaining / lo + LANDING_TOL, kind = i8)))
       ! Upper bound on the step, max_dt_increase_factor times the previous
       ! one or max_timestep, whichever is smaller: at least that many
       ! steps. Applied last, so it wins if the two bounds conflict
       n = max(n, ceiling(remaining / min(this%max_dt_increase_factor * &
            dt_prev, this%max_dt) - LANDING_TOL, kind = i8))
       ! The fitted step: n equal steps end exactly on the target
       time%dt = direction * remaining / real(n, dp)
    end if

    ! In an MPMD run the simulations advance in lockstep
    if (is_mpmd()) then
       call MPI_Allreduce(MPI_IN_PLACE, time%dt, 1, MPI_DOUBLE_PRECISION, &
            MPI_MIN, NEKO_GLOBAL_COMM, ierr)
    end if

    ! The projection spaces watch for a change of the step taken
    if (abs(abs(time%dt) - dt_prev) .gt. LANDING_TOL * dt_prev) then
       this%dt_last_change = 0
    end if

  end subroutine time_step_controller_land

end module time_step_controller
