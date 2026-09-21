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
  use num_types, only : dp
  use logger, only : neko_log, LOG_SIZE
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get_or_lookup_or_default
  use time_state, only : time_state_t
  use time_based_controller, only : TIME_TOL
  use comm, only : pe_size, global_pe_size, NEKO_GLOBAL_COMM
  use mpi_f08, only : MPI_MIN, MPI_MAX, MPI_IN_PLACE, MPI_Allreduce, &
       MPI_DOUBLE_PRECISION, MPI_INTEGER
  implicit none
  private

  !> Relative tolerance on the time step when landing on a scheduled time.
  !! The remaining time up to a scheduled time is split into steps that may
  !! exceed the step asked for by this fraction, so that round-off in the
  !! accumulation of the time does not split a remaining time that is a whole
  !! number of steps already into one step more. It is also the largest
  !! change of the step that is not registered as a change.
  real(kind=dp), public, parameter :: LANDING_TOL = 1.0e-6_dp

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
     !! it changed. -1 until the first time step, or as long as the time step
     !! cannot change, see `dt_may_change`.
     integer :: dt_last_change = -1
     real(kind=dp) :: alpha = 0.0_dp !< coefficient of running average
     real(kind=dp) :: max_dt_increase_factor = 0.0_dp
     real(kind=dp) :: min_dt_decrease_factor = 0.0_dp
     real(kind=dp) :: dev_tol = 0.0_dp
     !> Whether to shorten the time step so that the scheduled sampling and
     !! output times are reached exactly.
     logical :: exact_output_time = .false.
     !> Number of time steps ahead of a scheduled time over which the step is
     !! shortened to land on it.
     integer :: landing_steps = 10
     !> The time step asked for by the case or set by the CFL controller,
     !! before it is shortened to land on a scheduled time.
     real(kind=dp) :: dt_nominal = 0.0_dp
     !> Whether the time step about to be taken is shortened to land on a
     !! scheduled time.
     logical :: dt_is_landing = .false.
     !> The time step taken last, to tell whether the one about to be taken
     !! differs from it.
     real(kind=dp) :: dt_previous = 0.0_dp
     !> The scheduled time being landed on, only used for the log.
     real(kind=dp) :: landing_target = huge(0.0_dp)
   contains
     !> Initialize object.
     procedure, pass(this) :: init => time_step_controller_init
     !> Set time stepping
     procedure, pass(this) :: set_dt => time_step_controller_set_dt
     !> The first time step of a variable time step run.
     procedure, pass(this) :: first_dt => time_step_controller_first_dt
     !> Shorten the time step to land on the next scheduled time.
     procedure, pass(this) :: land => time_step_controller_land
     !> Whether the time step can change during the run.
     procedure, pass(this) :: dt_may_change => &
          time_step_controller_dt_may_change

  end type time_step_controller_t

contains

  !> Constructor
  !! @param order order of the interpolation
  subroutine time_step_controller_init(this, params)
    class(time_step_controller_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    integer :: exact_any, exact_all, ierr

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

    ! Landing on the scheduled times works for a fixed and a variable step.
    call json_get_or_default(params, 'exact_output_time', &
         this%exact_output_time, .false.)
    this%landing_steps = 10
    if (this%exact_output_time) then
       call json_get_or_lookup_or_default(params, 'output_landing_steps', &
            this%landing_steps, 10)
       if (this%landing_steps .lt. 1) then
          call neko_log%error('output_landing_steps must be at least 1')
       end if
    end if

    ! In an MPMD run the landing takes a collective at every step, so the
    ! coupled simulations have to agree on it, or the ones landing would
    ! wait forever for the ones that do not.
    if (pe_size .ne. global_pe_size) then
       exact_any = merge(1, 0, this%exact_output_time)
       exact_all = exact_any
       call MPI_Allreduce(MPI_IN_PLACE, exact_any, 1, MPI_INTEGER, &
            MPI_MAX, NEKO_GLOBAL_COMM, ierr)
       call MPI_Allreduce(MPI_IN_PLACE, exact_all, 1, MPI_INTEGER, &
            MPI_MIN, NEKO_GLOBAL_COMM, ierr)
       if (exact_any .ne. exact_all) then
          call neko_log%error('exact_output_time must be set in all the &
               &coupled cases of an MPMD run, or in none')
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

    dt = min(this%cfl_trg / cfl * time%dt, this%init_dt)
    dt = max(min(dt, this%max_dt), this%min_dt)

  end function time_step_controller_first_dt

  !> Whether the time step can change during the run, as it does with a
  !! variable time step and when landing on the scheduled times. If so,
  !! `dt_last_change` counts the steps since it last did.
  pure function time_step_controller_dt_may_change(this) result(may_change)
    class(time_step_controller_t), intent(in) :: this
    logical :: may_change

    may_change = this%is_variable_dt .or. this%exact_output_time

  end function time_step_controller_dt_may_change

  !> Set new dt based on cfl if requested
  !! @param time The time state, whose `dt` is the step taken last on entry
  !! and the step about to be taken on exit.
  !! @param cfl courant number of current iteration.
  !! @Algorithm:
  !! 1. Set the first time step such that cfl is the set one;
  !! 2. During time-stepping, adjust dt when cfl_avg is offset by 20%.
  !! A step shortened to land on a scheduled time (see `land`) is undone
  !! first, so that the controller always works on the step it asked for, and
  !! the shortening never carries over to the following steps.
  subroutine time_step_controller_set_dt(this, time, cfl)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=dp), intent(in) :: cfl
    real(kind=dp) :: dt_old, scaling_factor, global_min_dt
    real(kind=dp) :: cfl_nominal
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr

    this%dt_previous = time%dt
    cfl_nominal = cfl

    ! Undo the shortening applied to land on a scheduled time. The CFL number
    ! was computed with the shortened step and is proportional to it, so it is
    ! rescaled to the step the controller asked for.
    if (this%dt_is_landing) then
       if (abs(time%dt) .gt. 0.0_dp) then
          cfl_nominal = cfl * this%dt_nominal / time%dt
       end if
       time%dt = this%dt_nominal
       this%dt_is_landing = .false.
    end if

    if (this%is_variable_dt) then

       ! Reset the average cfl if it is the first time step since the last
       ! change
       if (this%dt_last_change .eq. 0) then
          this%cfl_avg = cfl_nominal
       end if

       if (this%dt_last_change .eq. -1) then

          ! Set the first dt for desired cfl, or use the provided initial dt if
          ! it is smaller. Then clamp between max and min dt if provided.
          time%dt = this%first_dt(time, cfl_nominal)
          this%dt_last_change = 0
          this%cfl_avg = cfl_nominal

       else
          ! Calculate the average of cfl over the desired interval
          this%cfl_avg = this%alpha * cfl_nominal &
               + (1 - this%alpha) * this%cfl_avg

          if (abs(this%cfl_avg - this%cfl_trg) .ge. this%dev_tol*this%cfl_trg &
               .and. this%dt_last_change .ge. this%max_update_frequency &
               .or. this%dt_last_change .ge. this%min_update_frequency) then

             if (this%cfl_trg/cfl_nominal .ge. 1) then
                ! increase of time step
                scaling_factor = min(this%max_dt_increase_factor, &
                     this%cfl_trg/cfl_nominal)
             else
                ! reduction of time step
                scaling_factor = max(this%min_dt_decrease_factor, &
                     this%cfl_trg/cfl_nominal)
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
       if (pe_size .ne. global_pe_size) then
          global_min_dt = time%dt
          call MPI_Allreduce(MPI_IN_PLACE, global_min_dt, 1, &
               MPI_DOUBLE_PRECISION, MPI_MIN, NEKO_GLOBAL_COMM, ierr)

          ! If my dt is larger that the global min, mark a change
          if (time%dt .gt. global_min_dt) then
             time%dt = global_min_dt
             this%dt_last_change = 0
          end if
       end if

    else if (this%exact_output_time) then
       ! The step is fixed, but landing on the scheduled times can shorten
       ! it, so keep counting the steps since it last changed.
       this%dt_last_change = this%dt_last_change + 1
    end if

    ! The step the controller settled on, before landing on a scheduled time.
    this%dt_nominal = time%dt

  end subroutine time_step_controller_set_dt

  !> Shorten the time step so that the next scheduled sampling or output
  !! time, and the end of the simulation, are reached exactly.
  !! @param time The time state, whose `dt` is the step about to be taken and
  !! is shortened in place.
  !! @param time_to_next The time until the next scheduled time, as
  !! `time_to_next` of the controllers gives it, `huge(0.0_dp)` if there is
  !! none. The end of the simulation is a target of its own and needs not
  !! be included.
  !! @details To be called after `set_dt`, whose step is the one shortened.
  !! Once the target is within `landing_steps` steps, the remaining time is
  !! split into the smallest whole number of equal steps that do not exceed
  !! the step asked for (by more than `LANDING_TOL`), so the step only ever
  !! gets shorter, never longer, and the steps up to the target are all
  !! equal. The step asked for is restored by the next call to `set_dt`.
  !!
  !! No step is shorter than `TIME_TOL` times the step asked for, or than
  !! `min_timestep`. A scheduled time closer than that to another one, which
  !! takes two controllers with unrelated intervals, or to `end_time`, is
  !! executed within that shortest step of its time instead of landed on.
  !! The last step of the simulation is exempt from the floor, so that a run
  !! started closer than that to `end_time` still ends exactly there. In an
  !! MPMD run the shortened step is the smallest one over the simulations, so
  !! that they keep advancing in lockstep.
  subroutine time_step_controller_land(this, time, time_to_next)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=dp), intent(in) :: time_to_next
    real(kind=dp) :: dt, dt_new, dt_min, direction, remaining, global_min_dt
    real(kind=dp) :: remaining_end, window
    character(len=LOG_SIZE) :: log_buf
    integer :: nsteps, ierr
    logical :: landing_on_end, adjusted

    if (.not. this%exact_output_time) return

    dt = abs(time%dt)
    direction = sign(1.0_dp, time%dt)
    dt_new = dt
    adjusted = .false.
    remaining = huge(0.0_dp)

    if (dt .gt. 0.0_dp) then
       ! Also the end of the simulation is landed on. Nothing follows the
       ! last step, so it may be as short as it needs to be, and the end is
       ! the target whenever it is within the shortest step allowed of the
       ! next scheduled time: that time is then executed at the end, rather
       ! than landed on and followed by a step too short to be worth taking.
       dt_min = min(dt, max(TIME_TOL * dt, this%min_dt))
       remaining_end = direction * (time%end_time - time%t)
       landing_on_end = remaining_end .le. time_to_next + dt_min
       if (landing_on_end) then
          remaining = remaining_end
       else
          remaining = time_to_next
       end if

       ! The landing is engaged at least (1 + TIME_TOL) steps ahead of the
       ! target: a full step from closer than that would leave less than the
       ! shortest step allowed, and the target would be passed by a fraction
       ! of a step instead of landed on. With a variable step the step can
       ! grow before the next one is taken, which the window accounts for.
       window = 1.0_dp + TIME_TOL
       if (this%is_variable_dt) then
          window = 1.0_dp + TIME_TOL * max(1.0_dp, this%max_dt_increase_factor)
       end if
       window = max(real(this%landing_steps, dp), window) * &
            (1.0_dp + LANDING_TOL)

       if (remaining .gt. 0.0_dp .and. remaining .le. window * dt) then
          nsteps = max(1, ceiling(remaining / dt - LANDING_TOL))
          dt_new = remaining / real(nsteps, dp)

          ! Leave the step alone while the change is round-off, so that a
          ! step that divides the remaining time already is not registered
          ! as changed at every step. The last step is always taken, as it
          ! is the one that has to be exact.
          if (nsteps .gt. 1 .and. abs(dt_new - dt) .le. LANDING_TOL * dt) then
             dt_new = dt
          else
             adjusted = .true.
          end if

          ! A scheduled time closer than the shortest step allowed is passed,
          ! except by the last step of the run.
          if (.not. (landing_on_end .and. nsteps .eq. 1)) then
             dt_new = max(dt_new, dt_min)
          end if
          adjusted = adjusted .and. abs(dt_new - dt) .gt. 0.0_dp
       end if
    end if

    ! If running in mpmd, the shortened step is the minimum across simulations
    if (pe_size .ne. global_pe_size) then
       global_min_dt = dt_new
       call MPI_Allreduce(MPI_IN_PLACE, global_min_dt, 1, &
            MPI_DOUBLE_PRECISION, MPI_MIN, NEKO_GLOBAL_COMM, ierr)
       if (global_min_dt .lt. dt_new) then
          dt_new = global_min_dt
          adjusted = .true.
       end if
    end if

    if (adjusted) then
       this%dt_is_landing = .true.
       time%dt = direction * dt_new
       ! Report each scheduled time landed on once
       if (remaining .lt. huge(0.0_dp) .and. &
            abs(time%t + direction * remaining - this%landing_target) .gt. &
            LANDING_TOL * dt) then
          this%landing_target = time%t + direction * remaining
          write(log_buf, '(A,E15.7,1x,A,E15.7)') &
               'Landing on output time:', this%landing_target, &
               'New dt:', time%dt
          call neko_log%message(log_buf)
       end if
    end if

    ! The projection spaces watch for a change of the step, whatever its
    ! reason. The equal steps up to a target differ by round-off only.
    if (abs(time%dt - this%dt_previous) .gt. LANDING_TOL * dt) then
       this%dt_last_change = 0
    end if

  end subroutine time_step_controller_land

end module time_step_controller
