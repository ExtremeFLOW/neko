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
  use time_based_controller, only : time_based_controller_next_scheduled_time
  use comm, only : pe_size, global_pe_size, NEKO_GLOBAL_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_MIN, MPI_IN_PLACE, MPI_Allreduce, MPI_DOUBLE_PRECISION
  implicit none
  private

  !> Relative tolerance used when splitting the remaining time up to an output
  !! into an integer number of steps, so that an interval that is already a
  !! whole multiple of dt is not split one step too finely.
  real(kind=dp), parameter :: LANDING_TOL = 1.0e-9_dp

  !> Provides a tool to set time step dt
  type, public :: time_step_controller_t
     logical :: is_variable_dt
     real(kind=dp) :: cfl_trg = 0.0_dp
     real(kind=dp) :: cfl_avg = 0.0_dp
     real(kind=dp) :: init_dt = huge(0.0_dp)
     real(kind=dp) :: max_dt = 0.0_dp
     real(kind=dp) :: min_dt = 0.0_dp
     integer :: max_update_frequency = 0
     integer :: min_update_frequency = 0
     integer :: dt_last_change = 0
     real(kind=dp) :: alpha = 0.0_dp !< coefficient of running average
     real(kind=dp) :: max_dt_increase_factor = 0.0_dp
     real(kind=dp) :: min_dt_decrease_factor = 0.0_dp
     real(kind=dp) :: dev_tol = 0.0_dp
     !> Whether to shrink dt so that sampling and output times are hit exactly.
     logical :: exact_output_time = .false.
     !> Number of steps ahead of an output time at which dt starts to be
     !! shrunk to land on it.
     integer :: landing_steps = 10
     !> The dt asked for by the user or the CFL controller, i.e. before any
     !! reduction applied to land on an output time.
     real(kind=dp) :: dt_nominal = 0.0_dp
     !> Whether dt is currently reduced to land on an output time.
     logical :: dt_is_landing = .false.
     !> The output time currently being landed on, only used for logging.
     real(kind=dp) :: landing_target = huge(0.0_dp)
   contains
     !> Initialize object.
     procedure, pass(this) :: init => time_step_controller_init
     !> Set time stepping
     procedure, pass(this) :: set_dt => time_step_controller_set_dt
     !> Reduce dt to land exactly on the next sampling or output time.
     procedure, pass(this) :: land_on_output => &
          time_step_controller_land_on_output

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

    ! Landing on the output times works both for a fixed and a variable dt.
    call json_get_or_default(params, 'exact_output_time', &
         this%exact_output_time, .false.)
    if (this%exact_output_time) then
       call json_get_or_lookup_or_default(params, 'output_landing_steps', &
            this%landing_steps, 10)
       if (this%landing_steps .lt. 1) then
          call neko_log%error('output_landing_steps must be at least 1')
       end if
    end if

  end subroutine time_step_controller_init

  !> Set new dt based on cfl if requested
  !! @param dt time step in case_t.
  !! @param cfl courant number of current iteration.
  !! @param tstep the current time step.
  !! @Algorithm:
  !! 1. Set the first time step such that cfl is the set one;
  !! 2. During time-stepping, adjust dt when cfl_avg is offset by 20%.
  !! 3. If requested, shrink dt so that the next sampling or output time is
  !!    reached exactly.
  subroutine time_step_controller_set_dt(this, time, cfl)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=dp), intent(in) :: cfl
    real(kind=dp) :: dt_old, scaling_factor, global_min_dt
    real(kind=dp) :: dt_entry, cfl_nominal
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr

    dt_entry = time%dt
    cfl_nominal = cfl

    ! Undo the reduction that was applied to land on an output time, so that
    ! the CFL controller always works with the dt it asked for. The cfl is
    ! proportional to dt, so it can simply be rescaled to the unreduced step.
    if (this%dt_is_landing) then
       if (time%dt .gt. 0.0_dp) cfl_nominal = cfl * this%dt_nominal / time%dt
       time%dt = this%dt_nominal
       this%dt_is_landing = .false.
    end if

    ! Check if variable dt is requested
    if (this%is_variable_dt) then

       ! Reset the average cfl if it is the first time step since the last
       ! change
       if (this%dt_last_change .eq. 0) then
          this%cfl_avg = cfl_nominal
       end if

       if (this%dt_last_change .eq. -1) then

          ! Set the first dt for desired cfl, or use the provided initial dt if
          ! it is smaller. Then clamp between max and min dt if provided.
          time%dt = min(this%cfl_trg / cfl_nominal * time%dt, this%init_dt)
          time%dt = max(min(time%dt, this%max_dt), this%min_dt)
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
    end if

    ! The dt the controller settled on, before landing on an output time.
    this%dt_nominal = time%dt

    if (this%exact_output_time) call this%land_on_output(time)

    ! Whatever the reason, a change of dt invalidates e.g. the projection space
    if (time%dt .ne. dt_entry) this%dt_last_change = 0

  end subroutine time_step_controller_set_dt

  !> Reduce dt such that the next sampling or output time is reached exactly.
  !! @details Once the target is within `landing_steps` time steps, the
  !! remaining time is split into the smallest integer number of equal steps
  !! that are no larger than the dt asked for. The step size therefore only
  !! ever decreases, which keeps the scheme stable, and the time value seen by
  !! the outputs is the requested one rather than the first one past it.
  !! The reduction is undone by `set_dt` on the following step.
  !! @param time Current time, whose `dt` is modified in place.
  subroutine time_step_controller_land_on_output(this, time)
    class(time_step_controller_t), intent(inout) :: this
    type(time_state_t), intent(inout) :: time
    real(kind=dp) :: t_target, remaining, dt_new
    character(len=LOG_SIZE) :: log_buf
    integer :: nsteps

    if (time%dt .le. 0.0_dp) return

    ! The first of the scheduled sampling and output times, but never past the
    ! end of the simulation, which is a target in its own right.
    t_target = time_based_controller_next_scheduled_time(time)
    if (time%end_time .gt. time%t) t_target = min(t_target, time%end_time)
    if (t_target .ge. huge(0.0_dp)) return

    remaining = t_target - time%t
    if (remaining .le. 0.0_dp) return

    ! Leave dt alone until the target is within reach.
    if (remaining .gt. this%landing_steps * time%dt) return

    nsteps = max(ceiling(remaining / time%dt - LANDING_TOL), 1)
    dt_new = remaining / nsteps

    ! Never go below the floor the user asked for.
    if (dt_new .lt. this%min_dt) return

    ! Leave dt alone while the adjustment is pure round-off, so that a dt that
    ! already divides the interval is not flagged as changed on every step.
    ! The last step is always taken, since it is the one that has to be exact.
    if (nsteps .gt. 1 .and. &
         abs(dt_new - time%dt) .le. LANDING_TOL * time%dt) return

    if (t_target .ne. this%landing_target) then
       write(log_buf, '(A,E15.7,1x,A,E15.7)') &
            'Landing on output time:', t_target, 'New dt:', dt_new
       call neko_log%message(log_buf)
       this%landing_target = t_target
    end if

    time%dt = dt_new
    this%dt_is_landing = .true.

  end subroutine time_step_controller_land_on_output




end module time_step_controller
