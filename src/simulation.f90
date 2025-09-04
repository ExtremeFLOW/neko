! Copyright (c) 2020-2021, The Neko Authors
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
!> Simulation driver
module simulation
  use mpi_f08, only: MPI_Wtime
  use case, only : case_t
  use checkpoint, only : chkp_t
  use num_types, only : rp, dp
  use time_scheme_controller, only : time_scheme_controller_t
  use file, only : file_t
  use logger, only : LOG_SIZE, neko_log
  use jobctrl, only : jobctrl_time_limit
  use profiler, only : profiler_start, profiler_stop, &
       profiler_start_region, profiler_end_region
  use simcomp_executor, only : neko_simcomps
  use json_utils, only : json_get, json_get_or_default
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  implicit none
  private

  interface simulation_restart
     module procedure case_restart_from_parameters, &
          case_restart_from_checkpoint
  end interface simulation_restart

  public :: simulation_init, simulation_step, simulation_finalize, &
       simulation_restart

contains

  !> Initialise a simulation of a case
  subroutine simulation_init(C, dt_controller)
    type(case_t), intent(inout) :: C
    type(time_step_controller_t), intent(inout) :: dt_controller
    character(len=LOG_SIZE) :: log_buf
    logical :: found
    character(len=:), allocatable :: restart_file

    ! Restart the case if needed
    call C%params%get('case.restart_file', restart_file, found)
    if (found .and. len_trim(restart_file) .gt. 0) then
       ! Restart the case
       call simulation_restart(C)
    end if

    ! Write the initial logging message
    call neko_log%section('Starting simulation')
    write(log_buf, '(A, E15.7,A,E15.7,A)') &
         'T  : [', C%time%t, ',', C%time%end_time, ']'
    call neko_log%message(log_buf)
    if (.not. dt_controller%is_variable_dt) then
       write(log_buf, '(A, E15.7)') 'dt :  ', C%time%dt
    else
       write(log_buf, '(A, E15.7)') 'CFL :  ', dt_controller%cfl_trg
    end if
    call neko_log%message(log_buf)

    ! Execute outputs and user-init before time loop
    call neko_log%section('Preprocessing')
    call C%output_controller%execute(C%time)

    call C%user%initialize(C%time)
    call neko_log%end_section()
    call neko_log%newline()

  end subroutine simulation_init

  !> Finalize a simulation of a case
  subroutine simulation_finalize(C)
    type(case_t), intent(inout) :: C
    logical :: output_at_end

    ! Run a final output if specified in the json
    call json_get_or_default(C%params, 'case.output_at_end', &
         output_at_end, .true.)
    call C%output_controller%execute(C%time, output_at_end)

    if (.not. (output_at_end) .and. C%time%t .lt. C%time%end_time) then
       call simulation_joblimit_chkp(C, C%time%t)
    end if

    ! Finalize the user modules
    call C%user%finalize(C%time)

    call neko_log%end_section('Normal end.')

  end subroutine simulation_finalize

  !> Compute a single time-step of a case
  subroutine simulation_step(C, dt_controller, tstep_loop_start_time)
    type(case_t), intent(inout) :: C
    type(time_step_controller_t), intent(inout) :: dt_controller
    real(kind=dp), optional, intent(in) :: tstep_loop_start_time
    real(kind=dp) :: start_time, end_time, tstep_start_time
    real(kind=rp) :: cfl
    character(len=LOG_SIZE) :: log_buf

    ! Setup the time step, and start time
    call profiler_start_region('Time-Step')
    start_time = MPI_WTIME()
    tstep_start_time = start_time

    ! Compute the next time step size
    cfl = C%fluid%compute_cfl(C%time%dt)
    call dt_controller%set_dt(C%time, cfl)
    if (dt_controller%is_variable_dt) cfl = C%fluid%compute_cfl(C%time%dt)

    ! Advance time step from t to t+dt and print the status
    call simulation_settime(C%time, C%fluid%ext_bdf)
    call C%time%status()
    call neko_log%begin()

    call neko_log%section('Preprocessing')

    write(log_buf, "(A,F8.4,2x,A,1P,E14.7)") 'CFL:', cfl, 'dt:', C%time%dt
    call neko_log%message(log_buf)

    ! Run the preprocessing
    call C%user%preprocess(C%time)
    call neko_simcomps%preprocess(C%time)
    call neko_log%end_section()

    ! Fluid step
    call neko_log%section('Fluid')
    start_time = MPI_WTIME()
    call C%fluid%step(C%time, dt_controller)
    end_time = MPI_WTIME()
    write(log_buf, '(A,3X,E15.7)') 'Fluid step time (s):', end_time - start_time
    call neko_log%end_section(log_buf)

    ! Scalar step
    if (allocated(C%scalars)) then
       start_time = MPI_WTIME()
       call neko_log%section('Scalar')
       call C%scalars%step(C%time, C%fluid%ext_bdf, dt_controller)
       end_time = MPI_WTIME()
       write(log_buf, '(A,6X,E15.7)') 'Scalar step time:', end_time - start_time
       call neko_log%end_section(log_buf)
    end if

    ! Postprocessing
    call neko_log%section('Postprocessing')

    ! Execute all simulation components
    call neko_simcomps%compute(C%time)

    ! Run the user compute routine
    call C%user%compute(C%time)

    ! Run any IO needed.
    call C%output_controller%execute(C%time)

    call neko_log%end_section()

    ! End the step and print summary
    end_time = MPI_WTIME()
    call neko_log%section('Step summary')
    write(log_buf, '(A20,I8,A6,E15.7)') &
         'Total time for step ', C%time%tstep, ' (s): ', &
         end_time-tstep_start_time
    call neko_log%message(log_buf)

    if (present(tstep_loop_start_time)) then
       write(log_buf, '(A34,E15.7)') &
            'Total elapsed time (s):           ', &
            end_time - tstep_loop_start_time
       call neko_log%message(log_buf)
    end if

    call neko_log%end_section()
    call neko_log%end()
    call profiler_end_region

  end subroutine simulation_step

  subroutine simulation_settime(time, ext_bdf)
    type(time_state_t), intent(inout) :: time
    type(time_scheme_controller_t), intent(inout), allocatable :: ext_bdf
    integer :: i

    if (allocated(ext_bdf)) then
       do i = 10, 2, -1
          time%tlag(i) = time%tlag(i-1)
          time%dtlag(i) = time%dtlag(i-1)
       end do

       time%dtlag(1) = time%dt
       time%tlag(1) = time%t
       if (ext_bdf%ndiff .eq. 0) then
          time%dtlag(2) = time%dt
          time%tlag(2) = time%t
       end if

       call ext_bdf%set_coeffs(time%dtlag)
    end if

    time%tstep = time%tstep + 1
    time%t = time%t + time%dt

  end subroutine simulation_settime

  !> Restart a case @a C from a given checkpoint
  subroutine case_restart_from_parameters(C)
    type(case_t), intent(inout) :: C
    type(file_t) :: chkpf, previous_meshf
    character(len=LOG_SIZE) :: log_buf
    character(len=:), allocatable :: restart_file
    character(len=:), allocatable :: restart_mesh_file
    real(kind=rp) :: tol
    logical :: found, check_cont
    integer :: i

    call json_get(C%params, 'case.restart_file', restart_file)
    call json_get_or_default(C%params, 'case.restart_mesh_file', &
         restart_mesh_file, "")

    if (restart_mesh_file .ne. "") then
       call previous_meshf%init(trim(restart_mesh_file))
       call previous_meshf%read(C%chkp%previous_mesh)

       call json_get_or_default(C%params, 'case.mesh2mesh_tolerance', &
            C%chkp%mesh2mesh_tol, 1e-6_rp)
    end if

    call neko_log%section('Restarting from checkpoint')
    write(log_buf, '(A,A)') 'File :   ', trim(restart_file)
    call neko_log%message(log_buf)

    call chkpf%init(trim(restart_file))
    call chkpf%read(C%chkp)

    call case_restart_from_checkpoint(C, C%chkp)

    ! Free the previous mesh
    call C%chkp%previous_mesh%free()

    write(log_buf, '(A,E15.7)') 'Time : ', C%time%t
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine case_restart_from_parameters

  !> Restart a case @a C from a given checkpoint
  subroutine case_restart_from_checkpoint(C, chkp)
    type(case_t), intent(inout) :: C
    type(chkp_t), intent(inout) :: chkp
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    ! Restart the time state and BDF coefficients
    call C%time%restart(chkp)
    do i = 1, size(C%time%dtlag)
       call C%fluid%ext_bdf%set_coeffs(C%time%dtlag)
    end do

    ! Restart the fluid and scalars
    call C%fluid%restart(chkp)
    if (allocated(C%scalars)) call C%scalars%restart(chkp)

    ! Restart the output controller
    call C%output_controller%set_counter(C%time)

    ! Restart the simulation components
    call neko_simcomps%restart(C%time)

  end subroutine case_restart_from_checkpoint

  !> Write a checkpoint at joblimit
  subroutine simulation_joblimit_chkp(C, t)
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    type(file_t) :: chkpf
    character(len=:), allocatable :: chkp_format
    character(len=LOG_SIZE) :: log_buf
    character(len=10) :: format_str
    logical :: found

    call json_get_or_default(C%params, 'case.checkpoint_format', chkp_format, &
         'default')
    call C%chkp%sync_host()
    if (chkp_format .eq. 'hdf5') then
       format_str = '.h5'
    else
       format_str = '.chkp'
    end if
    call chkpf%init(C%output_directory // 'joblimit'//trim(format_str))
    call chkpf%write(C%chkp, t)
    write(log_buf, '(A)') '! saving checkpoint >>>'
    call neko_log%message(log_buf)

  end subroutine simulation_joblimit_chkp

end module simulation
