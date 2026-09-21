!> @file simulation_adjoint.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Adjoint simulation driver
module simulation_adjoint
  use mpi_f08, only: MPI_WTIME
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp, dp
  use time_scheme_controller, only: time_scheme_controller_t
  use file, only: file_t
  use logger, only: LOG_SIZE, neko_log
  use profiler, only: profiler_start_region, profiler_end_region
  use json_utils, only: json_get_or_default
  use time_state, only : time_state_t
  use time_step_controller, only: time_step_controller_t
  use adjoint_case, only: adjoint_case_t
  use scratch_registry, only: neko_scratch_registry
  use device_math, only: device_glsc3
  use math, only: glsc3
  use vector, only: vector_t
  implicit none
  private

  public :: simulation_adjoint_init, simulation_adjoint_step, &
       simulation_adjoint_finalize, simulation_adjoint_restart

contains


  !> Initialise a simulation_adjoint of a case
  subroutine simulation_adjoint_init(C, dt_controller)
    type(adjoint_case_t), intent(inout) :: C
    type(time_step_controller_t), intent(inout) :: dt_controller
    character(len=LOG_SIZE) :: log_buf

    ! Write the initial logging message
    call neko_log%section('Adjoint Starting simulation')
    write(log_buf, '(A, E15.7,A,E15.7,A)') &
         'T  : [', C%time%t, ',', C%time%end_time, ']'
    call neko_log%message(log_buf)
    if (.not. dt_controller%is_variable_dt) then
       write(log_buf, '(A, E15.7)') 'dt :  ', C%time%dt
    else
       write(log_buf, '(A, E15.7)') 'CFL :  ', dt_controller%cfl_trg
    end if
    call neko_log%message(log_buf)

    ! Call stats, samplers and user-init before time loop
    call neko_log%section('Postprocessing')
    call C%output_controller%execute(C%time)
    call simulation_adjoint_norm_output(C, C%time)

    call C%case%user%initialize(C%time)
    call neko_log%end_section()
    call neko_log%newline()

  end subroutine simulation_adjoint_init

  !> Finalize a simulation of a case
  subroutine simulation_adjoint_finalize(C)
    type(adjoint_case_t), intent(inout) :: C
    logical :: output_at_end

    ! Run a final output if specified in the json
    call json_get_or_default(C%case%params, 'case.output_at_end', &
         output_at_end, .true.)
    call C%output_controller%execute(C%time, output_at_end)

    if (.not. (output_at_end) .and. C%time%t .lt. C%time%end_time) then
       call simulation_adjoint_joblimit_chkp(C, C%time%t)
    end if

    ! Finalize the user modules
    call C%case%user%finalize(C%time)

    call neko_log%end_section('Normal end.')

  end subroutine simulation_adjoint_finalize

  !> Compute a single time-step of an adjoint case
  subroutine simulation_adjoint_step(C, dt_controller, cfl, &
       tstep_loop_start_time, final_time)
    type(adjoint_case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: cfl
    type(time_step_controller_t), intent(inout) :: dt_controller
    real(kind=dp), intent(in) :: tstep_loop_start_time
    real(kind=rp), optional, intent(in) :: final_time
    real(kind=rp) :: t_bkp
    real(kind=dp) :: start_time, end_time, tstep_start_time
    character(len=LOG_SIZE) :: log_buf

    ! Setup the time step, and start time
    call profiler_start_region('Time-Step Adjoint')
    C%time%tstep = C%time%tstep + 1
    start_time = MPI_WTIME()
    tstep_start_time = start_time

    ! Compute the next time step
    ! NOTE. we should be wary here since CFL is based on the convective velocity
    ! not the adjoint velocity
    cfl = C%fluid_adj%compute_cfl(C%time%dt)
    call dt_controller%set_dt(C%time, cfl)
    if (dt_controller%is_variable_dt) cfl = C%fluid_adj%compute_cfl(C%time%dt)

    ! Calculate the cfl after the possibly varied dt
    ! cfl = C%fluid_adj%compute_cfl(C%time%dt)

    ! Advance time step from t to t+dt and print the status
    call simulation_settime(C%time, C%fluid_adj%ext_bdf)
    ! for cosmetic reasons we want the simulation to run backwards
    if (present(final_time)) then
       t_bkp = C%time%t
       C%time%t = final_time - t_bkp
    end if
    call C%time%status()

    call neko_log%begin()

    write(log_buf, '(A,E15.7,1x,A,E15.7)') 'CFL:', cfl, 'dt:', C%time%dt
    call neko_log%message(log_buf)

    ! Scalar step
    ! (Note that for the adjoint we should the adjoint_scalars first)
    if (allocated(C%adjoint_scalars)) then
       start_time = MPI_WTIME()
       call neko_log%section('Adjoint scalar')
       call C%adjoint_scalars%step(C%time, &
            C%case%fluid%ext_bdf, dt_controller)
       end_time = MPI_WTIME()
       write(log_buf, '(A,E15.7)') &
            'Scalar step time:      ', end_time-start_time
       call neko_log%end_section(log_buf)
    end if

    ! Fluid step
    call neko_log%section('Adjoint fluid')
    call C%fluid_adj%step(C%time, dt_controller)
    end_time = MPI_WTIME()
    write(log_buf, '(A,E15.7)') &
         'Fluid step time (s):   ', end_time-start_time
    call neko_log%end_section(log_buf)

    ! Postprocessing
    call neko_log%section('Postprocessing')

    ! Correct the time so the output fields are the same as the primal
    if (present(final_time)) then
       C%time%t = t_bkp
    end if

    ! Run any IO needed.
    call C%output_controller%execute(C%time)
    call simulation_adjoint_norm_output(C, C%time)

    call neko_log%end_section()

    ! End the step and print summary
    end_time = MPI_WTIME()
    call neko_log%section('Step summary')
    write(log_buf, '(A,I8,A,E15.7)') &
         'Total time for step ', C%time%tstep, ' (s): ', &
         end_time - tstep_start_time
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7)') &
         'Total elapsed time (s):           ', end_time-tstep_loop_start_time
    call neko_log%message(log_buf)

    call neko_log%end_section()
    call neko_log%end()
    call profiler_end_region


  end subroutine simulation_adjoint_step

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

    time%t = time%t + time%dt

  end subroutine simulation_settime

  !> Restart a case @a C from a given checkpoint
  subroutine simulation_adjoint_restart(C)
    type(adjoint_case_t), intent(inout) :: C
    integer :: i
    type(file_t) :: chkpf, previous_meshf
    character(len=LOG_SIZE) :: log_buf
    character(len=:), allocatable :: restart_file
    character(len=:), allocatable :: restart_mesh_file
    real(kind=rp) :: tol
    logical :: found

    call C%case%params%get('case.restart_file', restart_file, found)
    call C%case%params%get('case.restart_mesh_file', restart_mesh_file, found)

    if (found) then
       call previous_meshf%init(trim(restart_mesh_file))
       call previous_meshf%read(C%fluid_adj%chkp%previous_mesh)
    end if

    call C%case%params%get('case.mesh2mesh_tolerance', tol, found)

    if (found) C%case%fluid%chkp%mesh2mesh_tol = tol

    call chkpf%init(trim(restart_file))
    call chkpf%read(C%fluid_adj%chkp)
    C%time%dtlag = C%fluid_adj%chkp%dtlag
    C%time%tlag = C%fluid_adj%chkp%tlag

    ! Free the previous mesh, dont need it anymore
    do i = 1, size(C%time%dtlag)
       call C%case%fluid%ext_bdf%set_coeffs(C%time%dtlag)
    end do

    call C%fluid_adj%restart(C%case%chkp)
    call C%case%fluid%chkp%previous_mesh%free()
    if (allocated(C%adjoint_scalars)) then
       call C%adjoint_scalars%restart(C%case%chkp)
    end if

    C%time%t = real(C%case%fluid%chkp%restart_time(), kind=rp)
    call neko_log%section('Restarting from checkpoint')
    write(log_buf, '(A,A)') 'File :   ', trim(restart_file)
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7)') 'Time : ', C%time%t
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call C%output_controller%set_counter(C%time)
    if (C%norm_output_enabled) then
       call C%norm_output_ctrl%set_counter(C%time)
    end if
  end subroutine simulation_adjoint_restart

  subroutine simulation_adjoint_norm_output(C, time_output)
    type(adjoint_case_t), intent(inout) :: C
    type(time_state_t), intent(in) :: time_output
    type(vector_t), pointer :: data_line
    real(kind=rp) :: norm_l2
    integer :: n, idx

    if (.not. C%norm_output_enabled) return
    if (.not. C%norm_output_ctrl%check(time_output)) return

    n = C%fluid_adj%c_Xh%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       norm_l2 = device_glsc3(C%fluid_adj%u_adj%x_d, &
            C%fluid_adj%u_adj%x_d, C%fluid_adj%c_Xh%B_d, n) + &
            device_glsc3(C%fluid_adj%v_adj%x_d, &
            C%fluid_adj%v_adj%x_d, C%fluid_adj%c_Xh%B_d, n) + &
            device_glsc3(C%fluid_adj%w_adj%x_d, &
            C%fluid_adj%w_adj%x_d, C%fluid_adj%c_Xh%B_d, n)
    else
       norm_l2 = glsc3(C%fluid_adj%u_adj%x, C%fluid_adj%u_adj%x, &
            C%fluid_adj%c_Xh%B, n) + &
            glsc3(C%fluid_adj%v_adj%x, C%fluid_adj%v_adj%x, &
            C%fluid_adj%c_Xh%B, n) + &
            glsc3(C%fluid_adj%w_adj%x, C%fluid_adj%w_adj%x, &
            C%fluid_adj%c_Xh%B, n)
    end if

    norm_l2 = sqrt(norm_l2) / C%fluid_adj%c_Xh%volume

    call neko_scratch_registry%request(data_line, idx, 1, .false.)
    data_line%x = [norm_l2]
    call C%norm_output_file%write(data_line, time_output%t)
    call neko_scratch_registry%relinquish(idx)
    call C%norm_output_ctrl%register_execution()
  end subroutine simulation_adjoint_norm_output

  !> Write a checkpoint at joblimit
  subroutine simulation_adjoint_joblimit_chkp(C, t)
    type(adjoint_case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    type(file_t) :: chkpf
    character(len=:), allocatable :: chkp_format
    character(len=LOG_SIZE) :: log_buf
    character(len=10) :: format_str
    logical :: found

    call C%case%params%get('case.checkpoint_format', chkp_format, found)
    call C%case%fluid%chkp%sync_host()
    format_str = '.chkp'
    if (found) then
       if (chkp_format .eq. 'hdf5') then
          format_str = '.h5'
       end if
    end if
    call chkpf%init(C%case%output_directory // 'joblimit' // trim(format_str))
    call chkpf%write(C%case%fluid%chkp, t)
    write(log_buf, '(A)') '! saving checkpoint >>>'
    call neko_log%message(log_buf)

  end subroutine simulation_adjoint_joblimit_chkp

end module simulation_adjoint
