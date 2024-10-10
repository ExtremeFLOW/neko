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
  use mpi_f08
  use case, only : case_t
  use num_types, only : rp, dp
  use time_scheme_controller, only : time_scheme_controller_t
  use file, only : file_t
  use logger, only : LOG_SIZE, neko_log
  use jobctrl, only : jobctrl_time_limit
  use field, only : field_t
  use profiler, only : profiler_start, profiler_stop, &
       profiler_start_region, profiler_end_region
  use simcomp_executor, only : neko_simcomps
  use json_utils, only : json_get_or_default
  use time_step_controller, only : time_step_controller_t
  implicit none
  private

  public :: neko_solve

contains

  !> Main driver to solve a case @a C
  subroutine neko_solve(C)
    type(case_t), target, intent(inout) :: C
    real(kind=rp) :: t, cfl
    real(kind=dp) :: start_time_org, start_time, end_time, tstep_start_time
    character(len=LOG_SIZE) :: log_buf
    integer :: tstep
    character(len=:), allocatable :: restart_file
    logical :: output_at_end, found
    ! for variable_tsteping
    real(kind=rp) :: cfl_avrg = 0.0_rp
    type(time_step_controller_t) :: dt_controller

    t = 0d0
    tstep = 0
    call neko_log%section('Starting simulation')
    write(log_buf, '(A, E15.7,A,E15.7,A)') 'T  : [', 0d0, ',', C%end_time, ')'
    call neko_log%message(log_buf)
    call dt_controller%init(C%params)
    if (.not. dt_controller%if_variable_dt) then
       write(log_buf, '(A, E15.7)') 'dt :  ', C%dt
       call neko_log%message(log_buf)
    else
       write(log_buf, '(A, E15.7)') 'CFL :  ', dt_controller%set_cfl
       call neko_log%message(log_buf)
    end if

    call C%params%get('case.restart_file', restart_file, found)
    if (found .and. len_trim(restart_file) .gt. 0) then
       ! Restart the case
       call simulation_restart(C, t)

       ! Restart the simulation components
       call neko_simcomps%restart(t)
    end if

    !> Execute outputs and user-init before time loop
    call neko_log%section('Postprocessing')
    call C%output_controller%execute(t, tstep)

    call C%usr%user_init_modules(t, C%fluid%u, C%fluid%v, C%fluid%w,&
                                 C%fluid%p, C%fluid%c_Xh, C%params)
    call neko_log%end_section()
    call neko_log%newline()

    call profiler_start
    cfl = C%fluid%compute_cfl(C%dt)
    start_time_org = MPI_WTIME()

    do while (t .lt. C%end_time .and. (.not. jobctrl_time_limit()))
       call profiler_start_region('Time-Step')
       tstep = tstep + 1
       start_time = MPI_WTIME()
       tstep_start_time = start_time
       if (dt_controller%dt_last_change .eq. 0) then
          cfl_avrg = cfl
       end if
       call dt_controller%set_dt(C%dt, cfl, cfl_avrg, tstep)
       !calculate the cfl after the possibly varied dt
       cfl = C%fluid%compute_cfl(C%dt)

       call neko_log%status(t, C%end_time)
       write(log_buf, '(A,I6)') 'Time-step: ', tstep
       call neko_log%message(log_buf)
       call neko_log%begin()

       write(log_buf, '(A,E15.7,1x,A,E15.7)') 'CFL:', cfl, 'dt:', C%dt
       call neko_log%message(log_buf)
       call simulation_settime(t, C%dt, C%ext_bdf, C%tlag, C%dtlag, tstep)

       ! Run the preprocessing
       call neko_log%section('Preprocessing')
       call neko_simcomps%preprocess(t, tstep)
       call neko_log%end_section()

       call neko_log%section('Fluid')
       call C%fluid%step(t, tstep, C%dt, C%ext_bdf, dt_controller)
       end_time = MPI_WTIME()
       write(log_buf, '(A,E15.7)') &
            'Fluid step time (s):   ', end_time-start_time
       call neko_log%message(log_buf)
       write(log_buf, '(A,E15.7)') &
            'Total elapsed time (s):', end_time-start_time_org
       call neko_log%end_section(log_buf)

       ! Scalar step
       if (allocated(C%scalar)) then
          start_time = MPI_WTIME()
          call neko_log%section('Scalar')
          call C%scalar%step(t, tstep, C%dt, C%ext_bdf, dt_controller)
          end_time = MPI_WTIME()
          write(log_buf, '(A,E15.7)') &
            'Scalar step time:      ', end_time-start_time
          call neko_log%message(log_buf)
          write(log_buf, '(A,E15.7)') &
            'Total elapsed time (s):', end_time-start_time_org
          call neko_log%end_section(log_buf)
       end if

       call neko_log%section('Postprocessing')
       ! Execute all simulation components
       call neko_simcomps%compute(t, tstep)

       call C%output_controller%execute(t, tstep)

       ! Update material properties
       call C%usr%material_properties(t, tstep, C%fluid%rho,&
                                      C%fluid%mu, &
                                      C%scalar%cp, &
                                      C%scalar%lambda, &
                                      C%params)

       call C%usr%user_check(t, tstep, C%fluid%u, C%fluid%v, C%fluid%w, &
                             C%fluid%p, C%fluid%c_Xh, C%params)

       call neko_log%end_section()
       end_time = MPI_WTIME()
       call neko_log%section('Step summary')
       write(log_buf, '(A,I8,A,E15.7)') &
            'Total time for step ', tstep, ' (s): ', end_time-tstep_start_time
       call neko_log%message(log_buf)
       write(log_buf, '(A,E15.7)') &
            'Total elapsed time (s):           ', end_time-start_time_org
       call neko_log%message(log_buf)
       call neko_log%end_section()
       call neko_log%end()
       call profiler_end_region
    end do
    call profiler_stop

    call json_get_or_default(C%params, 'case.output_at_end',&
                             output_at_end, .true.)
    call C%output_controller%execute(t, tstep, output_at_end)

    if (.not. (output_at_end) .and. t .lt. C%end_time) then
       call simulation_joblimit_chkp(C, t)
    end if

    call C%usr%user_finalize_modules(t, C%params)

    call neko_log%end_section('Normal end.')

  end subroutine neko_solve

  subroutine simulation_settime(t, dt, ext_bdf, tlag, dtlag, step)
    real(kind=rp), intent(inout) :: t
    real(kind=rp), intent(in) :: dt
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    real(kind=rp), dimension(10) :: tlag
    real(kind=rp), dimension(10) :: dtlag
    integer, intent(in) :: step
    integer :: i


    do i = 10, 2, -1
       tlag(i) = tlag(i-1)
       dtlag(i) = dtlag(i-1)
    end do

    dtlag(1) = dt
    tlag(1) = t
    if (ext_bdf%ndiff .eq. 0) then
       dtlag(2) = dt
       tlag(2) = t
    end if

    t = t + dt

    call ext_bdf%set_coeffs(dtlag)

  end subroutine simulation_settime

  !> Restart a case @a C from a given checkpoint
  subroutine simulation_restart(C, t)
    implicit none
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    integer :: i
    type(file_t) :: chkpf, previous_meshf
    character(len=LOG_SIZE) :: log_buf
    character(len=:), allocatable :: restart_file
    character(len=:), allocatable :: restart_mesh_file
    real(kind=rp) :: tol
    logical :: found, check_cont

    call C%params%get('case.restart_file', restart_file, found)
    call C%params%get('case.restart_mesh_file', restart_mesh_file,&
                      found)

    if (found) then
       previous_meshf = file_t(trim(restart_mesh_file))
       call previous_meshf%read(C%fluid%chkp%previous_mesh)
    end if

    call C%params%get('case.mesh2mesh_tolerance', tol,&
                      found)

    if (found) C%fluid%chkp%mesh2mesh_tol = tol

    chkpf = file_t(trim(restart_file))
    call chkpf%read(C%fluid%chkp)
    C%dtlag = C%fluid%chkp%dtlag
    C%tlag = C%fluid%chkp%tlag

    !Free the previous mesh, dont need it anymore
    do i = 1, size(C%dtlag)
       call C%ext_bdf%set_coeffs(C%dtlag)
    end do

    call C%fluid%restart(C%dtlag, C%tlag)
    call C%fluid%chkp%previous_mesh%free()
    if (allocated(C%scalar)) &
       call C%scalar%restart( C%dtlag, C%tlag)

    t = C%fluid%chkp%restart_time()
    call neko_log%section('Restarting from checkpoint')
    write(log_buf, '(A,A)') 'File :   ', &
      trim(restart_file)
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7)') 'Time : ', t
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call C%output_controller%set_counter(t)
  end subroutine simulation_restart

  !> Write a checkpoint at joblimit
  subroutine simulation_joblimit_chkp(C, t)
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    type(file_t) :: chkpf
    character(len=:), allocatable :: chkp_format
    character(len=LOG_SIZE) :: log_buf
    character(len=10) :: format_str
    logical :: found

    call C%params%get('case.checkpoint_format', chkp_format, found)
    call C%fluid%chkp%sync_host()
    format_str = '.chkp'
    if (found) then
       if (chkp_format .eq. 'hdf5') then
          format_str = '.h5'
       end if
    end if
    chkpf = file_t('joblimit'//trim(format_str))
    call chkpf%write(C%fluid%chkp, t)
    write(log_buf, '(A)') '! saving checkpoint >>>'
    call neko_log%message(log_buf)

  end subroutine simulation_joblimit_chkp

end module simulation


