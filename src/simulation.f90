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
  use case
  use gather_scatter
  use time_scheme_controller
  use file
  use math
  use logger
  use device
  use device_math
  use jobctrl
  use field, only : field_t
  use profiler
  use math, only : col2
  use simulation_component_global, only : neko_simcomps
  use json_utils, only : json_get_or_default
  implicit none
  private

  public :: neko_solve

contains

  !> Main driver to solve a case @a C
  subroutine neko_solve(C)
    implicit none
    type(case_t), intent(inout) :: C
    real(kind=rp) :: t, cfl
    real(kind=dp) :: start_time_org, start_time, end_time
    character(len=LOG_SIZE) :: log_buf
    integer :: tstep, i
    character(len=:), allocatable :: restart_file
    logical :: output_at_end, found

    t = 0d0
    tstep = 0
    call neko_log%section('Starting simulation')
    write(log_buf,'(A, E15.7,A,E15.7,A)') 'T  : [', 0d0,',',C%end_time,')'
    call neko_log%message(log_buf)
    write(log_buf,'(A, E15.7)') 'dt :  ', C%dt
    call neko_log%message(log_buf)

    call C%params%get('case.restart_file', restart_file, found)
    if (found .and. len_trim(restart_file) .gt. 0) then
       ! Restart the case
       call simulation_restart(C, t)

       ! Restart the simulation components
       if (allocated(neko_simcomps)) then
         do i=1, size(neko_simcomps)
            call neko_simcomps(i)%simcomp%restart(t)
         end do
       end if
    end if

    !> Call stats, samplers and user-init before time loop
    call neko_log%section('Postprocessing')
    call C%q%eval(t, C%dt, tstep)
    call C%s%sample(t, tstep)

    call C%usr%user_init_modules(t, C%fluid%u, C%fluid%v, C%fluid%w,&
                                 C%fluid%p, C%fluid%c_Xh, C%params)
    call neko_log%end_section()
    call neko_log%newline()

    call profiler_start

    start_time_org = MPI_WTIME()
    do while (t .lt. C%end_time .and. (.not. jobctrl_time_limit()))
       call profiler_start_region('Time-Step')
       tstep = tstep + 1
       start_time = MPI_WTIME()
       cfl = C%fluid%compute_cfl(C%dt)
       call neko_log%status(t, C%end_time)
       write(log_buf, '(A,I6)') 'Time-step: ', tstep
       call neko_log%message(log_buf)
       call neko_log%begin()

       write(log_buf, '(A,E15.7,1x,A,E15.7)') 'CFL:', cfl, 'dt:', C%dt
       call neko_log%message(log_buf)

       call simulation_settime(t, C%dt, C%ext_bdf, C%tlag, C%dtlag, tstep)

       call neko_log%section('Fluid')
       call C%fluid%step(t, tstep, C%dt, C%ext_bdf)
       end_time = MPI_WTIME()
       write(log_buf, '(A,E15.7,A,E15.7)') &
            'Elapsed time (s):', end_time-start_time_org, ' Step time:', &
            end_time-start_time
       call neko_log%end_section(log_buf)

       ! Scalar step
       if (allocated(C%scalar)) then
          start_time = MPI_WTIME()
          call neko_log%section('Scalar')
          call C%scalar%step(t, tstep, C%dt, C%ext_bdf)
          end_time = MPI_WTIME()
          write(log_buf, '(A,E15.7,A,E15.7)') &
               'Elapsed time (s):', end_time-start_time_org, ' Step time:', &
               end_time-start_time
          call neko_log%end_section(log_buf)
       end if

       call neko_log%section('Postprocessing')
       ! Execute all simulation components
       if (allocated(neko_simcomps)) then
         do i=1, size(neko_simcomps)
            call neko_simcomps(i)%simcomp%compute(t, tstep)
         end do
       end if

       call C%q%eval(t, C%dt, tstep)
       call C%s%sample(t, tstep)

       ! Update material properties
       call C%usr%material_properties(t, tstep, C%material_properties%rho,&
                                      C%material_properties%mu, &
                                      C%material_properties%cp, &
                                      C%material_properties%lambda, &
                                      C%params)

       call C%usr%user_check(t, tstep,&
            C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, C%fluid%c_Xh, C%params)

       call neko_log%end_section()

       call neko_log%end()
       call profiler_end_region
    end do

    call profiler_stop

    call json_get_or_default(C%params, 'case.output_at_end',&
                             output_at_end, .true.)
    call C%s%sample(t, tstep, output_at_end)

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
    logical :: found

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
    
    C%dtlag(:) = C%dt 
    C%tlag(:) = t
    do i = 1, size(C%tlag)
       C%tlag(i) = t - i*C%dtlag(i)
    end do

    chkpf = file_t(trim(restart_file))
    call chkpf%read(C%fluid%chkp)
    !Free the previous mesh, dont need it anymore
    call C%fluid%chkp%previous_mesh%free()
    do i = 1, size(C%dtlag)
       call C%ext_bdf%set_coeffs(C%dtlag)
    end do
    
    call C%fluid%restart(C%dtlag, C%tlag)
    if (allocated(C%scalar)) call C%scalar%restart( C%dtlag, C%tlag)

    t = C%fluid%chkp%restart_time()
    call neko_log%section('Restarting from checkpoint')
    write(log_buf,'(A,A)') 'File :   ', &
         trim(restart_file)
    call neko_log%message(log_buf)
    write(log_buf,'(A,E15.7)') 'Time : ', t
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call C%s%set_counter(t)
  end subroutine simulation_restart

  !> Write a checkpoint at joblimit
  subroutine simulation_joblimit_chkp(C, t)
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    type(file_t) :: chkpf
    character(len=LOG_SIZE) :: log_buf

    call C%fluid%chkp%sync_host()
    chkpf = file_t('joblimit.chkp')
    call chkpf%write(C%fluid%chkp, t)
    write(log_buf, '(A)') '! saving checkpoint >>>'
    call neko_log%message(log_buf)

  end subroutine simulation_joblimit_chkp

end module simulation

