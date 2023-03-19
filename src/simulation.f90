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
  use ext_bdf_scheme
  use file
  use logger
  use jobctrl
  use profiler
  implicit none
  private

  public :: neko_solve
  
contains

  !> Main driver to solve a case @a C
  subroutine neko_solve(C)
    type(case_t), intent(inout) :: C
    real(kind=rp) :: t, cfl
    real(kind=dp) :: start_time_org, start_time, end_time
    character(len=LOG_SIZE) :: log_buf    
    integer :: tstep

    t = 0d0
    tstep = 0
    call neko_log%section('Starting simulation')
    write(log_buf,'(A, E15.7,A,E15.7,A)') 'T  : [', 0d0,',',C%params%T_end,')'
    call neko_log%message(log_buf)
    write(log_buf,'(A, E15.7)') 'dt :  ', C%params%dt
    call neko_log%message(log_buf)
    
    if (len_trim(C%params%restart_file) .gt. 0) then
       call simulation_restart(C, t)
    end if

    !> Call stats, samplers and user-init before time loop
    call neko_log%section('Postprocessing')       
    call C%q%eval(t, C%params%dt)
    call C%s%sample(t)
    call C%usr%user_init_modules(t, C%fluid%u, C%fluid%v, C%fluid%w,&
                                 C%fluid%p, C%fluid%c_Xh, C%params)
    call neko_log%end_section()
    call neko_log%newline()

    call profiler_start
    cfl = 0.0_rp
    start_time_org = MPI_WTIME()
    do while (t .lt. C%params%T_end .and. (.not. jobctrl_time_limit()) .and. cfl .lt. 0.7_rp)
       call profiler_start_region('Time-Step')
       tstep = tstep + 1
       start_time = MPI_WTIME()
       cfl = C%fluid%compute_cfl(C%params%dt)
       call neko_log%status(t, c%params%T_end)
       write(log_buf, '(A,I6)') 'Time-step: ', tstep
       call neko_log%message(log_buf)
       call neko_log%begin()

       write(log_buf, '(A,E15.7,1x,A,E15.7)') 'CFL:', cfl, 'dt:', C%params%dt
       call neko_log%message(log_buf)

       ! Fluid step 
       call simulation_settime(t, C%params%dt, C%ext_bdf, C%tlag, C%dtlag, tstep)

       call neko_log%section('Fluid')       
       call C%fluid%step(t, tstep, C%ext_bdf)
       end_time = MPI_WTIME()
       write(log_buf, '(A,E15.7,A,E15.7)') &
            'Elapsed time (s):', end_time-start_time_org, ' Step time:', &
            end_time-start_time
       call neko_log%end_section(log_buf)

       ! Scalar step
       if (allocated(C%scalar)) then
          start_time = MPI_WTIME()
          call neko_log%section('Scalar')       
          call C%scalar%step(t, tstep, C%ext_bdf)
          end_time = MPI_WTIME()
          write(log_buf, '(A,E15.7,A,E15.7)') &
               'Elapsed time (s):', end_time-start_time_org, ' Step time:', &
               end_time-start_time
          call neko_log%end_section(log_buf)
       end if                 

       call neko_log%section('Postprocessing')       
       call C%q%eval(t, C%params%dt)
       call C%s%sample(t)
       call C%usr%user_check(t, tstep,&
            C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, C%fluid%c_Xh, C%params)
       call neko_log%end_section()
       
       call neko_log%end()
       call profiler_end_region
    end do

    call profiler_stop

    if (t .lt. C%params%T_end) then
       call simulation_joblimit_chkp(C, t)
    end if
    
    call neko_log%end_section('Normal end.')
    
  end subroutine neko_solve

  subroutine simulation_settime(t, dt, ext_bdf, tlag, dtlag, step)
    real(kind=rp), intent(inout) :: t
    real(kind=rp), intent(in) :: dt
    type(ext_bdf_scheme_t), intent(inout) :: ext_bdf
    real(kind=rp), dimension(10) :: tlag
    real(kind=rp), dimension(10) :: dtlag
    integer, intent(in) :: step
    integer :: i
    

    do i = 10, 2, -1
       tlag(i) = tlag(i-1)
       dtlag(i) = dtlag(i-1)
    end do

    dtlag(1) = dt
    if (step .eq. 1) then
       dtlag(2) = dt
       tlag(2) = t
    end if

    t = t + dt

    call ext_bdf%set_bd(dtlag)
    call ext_bdf%set_abbd(dtlag)
    
  end subroutine simulation_settime

  !> Restart a case @a C from a given checkpoint
  subroutine simulation_restart(C, t)
    type(case_t), intent(inout) :: C
    real(kind=rp), intent(inout) :: t
    integer :: i
    type(file_t) :: chkpf
    character(len=LOG_SIZE) :: log_buf   


    chkpf = file_t(trim(C%params%restart_file))
    call chkpf%read(C%fluid%chkp)
    
    ! Make sure that continuity is maintained (important for interpolation) 
    call col2(C%fluid%u%x,C%fluid%c_Xh%mult,C%fluid%u%dof%size())
    call col2(C%fluid%v%x,C%fluid%c_Xh%mult,C%fluid%u%dof%size())
    call col2(C%fluid%w%x,C%fluid%c_Xh%mult,C%fluid%u%dof%size())
    call col2(C%fluid%p%x,C%fluid%c_Xh%mult,C%fluid%u%dof%size())
    select type (fld => C%fluid)
    type is(fluid_pnpn_t)
    do i = 1, fld%ulag%size()
       call col2(fld%ulag%lf(i)%x,fld%c_Xh%mult,fld%u%dof%size())
       call col2(fld%vlag%lf(i)%x,fld%c_Xh%mult,fld%u%dof%size())
       call col2(fld%wlag%lf(i)%x,fld%c_Xh%mult,fld%u%dof%size())
    end do
    end select
    if (allocated(C%scalar)) then
        call col2(C%scalar%s%x,C%scalar%c_Xh%mult, C%scalar%s%dof%size()) 
    end if

    call C%fluid%chkp%sync_device()
    call gs_op(C%fluid%gs_Xh,C%fluid%u,GS_OP_ADD)
    call gs_op(C%fluid%gs_Xh,C%fluid%v,GS_OP_ADD)
    call gs_op(C%fluid%gs_Xh,C%fluid%w,GS_OP_ADD)
    call gs_op(C%fluid%gs_Xh,C%fluid%p,GS_OP_ADD)
    select type (fld => C%fluid)
    type is(fluid_pnpn_t)
    do i = 1, fld%ulag%size()
       call gs_op(fld%gs_Xh,fld%ulag%lf(i),GS_OP_ADD)
       call gs_op(fld%gs_Xh,fld%vlag%lf(i),GS_OP_ADD)
       call gs_op(fld%gs_Xh,fld%wlag%lf(i),GS_OP_ADD)
    end do
    end select
 
    if (allocated(C%scalar)) then
       call gs_op(C%scalar%gs_Xh,C%scalar%s,GS_OP_ADD)
    end if
    t = C%fluid%chkp%restart_time()
    call neko_log%section('Restarting from checkpoint')
    write(log_buf,'(A,A)') 'File :   ', &
         trim(C%params%restart_file)
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

