!> Simulation driver
module simulation
  use case
  use abbdf
  use log
  implicit none
  private

  public :: neko_solve
  
contains

  !> Main driver to solve a case @a C
  subroutine neko_solve(C)
    type(case_t), intent(inout) :: C
    real(kind=dp) :: t, start_time_org, start_time, end_time, cfl
    character(len=LOG_SIZE) :: log_buf    
    integer :: tstep

    t = 0d0
    tstep = 0
    call neko_log%section('Simulation')
    write(log_buf,'(A, E15.7)') 'dt :  ', C%params%dt
    call neko_log%message(log_buf)
    write(log_buf,'(A, E15.7,A,E15.7,A)') 'T  : [', 0d0,',',C%params%T_end,')'
    call neko_log%message(log_buf)
    call neko_log%newline()
    
    start_time_org = MPI_WTIME()
    do while (t .lt. C%params%T_end)
       tstep = tstep + 1
       start_time = MPI_WTIME()
       cfl = C%fluid%compute_cfl(C%params%dt)
       call neko_log%status(t, c%params%T_end)
       write(log_buf, '(A,I)') 'Time-step: ', tstep
       call neko_log%message(log_buf)
       call neko_log%begin()

       write(log_buf, '(A,E15.7,1x,A,E15.7)') 'CFL:', cfl, 'dt:', C%params%dt
       call neko_log%message(log_buf)

       
       call simulation_settime(t, C%params%dt, C%ab_bdf, C%tlag, C%dtlag, tstep)
       call neko_log%section('Fluid')       
       call C%fluid%step(t, tstep, C%ab_bdf)
       end_time = MPI_WTIME()
       write(log_buf, '(A,E15.7,A,E15.7)') &
            'Elapsed time (s):', end_time-start_time_org, ' Step time:', &
            end_time-start_time
       call neko_log%end_section(log_buf)
       call C%usr%usr_chk(t, C%params%dt, tstep,&
            C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, C%fluid%c_Xh)
       call neko_log%end()
       call C%s%sample(t)
    end do
    call neko_log%end_section('normal end.')
  end subroutine neko_solve

  subroutine simulation_settime(t, dt, ab_bdf, tlag, dtlag, step)
    real(kind=dp), intent(inout) :: t
    real(kind=dp), intent(in) :: dt
    type(abbdf_t), intent(inout) :: ab_bdf
    real(kind=dp), dimension(10) :: tlag
    real(kind=dp), dimension(10) :: dtlag
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

    call ab_bdf%set_bd(dtlag)
    call ab_bdf%set_abbd(dtlag)
    
  end subroutine simulation_settime

end module simulation

