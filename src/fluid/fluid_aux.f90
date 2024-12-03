!> Auxiliary routines for fluid solvers
module fluid_aux
  use logger, only : neko_log, LOG_SIZE
  use num_types, only : rp
  use krylov, only : ksp_monitor_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use utils, only : neko_error
  implicit none
  private

  public :: fluid_step_info

contains

  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual
  subroutine fluid_step_info(step, t, dt, ksp_results)
    type(ksp_monitor_t), intent(in) :: ksp_results(4)
    integer, intent(in) :: step
    real(kind=rp), intent(in) :: t, dt
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    call neko_log%message('Pressure')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(1)%iter, &
         ksp_results(1)%res_start, ksp_results(1)%res_final
    call neko_log%message(log_buf)

    call neko_log%message('X-Velocity')
    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(2)%iter, &
         ksp_results(2)%res_start, ksp_results(2)%res_final
    call neko_log%message(log_buf)

    call neko_log%message('Y-Velocity')
    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(3)%iter, &
         ksp_results(3)%res_start, ksp_results(3)%res_final
    call neko_log%message(log_buf)

    call neko_log%message('Z-Velocity')
    write(log_buf, '(A,A,A)') 'Iterations:   ', &
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(4)%iter, &
         ksp_results(4)%res_start, ksp_results(4)%res_final
    call neko_log%message(log_buf)

    ! Check for convergence
    do i = 1, 4
       if (.not. ksp_results(i)%converged) then
          call neko_error("Fluid solver did not converge")
       end if
    end do

  end subroutine fluid_step_info

end module fluid_aux
