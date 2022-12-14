!> Auxiliary routines for fluid solvers
module scalar_aux
  use logger
  use num_types
  use krylov, only : ksp_monitor_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

contains

  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual 
  subroutine scalar_step_info(step, t, dt, ksp_results)
    type(ksp_monitor_t), intent(in) :: ksp_results(1)
    integer, intent(in) :: step
    real(kind=rp), intent(in) :: t, dt
    character(len=LOG_SIZE) :: log_buf
    integer :: i


    call neko_log%message('Scalar')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') ksp_results(1)%iter, &
         ksp_results(1)%res_start, ksp_results(1)%res_final
    call neko_log%message(log_buf)

    ! Check for divergence
    do i = 1, 4
       if (ieee_is_nan(ksp_results(i)%res_final)) then
          call neko_log%error("Scalar solver diverged")
          stop
       end if
    end do

  end subroutine scalar_step_info


end module scalar_aux
