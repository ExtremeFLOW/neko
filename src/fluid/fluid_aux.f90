!> Auxiliary routines for fluid solvers
module fluid_aux
  use logger, only : neko_log, LOG_SIZE
  use num_types, only : rp
  use krylov, only : ksp_monitor_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use utils, only : neko_error, neko_warning
  use comm, only : pe_rank
  implicit none
  private

  public :: fluid_step_info

contains

  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual
  subroutine fluid_step_info(step, t, dt, ksp_results, &
       full_stress_formulation, strict_convergence)
    type(ksp_monitor_t), dimension(:), intent(in) :: ksp_results
    integer, intent(in) :: step
    real(kind=rp), intent(in) :: t, dt
    logical, intent(in) :: full_stress_formulation
    logical, intent(in) :: strict_convergence
    character(len=LOG_SIZE) :: log_buf
    character(len=12) :: step_str
    integer :: i
    character(len=40) :: out_format

    ! Print the header
    write(log_buf, '((A5,7x),A3,(A5,5x),1x,A6,3x,A15,3x,A15)') &
         'Step:', ' | ', 'Field:', 'Iters:', &
         'Start residual:', 'Final residual:'
    call neko_log%message(log_buf)

    ! Define the output format
    out_format = '(A12,A3,A10,1x,I6,3x,E15.9,3x,E15.9)'
    write(step_str, '(I12)') step
    step_str = adjustl(step_str)

    ! Do the printing
    do i = 1, size(ksp_results)
       write(log_buf, out_format) &
            step_str, ' | ' , adjustl(ksp_results(i)%name), &
            ksp_results(i)%iter, &
            ksp_results(i)%res_start, ksp_results(i)%res_final
       call neko_log%message(log_buf)
    end do

    ! Check for convergence
    do i = 1, size(ksp_results)
       if (ieee_is_nan(ksp_results(i)%res_final)) then
          call neko_error("Fluid solver diverged")
       end if

       if ((.not. ksp_results(i)%converged) .and. (pe_rank .eq. 0)) then
          log_buf = 'Fluid solver did not converge for ' // &
               trim(ksp_results(i)%name)

          if (strict_convergence) then
             call neko_error(log_buf)
          else
             call neko_warning(log_buf)
          end if
       end if
    end do

  end subroutine fluid_step_info

end module fluid_aux
