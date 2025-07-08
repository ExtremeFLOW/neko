!> Auxiliary routines for fluid solvers
module scalar_aux
  use logger, only: LOG_SIZE
  use krylov, only : ksp_monitor_t
  use time_state, only : time_state_t
  use utils, only : neko_error, neko_warning
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

contains

  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual
  subroutine scalar_step_info(time, ksp_results, strict_convergence)
    type(ksp_monitor_t), dimension(:), intent(in) :: ksp_results
    type(time_state_t), intent(in) :: time
    logical, intent(in), optional :: strict_convergence
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    ! Do the printing
    call ksp_results(1)%print_header()
    do i = 1, size(ksp_results)
       call ksp_results(i)%print_result(time%tstep)
    end do

    ! Check for convergence
    do i = 1, size(ksp_results)
       if (ieee_is_nan(ksp_results(i)%res_final)) then
          call neko_error("Scalar solver diverged for " // &
               trim(ksp_results(i)%name))
       end if

       if (present(strict_convergence)) then
          if (.not. ksp_results(i)%converged) then
             log_buf = 'Scalar solver did not converge for ' &
                  // trim(ksp_results(i)%name)

             if (strict_convergence) then
                call neko_error(log_buf)
             else
                call neko_warning(log_buf)
             end if
          end if
       end if
    end do

  end subroutine scalar_step_info
end module scalar_aux
