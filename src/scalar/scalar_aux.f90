!> Auxiliary routines for fluid solvers
module scalar_aux
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use krylov, only : ksp_monitor_t
  use logger, only : LOG_SIZE
  use utils, only : neko_error, neko_warning
  use time_state, only : time_state_t
  implicit none
  private

  public :: scalar_step_info, scalar_step_info_reset_stabilized

  !> To track if the solver is stabilized
  logical :: stabilized = .false.

contains

  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual
  subroutine scalar_step_info(time, ksp_results, strict_convergence)
    type(ksp_monitor_t), dimension(:), intent(in) :: ksp_results
    type(time_state_t), intent(in) :: time
    logical, intent(in), optional :: strict_convergence
    logical :: converged, strict_conv
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    if (present(strict_convergence)) then
       strict_conv = strict_convergence
    else
       strict_conv = .false.
    end if

    ! Do the printing
    call ksp_results(1)%print_header()
    do i = 1, size(ksp_results)
       call ksp_results(i)%print_result(time%tstep)
    end do

    ! Check for convergence
    converged = .true.
    do i = 1, size(ksp_results)
       if (ieee_is_nan(ksp_results(i)%res_final)) then
          call neko_error("Scalar solver diverged for " // &
               trim(ksp_results(i)%name))
       end if

       if (.not. ksp_results(i)%converged) then
          converged = .false.
          log_buf = 'Scalar solver did not converge for ' &
               // trim(ksp_results(i)%name)

          if (.not. stabilized)then
             continue
          else if (strict_conv) then
             call neko_error(log_buf)
          else
             call neko_warning(log_buf)
          end if
       end if
    end do

    ! Update stabilized status
    if (.not. stabilized) stabilized = converged

  end subroutine scalar_step_info

  !> Resets the stabilized flag to false
  subroutine scalar_step_info_reset_stabilized()
    stabilized = .false.
  end subroutine scalar_step_info_reset_stabilized

end module scalar_aux
