! Copyright (c) 2025, The Neko Authors
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

!> Auxiliary routines for fluid solvers
module fluid_aux
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use krylov, only : ksp_monitor_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use time_state, only : time_state_t
  implicit none
  private

  public :: fluid_step_info, fluid_step_info_reset_stabilized

  !> To track if the solver is stabilized
  logical, public, protected :: stabilized = .false.

contains

  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual
  subroutine fluid_step_info(time, ksp_results, full_stress_formulation, &
       strict_convergence)
    type(ksp_monitor_t), dimension(:), intent(in) :: ksp_results
    type(time_state_t), intent(in) :: time
    logical, intent(in) :: full_stress_formulation
    logical, intent(in), optional :: strict_convergence
    logical :: converged, strict_conv
    character(len=LOG_SIZE) :: log_buf
    integer :: i, n

    n = size(ksp_results)
    if (full_stress_formulation) n = 2

    if (present(strict_convergence)) then
       strict_conv = strict_convergence
    else
       strict_conv = .false.
    end if

    ! Do the printing
    call ksp_results(1)%print_header()
    do i = 1, n
       call ksp_results(i)%print_result(time%tstep)
    end do

    ! Check for convergence
    converged = .true.
    do i = 1, n
       if (ieee_is_nan(ksp_results(i)%res_final)) then
          call neko_error("Fluid solver diverged for " // &
               trim(ksp_results(i)%name))
       end if

       if (.not. ksp_results(i)%converged) then
          converged = .false.
          log_buf = 'Fluid solver did not converge for ' &
               // trim(ksp_results(i)%name)

          if (.not. stabilized) then
             continue
          else if (strict_conv) then
             call neko_error(log_buf)
          else
             call neko_log%message(log_buf)
          end if
       end if
    end do

    ! Update stabilized status
    if (.not. stabilized) stabilized = converged

  end subroutine fluid_step_info

  !> Resets the stabilized flag to false
  subroutine fluid_step_info_reset_stabilized()
    stabilized = .false.
  end subroutine fluid_step_info_reset_stabilized

end module fluid_aux
