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
module scalar_aux
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use krylov, only : ksp_monitor_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error, neko_warning
  use time_state, only : time_state_t
  implicit none
  private

  public :: scalar_step_info, scalar_step_info_reset_stabilized

  !> To track if the solver is stabilized
  logical, public, protected :: stabilized = .false.

contains

  !> Prints for prs, velx, vely, velz the following:
  !! Number of iterations, start residual, end residual
  subroutine scalar_step_info(time, ksp_results, strict_convergence, &
       allow_stabilization)
    type(ksp_monitor_t), intent(in) :: ksp_results
    type(time_state_t), intent(in) :: time
    logical, intent(in), optional :: strict_convergence
    logical, intent(in), optional :: allow_stabilization
    logical :: converged, strict_conv, allow_stab
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    if (present(strict_convergence)) then
       strict_conv = strict_convergence
    else
       strict_conv = .false.
    end if

    if (present(allow_stabilization)) then
       allow_stab = allow_stabilization
    else
       allow_stab = .false.
    end if

    ! Do the printing
    call ksp_results%print_result(time%tstep)

    ! Check for convergence
    converged = .true.
    if (ieee_is_nan(ksp_results%res_final)) then
       call neko_error("Scalar solver diverged for " // trim(ksp_results%name))
    end if

    if (.not. ksp_results%converged) then
       converged = .false.
       log_buf = 'Scalar solver did not converge for ' // trim(ksp_results%name)

       if (.not. stabilized .and. allow_stab) then
          continue
       else if (strict_conv) then
          call neko_error(log_buf)
       else
          call neko_log%message(log_buf)
       end if
    end if

    ! Update stabilized status
    if (.not. stabilized) stabilized = converged

  end subroutine scalar_step_info

  !> Resets the stabilized flag to false
  subroutine scalar_step_info_reset_stabilized()
    stabilized = .false.
  end subroutine scalar_step_info_reset_stabilized

end module scalar_aux
