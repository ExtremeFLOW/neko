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
!
!> Implements the CPU kernel for the `spalding_t` type.
module spalding_cpu
  use num_types, only : rp
  use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: spalding_compute_cpu

contains
  !> Compute the wall shear stress on cpu using Spalding's model.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine spalding_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, B, tstep)
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, B
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess

    do i=1, n_nodes
       ! Sample the velocity
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       ! Project on tangential direction
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)

       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! Get initial guess for Newton solver
       if (tstep .eq. 1) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          guess = tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2
          guess = sqrt(sqrt(guess))
       end if

       utau = solve_cpu(magu, h(i), guess, nu(i), kappa, B)

       ! Distribute according to the velocity vector
       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu
    end do

  end subroutine spalding_compute_cpu

  !> Newton solver for the algebraic equation defined by the law on cpu.
  !! @param u The velocity value.
  !! @param y The wall-normal distance.
  !! @param guess Initial guess.
  !! @param nu The molecular kinematic viscosity.
  !! @param kappa The von Karman constant.
  !! @param B The log-law intercept.
  function solve_cpu(u, y, guess, nu, kappa, B) result(utau)
    real(kind=rp), intent(in) :: u
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: guess
    real(kind=rp), intent(in) :: nu, kappa, B
    real(kind=rp) :: yp, up, utau
    real(kind=rp) :: error, f, df, old
    integer :: niter, k, maxiter
    character(len=LOG_SIZE) :: log_msg

    utau = guess

    maxiter = 100

    do k=1, maxiter
       up = u / utau
       yp = y * utau / nu
       niter = k
       old = utau

       ! Evaluate function and its derivative
       f = (up + exp(-kappa*B)* &
            (exp(kappa*up) - 1.0_rp - kappa*up - 0.5_rp*(kappa*up)**2 - &
            1.0_rp/6*(kappa*up)**3) - yp)

       df = (-y / nu - u/utau**2 - kappa*up/utau*exp(-kappa*B) * &
            (exp(kappa*up) - 1 - kappa*up - 0.5*(kappa*up)**2))

       ! Update solution
       utau = utau - f / df

       error = abs((old - utau)/old)

       if (error < 1e-3) then
          exit
       endif

    enddo

    if (niter .eq. maxiter) then
       write(log_msg, *) "Newton not converged", error, f, utau, old, guess
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
    end if
  end function solve_cpu
end module spalding_cpu
