! Copyright (c) 2026, The Neko Authors
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
!> Implements the CPU kernel for the `reichardt_t` type.
module reichardt_cpu
  use num_types, only : rp
  use math, only : NEKO_EPS
  use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: reichardt_compute_cpu

contains
  !> Compute the wall shear stress on cpu using Reichardt's law.
  !! @param u The x component of the sampled velocity.
  !! @param v The y component of the sampled velocity.
  !! @param w The z component of the sampled velocity.
  !! @param n_x The x component of the wall normal.
  !! @param n_y The y component of the wall normal.
  !! @param n_z The z component of the wall normal.
  !! @param nu The kinematic viscosity at the wall.
  !! @param rho_w The density at the wall.
  !! @param h The wall-normal distance of the sampling point.
  !! @param tau_x The x component of the wall shear stress.
  !! @param tau_y The y component of the wall shear stress.
  !! @param tau_z The z component of the wall shear stress.
  !! @param n_nodes The number of wall nodes.
  !! @param kappa The von Karman coefficient.
  !! @param C The amplitude of the exponential correction.
  !! @param B1 The damping length scale, in wall units.
  !! @param B2 The decay length scale of the second exponential term, in
  !! wall units.
  !! @param tstep The current time-step.
  subroutine reichardt_compute_cpu(u, v, w, &
       n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, n_nodes, &
       kappa, C, B1, B2, tstep)
    integer, intent(in) :: n_nodes, tstep
    real(kind=rp), dimension(n_nodes), intent(in) :: u, v, w
    real(kind=rp), dimension(n_nodes), intent(in) :: rho_w
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, C, B1, B2
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess, rho

    !$omp parallel do private(i, ui, vi, wi, magu, utau, normu, guess, rho)
    do i = 1, n_nodes
       ! Load the sampled velocity
       ui = u(i)
       vi = v(i)
       wi = w(i)
       rho = rho_w(i)

       ! Project on tangential direction
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)

       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! No tangential velocity, no shear stress
       if (magu .le. NEKO_EPS) then
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
          cycle
       end if

       ! Get initial guess for Newton solver
       guess = tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2
       if (tstep .eq. 1 .or. guess .le. 0.0_rp) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          guess = sqrt(sqrt(guess) / rho)
       end if

       utau = solve_cpu(magu, h(i), guess, nu(i), kappa, C, B1, B2)

       ! Distribute according to the velocity vector
       tau_x(i) = -rho*utau**2 * ui / magu
       tau_y(i) = -rho*utau**2 * vi / magu
       tau_z(i) = -rho*utau**2 * wi / magu
    end do
    !$omp end parallel do

  end subroutine reichardt_compute_cpu

  !> Evaluate Reichardt's law of the wall, u+ as a function of y+.
  !! @param yp The wall-normal distance in wall units.
  !! @param kappa The von Karman coefficient.
  !! @param C The amplitude of the exponential correction.
  !! @param B1 The damping length scale, in wall units.
  !! @param B2 The decay length scale of the second exponential term.
  pure function reichardt_up(yp, kappa, C, B1, B2) result(up)
    real(kind=rp), intent(in) :: yp, kappa, C, B1, B2
    real(kind=rp) :: up

    up = log(1.0_rp + kappa*yp) / kappa + C * &
         (1.0_rp - exp(-yp/B1) - yp/B1 * exp(-yp/B2))
  end function reichardt_up

  !> Evaluate the derivative du+/dy+ of Reichardt's law of the wall.
  !! @param yp The wall-normal distance in wall units.
  !! @param kappa The von Karman coefficient.
  !! @param C The amplitude of the exponential correction.
  !! @param B1 The damping length scale, in wall units.
  !! @param B2 The decay length scale of the second exponential term.
  pure function reichardt_dup(yp, kappa, C, B1, B2) result(dup)
    real(kind=rp), intent(in) :: yp, kappa, C, B1, B2
    real(kind=rp) :: dup

    dup = 1.0_rp / (1.0_rp + kappa*yp) + C / B1 * &
         (exp(-yp/B1) - (1.0_rp - yp/B2) * exp(-yp/B2))
  end function reichardt_dup

  !> Newton solver for the algebraic equation defined by the law on cpu.
  !! Solves f(utau) = utau * u+(y utau / nu) - u = 0.
  !! @param u The velocity value.
  !! @param y The wall-normal distance.
  !! @param guess Initial guess.
  !! @param nu The molecular kinematic viscosity.
  !! @param kappa The von Karman constant.
  !! @param C The amplitude of the exponential correction.
  !! @param B1 The damping length scale, in wall units.
  !! @param B2 The decay length scale of the second exponential term.
  function solve_cpu(u, y, guess, nu, kappa, C, B1, B2) result(utau)
    real(kind=rp), intent(in) :: u
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: guess
    real(kind=rp), intent(in) :: nu, kappa, C, B1, B2
    real(kind=rp) :: yp, up, utau
    real(kind=rp) :: error, f, df, old, tol
    integer :: k, maxiter
    logical :: converged
    character(len=LOG_SIZE) :: log_msg

    utau = guess
    maxiter = 100
    tol = max(1e-8_rp, 10.0_rp * NEKO_EPS)
    converged = .false.
    error = 0.0_rp
    f = 0.0_rp

    do k = 1, maxiter
       old = utau
       yp = y * utau / nu
       up = reichardt_up(yp, kappa, C, B1, B2)

       ! Evaluate function and its derivative
       f = utau * up - u
       df = up + yp * reichardt_dup(yp, kappa, C, B1, B2)

       ! Update solution, keeping utau positive
       utau = utau - f / df
       if (utau .le. 0.0_rp) utau = 0.5_rp * old

       error = abs((old - utau)/old)

       if (error .lt. tol) then
          converged = .true.
          exit
       end if
    end do

    if (.not. converged) then
       ! Called from inside an OpenMP loop, so serialise the log write.
       !$omp critical
       write(log_msg, *) "Newton not converged", error, f, utau
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
       !$omp end critical
    end if
  end function solve_cpu
end module reichardt_cpu
