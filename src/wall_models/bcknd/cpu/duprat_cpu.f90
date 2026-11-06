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
!> Implements the CPU kernel for the `duprat_t` type.
module duprat_cpu
  use num_types, only : rp
  use math, only : NEKO_EPS
  use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: duprat_compute_cpu, duprat_update_dpds_cpu, duprat_u_star

  !> Number of sub-intervals of the composite quadrature for U*.
  integer, parameter :: N_SUB = 12
  !> Number of Gauss-Legendre points per sub-interval.
  integer, parameter :: N_GL = 5
  !> Gauss-Legendre nodes on [0, 1].
  real(kind=rp), parameter :: GL_X(N_GL) = [ &
       0.0469100770306680_rp, 0.2307653449471585_rp, 0.5_rp, &
       0.7692346550528415_rp, 0.9530899229693320_rp]
  !> Gauss-Legendre weights on [0, 1].
  real(kind=rp), parameter :: GL_W(N_GL) = [ &
       0.1184634425280945_rp, 0.2393143352496832_rp, &
       0.2844444444444444_rp, 0.2393143352496832_rp, &
       0.1184634425280945_rp]

contains
  !> Compute the wall shear stress on cpu using the model of Duprat et al.
  !! @param u The x component of the sampled velocity.
  !! @param v The y component of the sampled velocity.
  !! @param w The z component of the sampled velocity.
  !! @param n_x The x component of the wall normal.
  !! @param n_y The y component of the wall normal.
  !! @param n_z The z component of the wall normal.
  !! @param nu The kinematic viscosity at the wall.
  !! @param rho_w The density at the wall.
  !! @param h The wall-normal distance of the sampling point.
  !! @param dpds The wall-tangential pressure gradient, signed with respect
  !! to the local flow direction (positive for an adverse gradient).
  !! @param tau_x The x component of the wall shear stress.
  !! @param tau_y The y component of the wall shear stress.
  !! @param tau_z The z component of the wall shear stress.
  !! @param alpha The ratio u_tau^2 / u_tau_p^2, output for diagnostics.
  !! @param n_nodes The number of wall nodes.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term in the eddy
  !! viscosity.
  !! @param A The damping constant of the eddy viscosity.
  !! @param tstep The current time-step.
  subroutine duprat_compute_cpu(u, v, w, n_x, n_y, n_z, nu, rho_w, h, dpds, &
       tau_x, tau_y, tau_z, alpha, n_nodes, kappa, beta, A, tstep)
    integer, intent(in) :: n_nodes, tstep
    real(kind=rp), dimension(n_nodes), intent(in) :: u, v, w
    real(kind=rp), dimension(n_nodes), intent(in) :: rho_w, dpds
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), dimension(n_nodes), intent(inout) :: alpha
    real(kind=rp), intent(in) :: kappa, beta, A
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess, rho, up

    !$omp parallel do private(i, ui, vi, wi, magu, utau, normu, guess, rho, up)
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
          alpha(i) = 1.0_rp
          cycle
       end if

       ! Pressure velocity u_p = |(nu / rho) dp/ds|^(1/3), Duprat et al.
       ! (2011), Eq. (4). The pressure p in Neko is the dynamic pressure.
       up = (nu(i) * abs(dpds(i)) / rho)**(1.0_rp / 3.0_rp)

       ! Get initial guess for Newton solver
       guess = tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2
       if (tstep .eq. 1 .or. guess .le. 0.0_rp) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          guess = sqrt(sqrt(guess) / rho)
       end if

       utau = solve_cpu(magu, h(i), guess, nu(i), up, sign(1.0_rp, dpds(i)), &
            kappa, beta, A)

       ! Distribute according to the velocity vector
       tau_x(i) = -rho*utau**2 * ui / magu
       tau_y(i) = -rho*utau**2 * vi / magu
       tau_z(i) = -rho*utau**2 * wi / magu

       ! alpha = u_tau^2 / u_tau_p^2, 1 without pressure gradient
       alpha(i) = 1.0_rp
       if (up .gt. 0.0_rp) alpha(i) = utau**2 / (utau**2 + up**2)
    end do
    !$omp end parallel do

  end subroutine duprat_compute_cpu

  !> Update the filtered wall-tangential pressure gradient on cpu.
  !!
  !! At each wall node, the sampled pressure gradient is projected on the
  !! direction of the sampled wall-parallel velocity, which gives the
  !! streamwise gradient dp/ds of Duprat et al. (2011), positive when adverse
  !! to the local flow. It is limited to
  !! \f$ |\partial p / \partial s| \le (\rho / \nu) (u \nu / h)^{3/2} \f$,
  !! i.e. the pressure velocity u_p may not exceed the friction velocity of
  !! a laminar (Stokes) profile through the sampled velocity, which guards
  !! the model against unphysical pressure spikes.
  !!
  !! The limited gradient is then filtered in time with an exponential
  !! moving average, a first-order low-pass filter,
  !! \f$ (\partial p / \partial s)_{filt}^{n} = (1 - \epsilon)
  !! (\partial p / \partial s)_{filt}^{n-1} + \epsilon
  !! (\partial p / \partial s)^{n} \f$,
  !! where the weight \f$ \epsilon = 1 - e^{-\Delta t / T} \f$ follows from
  !! the filter time constant T and is exact for variable time steps.
  !! @param dpx The x component of the sampled pressure gradient.
  !! @param dpy The y component of the sampled pressure gradient.
  !! @param dpz The z component of the sampled pressure gradient.
  !! @param u The x component of the sampled velocity.
  !! @param v The y component of the sampled velocity.
  !! @param w The z component of the sampled velocity.
  !! @param n_x The x component of the wall normal.
  !! @param n_y The y component of the wall normal.
  !! @param n_z The z component of the wall normal.
  !! @param nu The kinematic viscosity at the wall.
  !! @param rho_w The density at the wall.
  !! @param h The wall-normal distance of the sampling point.
  !! @param eps The filter weight of the current gradient, in [0, 1].
  !! @param dpds The filtered pressure gradient, updated in place.
  !! @param n_nodes The number of wall nodes.
  subroutine duprat_update_dpds_cpu(dpx, dpy, dpz, u, v, w, n_x, n_y, n_z, &
       nu, rho_w, h, eps, dpds, n_nodes)
    integer, intent(in) :: n_nodes
    real(kind=rp), dimension(n_nodes), intent(in) :: dpx, dpy, dpz
    real(kind=rp), dimension(n_nodes), intent(in) :: u, v, w
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in) :: nu, rho_w, h
    real(kind=rp), intent(in) :: eps
    real(kind=rp), dimension(n_nodes), intent(inout) :: dpds
    integer :: i
    real(kind=rp) :: ui, vi, wi, normu, magu, dp, dp_max

    !$omp parallel do private(i, ui, vi, wi, normu, magu, dp, dp_max)
    do i = 1, n_nodes
       ! Wall-parallel part of the sampled velocity
       normu = u(i) * n_x(i) + v(i) * n_y(i) + w(i) * n_z(i)
       ui = u(i) - normu * n_x(i)
       vi = v(i) - normu * n_y(i)
       wi = w(i) - normu * n_z(i)
       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! Without a flow direction, keep the previous filtered value
       if (magu .le. NEKO_EPS) cycle

       ! Streamwise pressure gradient, positive when adverse
       dp = (dpx(i) * ui + dpy(i) * vi + dpz(i) * wi) / magu

       ! Limit to the Stokes bound, u_p <= sqrt(u nu / h)
       dp_max = rho_w(i) / nu(i) * (magu * nu(i) / h(i))**1.5_rp
       dp = sign(min(abs(dp), dp_max), dp)

       ! Exponential moving average
       dpds(i) = (1.0_rp - eps) * dpds(i) + eps * dp
    end do
    !$omp end parallel do

  end subroutine duprat_update_dpds_cpu

  !> Dimensionless velocity U*(y*) of the model of Duprat et al. (2011).
  !! Integrates Eq. (5),
  !! \f$ \frac{dU^*}{dy^*} = \frac{s_p (1 - \alpha)^{3/2} y^* + \alpha}
  !! {1 + \nu_t / \nu} \f$, with the eddy viscosity of Eq. (6),
  !! \f$ \frac{\nu_t}{\nu} = \kappa y^* [\alpha + y^* (1 - \alpha)^{3/2}]^\beta
  !! (1 - e^{-y^* / (1 + A \alpha^3)})^2 \f$.
  !! The integral is evaluated with a composite Gauss-Legendre rule in the
  !! stretched coordinate \f$ s = \ln(1 + y^*) \f$, which resolves both the
  !! near-wall damping and the logarithmic region.
  !! @param ys The upper integration limit y*.
  !! @param alpha The ratio u_tau^2 / u_tau_p^2.
  !! @param sp The sign of the pressure gradient relative to the flow.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  pure function duprat_u_star(ys, alpha, sp, kappa, beta, A) result(us)
    real(kind=rp), intent(in) :: ys, alpha, sp, kappa, beta, A
    real(kind=rp) :: us
    real(kind=rp) :: ds, s0, eta, apg, damp, nut
    integer :: j, k

    us = 0.0_rp
    if (ys .le. 0.0_rp) return

    apg = max(1.0_rp - alpha, 0.0_rp)**1.5_rp
    damp = 1.0_rp + A * alpha**3
    ds = log(1.0_rp + ys) / real(N_SUB, rp)

    do k = 1, N_SUB
       s0 = real(k - 1, rp) * ds
       do j = 1, N_GL
          eta = exp(s0 + GL_X(j) * ds) - 1.0_rp
          nut = kappa * eta * (alpha + eta * apg)**beta * &
               (1.0_rp - exp(-eta / damp))**2
          ! dU*/dy* times the Jacobian d(eta)/ds = 1 + eta
          us = us + GL_W(j) * ds * (sp * apg * eta + alpha) / &
               (1.0_rp + nut) * (1.0_rp + eta)
       end do
    end do
  end function duprat_u_star

  !> Residual f(utau) = u_tau_p U*(y*) - u of the model of Duprat et al.
  !! @param utau The friction velocity.
  !! @param u The tangential velocity magnitude.
  !! @param y The wall-normal distance.
  !! @param nu The kinematic viscosity.
  !! @param up The pressure velocity.
  !! @param sp The sign of the pressure gradient relative to the flow.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  pure function residual(utau, u, y, nu, up, sp, kappa, beta, A) result(f)
    real(kind=rp), intent(in) :: utau, u, y, nu, up, sp, kappa, beta, A
    real(kind=rp) :: f, utp, alpha

    ! Extended velocity scale u_tau_p = sqrt(u_tau^2 + u_p^2), Eq. (4).
    ! alpha is formed without the square root, so that it cannot exceed 1
    ! by rounding.
    utp = sqrt(utau**2 + up**2)
    alpha = utau**2 / (utau**2 + up**2)
    f = utp * duprat_u_star(y * utp / nu, alpha, sp, kappa, beta, A) - u
  end function residual

  !> Safeguarded Newton solver for the friction velocity.
  !! The residual is not monotone in utau for favourable pressure gradients,
  !! so plain Newton can stall. The root is therefore first bracketed by a
  !! sign change of the residual, and Newton steps with a finite-difference
  !! derivative are replaced by bisection whenever they leave the bracket.
  !! @param u The tangential velocity magnitude.
  !! @param y The wall-normal distance.
  !! @param guess Initial guess.
  !! @param nu The kinematic viscosity.
  !! @param up The pressure velocity.
  !! @param sp The sign of the pressure gradient relative to the flow.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  function solve_cpu(u, y, guess, nu, up, sp, kappa, beta, A) result(utau)
    real(kind=rp), intent(in) :: u, y, guess, nu, up, sp, kappa, beta, A
    real(kind=rp) :: utau
    real(kind=rp) :: f, df, delta, lo, hi, f_lo, f_hi, step, tol, tiny
    integer :: k
    integer, parameter :: maxiter = 100
    logical :: converged
    character(len=LOG_SIZE) :: log_msg

    tol = max(1e-8_rp, 10.0_rp * NEKO_EPS)
    tiny = sqrt(NEKO_EPS) * u
    converged = .false.

    ! Bracket the root, lo with f(lo) < 0 and hi with f(hi) > 0, starting
    ! from the guess. f grows without bound for large utau.
    lo = max(guess, tiny)
    hi = lo
    f_lo = residual(lo, u, y, nu, up, sp, kappa, beta, A)
    f_hi = f_lo
    do k = 1, maxiter
       if (f_hi .gt. 0.0_rp) exit
       hi = 2.0_rp * hi
       f_hi = residual(hi, u, y, nu, up, sp, kappa, beta, A)
    end do
    do k = 1, maxiter
       if (f_lo .lt. 0.0_rp .or. lo .le. tiny) exit
       lo = max(0.5_rp * lo, tiny)
       f_lo = residual(lo, u, y, nu, up, sp, kappa, beta, A)
    end do

    ! No sign change down to utau ~ 0: the shear stress vanishes.
    if (f_lo .ge. 0.0_rp) then
       utau = lo
       return
    end if

    utau = 0.5_rp * (lo + hi)
    if (guess .gt. lo .and. guess .lt. hi) utau = guess
    f = residual(utau, u, y, nu, up, sp, kappa, beta, A)

    do k = 1, maxiter
       ! Shrink the bracket
       if (f .lt. 0.0_rp) then
          lo = utau
       else
          hi = utau
       end if

       ! Newton step with a finite-difference derivative, or bisection if
       ! the step leaves the bracket.
       delta = sqrt(NEKO_EPS) * utau
       df = (residual(utau + delta, u, y, nu, up, sp, kappa, beta, A) - f) &
            / delta
       step = f / df
       if (df .le. 0.0_rp .or. utau - step .le. lo .or. &
            utau - step .ge. hi) then
          step = utau - 0.5_rp * (lo + hi)
       end if
       utau = utau - step

       if (abs(step) .lt. tol * utau .or. &
            (hi - lo) .lt. tol * utau) then
          converged = .true.
          exit
       end if
       f = residual(utau, u, y, nu, up, sp, kappa, beta, A)
    end do

    if (.not. converged) then
       ! Called from inside an OpenMP loop, so serialise the log write.
       !$omp critical
       write(log_msg, *) "Newton not converged", lo, hi, utau
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
       !$omp end critical
    end if
  end function solve_cpu
end module duprat_cpu
