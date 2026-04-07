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
!> Implements the CPU kernel for the `richardson_t` type.
module richardson_cpu
  use num_types, only : rp
  use utils, only : neko_error
  use logger, only : LOG_SIZE, neko_log
  use math, only : glsum, glmin, glmax
  implicit none
  private

  public :: richardson_compute_cpu

  abstract interface
     function tau_interface(magu, Ri_b, h, z0, l, kappa) result(tau)
       import rp
       real(kind=rp), intent(in) :: magu, Ri_b, h, z0, l, kappa
       real(kind=rp) :: tau
     end function tau_interface

     function heat_flux_interface(ti, ts, Ri_b, h, magu, z0h, Pr,&
          l, utau, kappa) result(heat_flux)
       import rp
       real(kind=rp), intent(in) :: ts, ti, Ri_b, h, magu
       real(kind=rp), intent(in) :: z0h, Pr, l, utau, kappa
       real(kind=rp) :: heat_flux
     end function heat_flux_interface

  end interface

  ! These will point to the correct functions
  ! depending on stability regime and bc_type.
  procedure(tau_interface), pointer :: tau_ptr => null()
  procedure(heat_flux_interface), pointer :: heat_flux_ptr => null()

contains

  !> Computes the Richardson number.
  subroutine compute_Ri_b(bc_type, g_dot_n, hi, ti, ts, magu, kappa, q, Ri_b)
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(in) :: hi, ti, ts
    real(kind=rp), intent(in) :: magu, kappa, g_dot_n
    real(kind=rp), intent(inout) :: q, Ri_b

    select case (bc_type)
    case ("neumann")
       Ri_b = - g_dot_n*hi / ti*q / (magu**3*kappa**2)
    case ("dirichlet")
       Ri_b = g_dot_n*hi/ti*(ti - ts)/magu**2
    case default
       call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
    end select
  end subroutine compute_Ri_b

  !> Initialises q when the temperature surface bc is dirichlet.
  subroutine assign_bc_value(bc_type,bc_value,q,ts,ti,kappa,utau,z0h,hi)
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(in) :: hi, ti, kappa, utau, z0h,bc_value
    real(kind=rp), intent(inout) :: q,ts

    select case (bc_type)
    case ("neumann")
       ! ts not used
       q = bc_value
    case ("dirichlet")
       ts = bc_value
       q = kappa*utau*(ts - ti)/log(hi/z0h)
    case default
       call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
    end select
  end subroutine assign_bc_value

  !> Sets the stability regime based on the Richardson number value (quite arbitrary).
  subroutine set_stability_regime(Ri_b,Ri_threshold)
    real(kind=rp), intent(in) :: Ri_b, Ri_threshold

    if (Ri_b > Ri_threshold) then
       tau_ptr => tau_stable
       heat_flux_ptr => heat_flux_stable
    elseif (Ri_b < -Ri_threshold) then
       tau_ptr => tau_convective
       heat_flux_ptr => heat_flux_convective
    else
       tau_ptr => tau_neutral
       heat_flux_ptr => heat_flux_neutral
    end if
  end subroutine set_stability_regime

  !> Main routine to compute the surface stresses based on richardson.
  !! @param tstep The current time-step
  subroutine richardson_compute_cpu(u, v, w, temp, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, mu, rho, g_vec, Pr, z0, z0h_in, bc_type, bc_value, tstep)
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w, temp
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), intent(in) :: kappa, z0, z0h_in, bc_value, mu, rho, Pr
    real(kind=rp), dimension(3), intent(in) :: g_vec
    real(kind=rp) :: g_dot_n
    character(len=*), intent(in) :: bc_type
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer :: i
    real(kind=rp) :: ui, vi, wi, hi
    real(kind=rp) :: normu, z0h
    real(kind=rp) :: l
    real(kind=rp), parameter :: tol = 0.001_rp
    real(kind=rp), parameter :: NR_step = 0.001_rp
    real(kind=rp), parameter :: Ri_threshold = 0.0001_rp
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp), dimension(n_nodes) :: utau, Ri_b, L_ob, magu, q, ti, ts

    do i=1, n_nodes
       ! Sample the variables
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       ti(i) = temp(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       hi = h(i)

       ! Project on horizontal directions
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       ! Compute velocity magnitude
       magu(i) = sqrt(ui**2 + vi**2 + wi**2)
       magu(i) = max(magu(i), 1.0e-6_rp)
       utau(i) = magu(i)*kappa / log(hi/z0)

       ! Compute thermal roughness length from Zilitinkevich, 1995
       if (z0h_in < 0) then
          ! z0h_in is interpreted as -C_Zil (Zilitinkevich constant) for z0h
          z0h = z0 * exp(z0h_in*sqrt((utau(i)*z0)/(mu/rho)))
       else
          z0h = z0h_in
       end if

       ! Get q, ts based on bc_type
       ! Maybe redundant, but needed to initialise Rib
       call assign_bc_value(bc_type,bc_value,q(i),ts(i),ti(i),kappa,utau(i),z0h,hi)

       ! Compute g along the normal (generalisation for hills and similar)
       g_dot_n = abs(g_vec(1)*n_x(i) + g_vec(2)*n_y(i) + g_vec(3)*n_z(i))

       ! Compute Richardson and set stability accordingly
       call compute_Ri_b(bc_type, g_dot_n, hi, ti(i), ts(i), magu(i), kappa, q(i), Ri_b(i))
       call set_stability_regime(Ri_b(i), Ri_threshold)

       ! Set length scale
       l = kappa * hi
       ! Compute u*
       utau(i) = sqrt(tau_ptr(magu(i), Ri_b(i), hi, z0, l, kappa))
       select case (bc_type)
       case ("neumann")
          !!! TEMPORARY: neutral log-law approximation
          ts(i) = ti(i) - (q(i) * Pr * log(hi/z0h)) / (max(utau(i), 1e-6_rp) * kappa)
          q(i) = q(i)
       case ("dirichlet")
          ! Compute q
          q(i) = heat_flux_ptr(ti(i), ts(i), Ri_b(i), hi, magu(i), z0h, Pr, l, utau(i), kappa)
       case default
          call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
       end select

       ! Distribute according to the velocity vector and bound magu to avoid 0 division
       magu(i) = max(magu(i), 1.0e-6_rp)
       tau_x(i) = -rho*utau(i)**2 * ui / magu(i)
       tau_y(i) = -rho*utau(i)**2 * vi / magu(i)
       tau_z(i) = -rho*utau(i)**2 * wi / magu(i)
       if (abs(Ri_b(i)) <= Ri_threshold) then
          ! Neutral (L_ob undefined)
          L_ob(i) = 1e10_rp
       else
          L_ob(i) = -(ts(i)*utau(i)**3)/(kappa*g_dot_n*q(i))
       end if
    end do

    ! Print diagnostics
    call neko_log%section('Wall model diagnostics')
    write(log_buf, '(A,E15.7)') 'mean min max'
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'Ri_b: ', &
         glsum(Ri_b, n_nodes) / n_nodes, &
         glmin(Ri_b, n_nodes), glmax(Ri_b, n_nodes)
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'L_ob: ', &
         glsum(L_ob, n_nodes) / n_nodes, &
         glmin(L_ob, n_nodes), glmax(L_ob, n_nodes)
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'utau: ', &
         glsum(utau, n_nodes) / n_nodes, &
         glmin(utau, n_nodes), glmax(utau, n_nodes)
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'magu: ', &
         glsum(magu, n_nodes) / n_nodes, &
         glmin(magu, n_nodes), glmax(magu, n_nodes)
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'ti: ', &
         glsum(ti, n_nodes) / n_nodes, &
         glmin(ti, n_nodes), glmax(ti, n_nodes)
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'q: ', &
         glsum(q, n_nodes) / n_nodes, &
         glmin(q, n_nodes), glmax(q, n_nodes)
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'ts: ', &
         glsum(ts, n_nodes) / n_nodes, &
         glmin(ts, n_nodes), glmax(ts, n_nodes)
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7,E15.7,E15.7)') 'hi: ', &
         glsum(h, n_nodes) / n_nodes, &
         glmin(h, n_nodes), glmax(h, n_nodes)
    call neko_log%message(trim(log_buf))
    call neko_log%end_section()

  end subroutine richardson_compute_cpu

  !> Similarity laws and corrections for the STABLE regime:
  !> Based on Mauritsen et al. 2007
  function tau_stable(magu, Ri_b, h, z0, l, kappa) result(tau)
    real(kind=rp), intent(in) :: magu, Ri_b, h, z0, l, kappa
    real(kind=rp) :: tau

    tau = magu**2/(log(h/z0)**2) * f_tau_stable(Ri_b)/ &
         f_tau_stable(0.0_rp) * (l/h)**2
  end function tau_stable

  function heat_flux_stable(ti, ts, Ri_b, h, magu, z0h, Pr,&
       l, utau, kappa) result(heat_flux)
    real(kind=rp), intent(in) :: ts, ti, Ri_b, h, magu
    real(kind=rp), intent(in) :: z0h, Pr, l, utau, kappa
    real(kind=rp) :: heat_flux

    heat_flux = (ti - ts)/(log(h/z0h)) * &
         f_theta_stable(Ri_b)/abs(f_theta_stable(0.0_rp)) * &
         (l/h) * utau/Pr
  end function heat_flux_stable

  function f_tau_stable(Ri_b) result(f_tau)
    real(kind=rp), intent(in) :: Ri_b
    real(kind=rp) :: f_tau

    f_tau = 0.17 * (0.25 + 0.75 / (1.0 + 4.0*Ri_b))
  end function f_tau_stable

  function f_theta_stable(Ri_b) result(f_theta)
    real(kind=rp), intent(in) :: Ri_b
    real(kind=rp) :: f_theta

    f_theta = -0.145 / (1.0 + 4.0 * Ri_b)
  end function f_theta_stable

  !> Similarity laws and corrections for the UNSTABLE (convective) regime:
  !> Based on Louis 1979
  function tau_convective(magu, Ri_b, h, z0, l, kappa) result(tau)
    real(kind=rp), intent(in) :: magu, Ri_b, h, z0, l, kappa
    real(kind=rp) :: tau
    real(kind=rp) :: a, b, c

    a = kappa / log(h/z0)
    b = 2.0
    c = 7.4 * a**2 * b * (h/z0)**0.5

    tau = a**2 * magu**2 * f_tau_convective(Ri_b, c)
  end function tau_convective

  function heat_flux_convective(ti, ts, Ri_b, h, magu, z0h, Pr,&
       l, utau, kappa) result(heat_flux)
    real(kind=rp), intent(in) :: ts, ti, Ri_b, h, magu
    real(kind=rp), intent(in) :: z0h, Pr, l, utau, kappa
    real(kind=rp) :: heat_flux
    real(kind=rp) :: a, b, c

    a = kappa / log(h/z0h)
    b = 2.0
    c = 5.3 * a**2 * b * (h/z0h)**0.5

    heat_flux = - a**2 / 0.74 * magu * &
         (ti - ts) * f_theta_convective(Ri_b, c)

  end function heat_flux_convective

  function f_tau_convective(Ri_b, c) result(f_tau)
    real(kind=rp), intent(in) :: Ri_b, c
    real(kind=rp) :: f_tau

    f_tau = 1.0 - 2*Ri_b / (1.0 + c * abs(Ri_b)**0.5)
  end function f_tau_convective

  function f_theta_convective(Ri_b, c) result(f_theta)
    real(kind=rp), intent(in) :: Ri_b, c
    real(kind=rp) :: f_theta

    f_theta = 1.0 - 2*Ri_b / (1.0 + c * abs(Ri_b)**0.5)
  end function f_theta_convective

  !> Similarity laws and corrections for the NEUTRAL regime:
  function tau_neutral(magu, Ri_b, h, z0, l, kappa) result(tau)
    real(kind=rp), intent(in) :: magu, Ri_b, h, z0, l, kappa
    real(kind=rp) :: tau

    tau = (kappa*magu/log(h/z0))**2
  end function tau_neutral

  function heat_flux_neutral(ti, ts, Ri_b, h, magu, z0h, Pr,&
       l, utau, kappa) result(heat_flux)
    real(kind=rp), intent(in) :: ts, ti, Ri_b, h, magu
    real(kind=rp), intent(in) :: z0h, Pr, l, utau, kappa
    real(kind=rp) :: heat_flux

    heat_flux = kappa*utau * (ti - ts)/log(h/z0h)
  end function heat_flux_neutral

end module richardson_cpu
