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
!> Implements the CPU kernel for the `richardson_t` type.
module richardson_cpu
  use num_types, only : rp
  use utils, only : neko_error
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  public :: richardson_compute_cpu

  abstract interface
     function tau_interface(magu, ri_b, h, z0, l, kappa) result(tau)
       import rp
       real(kind=rp), intent(in) :: magu, ri_b, h, z0, l, kappa
       real(kind=rp) :: tau
     end function tau_interface

     function heat_flux_interface(ti, ts, ri_b, h, magu, z1, pr,&
                                l, utau, kappa) result(heat_flux)
       import rp
       real(kind=rp), intent(in) :: ts, ti, ri_b, h, magu
       real(kind=rp), intent(in) :: z1, pr, l, utau, kappa
       real(kind=rp) :: heat_flux
     end function heat_flux_interface

  end interface


  ! These will point to the correct functions
  ! depending on stability regime and bc_type.
  procedure(tau_interface), pointer :: tau_ptr => null()
  procedure(heat_flux_interface), pointer :: heat_flux_ptr => null()


contains

  !> Computes the Richardson number.
  subroutine compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(in) :: g, hi, ti, ts, magu, kappa
    real(kind=rp), intent(inout) :: q, Ri_b

    select case (bc_type)
    case ("neumann")
       Ri_b = - g*hi / ti*q / (magu**3*kappa**2)
    case ("dirichlet")
       Ri_b = g*hi/ti*(ti - ts)/magu**2
    case default
       call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
    end select
  end subroutine compute_Ri_b

  !> Initialises q when the temperature surface bc is dirichlet.
  subroutine init_q(bc_type, hi, ti, ts, kappa, utau, z0h, q)
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(in) :: hi, ti, ts, kappa, utau, z0h
    real(kind=rp), intent(inout) :: q

    if (bc_type == "dirichlet") then
       q = kappa*utau*(ts - ti)/log(hi/z0h)
    end if
  end subroutine init_q

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
       kappa, z0, bc_type, zone_idx, h_idx, q, tstep)   ! maybe put q as a variable that can be either q o ts
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w, temp
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), intent(in) :: kappa, z0
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(inout) :: q ! only supports scalar at the moment
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer, intent(in) :: zone_idx ! only supports wall model on ONE boundary atm!
    integer :: ts_idx(3)
    integer, intent(in) :: h_idx
    integer :: i
    real(kind=rp) :: ui, vi, ti, ts, hi
    real(kind=rp) :: magu, utau, normu, z0h
    real(kind=rp) :: L_ob, Ri_b, l
    real(kind=rp), parameter :: g = 9.80665_rp
    real(kind=rp), parameter :: Ri_threshold = 0.00001_rp
    character(len=LOG_SIZE) :: log_buf

    ! Select the ts offset based on fid
    select case (zone_idx)
    case (1)
       ts_idx = [h_idx, 0, 0 ]
    case (2)
       ts_idx = [-h_idx, 0, 0]
    case (3)
       ts_idx = [0, h_idx, 0 ]
    case (4)
       ts_idx = [0, -h_idx, 0]
    case (5)
       ts_idx = [0, 0, h_idx ]
    case (6)
       ts_idx = [0, 0, -h_idx]
    case default
       call neko_error("The face index is not correct (richardson_cpu.f90)")
    end select

    ! debug only:
    ! ts  = 300.0_rp
    ! q = 0.05_rp
    ! call neko_registry%add_field(this%coef%dof, "sampling_height", &
    !      ignore_existing=.true.)
    ! h_field => neko_registry%get_field_by_name("sampling_height")

    do i=1, n_nodes
       ! Sample the variables
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       ti = temp(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       ts = temp(ind_r(i)-ts_idx(1), ind_s(i)-ts_idx(2), ind_t(i)-ts_idx(3), ind_e(i))
       hi = h(i)

       ! Project on horizontal directions
       normu = ui * n_x(i) + vi * n_y(i)
       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)

       ! Compute velocity magnitude
       magu = sqrt(ui**2 + vi**2)

       ! utau initialisation
       if (tstep < 1) then
          utau = sqrt( sqrt( tau_x(i)**2 + tau_y(i)**2 ) )
       else
          utau = magu*kappa / log(hi/z0)

       end if
       ! Compute thermal roughness length from Zilitinkevich, 1995
       z0h = z0 * exp(-0.1_rp*sqrt((utau*z0)/1.46e-5_rp))

       ! Get q, Ri_b, f_ptr, dfdl_ptr based on bc_type
       ! Maybe redundant, but needed to initialise Rib
       call init_q(bc_type, hi, ti, ts, kappa, utau, z0h, q)
       call compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)

       ! Compute Obukhov length
       if (tstep > 0) then
          ! Get q, Ri, f, dfdl based on bc_type
          call compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)
          call set_stability_regime(Ri_b,Ri_threshold)

          ! Set length scale
          l = kappa * hi
          ! Compute u*
          utau = sqrt(tau_ptr(magu, ri_b, hi, z0, l, kappa))
          select case (bc_type)
          case ("neumann")
            ! Todo: Compute ts from q here
            q = q
          case ("dirichlet")
             ! Compute q
             q = heat_flux_ptr(ti, ts, ri_b, hi, magu, z0h, 1.0_rp, l, utau, kappa)
          case default
             call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
          end select

       end if

       ! Distribute according to the velocity vector and bound magu to avoid 0 division
       magu = max(magu, 1.0e-6_rp)
       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = 0
       if (abs(Ri_b) <= Ri_threshold) then
             ! Neutral (L_ob undefined)
             L_ob = 0.0_rp
       else
             L_ob = -(ts*utau**3)/(kappa*g*q)
       end if
    end do

    ! Print some indicative quantities (these are just point quantities: don't trust 100%)
    call neko_log%section('Wall model quick look')
    write(log_buf, '(A,E15.7)') 'Ri_b: ', Ri_b
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'L_ob: ', L_ob
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'utau: ', utau
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'magu: ', magu
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'ts: ', ts
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'ti: ', ti
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'q: ', q
    call neko_log%message(trim(log_buf))
    write(log_buf, '(A,E15.7)') 'hi: ', hi
    call neko_log%message(trim(log_buf))
    call neko_log%end_section()

  end subroutine richardson_compute_cpu

  !> Similarity laws and corrections for the STABLE regime:
  !> Based on Mauritsen et al. 2007
  function tau_stable(magu, ri_b, h, z0, l, kappa) result(tau)
    real(kind=rp), intent(in) :: magu, ri_b, h, z0, l, kappa
    real(kind=rp) :: tau

    tau = magu**2/(log(h/z0)**2) * f_tau_stable(ri_b)/ &
          f_tau_stable(0.0_rp) * (l/h)**2
  end function tau_stable

  function heat_flux_stable(ti, ts, ri_b, h, magu, z1, pr,&
                            l, utau, kappa) result(heat_flux)
    real(kind=rp), intent(in) :: ts, ti, ri_b, h, magu
    real(kind=rp), intent(in) :: z1, pr, l, utau, kappa
    real(kind=rp) :: heat_flux

    heat_flux = (ti - ts)/(log(h/z1)) * &
                f_theta_stable(ri_b)/abs(f_theta_stable(0.0_rp)) * &
                (l/h) * utau/pr
  end function heat_flux_stable

  function f_tau_stable(ri_b) result(f_tau)
    real(kind=rp), intent(in) :: ri_b
    real(kind=rp) :: f_tau

    f_tau = 0.17 * (0.25 + 0.75 / (1 + 4*ri_b))
  end function f_tau_stable

  function f_theta_stable(ri_b) result(f_theta)
    real(kind=rp), intent(in) :: ri_b
    real(kind=rp) :: f_theta

    f_theta = -0.145 / (1 + 4 * ri_b)
  end function f_theta_stable

  !> Similarity laws and corrections for the UNSTABLE (convective) regime:
  !> Based on Louis 1979
  function tau_convective(magu, ri_b, h, z0, l, kappa) result(tau)
    real(kind=rp), intent(in) :: magu, ri_b, h, z0, l, kappa
    real(kind=rp) :: tau
    real(kind=rp) :: a, b, c

    a =  kappa / log(h/z0)
    b = 2
    c = 7.4 * a**2 * b * (h/z0)**0.5

    tau = a**2 * magu**2 * f_tau_convective(ri_b, c)
  end function tau_convective

  function heat_flux_convective(ti, ts, ri_b, h, magu, z1, pr,&
                                l, utau, kappa) result(heat_flux)
    real(kind=rp), intent(in) :: ts, ti, ri_b, h, magu
    real(kind=rp), intent(in) :: z1, pr, l, utau, kappa
    real(kind=rp) :: heat_flux
    real(kind=rp) :: a, b, c

    a = kappa / log(h/z1)
    b = 2
    c = 5.3 * a**2 * b * (h/z1)**0.5

    heat_flux = - a**2 / 0.74 * magu * &
                       (ti - ts) * f_theta_convective(ri_b, c)

  end function heat_flux_convective

  function f_tau_convective(ri_b, c) result(f_tau)
    real(kind=rp), intent(in) :: ri_b, c
    real(kind=rp) :: f_tau

    f_tau = 1 - 2*ri_b / (1 + c * abs(ri_b)**0.5)
  end function f_tau_convective

  function f_theta_convective(ri_b, c) result(f_theta)
    real(kind=rp), intent(in) :: ri_b, c
    real(kind=rp) :: f_theta

    f_theta = 1 - 2*ri_b / (1 + c * abs(ri_b)**0.5)
  end function f_theta_convective

  !> Similarity laws and corrections for the NEUTRAL regime:
  function tau_neutral(magu, ri_b, h, z0, l, kappa) result(tau)
    real(kind=rp), intent(in) :: magu, ri_b, h, z0, l, kappa
    real(kind=rp) :: tau

    tau = (kappa*magu/log(h/z0))**2
  end function tau_neutral

  function heat_flux_neutral(ti, ts, ri_b, h, magu, z1, pr,&
                            l, utau, kappa) result(heat_flux)
    real(kind=rp), intent(in) :: ts, ti, ri_b, h, magu
    real(kind=rp), intent(in) :: z1, pr, l, utau, kappa
    real(kind=rp) :: heat_flux

    heat_flux = kappa*utau *(ti - ts)/log(h/z1)
  end function heat_flux_neutral

end module richardson_cpu
