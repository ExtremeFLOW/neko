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
!> Implements the CPU kernel for the `most_t` type.
module most_cpu
  use num_types, only : rp
  use utils, only : neko_error, neko_warning
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  public :: most_compute_cpu

  abstract interface
     function slaw_m_interface(z,L_ob,z0) result(slaw)
       import rp
       real(kind=rp), intent(in) :: z, L_ob, z0
       real(kind=rp) :: slaw
     end function slaw_m_interface

     function slaw_h_interface(z,L_ob,z0h) result(slaw)
       import rp
       real(kind=rp), intent(in) :: z, L_ob, z0h
       real(kind=rp) :: slaw
     end function slaw_h_interface

     function corr_m_interface(z,L_ob) result(corr)
       import rp
       real(kind=rp), intent(in) :: z, L_ob
       real(kind=rp) :: corr
     end function corr_m_interface

     function corr_h_interface(z,L_ob) result(corr)
       import rp
       real(kind=rp), intent(in) :: z, L_ob
       real(kind=rp) :: corr
     end function corr_h_interface

     function f_interface(Ri_b, z, z0, z0h, L_ob, slaw_m, slaw_h) result(f)
       import rp, slaw_m_interface, slaw_h_interface
       real(kind=rp), intent(in) :: Ri_b, z, z0, z0h, L_ob
       real(kind=rp) :: f
       procedure(slaw_m_interface) :: slaw_m
       procedure(slaw_h_interface) :: slaw_h
     end function f_interface

     function dfdl_interface(l_upper, l_lower, z, z0, z0h, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
       import rp, slaw_m_interface, slaw_h_interface
       real(kind=rp), intent(in) :: l_upper, l_lower, z, z0, z0h, L_ob, fd_h
       real(kind=rp) :: dfdl
       procedure(slaw_m_interface) :: slaw_m
       procedure(slaw_h_interface) :: slaw_h
     end function dfdl_interface
  end interface


  ! These will point to the correct functions
  ! depending on stability regime and bc_type.
  procedure(slaw_m_interface), pointer :: slaw_m_ptr => null()
  procedure(slaw_h_interface), pointer :: slaw_h_ptr => null()
  procedure(corr_m_interface), pointer :: corr_m_ptr => null()
  procedure(corr_h_interface), pointer :: corr_h_ptr => null()
  procedure(f_interface), pointer :: f_ptr => null()
  procedure(dfdl_interface), pointer :: dfdl_ptr => null()

contains

  !> Selects different expressions for the similarity functions in  MOST
  !> based on the type of bottom boundary condition for temperature.
  subroutine select_bc_operators(bc_type,bc_value,q,ts,ti,kappa,utau,z0h,hi)
    character(len=*), intent(in) :: bc_type
    real(kind=rp), intent(in) :: hi, ti, kappa, utau, z0h, bc_value
    real(kind=rp), intent(inout) :: q,ts
    select case (bc_type)
    case ("neumann")
       q = bc_value
       f_ptr => f_neumann
       dfdl_ptr => dfdl_neumann
    case ("dirichlet")
       ts = bc_value
       q = kappa*utau*(ts - ti)/log(hi/z0h)
       f_ptr => f_dirichlet
       dfdl_ptr => dfdl_dirichlet
    case default
       call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
    end select
  end subroutine select_bc_operators

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

  !> Sets the stability regime based on the Richardson number value (quite arbitrary).
  subroutine set_stability_regime(Ri_b,Ri_threshold)
    real(kind=rp), intent(in) :: Ri_b, Ri_threshold

    if (Ri_b > Ri_threshold) then
       slaw_m_ptr => slaw_m_stable
       slaw_h_ptr => slaw_h_stable
       corr_m_ptr => corr_m_stable
       corr_h_ptr => corr_h_stable
    elseif (Ri_b < -Ri_threshold) then
       slaw_m_ptr => slaw_m_convective
       slaw_h_ptr => slaw_h_convective
       corr_m_ptr => corr_m_convective
       corr_h_ptr => corr_h_convective
    else
       slaw_m_ptr => slaw_m_neutral
       slaw_h_ptr => slaw_h_neutral
    end if
  end subroutine set_stability_regime

  !> Main routine to compute the surface stresses based on MOST.
  !! @param tstep The current time-step
  subroutine most_compute_cpu(u, v, w, temp, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, z0, z0h_in, bc_type, bc_value, tstep)   
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w, temp
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), intent(in) :: kappa, z0, z0h_in, bc_value
    character(len=*), intent(in) :: bc_type
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer :: ts_idx(3)
    integer :: i, count
    integer, parameter :: max_count = 20
    real(kind=rp) :: ui, vi, wi, ti, ts, q, hi
    real(kind=rp) :: magu, utau, normu, z0h
    real(kind=rp) :: L_ob, L_upper, L_lower, L_old
    real(kind=rp) :: Ri_b, f, dfdl, fd_h, L_new, L_sign
    real(kind=rp), parameter :: g = 9.80665_rp
    real(kind=rp), parameter :: tol = 0.001_rp
    real(kind=rp), parameter :: NR_step = 0.001_rp
    real(kind=rp), parameter :: Ri_threshold = 0.0001_rp
    character(len=LOG_SIZE) :: log_buf

    do i=1, n_nodes
       ! Sample the variables
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       ti = temp(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       hi = h(i)

       ! Project on horizontal directions
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       ! Compute velocity magnitude
       magu = sqrt(ui**2 + vi**2)
       utau = magu*kappa / log(hi/z0)

       ! Compute thermal roughness length from Zilitinkevich, 1995
       if (z0h_in < 0) then
         z0h = z0 * exp(-0.1_rp*sqrt((utau*z0)/1.46e-5_rp))
       else 
         z0h = z0h_in
       end if 

       ! Get q, Ri_b, f_ptr, dfdl_ptr based on bc_type
       ! Maybe redundant, but needed to initialise Rib
       call select_bc_operators(bc_type,bc_value,q,ts,ti,kappa,utau,z0h,hi)
       call compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)

       ! Get q, Ri, f, dfdl based on bc_type
       call compute_Ri_b(bc_type, g, hi, ti, ts, magu, kappa, q, Ri_b)
       call set_stability_regime(Ri_b,Ri_threshold)
       if (abs(Ri_b) <= Ri_threshold) then
          ! Neutral (L_ob undefined)
          L_ob = 0.0_rp
       else
          ! Determine target regime sign
          if (Ri_b > 0.0_rp) then
             L_ob = hi / max(Ri_b, Ri_threshold) ! Stable guess
             L_sign = 1.0_rp
          else
             L_ob = hi / min(Ri_b, -Ri_threshold) ! Convective guess
             L_sign = -1.0_rp
          end if
            L_old = 1.0e10_rp
          count = 0
          
          ! Find Obukhov length
          do while ((abs(L_old - L_ob)/abs(L_ob) > tol) .and. (count < max_count))
             ! Switch between stable and convective based on bulk Richardson (Ri_b)
             L_old = L_ob
             count = count + 1
               fd_h = NR_step*L_ob
             L_upper = L_ob + fd_h
             L_lower = L_ob - fd_h
               ! Compute L_ob based on stability and bc_type
             if (.not. associated(f_ptr) .or. .not. associated(dfdl_ptr)) then
                call neko_error("Unassociated pointer for f or dfdl")
             end if
               f = f_ptr(Ri_b, hi, z0, z0h, L_ob, slaw_m_ptr, slaw_h_ptr)
             dfdl = dfdl_ptr(l_upper, l_lower, hi, z0, z0h, L_ob, slaw_m_ptr, slaw_h_ptr, fd_h)
             if (abs(dfdl) < 1.0e-12_rp) call neko_error("Division by zero in dfdl")
             L_new = L_ob - f/dfdl
               ! Avoid regime crossing during Newton iter (otherwise crash)
             if (L_new*L_sign <= 0.0_rp) then
                ! "damp update" (stay on same side)
                L_new = 0.5_rp * L_ob
             end if
               ! Bound L_ob
             L_ob = sign(max(abs(L_new), 1.0e-6_rp), L_sign)
             L_ob = sign(min(abs(L_ob), 1.0e6_rp), L_sign)
          end do
            !!!! this might lead to more crashes
          if (abs(L_ob) > 5e4_rp .or. abs(L_ob) < 1e-5_rp) then
             count = max_count
             call neko_warning("Obukhov length did not converge (MOST wall model)")
          end if
       end if

      ! Based on stability and bc_type, compute utau/q
      select case (bc_type)
      case ("neumann")
         ! Compute u* with the new Obukhov length
         utau = kappa*magu/slaw_m_ptr(hi, L_ob, z0)
      case ("dirichlet")
         ! Compute u* with the new Obukhov length
         utau = kappa*magu/slaw_m_ptr(hi, L_ob, z0)
         ! and compute q from here
         q = kappa*utau*(ts - ti)/slaw_h_ptr(hi, L_ob, z0h)
      case default
         call neko_error("Invalid specified temperature b.c. type ('neumann' or 'dirichlet'?)")
      end select

       ! Distribute according to the velocity vector and bound magu to avoid 0 division
       magu = max(magu, 1.0e-6_rp)
       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu
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

  end subroutine most_compute_cpu

  !> Similarity laws and corrections for the STABLE regime:
  function slaw_m_stable(z,L_ob,z0) result(slaw)
    real(kind=rp), intent(in) :: z,L_ob,z0
    real(kind=rp) :: slaw

    slaw = log(z/z0)-corr_m_stable(z,L_ob)+corr_m_stable(z0,L_ob)
  end function slaw_m_stable

  function slaw_h_stable(z,L_ob,z0h) result(slaw)
    real(kind=rp), intent(in) :: z,L_ob,z0h
    real(kind=rp) :: slaw

    slaw = log(z/z0h)-corr_h_stable(z,L_ob)+corr_h_stable(z0h,L_ob)
  end function slaw_h_stable

  function corr_m_stable(z,L_ob) result(corr)
    real(kind=rp), intent(in) :: z,L_ob
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d
    real(kind=rp) :: zeta
    zeta = z/L_ob
    a = 1.0_rp
    b = 0.6666666_rp
    c = 5.0_rp
    d = 0.35_rp
    corr = - a*zeta - b*(zeta-c/d)*exp(-d*zeta) - b*c/d
  end function corr_m_stable

  function corr_h_stable(z,L_ob) result(corr)
    real(kind=rp), intent(in) :: z,L_ob
    real(kind=rp) :: corr
    real(kind=rp) :: a, b, c, d
    real(kind=rp) :: zeta

    zeta = z/L_ob
    a = 1.0_rp
    b = 0.6666666_rp
    c = 5.0_rp
    d = 0.35_rp
    corr = -b * (zeta-c/d)*exp(-d*zeta)-(1.0_rp+ 0.6666666_rp * a * zeta)**1.5_rp-b*c/d + 1.0_rp
  end function corr_h_stable

  !> Similarity laws and corrections for the UNSTABLE (convective) regime:
  function slaw_m_convective(z,L_ob,z0) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0
    real(kind=rp) :: slaw

    slaw = log(z/z0) - corr_m_convective(z, L_ob) + corr_m_convective(z0, L_ob)
  end function slaw_m_convective

  function slaw_h_convective(z,L_ob,z0h) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0h
    real(kind=rp) :: slaw

    slaw = log(z/z0h) - corr_h_convective(z, L_ob) + corr_h_convective(z0h, L_ob)
  end function slaw_h_convective

  function corr_m_convective(z,L_ob) result(corr)
    real(kind=rp), intent(in) :: z, L_ob
    real(kind=rp) :: xi, pi, zeta
    real(kind=rp) :: corr

    zeta = z/L_ob
    pi = 4*atan(1.0_rp)
    xi = (1.0_rp - 16.0_rp*zeta)**0.25_rp
    corr = 2*log(0.5_rp*(1 + xi)) + log(0.5_rp*(1 + xi**2)) - 2*atan(xi) + pi/2
  end function corr_m_convective

  function corr_h_convective(z,L_ob) result(corr)
    real(kind=rp), intent(in) :: z, L_ob
    real(kind=rp) :: zeta, pi, xi
    real(kind=rp) :: corr

    zeta = z/L_ob
    pi = 4*atan(1.0_rp)
    xi = (1.0_rp - 16.0_rp*zeta)**0.25_rp
    corr = 2*log(0.5_rp*(1 + xi**2))
  end function corr_h_convective

  !> Similarity laws and corrections for the NEUTRAL regime:
  function slaw_m_neutral(z,L_ob,z0) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0
    real(kind=rp) :: slaw

    slaw = log(z/z0)
  end function slaw_m_neutral

  function slaw_h_neutral(z,L_ob,z0h) result(slaw)
    real(kind=rp), intent(in) :: z, L_ob, z0h
    real(kind=rp) :: slaw

    slaw = log(z/z0h)
  end function slaw_h_neutral

  !> Simialrity laws (different for neumann and dirichlet bc's)
  function f_neumann(Ri_b, z, z0, z0h, L_ob, slaw_m, slaw_h) result(f)
    real(kind=rp), intent(in) :: Ri_b, z, z0, z0h, L_ob
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: f

    f = (Ri_b - z/L_ob/slaw_m(z, L_ob, z0)**3)
  end function f_neumann

  function dfdl_neumann(l_upper, l_lower, z, z0, z0h, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
    real(kind=rp), intent(in) :: l_upper, l_lower, z, z0, z0h, L_ob, fd_h
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: dfdl

    dfdl = (-z/l_upper/slaw_m(z, l_upper, z0)**3) ! conv
    dfdl = dfdl + (z/l_lower/slaw_m(z, l_lower, z0)**3) ! conv
    dfdl = dfdl/(2*fd_h)
  end function dfdl_neumann

  function f_dirichlet(Ri_b, z, z0, z0h, L_ob, slaw_m, slaw_h) result(f)
    real(kind=rp), intent(in) :: Ri_b, z, z0, z0h, L_ob
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: f

    f = (Ri_b - z/L_ob*slaw_h(z, L_ob, z0h)/slaw_m(z, L_ob, z0)**2) ! conv
  end function f_dirichlet

  function dfdl_dirichlet(l_upper, l_lower, z, z0, z0h, L_ob, slaw_m, slaw_h, fd_h) result(dfdl)
    real(kind=rp), intent(in) :: l_upper, l_lower, z, z0, z0h, L_ob, fd_h
    procedure(slaw_m_interface) :: slaw_m
    procedure(slaw_h_interface) :: slaw_h
    real(kind=rp) :: dfdl

    dfdl = (-z/l_upper*slaw_h(z, l_upper, z0h)/slaw_m(z, l_upper, z0)**2) ! conv
    dfdl = dfdl + (z/l_lower*slaw_h(z, l_lower, z0h)/slaw_m(z, l_lower, z0)**2) ! conv
    dfdl = dfdl/(2*fd_h)
  end function dfdl_dirichlet


end module most_cpu
