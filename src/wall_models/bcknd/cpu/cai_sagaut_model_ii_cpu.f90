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
!> CPU backend for `cai_sagaut_model_ii_t`.
module cai_sagaut_model_ii_cpu
  use num_types, only : rp
  use math, only : lambert_w0, NEKO_EPS
  implicit none
  private

  public :: cai_sagaut_model_ii_compute_cpu

contains
  !> Evaluate wall shear stresses with the CPU Model-II kernel.
  !! @param u The sampled x-velocity field.
  !! @param v The sampled y-velocity field.
  !! @param w The sampled z-velocity field.
  !! @param ind_r The r-index array for sampled GLL points.
  !! @param ind_s The s-index array for sampled GLL points.
  !! @param ind_t The t-index array for sampled GLL points.
  !! @param ind_e The element-index array for sampled GLL points.
  !! @param n_x The x-component of the wall normals.
  !! @param n_y The y-component of the wall normals.
  !! @param n_z The z-component of the wall normals.
  !! @param nu The sampled kinematic viscosity at wall points.
  !! @param rho_w The sampled density at wall points.
  !! @param h The wall-model sampling distances.
  !! @param tau_x The x-component of the wall shear stress.
  !! @param tau_y The y-component of the wall shear stress.
  !! @param tau_z The z-component of the wall shear stress.
  !! @param n_nodes The number of wall points.
  !! @param lx The number of GLL points per direction.
  !! @param nelv The number of velocity elements.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  !! @param p The blending exponent.
  !! @param s The blending scale.
  subroutine cai_sagaut_model_ii_compute_cpu(u, v, w, ind_r, ind_s, ind_t, &
       ind_e, n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, n_nodes, &
       lx, nelv, kappa, B, p, s)
    integer, intent(in) :: n_nodes, lx, nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    real(kind=rp), dimension(n_nodes), intent(in) :: rho_w
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    real(kind=rp), dimension(n_nodes), intent(out) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, B, p, s
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, rho
    real(kind=rp) :: rey, e_const, blend, up, warg

    e_const = exp(kappa * B)

    do i = 1, n_nodes
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       rho = rho_w(i)

       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)

       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       magu = sqrt(ui**2 + vi**2 + wi**2)

       if (magu < NEKO_EPS) then
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
          cycle
       end if

       rey = magu * h(i) / nu(i)
       blend = exp(-((rey / s) ** p))
       warg = kappa * e_const * rey
       up = blend * sqrt(rey) + (1.0_rp - blend) * lambert_w0(warg, 1) / kappa
       utau = magu / (up + NEKO_EPS)

       tau_x(i) = -rho * utau**2 * ui / (magu + NEKO_EPS)
       tau_y(i) = -rho * utau**2 * vi / (magu + NEKO_EPS)
       tau_z(i) = -rho * utau**2 * wi / (magu + NEKO_EPS)
    end do
  end subroutine cai_sagaut_model_ii_compute_cpu

end module cai_sagaut_model_ii_cpu
