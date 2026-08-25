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
!> Implements the CPU kernel for the `rough_log_law_t` type.
module rough_log_law_cpu
  use num_types, only : rp
  implicit none
  private

  public :: rough_log_law_compute_cpu

contains
  !> Compute the wall shear stress on CPU using the rough log-law model.
  !! @param tstep The current time-step.
  subroutine rough_log_law_compute_cpu(u, v, w, &
       n_x, n_y, n_z, h, tau_x, tau_y, tau_z, n_nodes, &
       kappa, rho_w, B, z0, tstep)
    integer, intent(in) :: n_nodes, tstep
    real(kind=rp), dimension(n_nodes), intent(in) :: u, v, w
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h
    real(kind=rp), dimension(n_nodes), intent(in) :: rho_w
    real(kind=rp), intent(in) :: kappa, B, z0
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, rho

    !$omp parallel do private(i, ui, vi, wi, magu, utau, normu, rho)
    do i=1, n_nodes
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

       ! Compute the wall shear stress using the rough log-law
       if (h(i) > z0) then
          utau = (magu - B) * kappa / log(h(i) / z0)
       else
          utau = 0.0_rp
       end if

       ! Distribute according to the velocity vector
       tau_x(i) = -rho*utau**2 * ui / magu
       tau_y(i) = -rho*utau**2 * vi / magu
       tau_z(i) = -rho*utau**2 * wi / magu
    end do
    !$omp end parallel do

  end subroutine rough_log_law_compute_cpu

end module rough_log_law_cpu
