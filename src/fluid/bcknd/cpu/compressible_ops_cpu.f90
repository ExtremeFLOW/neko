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
!> CPU implementation of compressible flow operations
module compressible_ops_cpu
  use num_types, only : rp
  implicit none
  private

  public :: compressible_ops_cpu_compute_max_wave_speed, &
            compressible_ops_cpu_compute_entropy

contains

  !> Compute maximum wave speed for compressible flows on CPU
  subroutine compressible_ops_cpu_compute_max_wave_speed(max_wave_speed, u, v, w, gamma, p, rho, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: gamma
    real(kind=rp), dimension(n), intent(in) :: u, v, w, p, rho
    real(kind=rp), dimension(n), intent(inout) :: max_wave_speed
    integer :: i
    real(kind=rp) :: vel_mag, sound_speed

    ! Compute maximum wave speed: |u| + c = sqrt(u^2 + v^2 + w^2) + sqrt(gamma * p / rho)
    do concurrent (i = 1:n)
       vel_mag = sqrt(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
       sound_speed = sqrt(gamma * p(i) / rho(i))
       max_wave_speed(i) = vel_mag + sound_speed
    end do

  end subroutine compressible_ops_cpu_compute_max_wave_speed

  !> Compute entropy field S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho)) on CPU
  subroutine compressible_ops_cpu_compute_entropy(S, p, rho, gamma, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: gamma
    real(kind=rp), dimension(n), intent(in) :: p, rho
    real(kind=rp), dimension(n), intent(inout) :: S
    integer :: i

    ! Compute entropy: S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
    do concurrent (i = 1:n)
       S(i) = (1.0_rp / (gamma - 1.0_rp)) * rho(i) * &
            (log(p(i)) - gamma * log(rho(i)))
    end do

  end subroutine compressible_ops_cpu_compute_entropy

end module compressible_ops_cpu
