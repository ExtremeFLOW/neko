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
       compressible_ops_cpu_compute_entropy, compressible_ops_cpu_update_uvw, &
       compressible_ops_cpu_update_mxyz_p_ruvw, compressible_ops_cpu_update_e


contains

  !> Compute maximum wave speed for compressible flows on CPU
  subroutine compressible_ops_cpu_compute_max_wave_speed(max_wave_speed, &
       u, v, w, gamma, p, rho, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: gamma
    real(kind=rp), dimension(n), intent(in) :: u, v, w, p, rho
    real(kind=rp), dimension(n), intent(inout) :: max_wave_speed
    integer :: i
    real(kind=rp) :: vel_mag, sound_speed

    ! Compute maximum wave speed:
    ! |u| + c = sqrt(u^2 + v^2 + w^2) + sqrt(gamma * p / rho)
    do concurrent (i = 1:n)
       vel_mag = sqrt(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
       sound_speed = sqrt(gamma * p(i) / rho(i))
       max_wave_speed(i) = vel_mag + sound_speed
    end do

  end subroutine compressible_ops_cpu_compute_max_wave_speed

  !> Compute entropy field
  !! S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho)) on CPU
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

  !> Update u,v,w fields
  subroutine compressible_ops_cpu_update_uvw(u, v, w, m_x, m_y, m_z, rho, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u, v, w
    real(kind=rp), dimension(n), intent(in) :: m_x, m_y, m_z, rho
    integer :: i

    do concurrent (i = 1:n)
       u(i) = m_x(i) / rho(i)
       v(i) = m_y(i) / rho(i)
       w(i) = m_z(i) / rho(i)
    end do

  end subroutine compressible_ops_cpu_update_uvw

  !> Update m_x, m_y, m_z, p, ruvw, fields
  subroutine compressible_ops_cpu_update_mxyz_p_ruvw(m_x, m_y, m_z, p, ruvw, &
       u, v, w, E, rho, gamma, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: m_x, m_y, m_z, p, ruvw
    real(kind=rp), dimension(n), intent(in) :: u, v, w, E, rho
    real(kind=rp), intent(in) :: gamma
    real(kind=rp) :: tmp
    integer :: i

    do concurrent (i = 1:n)
       m_x(i) = u(i) * rho(i)
       m_y(i) = v(i) * rho(i)
       m_z(i) = w(i) * rho(i)
    end do

    !Update p = (gamma - 1) * (E - 0.5 * rho * (u^2 + v^2 + w^2))
    do concurrent (i = 1:n)
       tmp = 0.5_rp * rho(i) * (u(i)**2 + v(i)**2 + w(i)**2)
       p(i) = (gamma - 1.0_rp) * (E(i) - tmp)
       ruvw(i) = tmp
    end do

  end subroutine compressible_ops_cpu_update_mxyz_p_ruvw

  !> Update E field
  subroutine compressible_ops_cpu_update_e(E, p, ruvw, gamma, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: E, p
    ! ruvw = 0.5 * rho * (u^2 + v^2 + w^2)
    real(kind=rp), dimension(n), intent(in) :: ruvw
    real(kind=rp), intent(in) :: gamma
    integer :: i


    do concurrent (i = 1:n)
       ! Ensure pressure is positive
       p(i) = max(p(i), 1.0e-12_rp)
       ! E = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2 + w^2)
       E(i) = p(i) * (1.0_rp / (gamma - 1.0_rp)) + ruvw(i)
    end do
  end subroutine compressible_ops_cpu_update_e

end module compressible_ops_cpu
