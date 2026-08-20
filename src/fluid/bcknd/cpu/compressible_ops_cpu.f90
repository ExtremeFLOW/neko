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
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, public, parameter :: EULER_STATE_OK = 0
  integer, public, parameter :: EULER_STATE_INVALID_GAMMA = 1
  integer, public, parameter :: EULER_STATE_INVALID_DENSITY = 2
  integer, public, parameter :: EULER_STATE_INVALID_INTERNAL_ENERGY = 3

  public :: compressible_ops_cpu_compute_max_wave_speed, &
       compressible_ops_cpu_compute_entropy, compressible_ops_cpu_update_uvw, &
       compressible_ops_cpu_update_mxyz_p_ruvw, compressible_ops_cpu_update_e, &
       compressible_ops_cpu_conserved_to_primitive


contains

  !> Convert conserved Euler states to primitives without modifying energy.
  !> @param rho Conserved density
  !> @param m_x,m_y,m_z Conserved momentum components
  !> @param E Conserved total-energy density
  !> @param gamma Ratio of specific heats
  !> @param internal_energy_floor Numerical internal-energy-density floor
  !> @param rho_primitive Primitive density
  !> @param u,v,w Primitive velocity components
  !> @param p Pressure
  !> @param sound_speed Sound speed
  !> @param internal_energy Internal-energy density
  !> @param status Admissibility status for each state
  !> @param n Number of states
  subroutine compressible_ops_cpu_conserved_to_primitive(rho, m_x, m_y, m_z, &
       E, gamma, internal_energy_floor, rho_primitive, &
       u, v, w, p, sound_speed, internal_energy, status, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: rho, m_x, m_y, m_z, E
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), dimension(n), intent(out) :: rho_primitive
    real(kind=rp), dimension(n), intent(out) :: u, v, w, p, sound_speed
    real(kind=rp), dimension(n), intent(out) :: internal_energy
    integer, dimension(n), intent(out) :: status
    real(kind=rp) :: kinetic_energy
    integer :: i

    if (.not. ieee_is_finite(gamma) .or. gamma .le. 1.0_rp) then
       rho_primitive = 0.0_rp
       u = 0.0_rp
       v = 0.0_rp
       w = 0.0_rp
       p = 0.0_rp
       sound_speed = 0.0_rp
       internal_energy = 0.0_rp
       status = EULER_STATE_INVALID_GAMMA
       return
    end if

    !$omp parallel do simd private(kinetic_energy)
    do i = 1, n
       rho_primitive(i) = rho(i)
       u(i) = 0.0_rp
       v(i) = 0.0_rp
       w(i) = 0.0_rp
       p(i) = 0.0_rp
       sound_speed(i) = 0.0_rp
       internal_energy(i) = 0.0_rp

       if (.not. ieee_is_finite(rho(i)) .or. rho(i) .le. 0.0_rp) then
          status(i) = EULER_STATE_INVALID_DENSITY
          cycle
       end if
       if (.not. ieee_is_finite(m_x(i)) .or. &
            .not. ieee_is_finite(m_y(i)) .or. &
            .not. ieee_is_finite(m_z(i)) .or. &
            .not. ieee_is_finite(E(i))) then
          status(i) = EULER_STATE_INVALID_INTERNAL_ENERGY
          cycle
       end if

       u(i) = m_x(i) / rho(i)
       v(i) = m_y(i) / rho(i)
       w(i) = m_z(i) / rho(i)
       kinetic_energy = 0.5_rp * &
            (m_x(i)**2 + m_y(i)**2 + m_z(i)**2) / rho(i)
       internal_energy(i) = E(i) - kinetic_energy
       if (.not. ieee_is_finite(internal_energy(i)) .or. &
            internal_energy(i) .lt. internal_energy_floor) then
          status(i) = EULER_STATE_INVALID_INTERNAL_ENERGY
          cycle
       end if

       p(i) = (gamma - 1.0_rp) * internal_energy(i)
       sound_speed(i) = sqrt(gamma * p(i) / rho(i))
       if (.not. ieee_is_finite(sound_speed(i))) then
          status(i) = EULER_STATE_INVALID_INTERNAL_ENERGY
          cycle
       end if
       status(i) = EULER_STATE_OK
    end do
    !$omp end parallel do simd
  end subroutine compressible_ops_cpu_conserved_to_primitive

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
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd private(vel_mag, sound_speed)
    do i = 1, n
       vel_mag = sqrt(u(i)*u(i) + v(i)*v(i) + &
            w(i)*w(i))
       sound_speed = sqrt(gamma * p(i) / rho(i))
       max_wave_speed(i) = vel_mag + sound_speed
    end do
    !$omp end parallel do simd

  end subroutine compressible_ops_cpu_compute_max_wave_speed

  !> Compute entropy field
  !! S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho)) on CPU
  subroutine compressible_ops_cpu_compute_entropy(S, p, rho, gamma, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: gamma
    real(kind=rp), dimension(n), intent(in) :: p, rho
    real(kind=rp), dimension(n), intent(inout) :: S
    integer :: i
    real(kind=rp) :: inv_gamma_m1

    inv_gamma_m1 = 1.0_rp / (gamma - 1.0_rp)

    ! Compute entropy: S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       S(i) = inv_gamma_m1 * rho(i) * &
            (log(p(i)) - gamma * log(rho(i)))
    end do
    !$omp end parallel do simd

  end subroutine compressible_ops_cpu_compute_entropy

  !> Update u,v,w fields
  subroutine compressible_ops_cpu_update_uvw(u, v, w, m_x, m_y, m_z, rho, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u, v, w
    real(kind=rp), dimension(n), intent(in) :: m_x, m_y, m_z, rho
    integer :: i

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       u(i) = m_x(i) / rho(i)
       v(i) = m_y(i) / rho(i)
       w(i) = m_z(i) / rho(i)
    end do
    !$omp end parallel do simd

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

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd private(tmp)
    do i = 1, n
       m_x(i) = u(i) * rho(i)
       m_y(i) = v(i) * rho(i)
       m_z(i) = w(i) * rho(i)
       tmp = 0.5_rp * rho(i) * (u(i)**2 + v(i)**2 + w(i)**2)
       p(i) = (gamma - 1.0_rp) * (E(i) - tmp)
       ruvw(i) = tmp
    end do
    !$omp end parallel do simd

  end subroutine compressible_ops_cpu_update_mxyz_p_ruvw

  !> Update E field
  subroutine compressible_ops_cpu_update_e(E, p, ruvw, gamma, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: E, p
    ! ruvw = 0.5 * rho * (u^2 + v^2 + w^2)
    real(kind=rp), dimension(n), intent(in) :: ruvw
    real(kind=rp), intent(in) :: gamma
    integer :: i
    real(kind=rp) :: inv_gamma_m1

    inv_gamma_m1 = 1.0_rp / (gamma - 1.0_rp)

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp parallel do simd
    do i = 1, n
       ! Ensure pressure is positive
       p(i) = max(p(i), 1.0e-12_rp)
       ! E = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2 + w^2)
       E(i) = p(i) * inv_gamma_m1 + ruvw(i)
    end do
    !$omp end parallel do simd
  end subroutine compressible_ops_cpu_update_e

end module compressible_ops_cpu
