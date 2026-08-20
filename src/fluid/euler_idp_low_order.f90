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
!     copyright notice, this list of conditions and the following disclaimer
!     in the documentation and/or other materials provided with the
!     distribution.
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
!> Algebraic operations for the low-order invariant-domain Euler scheme.
module euler_idp_low_order
  use num_types, only : rp
  use euler_idp, only : EULER_IDP_NCOMP
  implicit none
  private

  public :: euler_idp_maximum_wave_speed
  public :: euler_idp_flux_dot_vector
  public :: euler_idp_bar_state
  public :: euler_idp_internal_energy
  public :: euler_idp_internal_energy_timestep

contains

  !> Guaranteed upper bound for the two-state Euler Riemann fan speed.
  pure real(kind=rp) function euler_idp_maximum_wave_speed(left, right, &
       normal, gamma) result(speed)
    real(kind=rp), intent(in) :: left(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: normal(3), gamma
    real(kind=rp) :: reverse_normal(3)

    reverse_normal = -normal
    speed = max(euler_idp_ordered_wave_speed(left, right, normal, gamma), &
         euler_idp_ordered_wave_speed(right, left, reverse_normal, gamma))
  end function euler_idp_maximum_wave_speed

  !> Euler flux contracted with a Cartesian vector.
  pure subroutine euler_idp_flux_dot_vector(state, vector, gamma, flux)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: vector(3), gamma
    real(kind=rp), intent(out) :: flux(EULER_IDP_NCOMP)
    real(kind=rp) :: pressure, velocity(3), velocity_dot_vector

    velocity = state(2:4) / state(1)
    pressure = (gamma - 1.0_rp) * euler_idp_internal_energy(state)
    velocity_dot_vector = dot_product(velocity, vector)
    flux(1) = dot_product(state(2:4), vector)
    flux(2:4) = state(2:4) * velocity_dot_vector + pressure * vector
    flux(5) = (state(5) + pressure) * velocity_dot_vector
  end subroutine euler_idp_flux_dot_vector

  !> Construct the common bar state associated with one undirected edge.
  pure subroutine euler_idp_bar_state(left, right, flux_difference, &
       viscosity, bar_state)
    real(kind=rp), intent(in) :: left(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: flux_difference(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: viscosity
    real(kind=rp), intent(out) :: bar_state(EULER_IDP_NCOMP)

    bar_state = 0.5_rp * (left + right) - &
         flux_difference / (2.0_rp * viscosity)
  end subroutine euler_idp_bar_state

  !> Internal-energy density of one conserved Euler state.
  pure real(kind=rp) function euler_idp_internal_energy(state) result(value)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)

    value = state(5) - 0.5_rp * dot_product(state(2:4), state(2:4)) / &
         state(1)
  end function euler_idp_internal_energy

  !> Maximum timestep along one update direction respecting admissibility.
  pure real(kind=rp) function euler_idp_internal_energy_timestep(state, &
       residual, internal_energy_floor, upper_bound) result(dt_limit)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: residual(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp), intent(in) :: upper_bound
    real(kind=rp) :: lower, upper, midpoint
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration

    dt_limit = upper_bound
    if (residual(1) .gt. 0.0_rp) then
       dt_limit = min(dt_limit, state(1) / residual(1))
    end if
    dt_limit = max(0.0_rp, dt_limit)

    trial = state - dt_limit * residual
    if (trial(1) .gt. 0.0_rp) then
       if (euler_idp_internal_energy(trial) .ge. &
            internal_energy_floor) return
    end if

    lower = 0.0_rp
    upper = dt_limit
    do iteration = 1, 64
       midpoint = 0.5_rp * (lower + upper)
       trial = state - midpoint * residual
       if (trial(1) .gt. 0.0_rp) then
          if (euler_idp_internal_energy(trial) .ge. &
               internal_energy_floor) then
             lower = midpoint
          else
             upper = midpoint
          end if
       else
          upper = midpoint
       end if
    end do
    dt_limit = lower
  end function euler_idp_internal_energy_timestep

  !> Non-iterative Guermond-Popov upper bound for ordered Riemann data.
  pure real(kind=rp) function euler_idp_ordered_wave_speed(left, right, &
       normal, gamma) result(speed)
    real(kind=rp), intent(in) :: left(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: normal(3), gamma
    real(kind=rp) :: a_left, a_right, a_min, a_max
    real(kind=rp) :: pressure_left, pressure_right, pressure_min, pressure_max
    real(kind=rp) :: density_min
    real(kind=rp) :: velocity_left, velocity_right
    real(kind=rp) :: exponent, pressure_ratio, phi_min, phi_max
    real(kind=rp) :: pressure_two_rarefaction, pressure_upper
    real(kind=rp) :: coefficient_a, coefficient_b, numerator
    real(kind=rp) :: left_speed, right_speed

    pressure_left = (gamma - 1.0_rp) * &
         euler_idp_internal_energy(left)
    pressure_right = (gamma - 1.0_rp) * &
         euler_idp_internal_energy(right)
    a_left = sqrt(gamma * pressure_left / left(1))
    a_right = sqrt(gamma * pressure_right / right(1))
    velocity_left = dot_product(left(2:4), normal) / left(1)
    velocity_right = dot_product(right(2:4), normal) / right(1)

    if (pressure_left .le. pressure_right) then
       pressure_min = pressure_left
       pressure_max = pressure_right
       density_min = left(1)
       a_min = a_left
       a_max = a_right
    else
       pressure_min = pressure_right
       pressure_max = pressure_left
       density_min = right(1)
       a_min = a_right
       a_max = a_left
    end if

    exponent = (gamma - 1.0_rp) / (2.0_rp * gamma)
    pressure_ratio = (pressure_min / pressure_max)**exponent
    phi_min = 2.0_rp * a_max * (pressure_ratio - 1.0_rp) / &
         (gamma - 1.0_rp) + velocity_right - velocity_left
    if (phi_min .ge. 0.0_rp) then
       speed = max(max(-(velocity_left - a_left), 0.0_rp), &
            max(velocity_right + a_right, 0.0_rp))
       return
    end if

    coefficient_a = 2.0_rp / ((gamma + 1.0_rp) * density_min)
    coefficient_b = pressure_min * (gamma - 1.0_rp) / &
         (gamma + 1.0_rp)
    phi_max = (pressure_max - pressure_min) * &
         sqrt(coefficient_a / (pressure_max + coefficient_b)) + &
         velocity_right - velocity_left
    numerator = a_min + a_max - 0.5_rp * (gamma - 1.0_rp) * &
         (velocity_right - velocity_left)
    pressure_two_rarefaction = pressure_min * &
         (numerator / (a_min + a_max * pressure_ratio))** &
         (2.0_rp * gamma / (gamma - 1.0_rp))
    if (phi_max .lt. 0.0_rp) then
       pressure_upper = pressure_two_rarefaction
    else
       pressure_upper = min(pressure_max, pressure_two_rarefaction)
    end if

    left_speed = velocity_left - a_left * sqrt(1.0_rp + &
         max((pressure_upper - pressure_left) / pressure_left, 0.0_rp) * &
         (gamma + 1.0_rp) / (2.0_rp * gamma))
    right_speed = velocity_right + a_right * sqrt(1.0_rp + &
         max((pressure_upper - pressure_right) / pressure_right, 0.0_rp) * &
         (gamma + 1.0_rp) / (2.0_rp * gamma))
    speed = max(max(-left_speed, 0.0_rp), max(right_speed, 0.0_rp))
  end function euler_idp_ordered_wave_speed

end module euler_idp_low_order
