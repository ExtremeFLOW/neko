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
!> Convex limiting operations for conservative Euler edge corrections.
module euler_idp_limiter
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use num_types, only : rp
  use euler_idp, only : EULER_IDP_NCOMP
  implicit none
  private

  public :: euler_idp_state_is_admissible
  public :: euler_idp_specific_entropy
  public :: euler_idp_entropy_is_admissible
  public :: euler_idp_local_entropy_bounds
  public :: euler_idp_relax_density_bounds
  public :: euler_idp_limit_endpoint
  public :: euler_idp_limit_edge

contains

  !> Test positive density and the internal-energy-density threshold.
  pure logical function euler_idp_state_is_admissible(state, &
       internal_energy_floor) result(admissible)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp) :: energy_margin

    admissible = .false.
    if (.not. all(ieee_is_finite(state))) return
    if (state(1) .le. 0.0_rp) return
    energy_margin = state(1) * (state(5) - internal_energy_floor) - &
         0.5_rp * dot_product(state(2:4), state(2:4))
    admissible = ieee_is_finite(energy_margin) .and. &
         energy_margin .ge. 0.0_rp
  end function euler_idp_state_is_admissible

  !> Return log(p) - gamma*log(rho) for a gamma-law gas.
  pure real(kind=rp) function euler_idp_specific_entropy(state, gamma) &
       result(entropy)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: gamma
    real(kind=rp) :: internal_energy, pressure

    entropy = -huge(1.0_rp)
    if (.not. all(ieee_is_finite(state))) return
    if (state(1) .le. 0.0_rp) return
    internal_energy = state(5) - &
         0.5_rp * dot_product(state(2:4), state(2:4)) / state(1)
    pressure = (gamma - 1.0_rp) * internal_energy
    if (.not. ieee_is_finite(pressure) .or. pressure .le. 0.0_rp) return
    entropy = log(pressure) - gamma * log(state(1))
  end function euler_idp_specific_entropy

  !> Test the Euler admissible set and a local minimum entropy constraint.
  pure logical function euler_idp_entropy_is_admissible(state, gamma, &
       entropy_lower, internal_energy_floor) result(admissible)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: gamma, entropy_lower
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp) :: entropy, tolerance

    admissible = euler_idp_state_is_admissible(state, &
         internal_energy_floor)
    if (.not. admissible) return
    entropy = euler_idp_specific_entropy(state, gamma)
    tolerance = 64.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, abs(entropy_lower))
    admissible = ieee_is_finite(entropy) .and. &
         entropy .ge. entropy_lower - tolerance
  end function euler_idp_entropy_is_admissible

  !> Compute one-ring entropy minima from an immutable stage field.
  pure subroutine euler_idp_local_entropy_bounds(stage_entropy, left, right, &
       entropy_lower_bound)
    real(kind=rp), intent(in) :: stage_entropy(:,:,:,:)
    integer, intent(in) :: left(:,:), right(:,:)
    real(kind=rp), intent(out) :: entropy_lower_bound(:,:,:,:)
    integer :: edge

    entropy_lower_bound = stage_entropy
    do edge = 1, size(left, 2)
       associate(a => left(:,edge), b => right(:,edge))
         entropy_lower_bound(a(1),a(2),a(3),a(4)) = min( &
              entropy_lower_bound(a(1),a(2),a(3),a(4)), &
              stage_entropy(b(1),b(2),b(3),b(4)))
         entropy_lower_bound(b(1),b(2),b(3),b(4)) = min( &
              entropy_lower_bound(b(1),b(2),b(3),b(4)), &
              stage_entropy(a(1),a(2),a(3),a(4)))
       end associate
    end do
  end subroutine euler_idp_local_entropy_bounds

  !> Apply the averaging relaxation of Guermond et al. to density bounds.
  pure subroutine euler_idp_relax_density_bounds(strict_lower, strict_upper, &
       second_difference, nodal_mass, domain_volume, dimension, &
       relaxed_lower, relaxed_upper)
    real(kind=rp), intent(in) :: strict_lower, strict_upper
    real(kind=rp), intent(in) :: second_difference
    real(kind=rp), intent(in) :: nodal_mass, domain_volume
    integer, intent(in) :: dimension
    real(kind=rp), intent(out) :: relaxed_lower, relaxed_upper
    real(kind=rp) :: relaxation, r_h

    r_h = (nodal_mass / domain_volume)**(1.5_rp / real(dimension, rp))
    relaxation = abs(second_difference)
    relaxed_lower = max((1.0_rp - r_h) * strict_lower, &
         strict_lower - relaxation)
    relaxed_upper = strict_upper + relaxation
  end subroutine euler_idp_relax_density_bounds

  !> Limit one directed auxiliary correction to the Euler admissible set.
  pure subroutine euler_idp_limit_endpoint(base, correction, density_lower, &
       density_upper, entropy_lower, gamma, internal_energy_floor, limit, &
       density_limited, energy_limited, entropy_limited)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: density_lower, density_upper, entropy_lower
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(out) :: limit
    logical, intent(out) :: density_limited, energy_limited, entropy_limited
    real(kind=rp) :: entropy_tolerance
    real(kind=rp) :: trial(EULER_IDP_NCOMP)

    limit = 1.0_rp
    density_limited = .false.
    energy_limited = .false.
    entropy_limited = .false.
    if (.not. euler_idp_state_is_admissible(base, &
         internal_energy_floor) .or. base(1) .lt. density_lower .or. &
         base(1) .gt. density_upper) then
       limit = 0.0_rp
       density_limited = .true.
       energy_limited = .true.
       return
    end if
    if (.not. euler_idp_entropy_is_admissible(base, gamma, entropy_lower, &
         internal_energy_floor)) then
       limit = 0.0_rp
       entropy_limited = .true.
       return
    end if

    if (correction(1) .lt. 0.0_rp) then
       limit = min(limit, (base(1) - density_lower) / (-correction(1)))
    else if (correction(1) .gt. 0.0_rp) then
       limit = min(limit, (density_upper - base(1)) / correction(1))
    end if
    if (correction(1) .ne. 0.0_rp) then
       limit = max(0.0_rp, limit)
       density_limited = limit .lt. 1.0_rp
       if (density_limited .and. limit .gt. 0.0_rp) then
          limit = limit * (1.0_rp - 32.0_rp * epsilon(1.0_rp))
       end if
    end if

    trial = base + limit * correction
    energy_limited = .not. euler_idp_state_is_admissible(trial, &
         internal_energy_floor)
    entropy_limited = .not. euler_idp_entropy_is_admissible(trial, gamma, &
         entropy_lower, internal_energy_floor)
    if (.not. energy_limited .and. .not. entropy_limited) return

    if (energy_limited) then
       call limit_internal_energy(base, correction, internal_energy_floor, &
            limit)
    end if

    trial = base + limit * correction
    if (.not. euler_idp_entropy_is_admissible(trial, gamma, entropy_lower, &
         internal_energy_floor)) then
       entropy_tolerance = 64.0_rp * epsilon(1.0_rp) * &
            max(1.0_rp, abs(entropy_lower))
       call limit_entropy(base, correction, entropy_lower - &
            entropy_tolerance, gamma, limit)
    end if
    if (limit .gt. 0.0_rp) then
       limit = limit * (1.0_rp - 32.0_rp * epsilon(1.0_rp))
    end if
  end subroutine euler_idp_limit_endpoint

  !> Limit a segment with a safeguarded Newton--secant line search.
  pure subroutine limit_internal_energy(base, correction, energy_floor, &
       limit)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: energy_floor
    real(kind=rp), intent(inout) :: limit
    real(kind=rp) :: left, right, next_left, next_right
    real(kind=rp) :: value_left, value_right, value_next, value_tolerance
    real(kind=rp) :: slope
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration

    left = 0.0_rp
    right = limit
    value_left = internal_energy_margin(base, energy_floor)
    trial = base + right * correction
    value_right = internal_energy_margin(trial, energy_floor)
    value_tolerance = 256.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, abs(value_left), abs(value_right))
    if (value_left .le. 0.0_rp) then
       limit = 0.0_rp
       return
    end if

    do iteration = 1, 16
       if (right - left .le. 256.0_rp * epsilon(1.0_rp) * &
            max(1.0_rp, right)) exit
       if (value_left .le. value_right) exit

       slope = (value_right - value_left) / (right - left)
       next_left = left - value_left / slope
       if (.not. ieee_is_finite(next_left) .or. next_left .le. left .or. &
            next_left .ge. right) exit
       trial = base + next_left * correction
       value_next = internal_energy_margin(trial, energy_floor)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .lt. -value_tolerance) exit
       left = next_left
       if (abs(value_next) .le. value_tolerance) exit
       value_left = value_next

       trial = base + right * correction
       slope = internal_energy_margin_derivative(trial, correction)
       if (.not. ieee_is_finite(slope) .or. slope .ge. 0.0_rp) exit
       next_right = right - value_right / slope
       if (.not. ieee_is_finite(next_right) .or. next_right .le. left .or. &
            next_right .ge. right) exit
       trial = base + next_right * correction
       value_next = internal_energy_margin(trial, energy_floor)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .gt. value_tolerance) exit
       right = next_right
       if (abs(value_next) .le. value_tolerance) exit
       value_right = value_next
    end do
    limit = left
  end subroutine limit_internal_energy

  !> Limit the concave gamma-law entropy constraint by Newton--secant.
  pure subroutine limit_entropy(base, correction, entropy_lower, gamma, &
       limit)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower, gamma
    real(kind=rp), intent(inout) :: limit
    real(kind=rp) :: left, right, next_left, next_right
    real(kind=rp) :: value_left, value_right, value_next, value_tolerance
    real(kind=rp) :: slope
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration

    left = 0.0_rp
    right = limit
    value_left = entropy_margin(base, entropy_lower, gamma)
    trial = base + right * correction
    value_right = entropy_margin(trial, entropy_lower, gamma)
    value_tolerance = 256.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, abs(value_left), abs(value_right))
    if (value_left .le. 0.0_rp) then
       limit = 0.0_rp
       return
    end if

    do iteration = 1, 16
       if (right - left .le. 256.0_rp * epsilon(1.0_rp) * &
            max(1.0_rp, right)) exit
       if (value_left .le. value_right) exit

       slope = (value_right - value_left) / (right - left)
       next_left = left - value_left / slope
       if (.not. ieee_is_finite(next_left) .or. next_left .le. left .or. &
            next_left .ge. right) exit
       trial = base + next_left * correction
       value_next = entropy_margin(trial, entropy_lower, gamma)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .lt. -value_tolerance) exit
       left = next_left
       if (abs(value_next) .le. value_tolerance) exit
       value_left = value_next

       trial = base + right * correction
       slope = entropy_margin_derivative(trial, correction, entropy_lower, &
            gamma)
       if (.not. ieee_is_finite(slope) .or. slope .ge. 0.0_rp) exit
       next_right = right - value_right / slope
       if (.not. ieee_is_finite(next_right) .or. next_right .le. left .or. &
            next_right .ge. right) exit
       trial = base + next_right * correction
       value_next = entropy_margin(trial, entropy_lower, gamma)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .gt. value_tolerance) exit
       right = next_right
       if (abs(value_next) .le. value_tolerance) exit
       value_right = value_next
    end do
    limit = left
  end subroutine limit_entropy

  !> Internal-energy-density margin along a positive-density segment.
  pure real(kind=rp) function internal_energy_margin(state, energy_floor) &
       result(margin)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: energy_floor

    margin = state(5) - energy_floor - &
         0.5_rp * dot_product(state(2:4), state(2:4)) / state(1)
  end function internal_energy_margin

  !> Directional derivative of the internal-energy-density margin.
  pure real(kind=rp) function internal_energy_margin_derivative(state, &
       correction) result(derivative)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp) :: momentum_squared

    momentum_squared = dot_product(state(2:4), state(2:4))
    derivative = correction(5) - &
         dot_product(state(2:4), correction(2:4)) / state(1) + &
         0.5_rp * momentum_squared * correction(1) / state(1)**2
  end function internal_energy_margin_derivative

  !> Concave form of the gamma-law minimum-entropy constraint.
  pure real(kind=rp) function entropy_margin(state, entropy_lower, gamma) &
       result(margin)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower, gamma
    real(kind=rp) :: pressure

    pressure = (gamma - 1.0_rp) * internal_energy_margin(state, 0.0_rp)
    margin = pressure - exp(entropy_lower) * state(1)**gamma
  end function entropy_margin

  !> Directional derivative of the concave entropy margin.
  pure real(kind=rp) function entropy_margin_derivative(state, correction, &
       entropy_lower, gamma) result(derivative)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower, gamma

    derivative = (gamma - 1.0_rp) * &
         internal_energy_margin_derivative(state, correction) - &
         exp(entropy_lower) * gamma * state(1)**(gamma - 1.0_rp) * &
         correction(1)
  end function entropy_margin_derivative

  !> Return one symmetric coefficient for a complete vector edge correction.
  pure subroutine euler_idp_limit_edge(left_base, right_base, &
       left_correction, right_correction, left_density_lower, &
       left_density_upper, right_density_lower, right_density_upper, &
       left_entropy_lower, right_entropy_lower, gamma, &
       internal_energy_floor, limit, density_limited, energy_limited, &
       entropy_limited)
    real(kind=rp), intent(in) :: left_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_density_lower, left_density_upper
    real(kind=rp), intent(in) :: right_density_lower, right_density_upper
    real(kind=rp), intent(in) :: left_entropy_lower, right_entropy_lower
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(out) :: limit
    logical, intent(out) :: density_limited, energy_limited, entropy_limited
    real(kind=rp) :: left_limit, right_limit
    logical :: left_density_limited, right_density_limited
    logical :: left_energy_limited, right_energy_limited
    logical :: left_entropy_limited, right_entropy_limited

    call euler_idp_limit_endpoint(left_base, left_correction, &
         left_density_lower, left_density_upper, left_entropy_lower, gamma, &
         internal_energy_floor, left_limit, &
         left_density_limited, left_energy_limited, left_entropy_limited)
    call euler_idp_limit_endpoint(right_base, right_correction, &
         right_density_lower, right_density_upper, right_entropy_lower, &
         gamma, internal_energy_floor, right_limit, &
         right_density_limited, right_energy_limited, &
         right_entropy_limited)
    limit = min(left_limit, right_limit)
    density_limited = left_density_limited .or. right_density_limited
    energy_limited = left_energy_limited .or. right_energy_limited
    entropy_limited = left_entropy_limited .or. right_entropy_limited
  end subroutine euler_idp_limit_edge

end module euler_idp_limiter
