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
       density_upper, internal_energy_floor, limit, density_limited, &
       energy_limited)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: density_lower, density_upper
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp), intent(out) :: limit
    logical, intent(out) :: density_limited, energy_limited
    real(kind=rp) :: lower, upper, midpoint
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration

    limit = 1.0_rp
    density_limited = .false.
    energy_limited = .false.
    if (.not. euler_idp_state_is_admissible(base, &
         internal_energy_floor) .or. base(1) .lt. density_lower .or. &
         base(1) .gt. density_upper) then
       limit = 0.0_rp
       density_limited = .true.
       energy_limited = .true.
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
    if (euler_idp_state_is_admissible(trial, internal_energy_floor)) return

    energy_limited = .true.
    lower = 0.0_rp
    upper = limit
    do iteration = 1, 64
       midpoint = 0.5_rp * (lower + upper)
       trial = base + midpoint * correction
       if (euler_idp_state_is_admissible(trial, &
            internal_energy_floor)) then
          lower = midpoint
       else
          upper = midpoint
       end if
    end do
    limit = lower
    if (limit .gt. 0.0_rp) then
       limit = limit * (1.0_rp - 32.0_rp * epsilon(1.0_rp))
    end if
  end subroutine euler_idp_limit_endpoint

  !> Return one symmetric coefficient for a complete vector edge correction.
  pure subroutine euler_idp_limit_edge(left_base, right_base, &
       left_correction, right_correction, left_density_lower, &
       left_density_upper, right_density_lower, right_density_upper, &
       internal_energy_floor, limit, density_limited, energy_limited)
    real(kind=rp), intent(in) :: left_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_density_lower, left_density_upper
    real(kind=rp), intent(in) :: right_density_lower, right_density_upper
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp), intent(out) :: limit
    logical, intent(out) :: density_limited, energy_limited
    real(kind=rp) :: left_limit, right_limit
    logical :: left_density_limited, right_density_limited
    logical :: left_energy_limited, right_energy_limited

    call euler_idp_limit_endpoint(left_base, left_correction, &
         left_density_lower, left_density_upper, internal_energy_floor, &
         left_limit, left_density_limited, left_energy_limited)
    call euler_idp_limit_endpoint(right_base, right_correction, &
         right_density_lower, right_density_upper, internal_energy_floor, &
         right_limit, right_density_limited, right_energy_limited)
    limit = min(left_limit, right_limit)
    density_limited = left_density_limited .or. right_density_limited
    energy_limited = left_energy_limited .or. right_energy_limited
  end subroutine euler_idp_limit_edge

end module euler_idp_limiter
