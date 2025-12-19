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
!> CPU backend for entropy viscosity regularization
module entropy_viscosity_cpu
  use num_types, only : rp
  implicit none
  private

  public :: entropy_viscosity_compute_residual_cpu, &
            entropy_viscosity_compute_viscosity_cpu, &
            entropy_viscosity_apply_element_max_cpu, &
            entropy_viscosity_clamp_to_low_order_cpu, &
            entropy_viscosity_smooth_divide_cpu

contains

  !> Compute entropy residual on CPU
  !! @param entropy_residual Output entropy residual field
  !! @param S Current entropy field
  !! @param S_lag1 First lagged entropy field
  !! @param S_lag2 Second lagged entropy field
  !! @param S_lag3 Third lagged entropy field
  !! @param bdf_coeffs BDF time scheme coefficients
  !! @param dt Time step size
  !! @param n Number of points
  subroutine entropy_viscosity_compute_residual_cpu(entropy_residual, &
       S, S_lag1, S_lag2, S_lag3, bdf_coeffs, dt, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: entropy_residual
    real(kind=rp), dimension(n), intent(in) :: S, S_lag1, S_lag2, S_lag3
    real(kind=rp), intent(in) :: bdf_coeffs(4)
    real(kind=rp), intent(in) :: dt
    integer :: i

    do concurrent (i = 1:n)
       entropy_residual(i) = (bdf_coeffs(1) * S(i) &
                              - bdf_coeffs(2) * S_lag1(i) &
                              - bdf_coeffs(3) * S_lag2(i) &
                              - bdf_coeffs(4) * S_lag3(i)) / dt
    end do

  end subroutine entropy_viscosity_compute_residual_cpu

  !> Compute viscosity from entropy residual on CPU
  !! @param reg_coeff Output regularization coefficient field
  !! @param entropy_residual Entropy residual field
  !! @param h Mesh size field
  !! @param c_avisc_entropy Entropy viscosity constant
  !! @param n_S Normalization factor
  !! @param n Number of points
  subroutine entropy_viscosity_compute_viscosity_cpu(reg_coeff, &
       entropy_residual, h, c_avisc_entropy, n_S, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: reg_coeff
    real(kind=rp), dimension(n), intent(in) :: entropy_residual, h
    real(kind=rp), intent(in) :: c_avisc_entropy, n_S
    integer :: i

    do concurrent (i = 1:n)
       reg_coeff(i) = c_avisc_entropy * h(i) * h(i) * entropy_residual(i) / n_S
    end do

  end subroutine entropy_viscosity_compute_viscosity_cpu

  !> Apply element-wise maximum on CPU
  !! @param reg_coeff Regularization coefficient field (modified in-place)
  !! @param lx Polynomial order
  !! @param nelv Number of elements
  subroutine entropy_viscosity_apply_element_max_cpu(reg_coeff, lx, nelv)
    integer, intent(in) :: lx, nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: reg_coeff
    integer :: i, j, k, el
    real(kind=rp) :: max_visc_el

    do el = 1, nelv
       max_visc_el = 0.0_rp

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                max_visc_el = max(max_visc_el, reg_coeff(i,j,k,el))
             end do
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                reg_coeff(i,j,k,el) = max_visc_el
             end do
          end do
       end do
    end do

  end subroutine entropy_viscosity_apply_element_max_cpu

  !> Clamp regularization coefficient to low-order viscosity on CPU
  !! @param reg_coeff Regularization coefficient field (modified in-place)
  !! @param h Mesh size field
  !! @param max_wave_speed Maximum wave speed field
  !! @param c_avisc_low Low-order viscosity constant
  !! @param n Number of points
  subroutine entropy_viscosity_clamp_to_low_order_cpu(reg_coeff, &
       h, max_wave_speed, c_avisc_low, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: reg_coeff
    real(kind=rp), dimension(n), intent(in) :: h, max_wave_speed
    real(kind=rp), intent(in) :: c_avisc_low
    integer :: i
    real(kind=rp) :: low_order_visc

    do concurrent (i = 1:n)
       low_order_visc = c_avisc_low * h(i) * max_wave_speed(i)
       reg_coeff(i) = min(reg_coeff(i), low_order_visc)
    end do

  end subroutine entropy_viscosity_clamp_to_low_order_cpu

  !> Divide by multiplicity for smoothing on CPU
  !! @param reg_coeff Regularization coefficient field (modified in-place)
  !! @param temp_field Temporary field with summed values
  !! @param mult_field Multiplicity field
  !! @param n Number of points
  subroutine entropy_viscosity_smooth_divide_cpu(reg_coeff, &
       temp_field, mult_field, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(out) :: reg_coeff
    real(kind=rp), dimension(n), intent(in) :: temp_field, mult_field
    integer :: i

    do concurrent (i = 1:n)
       reg_coeff(i) = temp_field(i) / mult_field(i)
    end do

  end subroutine entropy_viscosity_smooth_divide_cpu

end module entropy_viscosity_cpu

