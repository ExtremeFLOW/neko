! Copyright (c) 2024, The Neko Authors
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
!> A module containing filter functions and subroutines. These functions
!! are used to modify fields in a way that is useful for various
!! simulations.
module filters
  use field, only: field_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use num_types, only: rp
  use utils, only: neko_error
  implicit none

  private
  public :: smooth_step_field, permeability_field, step_function_field

contains

  !> @brief Apply a smooth step function to a field.
  !! @details The smooth step function is defined as:
  !! \f[
  !! t = (x - edge0) / (edge1 - edge0)
  !!  f(t) = \begin{cases}
  !!            t^3 (t (6x - 15) + 10), & t \in [0, 1] \\
  !!              0, & t \leq 0 \\
  !!              1, & t \geq 1 \\
  !!          \end{cases}
  !! \f]
  !! @note The step can be inverted by swapping edge0 and edge1.
  !!
  !! @param[in,out] F Field to be modified.
  !! @param[in] edge0 Edge giving output 0.
  !! @param[in] edge1 Edge giving output 1.
  subroutine smooth_step_field(F, edge0, edge1)
    use filters_cpu, only: smooth_step_cpu

    type(field_t), intent(inout) :: F
    real(kind=rp), intent(in) :: edge0, edge1

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('smooth_step_field: not implemented for device')
    else
       F%x = smooth_step_cpu(F%x, edge0, edge1)
    end if
  end subroutine smooth_step_field

  !> @brief Apply a permeability function to a field.
  !! @details The permeability function is defined as:
  !! \f[ k(x) = k_0 + (k_1 - k_0) x \frac{q + 1}{q + x}} \f]
  !! @param[in,out] F Field to be modified.
  !! @param[in] k_0 Permeability at x=0.
  !! @param[in] k_1 Permeability at x=1.
  !! @param[in] q Penalty factor.
  subroutine permeability_field(F_out, x, k_0, k_1, q)
    use filters_cpu, only: permeability_cpu

    type(field_t), intent(inout) :: F_out
    type(field_t), intent(in) :: x
    real(kind=rp), intent(in) :: k_0, k_1
    real(kind=rp), intent(in) :: q

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('permeability_field: not implemented for device')
    else
       F_out%x = permeability_cpu(x%x, k_0, k_1, q)
    end if
  end subroutine permeability_field

  !> @brief Apply a step function to a field.
  !! @param[in,out] F Field to be modified.
  !! @param[in] x0 Position of the step.
  !! @param[in] value0 Value of the field before the step.
  !! @param[in] value1 Value of the field after the step.
  subroutine step_function_field(F, x0, value0, value1)
    use filters_cpu, only: step_function_cpu

    type(field_t), intent(inout) :: F
    real(kind=rp), intent(in) :: x0, value0, value1

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('step_function_field: not implemented for device')
    else
       F%x = step_function_cpu(F%x, x0, value0, value1)
    end if
  end subroutine step_function_field

end module filters
