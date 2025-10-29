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
!> Implements the device kernel for the `wale_t` type.
module wale_device
  use num_types, only : rp
  use math, only : NEKO_EPS
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use device_math, only : device_col2
  use device_wale_nut, only : device_wale_nut_compute
  implicit none
  private

  public :: wale_compute_device

contains

  !> Compute eddy viscosity on the device.
  !! @param if_ext If extrapolate the velocity field to evaluate
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c The Vreman model constant
  subroutine wale_compute_device(if_ext, t, tstep, coef, nut, delta, c)
    logical, intent(in) :: if_ext
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c
    ! This is the alpha tensor in the paper
    type(field_t), pointer :: g11, g12, g13, g21, g22, g23, g31, g32, g33
    type(field_t), pointer :: u, v, w

    integer :: temp_indices(9)
    integer :: e, i

    if (if_ext .eqv. .true.) then
       u => neko_field_registry%get_field_by_name("u_e")
       v => neko_field_registry%get_field_by_name("v_e")
       w => neko_field_registry%get_field_by_name("w_e")
    else
       u => neko_field_registry%get_field_by_name("u")
       v => neko_field_registry%get_field_by_name("v")
       w => neko_field_registry%get_field_by_name("w")
    end if

    call neko_scratch_registry%request_field(g11, temp_indices(1))
    call neko_scratch_registry%request_field(g12, temp_indices(2))
    call neko_scratch_registry%request_field(g13, temp_indices(3))
    call neko_scratch_registry%request_field(g21, temp_indices(4))
    call neko_scratch_registry%request_field(g22, temp_indices(5))
    call neko_scratch_registry%request_field(g23, temp_indices(6))
    call neko_scratch_registry%request_field(g31, temp_indices(7))
    call neko_scratch_registry%request_field(g32, temp_indices(8))
    call neko_scratch_registry%request_field(g33, temp_indices(9))


    ! Compute the derivatives of the velocity (the alpha tensor)
    call dudxyz (g11%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (g12%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (g13%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz (g21%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (g22%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (g23%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz (g31%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (g32%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (g33%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call coef%gs_h%op(g11, GS_OP_ADD)
    call coef%gs_h%op(g12, GS_OP_ADD)
    call coef%gs_h%op(g13, GS_OP_ADD)
    call coef%gs_h%op(g21, GS_OP_ADD)
    call coef%gs_h%op(g22, GS_OP_ADD)
    call coef%gs_h%op(g23, GS_OP_ADD)
    call coef%gs_h%op(g31, GS_OP_ADD)
    call coef%gs_h%op(g32, GS_OP_ADD)
    call coef%gs_h%op(g33, GS_OP_ADD)

    call device_wale_nut_compute(g11%x_d, g12%x_d, g13%x_d, &
         g21%x_d, g22%x_d, g23%x_d, &
         g31%x_d, g32%x_d, g33%x_d, &
         delta%x_d, nut%x_d, coef%mult_d, &
         c, NEKO_EPS, g11%dof%size())

    call coef%gs_h%op(nut, GS_OP_ADD)
    call device_col2(nut%x_d, coef%mult_d, nut%dof%size())

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine wale_compute_device

end module wale_device

