! Copyright (c) 2023-2024, The Neko Authors
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
!> Implements the device kernel for the `vreman_t` type.
module vreman_device
  use num_types, only : rp
  use math, only : NEKO_EPS
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use device_math, only : device_col2
  use device_vreman_nut, only : device_vreman_nut_compute
  implicit none
  private

  public :: vreman_compute_device

contains

  !> Compute eddy viscosity on the device.
  !! @param if_ext If extrapolate the velocity field to evaluate
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c The Vreman model constant
  subroutine vreman_compute_device(if_ext, t, tstep, coef, nut, delta, c)
    logical, intent(in) :: if_ext
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c
    ! This is the alpha tensor in the paper
    type(field_t), pointer :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    type(field_t), pointer :: u, v, w

    type(field_t), pointer :: beta11
    type(field_t), pointer :: beta12
    type(field_t), pointer :: beta13
    type(field_t), pointer :: beta22
    type(field_t), pointer :: beta23
    type(field_t), pointer :: beta33
    type(field_t), pointer :: b_beta
    type(field_t), pointer :: aijaij
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

    call neko_scratch_registry%request_field(a11, temp_indices(1))
    call neko_scratch_registry%request_field(a12, temp_indices(2))
    call neko_scratch_registry%request_field(a13, temp_indices(3))
    call neko_scratch_registry%request_field(a21, temp_indices(4))
    call neko_scratch_registry%request_field(a22, temp_indices(5))
    call neko_scratch_registry%request_field(a23, temp_indices(6))
    call neko_scratch_registry%request_field(a31, temp_indices(7))
    call neko_scratch_registry%request_field(a32, temp_indices(8))
    call neko_scratch_registry%request_field(a33, temp_indices(9))

    ! Compute the derivatives of the velocity (the alpha tensor)
    call dudxyz (a11%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (a12%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (a13%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz (a21%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (a22%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (a23%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz (a31%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (a32%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (a33%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call coef%gs_h%op(a11, GS_OP_ADD)
    call coef%gs_h%op(a12, GS_OP_ADD)
    call coef%gs_h%op(a13, GS_OP_ADD)
    call coef%gs_h%op(a21, GS_OP_ADD)
    call coef%gs_h%op(a22, GS_OP_ADD)
    call coef%gs_h%op(a23, GS_OP_ADD)
    call coef%gs_h%op(a31, GS_OP_ADD)
    call coef%gs_h%op(a32, GS_OP_ADD)
    call coef%gs_h%op(a33, GS_OP_ADD)

    call device_vreman_nut_compute(a11%x_d, a12%x_d, a13%x_d, &
         a21%x_d, a22%x_d, a23%x_d, &
         a31%x_d, a32%x_d, a33%x_d, &
         delta%x_d, nut%x_d, coef%mult_d, &
         c, NEKO_EPS, a11%dof%size())

    call coef%gs_h%op(nut, GS_OP_ADD)
    call device_col2(nut%x_d, coef%mult_d, nut%dof%size())

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine vreman_compute_device

end module vreman_device

