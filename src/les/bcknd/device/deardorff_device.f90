! Copyright (c) 2025-2026, The Neko Authors
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
!> Implements the device kernel for the `deardorff_t` type.
module deardorff_device
  use num_types, only : rp
  use utils, only : neko_error
  use math, only : NEKO_EPS
  use scratch_registry, only : neko_scratch_registry
  use registry, only : neko_registry
  use field, only : field_t
  use operators, only : grad
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use device_math, only : device_col2
  use device_deardorff_nut, only : device_deardorff_nut_compute
  implicit none
  private

  public :: deardorff_compute_device

contains

  !> Compute eddy viscosity on the device.
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param temperature_field_name The name of the temperature field.
  !! @param TKE_field_name The name of the TKE field.
  !! @param nut The eddy viscosity field.
  !! @param temperature_alphat The eddy diffusivity field for temperature.
  !! @param TKE_alphat The eddy diffusivity field for TKE.
  !! @param TKE_source The source terms for TKE equation.
  !! @param delta The LES lengthscale.
  !! @param c_k The deardorff model constant
  !! @param T0 The reference temperature.
  !! @param g The gravitational acceleration vector.
  subroutine deardorff_compute_device(t, tstep, coef, &
       temperature_field_name, TKE_field_name, &
       nut, temperature_alphat, &
       TKE_alphat, TKE_source, &
       delta, c_k, T0, g)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: temperature_field_name
    character(len=*), intent(in) :: TKE_field_name
    type(field_t), intent(inout) :: nut, temperature_alphat
    type(field_t), intent(inout) :: TKE_alphat, TKE_source
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c_k, T0, g(3)

    type(field_t), pointer :: TKE, temperature
    type(field_t), pointer :: dTdx, dTdy, dTdz
    type(field_t), pointer :: u, v, w
    type(field_t), pointer :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    integer :: temp_indices(12)

    TKE => neko_registry%get_field_by_name(TKE_field_name)
    temperature => neko_registry%get_field_by_name(temperature_field_name)

    u => neko_registry%get_field_by_name("u")
    v => neko_registry%get_field_by_name("v")
    w => neko_registry%get_field_by_name("w")

    call neko_scratch_registry%request_field(dTdx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(dTdy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(dTdz, temp_indices(3), .false.)

    ! Calculate vertical temperature gradients
    call grad(dTdx%x, dTdy%x, dTdz%x, temperature%x, coef)

    call coef%gs_h%op(dTdx, GS_OP_ADD)
    call coef%gs_h%op(dTdy, GS_OP_ADD)
    call coef%gs_h%op(dTdz, GS_OP_ADD)
    call device_col2(dTdx%x_d, coef%mult_d, nut%dof%size())
    call device_col2(dTdy%x_d, coef%mult_d, nut%dof%size())
    call device_col2(dTdz%x_d, coef%mult_d, nut%dof%size())

    ! Compute velocity gradients
    call neko_scratch_registry%request_field(a11, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(a12, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(a13, temp_indices(6), .false.)
    call neko_scratch_registry%request_field(a21, temp_indices(7), .false.)
    call neko_scratch_registry%request_field(a22, temp_indices(8), .false.)
    call neko_scratch_registry%request_field(a23, temp_indices(9), .false.)
    call neko_scratch_registry%request_field(a31, temp_indices(10), .false.)
    call neko_scratch_registry%request_field(a32, temp_indices(11), .false.)
    call neko_scratch_registry%request_field(a33, temp_indices(12), .false.)

    call grad(a11%x, a12%x, a13%x, u%x, coef)
    call grad(a21%x, a22%x, a23%x, v%x, coef)
    call grad(a31%x, a32%x, a33%x, w%x, coef)

    call coef%gs_h%op(a11, GS_OP_ADD)
    call coef%gs_h%op(a12, GS_OP_ADD)
    call coef%gs_h%op(a13, GS_OP_ADD)
    call coef%gs_h%op(a21, GS_OP_ADD)
    call coef%gs_h%op(a22, GS_OP_ADD)
    call coef%gs_h%op(a23, GS_OP_ADD)
    call coef%gs_h%op(a31, GS_OP_ADD)
    call coef%gs_h%op(a32, GS_OP_ADD)
    call coef%gs_h%op(a33, GS_OP_ADD)

    call device_col2(a11%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a12%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a13%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a21%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a22%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a23%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a31%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a32%x_d, coef%mult_d, nut%dof%size())
    call device_col2(a33%x_d, coef%mult_d, nut%dof%size())

    call device_deardorff_nut_compute(TKE%x_d, &
         dTdx%x_d, dTdy%x_d, dTdz%x_d, &
         a11%x_d, a12%x_d, a13%x_d, &
         a21%x_d, a22%x_d, a23%x_d, &
         a31%x_d, a32%x_d, a33%x_d, &
         delta%x_d, nut%x_d, temperature_alphat%x_d, &
         TKE_alphat%x_d, TKE_source%x_d, &
         c_k, T0, g, NEKO_EPS, a11%dof%size())

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine deardorff_compute_device

end module deardorff_device

