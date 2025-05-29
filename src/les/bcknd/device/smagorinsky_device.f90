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
!> Implements the device kernel for the `smagorinsky_t` type.
module smagorinsky_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : strain_rate
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use device_math, only : device_col2
  use device_smagorinsky_nut, only : device_smagorinsky_nut_compute
  implicit none
  private

  public :: smagorinsky_compute_device

contains

  !> Compute eddy viscosity on the device.
  !! @param if_ext If extrapolate the velocity field to evaluate
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c_s The smagorinsky model constant
  subroutine smagorinsky_compute_device(if_ext, t, tstep, coef, nut, delta, c_s)
    logical, intent(in) :: if_ext
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c_s
    type(field_t), pointer :: u, v, w
    ! double of the strain rate tensor
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    integer :: temp_indices(6)
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

    call neko_scratch_registry%request_field(s11, temp_indices(1))
    call neko_scratch_registry%request_field(s22, temp_indices(2))
    call neko_scratch_registry%request_field(s33, temp_indices(3))
    call neko_scratch_registry%request_field(s12, temp_indices(4))
    call neko_scratch_registry%request_field(s13, temp_indices(5))
    call neko_scratch_registry%request_field(s23, temp_indices(6))

    ! Compute the strain rate tensor
    call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, u, v, w, coef)

    call coef%gs_h%op(s11, GS_OP_ADD)
    call coef%gs_h%op(s22, GS_OP_ADD)
    call coef%gs_h%op(s33, GS_OP_ADD)
    call coef%gs_h%op(s12, GS_OP_ADD)
    call coef%gs_h%op(s13, GS_OP_ADD)
    call coef%gs_h%op(s23, GS_OP_ADD)

    call device_smagorinsky_nut_compute(s11%x_d, s22%x_d, s33%x_d, &
         s12%x_d, s13%x_d, s23%x_d, &
         delta%x_d, nut%x_d, coef%mult_d, &
         c_s, s11%dof%size())

    call coef%gs_h%op(nut, GS_OP_ADD)
    call device_col2(nut%x_d, coef%mult_d, nut%dof%size())

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine smagorinsky_compute_device

end module smagorinsky_device

