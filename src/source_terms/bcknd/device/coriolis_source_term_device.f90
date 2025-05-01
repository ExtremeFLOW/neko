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
!> Implements the device kernel for the `coriolis_source_term_t` type.
module coriolis_source_term_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use device_math, only : device_copy, device_add3s2, device_add2, device_cadd
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  public :: coriolis_source_term_compute_device

contains

  !> Computes the Coriolis source term on the device.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param fields The right-hand side, which should be the velocity components.
  !! @param omega The rotation vector.
  !! @param omega The geostrophic wind.
  subroutine coriolis_source_term_compute_device(u, v, w, fields, omega, u_geo)
    type(field_t), intent(in) :: u, v, w
    type(field_list_t), intent(inout) :: fields
    real(kind=rp), intent(in) :: omega(3)
    real(kind=rp), intent(in) :: u_geo(3)
    integer :: n
    type(field_t), pointer :: fu, fv, fw
    real(kind=rp) :: ui, vi, wi
    integer :: tmp_index(6)
    type(field_t), pointer :: tmp_u, tmp_v, tmp_w
    type(field_t), pointer :: tmp_fu, tmp_fv, tmp_fw

    call neko_scratch_registry%request_field(tmp_u, tmp_index(1))
    call neko_scratch_registry%request_field(tmp_v, tmp_index(2))
    call neko_scratch_registry%request_field(tmp_w, tmp_index(3))
    call neko_scratch_registry%request_field(tmp_fu, tmp_index(4))
    call neko_scratch_registry%request_field(tmp_fv, tmp_index(5))
    call neko_scratch_registry%request_field(tmp_fw, tmp_index(6))

    ! The RHS components
    fu => fields%get_by_index(1)
    fv => fields%get_by_index(2)
    fw => fields%get_by_index(3)

    n = fu%size()

    ! Copy over velocity arrays for manipulation
    call device_copy(tmp_u%x_d, u%x_d, n)
    call device_copy(tmp_v%x_d, v%x_d, n)
    call device_copy(tmp_w%x_d, w%x_d, n)

    ! velocity minus u_geo
    call device_cadd(tmp_u%x_d, -u_geo(1), n)
    call device_cadd(tmp_v%x_d, -u_geo(2), n)
    call device_cadd(tmp_w%x_d, -u_geo(3), n)

    ! The Coriolis term to be added to the RHS
    call device_add3s2(tmp_fu%x_d, tmp_w%x_d, tmp_v%x_d, -2.0_rp * omega(2), &
         2.0_rp * omega(3), n)
    call device_add3s2(tmp_fv%x_d, tmp_u%x_d, tmp_w%x_d, -2.0_rp * omega(3), &
         2.0_rp * omega(1), n)
    call device_add3s2(tmp_fw%x_d, tmp_v%x_d, tmp_u%x_d, -2.0_rp * omega(1), &
         2.0_rp * omega(2), n)

    ! Add the Coriolis term to the RHS
    call device_add2(fu%x_d, tmp_fu%x_d, n)
    call device_add2(fv%x_d, tmp_fv%x_d, n)
    call device_add2(fw%x_d, tmp_fw%x_d, n)

    call neko_scratch_registry%relinquish_field(tmp_index)

  end subroutine coriolis_source_term_compute_device

end module coriolis_source_term_device
