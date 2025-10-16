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
!> Implements the device kernel for the `centrifugal_source_term_t` type.

module centrifugal_source_term_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use dofmap, only : dofmap_t
  use device_math, only : device_add3s2, device_add2, device_cadd2
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  public :: centrifugal_source_term_compute_device

contains

  !> Computes the centrifugal source term on the device.
  !! @param omega The rotation vector.
  !! @param ref_point The reference point on the rotation axis.
  !! @param fields The right-hand side, which should be the velocity components.
  subroutine centrifugal_source_term_compute_device(omega, ref_point, fields)
    real(kind=rp), intent(in) :: omega(3)
    real(kind=rp), intent(in) :: ref_point(3)
    type(field_list_t), intent(inout) :: fields
    type(dofmap_t), pointer :: dof
    integer :: n
    type(field_t), pointer :: fu, fv, fw
    integer :: tmp_index(6)
    type(field_t), pointer :: tmp_rx, tmp_ry, tmp_rz
    type(field_t), pointer :: tmp_cx, tmp_cy, tmp_cz

    dof => fields%dof(1)
    n = fields%item_size(1)

    fu => fields%get_by_index(1)
    fv => fields%get_by_index(2)
    fw => fields%get_by_index(3)

    call neko_scratch_registry%request_field(tmp_rx, tmp_index(1))
    call neko_scratch_registry%request_field(tmp_ry, tmp_index(2))
    call neko_scratch_registry%request_field(tmp_rz, tmp_index(3))
    call neko_scratch_registry%request_field(tmp_cx, tmp_index(4))
    call neko_scratch_registry%request_field(tmp_cy, tmp_index(5))
    call neko_scratch_registry%request_field(tmp_cz, tmp_index(6))

    ! displacement with respect to reference point
    call device_cadd2(tmp_rx%x_d, dof%x_d, -ref_point(1), n)
    call device_cadd2(tmp_ry%x_d, dof%y_d, -ref_point(2), n)
    call device_cadd2(tmp_rz%x_d, dof%z_d, -ref_point(3), n)

    ! Omega x r
    call device_add3s2(tmp_cx%x_d, tmp_rz%x_d, tmp_ry%x_d, omega(2), &
         -omega(3), n)
    call device_add3s2(tmp_cy%x_d, tmp_rx%x_d, tmp_rz%x_d, omega(3), &
         -omega(1), n)
    call device_add3s2(tmp_cz%x_d, tmp_ry%x_d, tmp_rx%x_d, omega(1), &
         -omega(2), n)

    ! - Omega x (Omega x r)
    call device_add3s2(tmp_rx%x_d, tmp_cz%x_d, tmp_cy%x_d, -omega(2), &
         omega(3), n)
    call device_add3s2(tmp_ry%x_d, tmp_cx%x_d, tmp_cz%x_d, -omega(3), &
         omega(1), n)
    call device_add3s2(tmp_rz%x_d, tmp_cy%x_d, tmp_cx%x_d, -omega(1), &
         omega(2), n)

    ! Add the centrifugal term to the RHS
    call device_add2(fu%x_d, tmp_rx%x_d, n)
    call device_add2(fv%x_d, tmp_ry%x_d, n)
    call device_add2(fw%x_d, tmp_rz%x_d, n)

    call neko_scratch_registry%relinquish_field(tmp_index)

  end subroutine centrifugal_source_term_compute_device

end module centrifugal_source_term_device
