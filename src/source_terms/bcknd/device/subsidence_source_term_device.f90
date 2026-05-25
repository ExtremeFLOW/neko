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
!> Implements the device interface for the `subsidence_source_term_t` type.
module subsidence_source_term_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use scratch_registry, only : neko_scratch_registry
  use coefs, only : coef_t
  use operators, only : grad
  use gs_ops, only : GS_OP_ADD
  use device_math, only : device_copy, device_cmult, device_add2s2, device_addcol3, &
       device_add2, device_col2
  implicit none
  private

  public :: subsidence_source_term_compute_device

contains

  !> Prepares device arrays and performs device operations to add the subsidence
  !! source: fu += w_sub * (dir . grad(u)), fv += w_sub * (dir . grad(v)).
  !! This routine computes gradients on the host (via `grad`), performs the
  !! gather-scatter and normalization, moves the gradient arrays to device
  !! scaling (via `device_col2`), and uses device elementwise ops to form and
  !! add the products. The low-level device kernel is intentionally not
  !! implemented here.
  subroutine subsidence_source_term_compute_device(u, v, fields, coef, w_sub, vertical_dir)
    type(field_t), intent(in) :: u, v, w_sub
    type(field_list_t), intent(inout) :: fields
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: vertical_dir(3)

    integer :: n
    integer :: tmp_idx(5)
    type(field_t), pointer :: fu, fv
    type(field_t), pointer :: dx, dy, dz
    type(field_t), pointer :: deriv, prod

    ! RHS fields
    fu => fields%get_by_index(1)
    fv => fields%get_by_index(2)

    n = fu%size()

    ! Request scratch fields: dx,dy,dz and a temporary for device work
    call neko_scratch_registry%request_field(dx, tmp_idx(1), .false.)
    call neko_scratch_registry%request_field(dy, tmp_idx(2), .false.)
    call neko_scratch_registry%request_field(dz, tmp_idx(3), .false.)
    call neko_scratch_registry%request_field(deriv, tmp_idx(4), .false.)

    ! Compute host gradients for u
    call grad(dx%x, dy%x, dz%x, u%x, coef)

    ! Gather-scatter and normalize on host
    call coef%gs_h%op(dx, GS_OP_ADD)
    call coef%gs_h%op(dy, GS_OP_ADD)
    call coef%gs_h%op(dz, GS_OP_ADD)

    call device_col2(dx%x_d, coef%mult_d, n)
    call device_col2(dy%x_d, coef%mult_d, n)
    call device_col2(dz%x_d, coef%mult_d, n)

    ! deriv = vertical_dir(1)*dx + vertical_dir(2)*dy + vertical_dir(3)*dz
    call device_copy(deriv%x_d, dx%x_d, n)
    call device_cmult(deriv%x_d, vertical_dir(1), n)
    call device_add2s2(deriv%x_d, dy%x_d, vertical_dir(2), n)
    call device_add2s2(deriv%x_d, dz%x_d, vertical_dir(3), n)

    ! prod = w_sub * deriv
    call neko_scratch_registry%request_field(prod, tmp_idx(5), .false.)
    call device_cmult(prod%x_d, 0.0_rp, n)
    call device_addcol3(prod%x_d, w_sub%x_d, deriv%x_d, n)

    ! Add to RHS: fu += prod
    call device_add2(fu%x_d, prod%x_d, n)

    ! Compute host gradients for v
    call grad(dx%x, dy%x, dz%x, v%x, coef)

    ! Gather-scatter and normalize on host
    call coef%gs_h%op(dx, GS_OP_ADD)
    call coef%gs_h%op(dy, GS_OP_ADD)
    call coef%gs_h%op(dz, GS_OP_ADD)

    call device_col2(dx%x_d, coef%mult_d, n)
    call device_col2(dy%x_d, coef%mult_d, n)
    call device_col2(dz%x_d, coef%mult_d, n)

    ! deriv = vertical_dir(1)*dx + vertical_dir(2)*dy + vertical_dir(3)*dz
    call device_copy(deriv%x_d, dx%x_d, n)
    call device_cmult(deriv%x_d, vertical_dir(1), n)
    call device_add2s2(deriv%x_d, dy%x_d, vertical_dir(2), n)
    call device_add2s2(deriv%x_d, dz%x_d, vertical_dir(3), n)

    ! prod = w_sub * deriv (reuse same prod scratch)
    call device_cmult(prod%x_d, 0.0_rp, n)
    call device_addcol3(prod%x_d, w_sub%x_d, deriv%x_d, n)

    ! Add to RHS: fv += prod
    call device_add2(fv%x_d, prod%x_d, n)

    ! Release scratch
    call neko_scratch_registry%relinquish_field(tmp_idx)

  end subroutine subsidence_source_term_compute_device

end module subsidence_source_term_device
