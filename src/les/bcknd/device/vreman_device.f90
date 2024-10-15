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
  use field_list, only : field_list_t
  use device_math, only : device_cadd, device_col2, device_col3, &
                          device_addcol3, device_subcol3, device_add4, &
                          device_invcol2, device_vecsqrt1, device_cmult, &
                          device_rmneg
  use math, only : NEKO_EPS
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  implicit none
  private

  public :: vreman_compute_device

contains

  !> Compute eddy viscosity on the device.
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c The Vreman model constant
  subroutine vreman_compute_device(t, tstep, coef, nut, delta, c)
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
    integer :: temp_indices(17)
    integer :: e, i

    u => neko_field_registry%get_field_by_name("u")
    v => neko_field_registry%get_field_by_name("v")
    w => neko_field_registry%get_field_by_name("w")

    call neko_scratch_registry%request_field(a11, temp_indices(1))
    call neko_scratch_registry%request_field(a12, temp_indices(2))
    call neko_scratch_registry%request_field(a13, temp_indices(3))
    call neko_scratch_registry%request_field(a21, temp_indices(4))
    call neko_scratch_registry%request_field(a22, temp_indices(5))
    call neko_scratch_registry%request_field(a23, temp_indices(6))
    call neko_scratch_registry%request_field(a31, temp_indices(7))
    call neko_scratch_registry%request_field(a32, temp_indices(8))
    call neko_scratch_registry%request_field(a33, temp_indices(9))
    call neko_scratch_registry%request_field(beta11, temp_indices(10))
    call neko_scratch_registry%request_field(beta12, temp_indices(11))
    call neko_scratch_registry%request_field(beta13, temp_indices(12))
    call neko_scratch_registry%request_field(beta22, temp_indices(13))
    call neko_scratch_registry%request_field(beta23, temp_indices(14))
    call neko_scratch_registry%request_field(beta33, temp_indices(15))
    call neko_scratch_registry%request_field(b_beta, temp_indices(16))
    call neko_scratch_registry%request_field(aijaij, temp_indices(17))


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

    call coef%gs_h%op(a11%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a12%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a13%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a21%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a22%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a23%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a31%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a32%x, a11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(a33%x, a11%dof%size(), GS_OP_ADD)
    
    call device_col2(a11%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a12%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a13%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a21%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a22%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a23%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a31%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a32%x_d, coef%mult_d, a11%dof%size())
    call device_col2(a33%x_d, coef%mult_d, a11%dof%size())

    ! beta_ij = alpha_mi alpha_mj
    call device_col3(beta11%x_d, a11%x_d, a11%x_d, a11%dof%size())
    call device_addcol3(beta11%x_d, a21%x_d, a21%x_d, a11%dof%size())
    call device_addcol3(beta11%x_d, a31%x_d, a31%x_d, a11%dof%size())

    call device_col3(beta22%x_d, a12%x_d, a12%x_d, a11%dof%size())
    call device_addcol3(beta22%x_d, a22%x_d, a22%x_d, a11%dof%size())
    call device_addcol3(beta22%x_d, a32%x_d, a32%x_d, a11%dof%size())

    call device_col3(beta33%x_d, a13%x_d, a13%x_d, a11%dof%size())
    call device_addcol3(beta33%x_d, a23%x_d, a23%x_d, a11%dof%size())
    call device_addcol3(beta33%x_d, a33%x_d, a33%x_d, a11%dof%size())

    call device_col3(beta12%x_d, a11%x_d, a12%x_d, a11%dof%size())
    call device_addcol3(beta12%x_d, a21%x_d, a22%x_d, a11%dof%size())
    call device_addcol3(beta12%x_d, a31%x_d, a32%x_d, a11%dof%size())

    call device_col3(beta13%x_d, a11%x_d, a13%x_d, a11%dof%size())
    call device_addcol3(beta13%x_d, a21%x_d, a23%x_d, a11%dof%size())
    call device_addcol3(beta13%x_d, a31%x_d, a33%x_d, a11%dof%size())

    call device_col3(beta23%x_d, a12%x_d, a13%x_d, a11%dof%size())
    call device_addcol3(beta23%x_d, a22%x_d, a23%x_d, a11%dof%size())
    call device_addcol3(beta23%x_d, a32%x_d, a33%x_d, a11%dof%size())

    call device_col3(b_beta%x_d, beta11%x_d, beta22%x_d, a11%dof%size())
    call device_subcol3(b_beta%x_d, beta12%x_d, beta12%x_d, a11%dof%size())
    call device_addcol3(b_beta%x_d, beta11%x_d, beta33%x_d, a11%dof%size())
    call device_subcol3(b_beta%x_d, beta13%x_d, beta13%x_d, a11%dof%size())
    call device_addcol3(b_beta%x_d, beta22%x_d, beta33%x_d, a11%dof%size())
    call device_subcol3(b_beta%x_d, beta23%x_d, beta23%x_d, a11%dof%size())

    ! Zero the negative element of b_beta
    call device_rmneg(b_beta%x_d, a11%dof%size())
    
    ! alpha_ij alpha_ij
    call device_add4(aijaij%x_d, beta11%x_d, beta22%x_d, &
                     beta33%x_d, a11%dof%size())
    
    ! Get the eddy viscosity
    call device_cadd(aijaij%x_d, NEKO_EPS, a11%dof%size())
    call device_invcol2(b_beta%x_d, aijaij%x_d, a11%dof%size())
    call device_vecsqrt1(b_beta%x_d, a11%dof%size())
    call device_col3(nut%x_d, delta%x_d, b_beta%x_d, a11%dof%size())
    call device_col2(nut%x_d, delta%x_d, a11%dof%size())
    call device_cmult(nut%x_d, c, a11%dof%size())

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine vreman_compute_device

end module vreman_device

