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
!> Implements the device kernel for the `sponge_source_term_t` type.
module sponge_source_term_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use device_math, only : device_sub3, device_col2, device_add2s2
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  public :: sponge_source_term_compute_device

contains

  !> Computes the sponge source term on the device.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param fields The right-hand side, which should be the velocity components.
  !! @param omega The rotation vector.
  !! @param omega The geostrophic wind.
  subroutine sponge_source_term_compute_device(fields, u, v, w, &
       u_bf, v_bf, w_bf, fringe, a_x, a_y, a_z)
    type(field_list_t), intent(inout) :: fields
    type(field_t), intent(in) :: u, v, w, fringe, u_bf, v_bf, w_bf
    real(kind=rp), intent(in) :: a_x, a_y, a_z
    integer :: n
    type(field_t), pointer :: fu, fv, fw
    integer :: tmp_index
    type(field_t), pointer :: wk

    call neko_scratch_registry%request_field(wk, tmp_index, .false.)

    ! The RHS components
    fu => fields%get_by_index(1)
    fv => fields%get_by_index(2)
    fw => fields%get_by_index(3)

    n = fu%size()

    ! wk = u_bf - u
    call device_sub3(wk%x_d, u_bf%x_d, u%x_d, n)
    ! wk = fringe * wk = fringe * (u_bf - u)
    call device_col2(wk%x_d, fringe%x_d, n)
    ! fu = fu + amplitude(1)*wk = fu + amplitude(1)*fringe*(u_bf - u)
    call device_add2s2(fu%x_d, wk%x_d, a_x, n)

    call device_sub3(wk%x_d, v_bf%x_d, v%x_d, n)
    call device_col2(wk%x_d, fringe%x_d, n)
    call device_add2s2(fv%x_d, wk%x_d, a_y, n)

    call device_sub3(wk%x_d, w_bf%x_d, w%x_d, n)
    call device_col2(wk%x_d, fringe%x_d, n)
    call device_add2s2(fw%x_d, wk%x_d, a_z, n)

    call neko_scratch_registry%relinquish_field(tmp_index)

  end subroutine sponge_source_term_compute_device

end module sponge_source_term_device
