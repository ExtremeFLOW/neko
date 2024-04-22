! Copyright (c) 2023, The Neko Authors
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
!> Implements the CPU kernel for the `smagorinsky_t` type.
module smagorinsky_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use math, only : cadd, NEKO_EPS
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use math, only : col2
  implicit none
  private

  public :: smagorinsky_compute_cpu

contains

  !> Compute eddy viscosity on the CPU.
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c_s The smagorinsky model constant
  subroutine smagorinsky_compute_cpu(t, tstep, coef, nut, delta, c_s)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c_s
    ! velocity gradient fields
    type(field_t), pointer :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    type(field_t), pointer :: u, v, w
    ! rate-of-strain tensor
    real(kind=rp) :: s11, s22, s33, s12, s13, s23
    real(kind=rp) :: s_abs
    integer :: temp_indices(9)
    integer :: e, i
    
    u => neko_field_registry%get_field_by_name("u")
    v => neko_field_registry%get_field_by_name("v")
    w => neko_field_registry%get_field_by_name("u")

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

    do e=1, coef%msh%nelv
       do i=1, coef%Xh%lxyz
          ! alpha_ij alpha_ij
          s11 = 0.5_rp * (a11%x(i,1,1,e) + a11%x(i,1,1,e))
          s22 = 0.5_rp * (a22%x(i,1,1,e) + a22%x(i,1,1,e))
          s33 = 0.5_rp * (a33%x(i,1,1,e) + a33%x(i,1,1,e))
          s12 = 0.5_rp * (a12%x(i,1,1,e) + a12%x(i,1,1,e))
          s13 = 0.5_rp * (a13%x(i,1,1,e) + a13%x(i,1,1,e))
          s23 = 0.5_rp * (a23%x(i,1,1,e) + a23%x(i,1,1,e))

          s_abs = sqrt(2 * (s11*s11 + s22*s22 + s33*s33) + &
                  4 * (s12*s12 + s13*s13 + s23*s23))

          nut%x(i,1,1,e) = c_s**2 * delta%x(i,1,1,e)**2 * s_abs
       end do
    end do

    call coef%gs_h%op(nut%x, nut%dof%size(), GS_OP_ADD)
    call col2(nut%x, coef%mult, nut%dof%size())

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine smagorinsky_compute_cpu

end module smagorinsky_cpu

