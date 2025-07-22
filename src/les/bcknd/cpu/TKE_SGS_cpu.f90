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
!> Implements the CPU kernel for the `TKE_SGS_t` type.
module TKE_SGS_cpu
  use utils, only : neko_error
  use num_types, only : rp
  use field_list, only : field_list_t
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : strain_rate
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use math, only : col2
  implicit none
  private

  public :: TKE_SGS_compute_cpu

contains

  !> Compute eddy viscosity on the CPU.
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c_k The TKE_SGS model constant
  subroutine TKE_SGS_compute_cpu(t, tstep, coef, nut, delta, c_k, &
                                 T0, g, vert_dir)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c_k, T0, g
    character(len=*), intent(in) :: vert_dir
    type(field_t), pointer :: TKE, T, dTdz
    integer :: temp_indices
    real(kind=rp) :: l, N2
    integer :: e, i

    TKE => neko_field_registry%get_field_by_name("TKE")
    T => neko_field_registry%get_field_by_name("Temperature")

    call neko_scratch_registry%request_field(dTdz, temp_indices)

    select case (vert_dir)
    case ("x")
       call dudxyz(dTdz%x, T%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    case ("y")
       call dudxyz(dTdz%x, T%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    case ("z")
       call dudxyz(dTdz%x, T%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    case default
       call neko_error("Invalid specified vertical direction.")
    end select

    call coef%gs_h%op(dTdz, GS_OP_ADD)
    call col2(dTdz%x, coef%mult, nut%dof%size())

    do concurrent (i = 1:coef%dof%size())
       N2 = dTdz%x(i,1,1,1) * g / T0
       if (N2 .gt. 0.0_rp) then
          l = 0.76_rp * sqrt(TKE%x(i,1,1,1) / abs(N2))
          l = min(l, delta%x(i,1,1,1))
       else
          l = delta%x(i,1,1,1)
       end if

       nut%x(i,1,1,e) = c_k * l * sqrt(TKE%x(i,1,1,1))
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine TKE_SGS_compute_cpu

end module TKE_SGS_cpu

