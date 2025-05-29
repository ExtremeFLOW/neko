! Copyright (c) 2023-2025, The Neko Authors
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
!> Implements the CPU kernel for the `wale_t` type.
module wale_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz, strain_rate
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  use math, only : col2, NEKO_EPS
  implicit none
  private

  public :: wale_compute_cpu

contains

  !> Compute eddy viscosity on the CPU.
  !! @param if_ext If extrapolate the velocity field to evaluate
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c_w The wale model constant
  subroutine wale_compute_cpu(if_ext, t, tstep, coef, nut, delta, c_w)
    logical, intent(in) :: if_ext
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c_w
    type(field_t), pointer :: u, v, w
    ! the velocity gradient tensor
    type(field_t), pointer :: g11, g12, g13, g21, g22, g23, g31, g32, g33
    ! strain rate tensor
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23

    real(kind=rp) :: gsqr_11, gsqr_12, gsqr_13, gsqr_21, gsqr_22, gsqr_23, gsqr_31, gsqr_32, gsqr_33
    real(kind=rp) :: sd11, sd22, sd33, sd12, sd13, sd23
    real(kind=rp) :: Sdij_Sdij
    real(kind=rp) :: Sij_Sij
    real(kind=rp) :: OP_wale
    integer :: temp_indices(15)
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

    call neko_scratch_registry%request_field(g11, temp_indices(1))
    call neko_scratch_registry%request_field(g12, temp_indices(2))
    call neko_scratch_registry%request_field(g13, temp_indices(3))
    call neko_scratch_registry%request_field(g21, temp_indices(4))
    call neko_scratch_registry%request_field(g22, temp_indices(5))
    call neko_scratch_registry%request_field(g23, temp_indices(6))
    call neko_scratch_registry%request_field(g31, temp_indices(7))
    call neko_scratch_registry%request_field(g32, temp_indices(8))
    call neko_scratch_registry%request_field(g33, temp_indices(9))

    call neko_scratch_registry%request_field(s11, temp_indices(10))
    call neko_scratch_registry%request_field(s22, temp_indices(11))
    call neko_scratch_registry%request_field(s33, temp_indices(12))
    call neko_scratch_registry%request_field(s12, temp_indices(13))
    call neko_scratch_registry%request_field(s13, temp_indices(14))
    call neko_scratch_registry%request_field(s23, temp_indices(15))


    ! Compute the velocity gradient tensor g
    call dudxyz(g11%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(g12%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(g13%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz(g21%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(g22%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(g23%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz(g31%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(g32%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(g33%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call coef%gs_h%op(g11, GS_OP_ADD)
    call coef%gs_h%op(g12, GS_OP_ADD)
    call coef%gs_h%op(g13, GS_OP_ADD)
    call coef%gs_h%op(g21, GS_OP_ADD)
    call coef%gs_h%op(g22, GS_OP_ADD)
    call coef%gs_h%op(g23, GS_OP_ADD)
    call coef%gs_h%op(g31, GS_OP_ADD)
    call coef%gs_h%op(g32, GS_OP_ADD)
    call coef%gs_h%op(g33, GS_OP_ADD)

    ! Compute components of the strain rate tensor
    call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, u, v, w, coef)

    call coef%gs_h%op(s11, GS_OP_ADD)
    call coef%gs_h%op(s22, GS_OP_ADD)
    call coef%gs_h%op(s33, GS_OP_ADD)
    call coef%gs_h%op(s12, GS_OP_ADD)
    call coef%gs_h%op(s13, GS_OP_ADD)
    call coef%gs_h%op(s23, GS_OP_ADD)


    do concurrent(e = 1:coef%msh%nelv)
       do concurrent(i = 1:coef%Xh%lxyz)
          ! gij^2 = g_ik * g_kj
          gsqr_11 = g11%x(i,1,1,e)*g11%x(i,1,1,e) + &
                g12%x(i,1,1,e)*g21%x(i,1,1,e) + &
                g13%x(i,1,1,e)*g31%x(i,1,1,e)
          gsqr_12 = g11%x(i,1,1,e)*g12%x(i,1,1,e) + &
                g12%x(i,1,1,e)*g22%x(i,1,1,e) + &
                g13%x(i,1,1,e)*g32%x(i,1,1,e)
          gsqr_13 = g11%x(i,1,1,e)*g13%x(i,1,1,e) + &
                g12%x(i,1,1,e)*g23%x(i,1,1,e) + &
                g13%x(i,1,1,e)*g33%x(i,1,1,e)
          gsqr_21 = g21%x(i,1,1,e)*g11%x(i,1,1,e) + &
                g22%x(i,1,1,e)*g21%x(i,1,1,e) + &
                g23%x(i,1,1,e)*g31%x(i,1,1,e)
          gsqr_22 = g21%x(i,1,1,e)*g12%x(i,1,1,e) + &
                g22%x(i,1,1,e)*g22%x(i,1,1,e) + &
                g23%x(i,1,1,e)*g32%x(i,1,1,e)
          gsqr_23 = g21%x(i,1,1,e)*g13%x(i,1,1,e) + &
                g22%x(i,1,1,e)*g23%x(i,1,1,e) + &
                g23%x(i,1,1,e)*g33%x(i,1,1,e)
          gsqr_31 = g31%x(i,1,1,e)*g11%x(i,1,1,e) + &
                g32%x(i,1,1,e)*g21%x(i,1,1,e) + &
                g33%x(i,1,1,e)*g31%x(i,1,1,e)
          gsqr_32 = g31%x(i,1,1,e)*g12%x(i,1,1,e) + &
                g32%x(i,1,1,e)*g22%x(i,1,1,e) + &
                g33%x(i,1,1,e)*g32%x(i,1,1,e)
          gsqr_33 = g31%x(i,1,1,e)*g13%x(i,1,1,e) + &
                g32%x(i,1,1,e)*g23%x(i,1,1,e) + &
                g33%x(i,1,1,e)*g33%x(i,1,1,e)

          ! sdij components
          sd11 = gsqr_11 - ( (gsqr_11 + gsqr_22 + gsqr_33) / 3.0_rp)
          sd22 = gsqr_22 - ( (gsqr_11 + gsqr_22 + gsqr_33) / 3.0_rp)
          sd33 = gsqr_33 - ( (gsqr_11 + gsqr_22 + gsqr_33) / 3.0_rp)
          sd12 = 0.5_rp * (gsqr_12 + gsqr_21)
          sd13 = 0.5_rp * (gsqr_13 + gsqr_31)
          sd23 = 0.5_rp * (gsqr_23 + gsqr_32)

          ! Sdij*Sdij
          Sdij_Sdij = sd11*sd11 + sd22*sd22 + sd33*sd33 + &
                            2.0_rp * (sd12*sd12 + sd13*sd13 + sd23*sd23)
          ! Sij*Sij
          Sij_Sij = s11%x(i,1,1,e)*s11%x(i,1,1,e) + s22%x(i,1,1,e)*s22%x(i,1,1,e) + &
                    s33%x(i,1,1,e)*s33%x(i,1,1,e) + 2.0_rp * (s12%x(i,1,1,e)*s12%x(i,1,1,e) + &
                    s13%x(i,1,1,e)*s13%x(i,1,1,e) + s23%x(i,1,1,e)*s23%x(i,1,1,e))

          ! Wale operator
          OP_wale = Sdij_Sdij**(3.0_rp / 2.0_rp) / &
                    max((Sij_Sij**(5.0_rp / 2.0_rp) + Sdij_Sdij**(5.0_rp / 4.0_rp)), NEKO_EPS)

          ! turbulent viscosity
          nut%x(i,1,1,e) = c_w**2 * delta%x(i,1,1,e)**2 * OP_wale * coef%mult(i,1,1,e)
       end do
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine wale_compute_cpu

end module wale_cpu

