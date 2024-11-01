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
!> Implements the CPU kernel for the `sigma_t` type.
!! Following Nicoud et al. "Using singular values to build a
!! subgrid-scale model for large-eddy simulations"
!! https://doi.org/10.1063/1.3623274

module sigma_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
   use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz
  use coefs, only : coef_t
  use gs_ops, only : GS_OP_ADD
  implicit none
  private

  public :: sigma_compute_cpu

contains

  !> Compute eddy viscosity on the CPU.
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c The Sigma model constant
  subroutine sigma_compute_cpu(t, tstep, coef, nut, delta, c)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: c
    ! This is the velocity gradient tensor
    type(field_t), pointer :: g11, g12, g13, g21, g22, g23, g31, g32, g33
    type(field_t), pointer :: u, v, w

    real(kind=rp) :: sigG11, sigG12, sigG13, sigG22, sigG23, sigG33
    real(kind=rp) :: sigma1, sigma2, sigma3
    real(kind=rp) :: Invariant1, Invariant2, Invariant3
    real(kind=rp) :: alpha1, alpha2, alpha3
    real(kind=rp) :: Dsigma
    real(kind=rp) :: pi_3 = 4.0_rp/3.0_rp*atan(1.0_rp)
    real(kind=rp) :: tmp1
    real(kind=rp) :: eps

    integer :: temp_indices(9)
    integer :: e, i

    ! some constant
    eps = 1.d-14


    ! get fields from registry
    u => neko_field_registry%get_field_by_name("u")
    v => neko_field_registry%get_field_by_name("v")
    w => neko_field_registry%get_field_by_name("u")

    call neko_scratch_registry%request_field(g11, temp_indices(1))
    call neko_scratch_registry%request_field(g12, temp_indices(2))
    call neko_scratch_registry%request_field(g13, temp_indices(3))
    call neko_scratch_registry%request_field(g21, temp_indices(4))
    call neko_scratch_registry%request_field(g22, temp_indices(5))
    call neko_scratch_registry%request_field(g23, temp_indices(6))
    call neko_scratch_registry%request_field(g31, temp_indices(7))
    call neko_scratch_registry%request_field(g32, temp_indices(8))
    call neko_scratch_registry%request_field(g33, temp_indices(9))


    ! Compute the derivatives of the velocity (the components of the g tensor)
    call dudxyz (g11%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (g12%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (g13%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz (g21%x, v%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (g22%x, v%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (g23%x, v%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call dudxyz (g31%x, w%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz (g32%x, w%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz (g33%x, w%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call coef%gs_h%op(g11%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g12%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g13%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g21%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g22%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g23%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g31%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g32%x, g11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(g33%x, g11%dof%size(), GS_OP_ADD)

    do concurrent (i = 1:g11%dof%size())
       g11%x(i,1,1,1) = g11%x(i,1,1,1) * coef%mult(i,1,1,1)
       g12%x(i,1,1,1) = g12%x(i,1,1,1) * coef%mult(i,1,1,1)
       g13%x(i,1,1,1) = g13%x(i,1,1,1) * coef%mult(i,1,1,1)
       g21%x(i,1,1,1) = g21%x(i,1,1,1) * coef%mult(i,1,1,1)
       g22%x(i,1,1,1) = g22%x(i,1,1,1) * coef%mult(i,1,1,1)
       g23%x(i,1,1,1) = g23%x(i,1,1,1) * coef%mult(i,1,1,1)
       g31%x(i,1,1,1) = g31%x(i,1,1,1) * coef%mult(i,1,1,1)
       g32%x(i,1,1,1) = g32%x(i,1,1,1) * coef%mult(i,1,1,1)
       g33%x(i,1,1,1) = g33%x(i,1,1,1) * coef%mult(i,1,1,1)
    end do

    do concurrent (e = 1:coef%msh%nelv)
       do concurrent (i = 1:coef%Xh%lxyz)
          ! G_ij = g^t g = g_mi g_mj
          sigG11 = g11%x(i,1,1,e)**2 + g21%x(i,1,1,e)**2 + g31%x(i,1,1,e)**2
          sigG22 = g12%x(i,1,1,e)**2 + g22%x(i,1,1,e)**2 + g32%x(i,1,1,e)**2
          sigG33 = g13%x(i,1,1,e)**2 + g23%x(i,1,1,e)**2 + g33%x(i,1,1,e)**2
          sigG12 = g11%x(i,1,1,e)*g12%x(i,1,1,e) + &
                   g21%x(i,1,1,e)*g22%x(i,1,1,e) + &
                   g31%x(i,1,1,e)*g32%x(i,1,1,e)
          sigG13 = g11%x(i,1,1,e)*g13%x(i,1,1,e) + &
                   g21%x(i,1,1,e)*g23%x(i,1,1,e) + &
                   g31%x(i,1,1,e)*g33%x(i,1,1,e)
          sigG23 = g12%x(i,1,1,e)*g13%x(i,1,1,e) + &
                   g22%x(i,1,1,e)*g23%x(i,1,1,e) + &
                   g32%x(i,1,1,e)*g33%x(i,1,1,e)

          !        If LAPACK compute eigenvalues of the semi-definite positive matrix G
          !        ..........to be done later on......
          !        ELSE use the analytical method as done in the following

          ! eigenvalues with the analytical method of Hasan et al. (2001)
          ! doi:10.1006/jmre.2001.2400
              if (abs(sigG11) .lt. eps) then
                 sigG11 = 0.0_rp
              end if
              if (abs(sigG12) .lt. eps) then
                 sigG12 = 0.0_rp
              end if
              if (abs(sigG13) .lt. eps) then
                 sigG13 = 0.0_rp
              end if
              if (abs(sigG22) .lt. eps) then
                 sigG22 = 0.0_rp
              end if
              if (abs(sigG23) .lt. eps) then
                 sigG23 = 0.0_rp
              end if
              if (abs(sigG33) .lt. eps) then
                 sigG33 = 0.0_rp
              end if

              if (abs(sigG12*sigG12 + &
                      sigG13*sigG13 + sigG23*sigG23) .lt. eps) then
                 !             G is diagonal
                 ! estimate the singular values according to:
                 sigma1 = sqrt(max(max(max(sigG11, sigG22), sigG33), 0.0_rp))
                 sigma3 = sqrt(max(min(min(sigG11, sigG22), sigG33), 0.0_rp))
                 Invariant1 = sigG11 + sigG22 + sigG33
                 sigma2 = sqrt(abs(Invariant1 - sigma1*sigma1 - sigma3*sigma3))
              else

                 !  estimation of invariants
                 Invariant1 = sigG11 + sigG22 + sigG33
                 Invariant2 = sigG11*sigG22 + sigG11*sigG33 + sigG22*sigG33 - &
                           (sigG12*sigG12 + sigG13*sigG13 + sigG23*sigG23)
                 Invariant3 = sigG11*sigG22*sigG33 + &
                              2.0_rp*sigG12*sigG13*sigG23 - &
                           (sigG11*sigG23*sigG23 + sigG22*sigG13*sigG13 + &
                            sigG33*sigG12*sigG12)

                 ! G is symmetric semi-definite positive matrix:
                 ! the invariants have to be larger-equal zero
                 ! which is obtained via forcing
                 Invariant1 = max(Invariant1, 0.0_rp)
                 Invariant2 = max(Invariant2, 0.0_rp)
                 Invariant3 = max(Invariant3, 0.0_rp)

                 ! compute the following angles from the invariants
                 alpha1 = Invariant1*Invariant1/9.0_rp - Invariant2/3.0_rp

                 ! since alpha1 is always positive (see Hasan et al. (2001))
                 ! forcing is applied
                 alpha1 = max(alpha1, 0.0_rp)

                 alpha2 = Invariant1*Invariant1*Invariant1/27.0_rp - &
                        Invariant1*Invariant2/6.0_rp + Invariant3/2.0_rp

                 ! since acos(alpha2/(alpha1^(3/2)))/3.0_rp only valid for
                 ! alpha2^2 < alpha1^3.0_rp and arccos(x) only valid for -1<=x<=1
                 !  alpha3 is between 0 and pi/3
                 tmp1 = alpha2/(alpha1**(3.0_rp/2.0_rp))

                 if (tmp1 .le. -1.0_rp) then
                    ! alpha3=pi/3 -> cos(alpha3)=0.5
                    ! compute the singular values
                    sigma1 = sqrt(max(Invariant1/3.0_rp + sqrt(alpha1), 0.0_rp))
                    sigma2 = sigma1
                    sigma3 = sqrt(Invariant1/3.0_rp - 2.0_rp*sqrt(alpha1))

                elseif (tmp1 .ge. 1.0_rp) then
                    ! alpha3=0.0_rp -> cos(alpha3)=1.0
                    sigma1 = sqrt(max(Invariant1/3.0_rp + 2.0_rp*sqrt(alpha1), &
                                      0.0_rp))
                    sigma2 = sqrt(Invariant1/3.0_rp - sqrt(alpha1))
                    sigma3 = sigma2
                else
                    alpha3 = acos(tmp1)/3.0_rp

                  if (abs(Invariant3) .lt. eps) then
                     ! In case of Invariant3=0, one or more eigenvalues are equal to zero
                     ! Therefore force sigma3 to 0 and compute sigma1 and sigma2
                     sigma1 = sqrt(max(Invariant1/3.0_rp + &
                                   2.0_rp*sqrt(alpha1)*cos(alpha3), 0.0_rp))
                     sigma2 = sqrt(abs(Invariant1 - sigma1*sigma1))
                     sigma3 = 0.0_rp
                  else
                     sigma1 = sqrt(max(Invariant1/3.0_rp + &
                                   2.0_rp*sqrt(alpha1)*cos(alpha3), 0.0_rp))
                     sigma2 = sqrt(Invariant1/3.0_rp - &
                                   2.0_rp*sqrt(alpha1)*cos(pi_3 + alpha3))
                     sigma3 = sqrt(abs(Invariant1 - &
                                       sigma1*sigma1-sigma2*sigma2))
                  end if ! Invariant3=0 ?
                end if ! tmp1
              end if ! G diagonal ?

              ! Estimate Dsigma
              if (sigma1 .gt. 0.0_rp) then
                 Dsigma = &
                   sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1*sigma1)
              else
                 Dsigma = 0.0_rp
              end if

              !clipping to avoid negative values
              Dsigma = max(Dsigma, 0.0_rp)

              ! estimate turbulent viscosity

              nut%x(i,1,1,e) = (c*delta%x(i,1,1,e))**2 * Dsigma


       end do
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine sigma_compute_cpu

end module sigma_cpu

