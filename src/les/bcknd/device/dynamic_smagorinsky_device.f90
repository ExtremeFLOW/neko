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
module dynamic_smagorinsky_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use math, only : cadd, col2, sub2, col3, cmult
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : strain_rate
  use coefs, only : coef_t
  use elementwise_filter, only : elementwise_filter_t
  use gs_ops, only : GS_OP_ADD
  implicit none
  private

  public :: dynamic_smagorinsky_compute_device

contains

  !> Compute eddy viscosity on the device.
  !! @param t The time value.
  !! @param tstep The current time-step.
  !! @param coef SEM coefficients.
  !! @param nut The SGS viscosity array.
  !! @param delta The LES lengthscale.
  !! @param c_dyn The DS model coefficient
  !! @param test_filter
  !! @param mij
  !! @param lij The Germano identity.
  !! @param num The numerator in the expression of c_dyn, i.e. <mij*lij>
  !! @param den The denominator in the expression of c_dyn, i.e. <mij*mij>
  subroutine dynamic_smagorinsky_compute_device(t, tstep, coef, nut, delta, &
                                             c_dyn, test_filter, mij, lij, num, den)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: nut
    type(field_t), intent(in) :: delta
    type(field_t), intent(inout) :: c_dyn
    type(elementwise_filter_t), intent(inout) :: test_filter
    type(field_t), intent(inout) :: mij(6), lij(6)
    type(field_t), intent(inout) :: num, den

    type(field_t), pointer :: u, v, w
    type(field_t) :: c_dyn_curr
    ! the strain rate tensor
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23, s_abs
    real(kind=rp) :: alpha ! running averaging coefficient
    integer :: temp_indices(7)
    integer :: i

    if (tstep .eq. 1) then
       alpha = 1.0_rp
    else
       alpha = 0.9_rp
    end if

    u => neko_field_registry%get_field_by_name("u")
    v => neko_field_registry%get_field_by_name("v")
    w => neko_field_registry%get_field_by_name("w")

    call neko_scratch_registry%request_field(s11, temp_indices(1))
    call neko_scratch_registry%request_field(s22, temp_indices(2))
    call neko_scratch_registry%request_field(s33, temp_indices(3))
    call neko_scratch_registry%request_field(s12, temp_indices(4))
    call neko_scratch_registry%request_field(s13, temp_indices(5))
    call neko_scratch_registry%request_field(s23, temp_indices(6))
    call neko_scratch_registry%request_field(s_abs, temp_indices(7))

    ! Compute the strain rate tensor
    call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, u, v, w, coef)

    call coef%gs_h%op(s11%x, s11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(s22%x, s11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(s33%x, s11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(s12%x, s11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(s13%x, s11%dof%size(), GS_OP_ADD)
    call coef%gs_h%op(s23%x, s11%dof%size(), GS_OP_ADD)

    call device_s_abs_compute(s_abs%x_d, s11%x_d, s22%x_d, s33%x_d, &
                                         s12%x_d, s13%x_d, s23%x_d, &
                                         coef%mult_d, s11%dof%size())

    call compute_lij_device(lij, u, v, w, test_filter, u%dof%size(), u%msh%nelv)
    call compute_nut_device(num, den, lij, mij, &
                            s11, s22, s33, s12, s13, s23, &
                            s_abs, test_filter, delta, alpha, &
                            u%dof%size(), u%msh%nelv)
    
    do concurrent (i =1:u%dof%size())
       if (den%x(i,1,1,1) .gt. 0.0_rp) then
          c_dyn%x(i,1,1,1) = 0.5_rp * (num%x(i,1,1,1)/den%x(i,1,1,1))
       else
          c_dyn%x(i,1,1,1) = 0.0_rp
       end if
       c_dyn%x(i,1,1,1) = max(c_dyn%x(i,1,1,1),0.0_rp)
       nut%x(i,1,1,1) = c_dyn%x(i,1,1,1) * delta%x(i,1,1,1)**2 * s_abs%x(i,1,1,1)
    end do

    call coef%gs_h%op(nut, GS_OP_ADD)
    call device_col2(nut%x_d, coef%mult_d, nut%dof%size())

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine dynamic_smagorinsky_compute_device

  !> Compute Germano Identity on the device.
  !!        ___ ___   _______
  !! L_ij = u_i*u_j - u_i*u_j
  !! @param lij The Germano identity.
  !! @param u x-velocity resolved (only filtered once)
  !! @param v y-velocity resolved (only filtered once)
  !! @param w z-velocity resolved (only filtered once)
  !! @param test_filter
  subroutine compute_lij_device(lij, u, v, w, test_filter, n, nelv)
    type(field_t), intent(inout) :: lij(6)
    type(field_t), pointer, intent(in) :: u, v, w
    type(elementwise_filter_t), intent(inout) :: test_filter
    integer, intent(in) :: n
    integer, intent(inout) :: nelv
    integer :: i
    !> filtered u,v,w by the test filter
    type(field_t) :: fu, fv, fw, fuu, fvv, fww, fuv, fuw, fvw
    integer :: temp_indices(9)

    call neko_scratch_registry%request_field(fu, temp_indices(1))
    call neko_scratch_registry%request_field(fv, temp_indices(2))
    call neko_scratch_registry%request_field(fw, temp_indices(3))
    call neko_scratch_registry%request_field(fuu, temp_indices(4))
    call neko_scratch_registry%request_field(fvv, temp_indices(5))
    call neko_scratch_registry%request_field(fww, temp_indices(6))
    call neko_scratch_registry%request_field(fuv, temp_indices(7))
    call neko_scratch_registry%request_field(fuw, temp_indices(8))
    call neko_scratch_registry%request_field(fvw, temp_indices(9))

    ! Use test filter for the velocity fields
    call test_filter%filter_3d(fu%x, u%x, nelv)
    call test_filter%filter_3d(fv%x, v%x, nelv)
    call test_filter%filter_3d(fw%x, w%x, nelv)

    !!!!!!
    !!!!!! PART I START
    call device_lij_compute_part1(lij(1)%x_d, lij(2)%x_d, lij(3)%x_d, &
                                  lij(4)%x_d, lij(5)%x_d, lij(6)%x_d, &
                                  u%x_d, v%x_d, w%x_d, &
                                  fu%x_d, fv%x_d, fw%x_d, &
                                  fuu%x_d, fvv%x_d, fww%x_d, &
                                  fuv%x_d, fuw%x_d, fvw%x_d, n)

    !! The first term
    do concurrent (i = 1:n)
       lij(1)%x(i,1,1,1) = fu(i) * fu(i)
       lij(2)%x(i,1,1,1) = fv(i) * fv(i)
       lij(3)%x(i,1,1,1) = fw(i) * fw(i)
       lij(4)%x(i,1,1,1) = fu(i) * fv(i)
       lij(5)%x(i,1,1,1) = fu(i) * fw(i)
       lij(6)%x(i,1,1,1) = fv(i) * fw(i)
    end do

    !! Subtract the second term:
    !! use test filter for the cross terms
    !! fu and fv are used as work array
    call col3(fuu, u%x, u%x, n)
    call col3(fvv, v%x, v%x, n)
    call col3(fww, w%x, w%x, n)
    call col3(fuv, u%x, v%x, n)
    call col3(fuw, u%x, w%x, n)
    call col3(fvw, v%x, w%x, n)
    !!!!!!
    !!!!!! PART I FINISH 

    call test_filter%filter_3d(fuu%x, fuu%x, nelv)
    call test_filter%filter_3d(fvv%x, fvv%x, nelv)
    call test_filter%filter_3d(fww%x, fww%x, nelv)
    call test_filter%filter_3d(fuv%x, fuv%x, nelv)
    call test_filter%filter_3d(fuw%x, fuw%x, nelv)
    call test_filter%filter_3d(fvw%x, fvw%x, nelv)

    !!!!!!
    !!!!!! PART II START
    call device_lij_compute_part2(lij(1)%x_d, lij(2)%x_d, lij(3)%x_d, &
                                  lij(4)%x_d, lij(5)%x_d, lij(6)%x_d, &
                                  fuu%x_d, fvv%x_d, fww%x_d, &
                                  fuv%x_d, fuw%x_d, fvw%x_d, n)
    call sub2(lij(1)%x, fuu, n)
    call sub2(lij(2)%x, fvv, n)
    call sub2(lij(3)%x, fww, n)
    call sub2(lij(4)%x, fuv, n)
    call sub2(lij(5)%x, fuw, n)
    call sub2(lij(6)%x, fvw, n)
    !!!!!!
    !!!!!! PART II FINISH

  end subroutine compute_lij_device

  !> Compute M_ij on the CPU.
  !!                              _____ ____   __________
  !! M_ij = ((delta_test/delta)^2 s_abs*s_ij - s_abs*s_ij)*(delta^2)
  !! @param Mij
  !! @param u x-velocity resolved (only filtered once)
  !! @param v y-velocity resolved (only filtered once)
  !! @param w z-velocity resolved (only filtered once)
  !! @param test_filter
  subroutine compute_mij_cpu(mij, s11, s22, s33, s12, s13, s23, &
                             s_abs, test_filter, delta, n, nelv)
    type(field_t), intent(inout) :: mij(6)
    type(field_t), intent(inout) :: s11, s22, s33, s12, s13, s23, s_abs
    type(elementwise_filter_t), intent(inout) :: test_filter
    type(field_t), intent(in) :: delta
    integer, intent(in) :: n
    integer, intent(inout) :: nelv

    real(kind=rp), dimension(n) :: fs11, fs22, fs33, fs12, fs13, fs23, fs_abs
    real(kind=rp) :: delta_ratio2 !test- to grid- filter ratio, squared
    integer :: i
    real(kind=rp) :: delta2

    delta_ratio2 = ((test_filter%nx-1)/(test_filter%nt-1))**2

    !! The first term:
    !!                      _____ ____
    !! (delta_test/delta)^2 s_abs*s_ij
    call test_filter%filter_3d(fs_abs, s_abs%x, nelv)

    call test_filter%filter_3d(fs11, s11%x, nelv)
    call col3(mij(1)%x, fs_abs, fs11, n)
    call cmult(mij(1)%x, delta_ratio2, n)

    call test_filter%filter_3d(fs22, s22%x, nelv)
    call col3(mij(2)%x, fs_abs, fs11, n)
    call cmult(mij(2)%x, delta_ratio2, n)

    call test_filter%filter_3d(fs33, s33%x, nelv)
    call col3(mij(3)%x, fs_abs, fs11, n)
    call cmult(mij(3)%x, delta_ratio2, n)

    call test_filter%filter_3d(fs12, s12%x, nelv)
    call col3(mij(4)%x, fs_abs, fs11, n)
    call cmult(mij(4)%x, delta_ratio2, n)

    call test_filter%filter_3d(fs13, s13%x, nelv)
    call col3(mij(5)%x, fs_abs, fs11, n)
    call cmult(mij(5)%x, delta_ratio2, n)

    call test_filter%filter_3d(fs23, s23%x, nelv)
    call col3(mij(6)%x, fs_abs, fs11, n)
    call cmult(mij(6)%x, delta_ratio2, n)

    !! Substract the second term:
    !!                      _____ ____   __________
    !! (delta_test/delta)^2 s_abs*s_ij - s_abs*s_ij
    !! fs11 and fs22 are used as work array
    call col3(fs11, s_abs%x, s11%x, n)
    call test_filter%filter_3d(fs22, fs11, nelv)
    call sub2(mij(1)%x, fs22, n)

    call col3(fs11, s_abs%x, s22%x, n)
    call test_filter%filter_3d(fs22, fs11, nelv)
    call sub2(mij(2)%x, fs22, n)

    call col3(fs11, s_abs%x, s33%x, n)
    call test_filter%filter_3d(fs22, fs11, nelv)
    call sub2(mij(3)%x, fs22, n)

    call col3(fs11, s_abs%x, s12%x, n)
    call test_filter%filter_3d(fs22, fs11, nelv)
    call sub2(mij(4)%x, fs22, n)

    call col3(fs11, s_abs%x, s13%x, n)
    call test_filter%filter_3d(fs22, fs11, nelv)
    call sub2(mij(5)%x, fs22, n)

    call col3(fs11, s_abs%x, s23%x, n)
    call test_filter%filter_3d(fs22, fs11, nelv)
    call sub2(mij(6)%x, fs22, n)

    !! Lastly multiplied by delta^2
    do concurrent (i = 1:n)
       delta2 = delta%x(i,1,1,1)**2
       mij(1)%x(i,1,1,1) = mij(1)%x(i,1,1,1) * delta2
       mij(2)%x(i,1,1,1) = mij(2)%x(i,1,1,1) * delta2
       mij(3)%x(i,1,1,1) = mij(3)%x(i,1,1,1) * delta2
       mij(4)%x(i,1,1,1) = mij(4)%x(i,1,1,1) * delta2
       mij(5)%x(i,1,1,1) = mij(5)%x(i,1,1,1) * delta2
       mij(6)%x(i,1,1,1) = mij(6)%x(i,1,1,1) * delta2
    end do

  end subroutine compute_mij_cpu

  !> Compute numerator and denominator for c_dyn on the CPU.
  !! @param num The numerator in the expression of c_dyn, i.e. <mij*lij>
  !! @param den The denominator in the expression of c_dyn, i.e. <mij*mij>
  !! @param mij
  !! @param lij The Germano identity.
  !! @param alpha The moving average coefficient 
  subroutine compute_num_den_cpu(num, den, lij, mij, alpha, n)
    type(field_t), intent(inout) :: num, den
    type(field_t), intent(in) :: lij(6), mij(6)
    real(kind=rp), intent(in) :: alpha
    integer, intent(in) :: n

    real(kind=rp), dimension(n) :: num_curr, den_curr
    integer :: i

    do concurrent (i = 1:n)
       num_curr(i) = mij(1)%x(i,1,1,1)*lij(1)%x(i,1,1,1) + &
                     mij(2)%x(i,1,1,1)*lij(2)%x(i,1,1,1) + &
                     mij(3)%x(i,1,1,1)*lij(3)%x(i,1,1,1) + &
                     2.0_rp*(mij(4)%x(i,1,1,1)*lij(4)%x(i,1,1,1) + &
                     mij(5)%x(i,1,1,1)*lij(5)%x(i,1,1,1) + &
                     mij(6)%x(i,1,1,1)*lij(6)%x(i,1,1,1))
       den_curr(i) = mij(1)%x(i,1,1,1)*mij(1)%x(i,1,1,1) + &
                     mij(2)%x(i,1,1,1)*mij(2)%x(i,1,1,1) + &
                     mij(3)%x(i,1,1,1)*mij(3)%x(i,1,1,1) + &
                     2.0_rp*(mij(4)%x(i,1,1,1)*mij(4)%x(i,1,1,1) + &
                     mij(5)%x(i,1,1,1)*mij(5)%x(i,1,1,1) + &
                     mij(6)%x(i,1,1,1)*mij(6)%x(i,1,1,1))
    end do

    ! running average over time
    do concurrent (i = 1:n)
       num%x(i,1,1,1) = alpha * num%x(i,1,1,1) + (1.0_rp - alpha) * num_curr(i)
       den%x(i,1,1,1) = alpha * den%x(i,1,1,1) + (1.0_rp - alpha) * den_curr(i)
    end do

  end subroutine compute_num_den_cpu

end module dynamic_smagorinsky_device

