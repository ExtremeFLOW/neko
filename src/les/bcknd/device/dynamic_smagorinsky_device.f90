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

    ! Use test filter for the velocity fields
    call test_filter%filter_3d(fu, u%x, nelv)
    call test_filter%filter_3d(fv, v%x, nelv)
    call test_filter%filter_3d(fw, w%x, nelv)

    !!!!!!
    !!!!!! PART I START

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

    call test_filter%filter_3d(fuu, fuu, nelv)
    call test_filter%filter_3d(fvv, fvv, nelv)
    call test_filter%filter_3d(fww, fww, nelv)
    call test_filter%filter_3d(fuv, fuv, nelv)
    call test_filter%filter_3d(fuw, fuw, nelv)
    call test_filter%filter_3d(fvw, fvw, nelv)

    !!!!!!
    !!!!!! PART II START
    call sub2(lij(1)%x, fv, n)
    call sub2(lij(2)%x, fv, n)
    call sub2(lij(3)%x, fv, n)
    call sub2(lij(4)%x, fv, n)
    call sub2(lij(5)%x, fv, n)
    call sub2(lij(6)%x, fv, n)
    !!!!!!
    !!!!!! PART II FINISH

  end subroutine compute_lij_device

end module dynamic_smagorinsky_device

