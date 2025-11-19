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
!> Implements the device kernel for the `dynamic_smagorinsky_t` type.
module dynamic_smagorinsky_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : strain_rate
  use coefs, only : coef_t
  use elementwise_filter, only : elementwise_filter_t
  use gs_ops, only : GS_OP_ADD
  use device_math, only : device_col2
  use device_dynamic_smagorinsky_nut, only : device_s_abs_compute, &
       device_lij_compute_part1, device_lij_compute_part2, &
       device_mij_compute_part1, device_mij_nut_compute_part2
  implicit none
  private

  public :: dynamic_smagorinsky_compute_device

contains

  !> Compute eddy viscosity on the device.
  !! @param if_ext If extrapolate the velocity field to evaluate
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
  subroutine dynamic_smagorinsky_compute_device(if_ext, t, tstep, coef, nut, &
       delta, c_dyn, test_filter, mij, lij, num, den)
    logical, intent(in) :: if_ext
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

    if (if_ext .eqv. .true.) then
       u => neko_field_registry%get_field_by_name("u_e")
       v => neko_field_registry%get_field_by_name("v_e")
       w => neko_field_registry%get_field_by_name("w_e")
    else
       u => neko_field_registry%get_field_by_name("u")
       v => neko_field_registry%get_field_by_name("v")
       w => neko_field_registry%get_field_by_name("w")
    end if

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

    call device_col2(s11%x_d, coef%mult_d, s11%dof%size())
    call device_col2(s22%x_d, coef%mult_d, s22%dof%size())
    call device_col2(s33%x_d, coef%mult_d, s33%dof%size())
    call device_col2(s12%x_d, coef%mult_d, s12%dof%size())
    call device_col2(s13%x_d, coef%mult_d, s13%dof%size())
    call device_col2(s23%x_d, coef%mult_d, s23%dof%size())

    call device_s_abs_compute(s_abs%x_d, s11%x_d, s22%x_d, s33%x_d, &
         s12%x_d, s13%x_d, s23%x_d, &
         s11%dof%size())

    call compute_lij_device(lij, u, v, w, test_filter, u%dof%size())
    call compute_nut_device(nut, c_dyn, num, den, lij, mij, &
         s11, s22, s33, s12, s13, s23, &
         s_abs, test_filter, delta, alpha, &
         coef, u%dof%size())

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
  !! @param n
  subroutine compute_lij_device(lij, u, v, w, test_filter, n)
    type(field_t), intent(inout) :: lij(6)
    type(field_t), pointer, intent(in) :: u, v, w
    type(elementwise_filter_t), intent(inout) :: test_filter
    integer, intent(in) :: n
    integer :: i
    !> filtered u,v,w by the test filter
    type(field_t), pointer :: fu, fv, fw, fuu, fvv, fww, fuv, fuw, fvw
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
    call test_filter%apply(fu, u)
    call test_filter%apply(fv, v)
    call test_filter%apply(fw, w)

    !!         ___ ___
    !! Compute u_i*u_j and u_i*u_j
    call device_lij_compute_part1(lij(1)%x_d, lij(2)%x_d, lij(3)%x_d, &
         lij(4)%x_d, lij(5)%x_d, lij(6)%x_d, &
         u%x_d, v%x_d, w%x_d, &
         fu%x_d, fv%x_d, fw%x_d, &
         fuu%x_d, fvv%x_d, fww%x_d, &
         fuv%x_d, fuw%x_d, fvw%x_d, n)

    !! Filter u_i*u_j by the test filter
    call test_filter%apply(fuu, fuu)
    call test_filter%apply(fvv, fvv)
    call test_filter%apply(fww, fww)
    call test_filter%apply(fuv, fuv)
    call test_filter%apply(fuw, fuw)
    call test_filter%apply(fvw, fvw)

    !! Assember the final form
    !!        ___ ___   _______
    !! L_ij = u_i*u_j - u_i*u_j
    call device_lij_compute_part2(lij(1)%x_d, lij(2)%x_d, lij(3)%x_d, &
         lij(4)%x_d, lij(5)%x_d, lij(6)%x_d, &
         fuu%x_d, fvv%x_d, fww%x_d, &
         fuv%x_d, fuw%x_d, fvw%x_d, n)
    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine compute_lij_device

  !> Compute M_ij and nut on the device.
  !!                              _____ ____   __________
  !! M_ij = ((delta_test/delta)^2 s_abs*s_ij - s_abs*s_ij)*(delta^2)
  !! @param nut The SGS viscosity array
  !! @param c_dyn The DS model coefficient
  !! @param num The numerator in the expression of c_dyn, i.e. <mij*lij>
  !! @param den The denominator in the expression of c_dyn, i.e. <mij*mij>
  !! @param lij The Germano identity
  !! @param Mij
  !! @param s11
  !! @param s22
  !! @param s33
  !! @param s12
  !! @param s13
  !! @param s23
  !! @param s_abs
  !! @param test_filter
  !! @param delta The filter size
  !! @param alpha The moving average coefficient
  !! @param n
  subroutine compute_nut_device(nut, c_dyn, num, den, lij, mij, &
       s11, s22, s33, s12, s13, s23, &
       s_abs, test_filter, delta, alpha, &
       coef, n)
    type(field_t), intent(inout) :: nut, c_dyn
    type(field_t), intent(inout) :: num, den
    type(field_t), intent(in) :: lij(6)
    type(field_t), intent(inout) :: mij(6)
    type(field_t), intent(inout) :: s11, s22, s33, s12, s13, s23, s_abs
    type(elementwise_filter_t), intent(inout) :: test_filter
    type(field_t), intent(in) :: delta
    real(kind=rp), intent(in) :: alpha
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n

    real(kind=rp) :: delta_ratio2
    integer :: temp_indices(13)
    type(field_t), pointer :: fs11, fs22, fs33, fs12, fs13, fs23, fs_abs, &
         fsabss11, fsabss22, fsabss33, &
         fsabss12, fsabss13, fsabss23

    delta_ratio2 = ((test_filter%nx-1.0_rp)/(test_filter%nt-1.0_rp))**2

    call neko_scratch_registry%request_field(fs11, temp_indices(1))
    call neko_scratch_registry%request_field(fs22, temp_indices(2))
    call neko_scratch_registry%request_field(fs33, temp_indices(3))
    call neko_scratch_registry%request_field(fs12, temp_indices(4))
    call neko_scratch_registry%request_field(fs13, temp_indices(5))
    call neko_scratch_registry%request_field(fs23, temp_indices(6))
    call neko_scratch_registry%request_field(fsabss11, temp_indices(7))
    call neko_scratch_registry%request_field(fsabss22, temp_indices(8))
    call neko_scratch_registry%request_field(fsabss33, temp_indices(9))
    call neko_scratch_registry%request_field(fsabss12, temp_indices(10))
    call neko_scratch_registry%request_field(fsabss13, temp_indices(11))
    call neko_scratch_registry%request_field(fsabss23, temp_indices(12))
    call neko_scratch_registry%request_field(fs_abs, temp_indices(13))

    !! Compute M_ij
    !!         _____     ____
    !! Compute s_abs and s_ij
    call test_filter%apply(fs_abs, s_abs)
    call test_filter%apply(fs11, s11)
    call test_filter%apply(fs22, s22)
    call test_filter%apply(fs33, s33)
    call test_filter%apply(fs12, s12)
    call test_filter%apply(fs13, s13)
    call test_filter%apply(fs23, s23)

    !!                              _____ ____
    !! Compute (delta_test/delta)^2 s_abs*s_ij and s_abs*s_ij
    call device_mij_compute_part1(mij(1)%x_d, mij(2)%x_d, mij(3)%x_d, &
         mij(4)%x_d, mij(5)%x_d, mij(6)%x_d, &
         s_abs%x_d, s11%x_d, s22%x_d, s33%x_d, &
         s12%x_d, s13%x_d, s23%x_d, &
         fs_abs%x_d, fs11%x_d, fs22%x_d, fs33%x_d, &
         fs12%x_d, fs13%x_d, fs23%x_d, &
         fsabss11%x_d, fsabss22%x_d, fsabss33%x_d, &
         fsabss12%x_d, fsabss13%x_d, fsabss23%x_d, &
         delta_ratio2, n)

    !! Filter s_abs*s_ij by the test filter
    call test_filter%apply(fsabss11, fsabss11)
    call test_filter%apply(fsabss22, fsabss22)
    call test_filter%apply(fsabss33, fsabss33)
    call test_filter%apply(fsabss12, fsabss12)
    call test_filter%apply(fsabss13, fsabss13)
    call test_filter%apply(fsabss23, fsabss23)

    !! Finalise the compute of Mij and nut
    call device_mij_nut_compute_part2(mij(1)%x_d, mij(2)%x_d, mij(3)%x_d, &
         mij(4)%x_d, mij(5)%x_d, mij(6)%x_d, &
         lij(1)%x_d, lij(2)%x_d, lij(3)%x_d, &
         lij(4)%x_d, lij(5)%x_d, lij(6)%x_d, &
         fsabss11%x_d, fsabss22%x_d, fsabss33%x_d, &
         fsabss12%x_d, fsabss13%x_d, fsabss23%x_d, &
         num%x_d, den%x_d, c_dyn%x_d, delta%x_d, &
         s_abs%x_d, nut%x_d, alpha, n)
    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine compute_nut_device

end module dynamic_smagorinsky_device

