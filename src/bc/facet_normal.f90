! Copyright (c) 2020-2023, The Neko Authors
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
!> Dirichlet condition applied in the facet normal direction
module facet_normal
  use device_facet_normal
  use num_types
  use math
  use coefs, only : coef_t
  use bc, only : bc_t
  use utils
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Dirichlet condition in facet normal direction
  type, public, extends(bc_t) :: facet_normal_t
   contains
     procedure, pass(this) :: apply_scalar => facet_normal_apply_scalar
     procedure, pass(this) :: apply_scalar_dev => facet_normal_apply_scalar_dev
     procedure, pass(this) :: apply_vector => facet_normal_apply_vector
     procedure, pass(this) :: apply_vector_dev => facet_normal_apply_vector_dev
     procedure, pass(this) :: apply_surfvec => facet_normal_apply_surfvec
     procedure, pass(this) :: apply_surfvec_dev => facet_normal_apply_surfvec_dev
     !> Destructor.
     procedure, pass(this) :: free => facet_normal_free
  end type facet_normal_t

contains

  !> No-op scalar apply
  subroutine facet_normal_apply_scalar(this, x, n, t, tstep)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine facet_normal_apply_scalar

  !> No-op scalar apply on device
  subroutine facet_normal_apply_scalar_dev(this, x_d, t, tstep)
    class(facet_normal_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

  end subroutine facet_normal_apply_scalar_dev

  !> No-op vector apply on device
  subroutine facet_normal_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(facet_normal_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

  end subroutine facet_normal_apply_vector_dev

  !> No-op vector apply
  subroutine facet_normal_apply_vector(this, x, y, z, n, t, tstep)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine facet_normal_apply_vector

  !> Apply in facet normal direction (vector valued)
  subroutine facet_normal_apply_surfvec(this, x, y, z, u, v, w, n, t, tstep)
    class(facet_normal_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    real(kind=rp), intent(inout), dimension(n) :: u
    real(kind=rp), intent(inout), dimension(n) :: v
    real(kind=rp), intent(inout), dimension(n) :: w
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k, idx(4), facet

    associate(c => this%coef)
      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, c%Xh%lx, c%Xh%lx, c%Xh%lx)
         select case(facet)
         case(1,2)
            x(k) = u(k) * c%nx(idx(2), idx(3), facet, idx(4)) &
                 * c%area(idx(2), idx(3), facet, idx(4))
            y(k) = v(k) * c%ny(idx(2), idx(3), facet, idx(4)) &
                 * c%area(idx(2), idx(3), facet, idx(4))
            z(k) = w(k) * c%nz(idx(2), idx(3), facet, idx(4)) &
                 * c%area(idx(2), idx(3), facet, idx(4))
         case(3,4)
            x(k) = u(k) * c%nx(idx(1), idx(3), facet, idx(4)) &
                 * c%area(idx(1), idx(3), facet, idx(4))
            y(k) = v(k) * c%ny(idx(1), idx(3), facet, idx(4)) &
                 * c%area(idx(1), idx(3), facet, idx(4))
            z(k) = w(k) * c%nz(idx(1), idx(3), facet, idx(4)) &
                 * c%area(idx(1), idx(3), facet, idx(4))
         case(5,6)
            x(k) = u(k) * c%nx(idx(1), idx(2), facet, idx(4)) &
                 * c%area(idx(1), idx(2), facet, idx(4))
            y(k) = v(k) * c%ny(idx(1), idx(2), facet, idx(4)) &
                 * c%area(idx(1), idx(2), facet, idx(4))
            z(k) = w(k) * c%nz(idx(1), idx(2), facet, idx(4)) &
                 * c%area(idx(1), idx(2), facet, idx(4))
         end select
      end do
    end associate

  end subroutine facet_normal_apply_surfvec

  !> Apply in facet normal direction (vector valued, device version)
  subroutine facet_normal_apply_surfvec_dev(this, x_d, y_d, z_d, &
                                            u_d, v_d, w_d, t, tstep)
    class(facet_normal_t), intent(inout), target :: this
    type(c_ptr) :: x_d, y_d, z_d, u_d, v_d, w_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    associate(c => this%coef)
      call device_facet_normal_apply_surfvec(this%msk_d, this%facet_d, &
                                             x_d, y_d, z_d, u_d, v_d, w_d, &
                                             c%nx_d, c%ny_d, c%nz_d, c%area_d, &
                                             c%Xh%lx, size(this%msk))
    end associate

  end subroutine facet_normal_apply_surfvec_dev

  !> Destructor
  subroutine facet_normal_free(this)
    class(facet_normal_t), target, intent(inout) :: this

    call this%free_base()

  end subroutine facet_normal_free

end module facet_normal
