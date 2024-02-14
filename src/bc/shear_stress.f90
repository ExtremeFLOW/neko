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
!> Defines a shear stress boundary condition for a vector field.
module shear_stress
    use num_types
    use bc, only : bc_t
    use, intrinsic :: iso_c_binding, only : c_ptr
    use utils, only : neko_error, nonlinear_index
    use coefs, only : coef_t
    implicit none
    private

    !> A shear stress boundary condition.
    !! This sets a given flux of the boundary-parallel components of the vector
    !! field and a homogenoeus Dirichlet condition on the wall-normal component.
    !! @note The condition is imposed weekly by adding an appropriate source
    !! term to the right-hand-side.
    type, public, extends(bc_t) :: shear_stress_t
       real(kind=rp), private :: flux_
       !> SEM coeffs.
       type(coef_t), pointer :: coef
     contains
       procedure, pass(this) :: apply_scalar => shear_stress_apply_scalar
       procedure, pass(this) :: apply_vector => shear_stress_apply_vector
       procedure, pass(this) :: apply_scalar_dev => shear_stress_apply_scalar_dev
       procedure, pass(this) :: apply_vector_dev => shear_stress_apply_vector_dev
       procedure, pass(this) :: init_shear_stress => shear_stress_init_shear_stress
       procedure, pass(this) :: flux => shear_stress_flux
    end type shear_stress_t

  contains

    !> Boundary condition apply for a generic shear_stress condition
    !! to a vector @a x
    subroutine shear_stress_apply_scalar(this, x, n, t, tstep)
      class(shear_stress_t), intent(inout) :: this
      integer, intent(in) :: n
      real(kind=rp), intent(inout),  dimension(n) :: x
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep
      integer :: i, m, k, facet
      ! Store non-linear index
      integer :: idx(4)

      m = this%msk(0)
      do i = 1, m
         k = this%msk(i)
         facet = this%facet(i)
         idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%lx,&
                               this%coef%Xh%lx)
         select case(facet)
         case(1,2)
            x(k) = x(k) + this%flux_*this%coef%area(idx(2), idx(3), facet, idx(4))
         case(3,4)
            x(k) = x(k) + this%flux_*this%coef%area(idx(1), idx(3), facet, idx(4))
         case(5,6)
            x(k) = x(k) + this%flux_*this%coef%area(idx(1), idx(2), facet, idx(4))
         end select
      end do
    end subroutine shear_stress_apply_scalar

    !> Boundary condition apply for a generic shear_stress condition
    !! to vectors @a x, @a y and @a z
    subroutine shear_stress_apply_vector(this, x, y, z, n, t, tstep)
      class(shear_stress_t), intent(inout) :: this
      integer, intent(in) :: n
      real(kind=rp), intent(inout),  dimension(n) :: x
      real(kind=rp), intent(inout),  dimension(n) :: y
      real(kind=rp), intent(inout),  dimension(n) :: z
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep
      integer :: i, m, k

      call neko_error("shear_stress bc not implemented for vectors")

    end subroutine shear_stress_apply_vector

    !> Boundary condition apply for a generic shear_stress condition
    !! to a vector @a x (device version)
    subroutine shear_stress_apply_scalar_dev(this, x_d, t, tstep)
      class(shear_stress_t), intent(inout), target :: this
      type(c_ptr) :: x_d
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep

      call neko_error("shear_stress bc not implemented on the device")

    end subroutine shear_stress_apply_scalar_dev

    !> Boundary condition apply for a generic shear_stress condition
    !! to vectors @a x, @a y and @a z (device version)
    subroutine shear_stress_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
      class(shear_stress_t), intent(inout), target :: this
      type(c_ptr) :: x_d
      type(c_ptr) :: y_d
      type(c_ptr) :: z_d
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep

      call neko_error("shear_stress bc not implemented on the device")

    end subroutine shear_stress_apply_vector_dev

    !> Constructor
    !> @param flux The desired flux.
    !> @param coef The SEM coefficients.
    subroutine shear_stress_init_shear_stress(this, flux, coef)
      class(shear_stress_t), intent(inout) :: this
      real(kind=rp), intent(in) :: flux
      type(coef_t), target, intent(in) :: coef

      this%flux_ = flux
      this%coef => coef
    end subroutine shear_stress_init_shear_stress

    !> Get the set flux.
    pure function shear_stress_flux(this) result(flux)
      class(shear_stress_t), intent(in) :: this
      real(kind=rp) :: flux

      flux = this%flux_
    end function shear_stress_flux
  end module shear_stress