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
!> Defines the `wall_model_bc_t` type.
module wall_model_bc
    use num_types
    use bc, only : bc_t
    use, intrinsic :: iso_c_binding, only : c_ptr
    use utils, only : neko_error, nonlinear_index
    use coefs, only : coef_t
    use wall_model, only : wall_model_t
    use rough_log_law, only : rough_log_law_t
    use spalding, only : spalding_t
    use shear_stress, only : shear_stress_t
    implicit none
    private

    !> A shear stress boundary condition, computing the stress values using a
    !! wall model.
    type, public, extends(shear_stress_t) :: wall_model_bc_t
       !> The wall model to compute the stress.
       class(wall_model_t), allocatable :: wall_model
     contains
       procedure, pass(this) :: apply_scalar => wall_model_bc_apply_scalar
       procedure, pass(this) :: apply_vector => wall_model_bc_apply_vector
       procedure, pass(this) :: apply_scalar_dev => wall_model_bc_apply_scalar_dev
       procedure, pass(this) :: apply_vector_dev => wall_model_bc_apply_vector_dev
       procedure, pass(this) :: init_wall_model_bc => &
         wall_model_bc_init_wall_model_bc
    end type wall_model_bc_t

  contains

    !> Apply shear stress for a scalar field @a x.
    subroutine wall_model_bc_apply_scalar(this, x, n, t, tstep)
      class(wall_model_bc_t), intent(inout) :: this
      integer, intent(in) :: n
      real(kind=rp), intent(inout),  dimension(n) :: x
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep

      call neko_error("The wall model bc is not applicable to scalar fields.")

    end subroutine wall_model_bc_apply_scalar

    !> Boundary condition apply for a generic wall_model_bc condition
    !! to vectors @a x, @a y and @a z
    subroutine wall_model_bc_apply_vector(this, x, y, z, n, t, tstep)
      class(wall_model_bc_t), intent(inout) :: this
      integer, intent(in) :: n
      real(kind=rp), intent(inout),  dimension(n) :: x
      real(kind=rp), intent(inout),  dimension(n) :: y
      real(kind=rp), intent(inout),  dimension(n) :: z
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep
      integer :: i, m, k, fid
      real(kind=rp) :: magtau

      call this%wall_model%compute(t, tstep)

      do i=1, this%msk(0)
        magtau = sqrt(this%wall_model%tau_x(i)**2 + this%wall_model%tau_y(i)**2&
                      + this%wall_model%tau_z(i)**2)

        ! Mark sampling nodes with a -1 for debugging
        this%wall_model%tau_field%x(this%wall_model%ind_r(i), &
                                    this%wall_model%ind_s(i), &
                                    this%wall_model%ind_t(i), &
                                    this%wall_model%ind_e(i)) = -1.0_rp
        this%wall_model%tau_field%x(this%msk(i),1,1,1) = magtau
      end do

      call this%shear_stress_t%set_stress(this%wall_model%tau_x, &
                                          this%wall_model%tau_z)
      call this%shear_stress_t%apply_vector(x, y, z, n, t, tstep)

    end subroutine wall_model_bc_apply_vector

    !> Boundary condition apply for a generic wall_model_bc condition
    !! to a vector @a x (device version)
    subroutine wall_model_bc_apply_scalar_dev(this, x_d, t, tstep)
      class(wall_model_bc_t), intent(inout), target :: this
      type(c_ptr) :: x_d
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep

      call neko_error("wall_model_bc bc not implemented on the device")

    end subroutine wall_model_bc_apply_scalar_dev

    !> Boundary condition apply for a generic wall_model_bc condition
    !! to vectors @a x, @a y and @a z (device version)
    subroutine wall_model_bc_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
      class(wall_model_bc_t), intent(inout), target :: this
      type(c_ptr) :: x_d
      type(c_ptr) :: y_d
      type(c_ptr) :: z_d
      real(kind=rp), intent(in), optional :: t
      integer, intent(in), optional :: tstep

      call neko_error("wall_model_bc bc not implemented on the device")

    end subroutine wall_model_bc_apply_vector_dev

    !> Constructor.
    !> @param coef The SEM coefficients.
    subroutine wall_model_bc_init_wall_model_bc(this, coef)
      class(wall_model_bc_t), intent(inout) :: this
      type(coef_t), target, intent(in) :: coef

      call this%shear_stress_t%init_shear_stress(coef)

    end subroutine wall_model_bc_init_wall_model_bc

  end module wall_model_bc