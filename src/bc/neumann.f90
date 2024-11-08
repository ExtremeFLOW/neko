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
!> Defines a Neumann boundary condition.
module neumann
  use num_types, only : rp
  use bc, only : bc_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error, nonlinear_index
  use coefs, only : coef_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use math, only : cfill
  use math, only : cfill, copy
  implicit none
  private

  !> A Neumann boundary condition for scalar fields.
  !! Sets the flux of the field to the chosen value.
  !! @note The condition is imposed weekly by adding an appropriate source term
  !! to the right-hand-side.
  type, public, extends(bc_t) :: neumann_t
     real(kind=rp), allocatable, private :: flux_(:)
     real(kind=rp), private ::  init_flux_
   contains
     procedure, pass(this) :: apply_scalar => neumann_apply_scalar
     procedure, pass(this) :: apply_vector => neumann_apply_vector
     procedure, pass(this) :: apply_scalar_dev => neumann_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => neumann_apply_vector_dev
     !> Constructor
     procedure, pass(this) :: init => neumann_init
     !> Constructor from components
     procedure, pass(this) :: init_from_components => &
        neumann_init_from_components
     procedure, pass(this) :: flux => neumann_flux
     !> Set the flux using a scalar.
     procedure, pass(this) :: set_flux_scalar => neumann_set_flux_scalar
     !> Set the flux using an array.
     procedure, pass(this) :: set_flux_array => neumann_set_flux_array
     !> Generic interface for setting the flux.
     generic :: set_flux => set_flux_scalar, set_flux_array
     !> Destructor.
     procedure, pass(this) :: free => neumann_free
     !> Finalize.
     procedure, pass(this) :: finalize => neumann_finalize
  end type neumann_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine neumann_init(this, coef, json)
    class(neumann_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    real(kind=rp) :: flux

    call this%init_base(coef)
    this%strong = .false.

    call json_get(json, "value", flux)
    this%init_flux_ = flux
  end subroutine neumann_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] g The value to apply at the boundary.
  subroutine neumann_init_from_components(this, coef, flux)
    class(neumann_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: flux

    call this%init_base(coef)
    this%init_flux_ = flux
  end subroutine neumann_init_from_components

  !> Boundary condition apply for a generic Neumann condition
  !! to a vector @a x
  subroutine neumann_apply_scalar(this, x, n, t, tstep)
    class(neumann_t), intent(inout) :: this
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
       select case (facet)
       case (1,2)
          x(k) = x(k) + this%flux_(i)*this%coef%area(idx(2), idx(3), facet, &
               idx(4))
       case (3,4)
          x(k) = x(k) + this%flux_(i)*this%coef%area(idx(1), idx(3), facet, &
               idx(4))
       case (5,6)
          x(k) = x(k) + this%flux_(i)*this%coef%area(idx(1), idx(2), facet, &
               idx(4))
       end select
    end do
  end subroutine neumann_apply_scalar

  !> Boundary condition apply for a generic Neumann condition
  !! to vectors @a x, @a y and @a z
  subroutine neumann_apply_vector(this, x, y, z, n, t, tstep)
    class(neumann_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("Neumann bc not implemented for vectors")

  end subroutine neumann_apply_vector

  !> Boundary condition apply for a generic Neumann condition
  !! to a vector @a x (device version)
  subroutine neumann_apply_scalar_dev(this, x_d, t, tstep)
    class(neumann_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("Neumann bc not implemented on the device")

  end subroutine neumann_apply_scalar_dev

  !> Boundary condition apply for a generic Neumann condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine neumann_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(neumann_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call neko_error("Neumann bc not implemented on the device")

  end subroutine neumann_apply_vector_dev

  !> Destructor
  subroutine neumann_free(this)
    class(neumann_t), target, intent(inout) :: this

    call this%free_base()

  end subroutine neumann_free

  !> Finalize by setting the flux.
  subroutine neumann_finalize(this)
    class(neumann_t), target, intent(inout) :: this

    call this%finalize_base()
    allocate(this%flux_(this%msk(0)))

    call cfill(this%flux_, this%init_flux_, this%msk(0))
  end subroutine neumann_finalize

  !> Get the flux.
  pure function neumann_flux(this) result(flux)
    class(neumann_t), intent(in) :: this
    real(kind=rp) :: flux(this%msk(0))

    flux = this%flux_
  end function neumann_flux

  !> Set the flux using a scalar.
  subroutine neumann_set_flux_scalar(this, flux)
    class(neumann_t), intent(inout) :: this
    real(kind=rp), intent(in) :: flux

    this%flux_ = flux
  end subroutine neumann_set_flux_scalar

  !> Set the flux using an array of values.
  !> @param flux The desired flux.
  subroutine neumann_set_flux_array(this, flux)
    class(neumann_t), intent(inout) :: this
    real(kind=rp), intent(in) :: flux(this%msk(0))

    call copy(this%flux_, flux, this%msk(0))
  end subroutine neumann_set_flux_array

end module neumann
