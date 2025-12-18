! Copyright (c) 2024-2025, The Neko Authors
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
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr
  use utils, only : neko_error, nonlinear_index
  use coefs, only : coef_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_lookup
  use math, only : cfill, copy, abscmp
  use vector, only: vector_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device_math, only : device_cfill, device_copy
  use device, only : device_memcpy, DEVICE_TO_HOST
  use device_neumann, only : device_neumann_apply_scalar, &
       device_neumann_apply_vector
  use time_state, only : time_state_t
  implicit none
  private

  !> A Neumann boundary condition.
  !! Sets the flux of the field to the chosen values.
  !! @details The condition is imposed weakly by adding an appropriate source
  !! term to the right-hand-side. At construction time, we check if the
  !! prescribed flux is zero and if so, the condition just does nothing. Setting
  !! the flux using the `set_flux` routine, automatically removes this
  !! assumption.
  type, public, extends(bc_t) :: neumann_t
     !> The flux values at the boundary. Each vector in the array corresponds to
     !> a component of the flux.
     type(vector_t), allocatable :: flux(:)
     !> An initial flux value set at construction. A constant for each
     !> component. Copied to `flux` at finalization. This is needed because we
     !> do not know the size of the boundary at construction time, only at
     !> finalization.
     real(kind=rp), allocatable, private :: init_flux_(:)
     !> Flag for whether we have a homogeneous Neumann condition, i.e.
     !! "do nothing".
     logical :: uniform_0 = .false.
   contains
     procedure, pass(this) :: apply_scalar => neumann_apply_scalar
     procedure, pass(this) :: apply_vector => neumann_apply_vector
     procedure, pass(this) :: apply_scalar_dev => neumann_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => neumann_apply_vector_dev
     !> Constructor
     procedure, pass(this) :: init => neumann_init
     !> Constructor from components, one flux value per component.
     procedure, pass(this) :: neumann_init_from_components_array
     !> Constructor from components, single flux for the single-component case.
     procedure, pass(this) :: neumann_init_from_components_single
     generic :: init_from_components => &
          neumann_init_from_components_array, &
          neumann_init_from_components_single
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
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    real(kind=rp) :: flux
    logical :: found

    call this%init_base(coef)
    this%strong = .false.

    ! Try to read array from json
    call json%get("flux", this%init_flux_, found)

    ! If we haven't found an array, try to read a single value
    if (.not. found) then
       call json_get_or_lookup(json, "flux", flux)
       allocate(this%init_flux_(1))
       this%init_flux_(1) = flux
    end if

    if ((size(this%init_flux_) .ne. 1) &
         .and. (size(this%init_flux_) .ne. 3)) then
       call neko_error("Neumann BC flux must be a scalar or a 3-component" // &
            " vector.")
    end if

    allocate(this%flux(size(this%init_flux_)))
  end subroutine neumann_init

  !> Constructor from components, using a flux array for vector components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] flux The value of the flux at the boundary.
  subroutine neumann_init_from_components_array(this, coef, flux)
    class(neumann_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: flux(3)

    call this%init_base(coef)
    this%init_flux_ = flux

    if ((size(this%init_flux_) .ne. 3)) then
       call neko_error("Neumann BC flux must be a scalar or a 3-component" // &
            " vector.")
    end if
  end subroutine neumann_init_from_components_array

  !> Constructor from components, using an signle flux.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] flux The value of the flux at the boundary.
  subroutine neumann_init_from_components_single(this, coef, flux)
    class(neumann_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: flux

    call this%init_base(coef)
    allocate(this%init_flux_(1))
    this%init_flux_(1) = flux

  end subroutine neumann_init_from_components_single

  !> Boundary condition apply for a generic Neumann condition
  !! to a vector @a x
  subroutine neumann_apply_scalar(this, x, n, time, strong)
    class(neumann_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m, k, facet
    ! Store non-linear index
    integer :: idx(4)
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    m = this%msk(0)
    if (.not. strong_) then
       do i = 1, m
          k = this%msk(i)
          facet = this%facet(i)
          idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%lx,&
               this%coef%Xh%lx)
          select case (facet)
          case (1,2)
             x(k) = x(k) + &
                  this%flux(1)%x(i)*this%coef%area(idx(2), idx(3), facet, idx(4))
          case (3,4)
             x(k) = x(k) + &
                  this%flux(1)%x(i)*this%coef%area(idx(1), idx(3), facet, idx(4))
          case (5,6)
             x(k) = x(k) + &
                  this%flux(1)%x(i)*this%coef%area(idx(1), idx(2), facet, idx(4))
          end select
       end do
    end if
  end subroutine neumann_apply_scalar

  !> Boundary condition apply for a generic Neumann condition
  !! to vectors @a x, @a y and @a z
  subroutine neumann_apply_vector(this, x, y, z, n, time, strong)
    class(neumann_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m, k, facet
    ! Store non-linear index
    integer :: idx(4)
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    m = this%msk(0)
    if (.not. strong_) then
       do i = 1, m
          k = this%msk(i)
          facet = this%facet(i)
          idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%lx,&
               this%coef%Xh%lx)
          select case (facet)
          case (1,2)
             x(k) = x(k) + &
                  this%flux(1)%x(i)*this%coef%area(idx(2), idx(3), facet, idx(4))
             y(k) = y(k) + &
                  this%flux(2)%x(i)*this%coef%area(idx(2), idx(3), facet, idx(4))
             z(k) = z(k) + &
                  this%flux(3)%x(i)*this%coef%area(idx(2), idx(3), facet, idx(4))
          case (3,4)
             x(k) = x(k) + &
                  this%flux(1)%x(i)*this%coef%area(idx(1), idx(3), facet, idx(4))
             y(k) = y(k) + &
                  this%flux(2)%x(i)*this%coef%area(idx(1), idx(3), facet, idx(4))
             z(k) = z(k) + &
                  this%flux(3)%x(i)*this%coef%area(idx(1), idx(3), facet, idx(4))
          case (5,6)
             x(k) = x(k) + &
                  this%flux(1)%x(i)*this%coef%area(idx(1), idx(2), facet, idx(4))
             y(k) = y(k) + &
                  this%flux(2)%x(i)*this%coef%area(idx(1), idx(2), facet, idx(4))
             z(k) = z(k) + &
                  this%flux(3)%x(i)*this%coef%area(idx(1), idx(2), facet, idx(4))
          end select
       end do
    end if
  end subroutine neumann_apply_vector

  !> Boundary condition apply for a generic Neumann condition
  !! to a vector @a x (device version)
  subroutine neumann_apply_scalar_dev(this, x_d, time, strong, strm)
    class(neumann_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (.not. this%uniform_0 .and. this%msk(0) .gt. 0 .and. &
         .not. strong_) then
       call device_neumann_apply_scalar(this%msk_d, this%facet_d, x_d, &
            this%flux(1)%x_d, this%coef%area_d, this%coef%Xh%lx, &
            size(this%msk), strm)
    end if
  end subroutine neumann_apply_scalar_dev

  !> Boundary condition apply for a generic Neumann condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine neumann_apply_vector_dev(this, x_d, y_d, z_d, &
       time, strong, strm)
    class(neumann_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (.not. this%uniform_0 .and. this%msk(0) .gt. 0 .and. &
         .not. strong_) then
       call device_neumann_apply_vector(this%msk_d, this%facet_d, &
            x_d, y_d, z_d, &
            this%flux(1)%x_d, this%flux(2)%x_d, this%flux(3)%x_d, &
            this%coef%area_d, this%coef%Xh%lx, &
            size(this%msk), strm)
    end if

  end subroutine neumann_apply_vector_dev

  !> Destructor
  subroutine neumann_free(this)
    class(neumann_t), target, intent(inout) :: this

    call this%free_base()

  end subroutine neumann_free

  !> Finalize by setting the flux.
  subroutine neumann_finalize(this, only_facets)
    class(neumann_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    integer :: i, j

    if (present(only_facets)) then
       if (only_facets .eqv. .false.) then
          call neko_error("For neumann_t, only_facets has to be true.")
       end if
    end if

    call this%finalize_base(.true.)

    ! Allocate flux vectors and assign to initial constant values
    do i = 1,size(this%init_flux_)
       call this%flux(i)%init(this%msk(0))
       this%flux(i) = this%init_flux_(i)
    end do

    this%uniform_0 = .true.

    do i = 1, 3
       this%uniform_0 = abscmp(this%init_flux_(i), 0.0_rp) .and. this%uniform_0
    end do
  end subroutine neumann_finalize

  !> Set the flux using a scalar.
  !! @param flux The desired flux.
  !! @param comp The component to set.
  subroutine neumann_set_flux_scalar(this, flux, comp)
    class(neumann_t), intent(inout) :: this
    real(kind=rp), intent(in) :: flux
    integer, intent(in) :: comp

    if (size(this%flux) .lt. comp) then
       call neko_error("Component index out of bounds in " // &
            "neumann_set_flux_scalar")
    end if

    this%flux(comp) = flux
    ! If we were uniform zero before, and this comp is set to zero, we are still
    ! uniform zero
    this%uniform_0 = abscmp(flux, 0.0_rp) .and. this%uniform_0

  end subroutine neumann_set_flux_scalar

  !> Set a flux component using a vector_t of values.
  !> @param flux The desired flux.
  !> @param comp The component to set.
  subroutine neumann_set_flux_array(this, flux, comp)
    class(neumann_t), intent(inout) :: this
    type(vector_t), intent(in) :: flux
    integer, intent(in) :: comp
    integer :: i

    if (size(this%flux) .lt. comp) then
       call neko_error("Component index out of bounds in " // &
            "neuman_set_flux_array")
    end if

    this%flux(comp) = flux

    ! Once a flux is set explicitly, we no longer assume it is uniform zero.
    this%uniform_0 = .false.

  end subroutine neumann_set_flux_array
end module neumann
