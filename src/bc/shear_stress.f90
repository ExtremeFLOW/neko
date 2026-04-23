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
!> Defines a shear stress boundary condition for a vector field.
!! Maintainer: Timofey Mukha.
module shear_stress
  use num_types, only : rp
  use bc, only : BC_TYPES
  use mixed_bc, only : mixed_bc_t
  use device_constrain_mixed_bc, only : device_constrain_mixed_bc_zero
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error, nonlinear_index
  use coefs, only : coef_t
  use neumann, only : neumann_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_lookup
  use vector, only : vector_t
  use time_state, only : time_state_t
  implicit none
  private

  !> A shear stress boundary condition.
  type, public, extends(mixed_bc_t) :: shear_stress_t
     !> Neumann condition for the x direction.
     type(neumann_t) :: neumann_x
     !> Neumann condition for the y direction.
     type(neumann_t) :: neumann_y
     !> Neumann condition for the z direction.
     type(neumann_t) :: neumann_z
   contains
     procedure, pass(this) :: apply_scalar => shear_stress_apply_scalar
     procedure, pass(this) :: apply_vector => shear_stress_apply_vector
     procedure, pass(this) :: apply_scalar_dev => shear_stress_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => shear_stress_apply_vector_dev
     !> Constructor.
     procedure, pass(this) :: init => shear_stress_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          shear_stress_init_from_components
     procedure, pass(this) :: set_stress_scalar => &
          shear_stress_set_stress_scalar
     procedure, pass(this) :: set_stress_array => &
          shear_stress_set_stress_array
     !> Set the shear stress to apply.
     generic :: set_stress => set_stress_scalar, set_stress_array
     !> Destructor.
     procedure, pass(this) :: free => shear_stress_free
     !> Finalize the construction.
     procedure, pass(this) :: finalize => shear_stress_finalize
  end type shear_stress_t

contains

  !> Apply shear stress for a scalar field @a x.
  subroutine shear_stress_apply_scalar(this, x, n, time, strong)
    class(shear_stress_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m, k, facet
    ! Store non-linear index
    integer :: idx(4)

    call neko_error("The shear stress bc is not applicable to scalar fields.")

  end subroutine shear_stress_apply_scalar

  !> Boundary condition apply for a generic shear_stress condition
  !! to vectors @a x, @a y and @a z
  subroutine shear_stress_apply_vector(this, x, y, z, n, time, strong)
    class(shear_stress_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_
    integer :: i, m, k
    real(kind=rp) :: normal(3), u_n

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       m = this%resolved_msk%size()

       do i = 1, m
          k = this%resolved_msk%get(i)
          normal = this%n%x(:,i)
          u_n = x(k) * normal(1) + y(k) * normal(2) + z(k) * normal(3)

          x(k) = x(k) - u_n * normal(1)
          y(k) = y(k) - u_n * normal(2)
          z(k) = z(k) - u_n * normal(3)
       end do
    else
       call this%neumann_x%apply_scalar(x, n, strong = .false.)
       call this%neumann_y%apply_scalar(y, n, strong = .false.)
       call this%neumann_z%apply_scalar(z, n, strong = .false.)
    end if

  end subroutine shear_stress_apply_vector

  !> Boundary condition apply for a generic shear_stress condition
  !! to a vector @a x (device version)
  subroutine shear_stress_apply_scalar_dev(this, x_d, time, strong, strm)
    class(shear_stress_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call neko_error("The shear stress bc is not applicable to scalar fields.")

  end subroutine shear_stress_apply_scalar_dev

  !> Boundary condition apply for a generic shear_stress condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine shear_stress_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(shear_stress_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_
    integer :: m

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       m = this%resolved_msk%size()
       if (m .gt. 0) then
          call device_constrain_mixed_bc_zero(this%resolved_msk%get_d(), &
               x_d, y_d, z_d, 1, 0, 0, this%n%x_d, this%t1%x_d, &
               this%t2%x_d, m, strm)
       end if
    else
       call this%neumann_x%apply_scalar_dev(x_d, strong = .false., strm = strm)
       call this%neumann_y%apply_scalar_dev(y_d, strong = .false., strm = strm)
       call this%neumann_z%apply_scalar_dev(z_d, strong = .false., strm = strm)
    end if

  end subroutine shear_stress_apply_vector_dev

  !> Constructor.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine shear_stress_init(this, coef, json)
    class(shear_stress_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) ::json
    real(kind=rp), allocatable :: value(:)

    call json_get_or_lookup(json, 'value', value)

    if (size(value) .ne. 3) then
       call neko_error ("The shear stress vector provided for the shear stress &
       & boundary condition should have 3 components.")
    end if

    call this%init_from_components(coef, value)
  end subroutine shear_stress_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] value The value of the shear stress to apply.
  subroutine shear_stress_init_from_components(this, coef, value)
    class(shear_stress_t), target, intent(inout) :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: value(3)

    call this%init_base(coef)
    this%bc_type = BC_TYPES%MIXED_CONSTRAINS_NORMAL

    call this%neumann_x%free()
    call this%neumann_y%free()
    call this%neumann_z%free()

    call this%neumann_x%init_from_components(this%coef, value(1))
    call this%neumann_y%init_from_components(this%coef, value(2))
    call this%neumann_z%init_from_components(this%coef, value(3))

  end subroutine shear_stress_init_from_components

  subroutine shear_stress_finalize(this)
    class(shear_stress_t), target, intent(inout) :: this
    call this%finalize_base()

    call this%neumann_x%mark_facets(this%marked_facet)
    call this%neumann_y%mark_facets(this%marked_facet)
    call this%neumann_z%mark_facets(this%marked_facet)

    call this%neumann_x%finalize()
    call this%neumann_y%finalize()
    call this%neumann_z%finalize()

  end subroutine shear_stress_finalize

  !> Set the value of the shear stress vector using 3 scalars.
  subroutine shear_stress_set_stress_scalar(this, tau_x, tau_y, tau_z)
    class(shear_stress_t), intent(inout) :: this
    real(kind=rp), intent(in) :: tau_x
    real(kind=rp), intent(in) :: tau_y
    real(kind=rp), intent(in) :: tau_z

    ! Calls finalize and allocates the flux arrays
    call this%neumann_x%set_flux(tau_x, 1)
    call this%neumann_y%set_flux(tau_y, 1)
    call this%neumann_z%set_flux(tau_z, 1)
  end subroutine shear_stress_set_stress_scalar

  !> Set the shear stress components.
  !! @param tau_x The x component of the stress.
  !! @param tau_y The y component of the stress.
  !! @param tau_z The z component of the stress.
  subroutine shear_stress_set_stress_array(this, tau_x, tau_y, tau_z)
    class(shear_stress_t), intent(inout) :: this
    type(vector_t), intent(in) :: tau_x
    type(vector_t), intent(in) :: tau_y
    type(vector_t), intent(in) :: tau_z

    call this%neumann_x%set_flux(tau_x, 1)
    call this%neumann_y%set_flux(tau_y, 1)
    call this%neumann_z%set_flux(tau_z, 1)

  end subroutine shear_stress_set_stress_array

  !> Destructor.
  subroutine shear_stress_free(this)
    class(shear_stress_t), target, intent(inout) :: this
    call this%free_mixed()

    call this%neumann_x%free
    call this%neumann_y%free
    call this%neumann_z%free

  end subroutine shear_stress_free
end module shear_stress
