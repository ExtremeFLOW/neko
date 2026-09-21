! Copyright (c) 2026, The Neko Authors
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
!> Template for a user-defined scalar Pn/Pn boundary condition.
module scalar_pnpn_bc_template
  use bc, only : bc_t
  use coefs, only : coef_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  use json_module, only : json_file
  use num_types, only : rp
  use scalar_pnpn, only : scalar_pnpn_bc_allocate, &
       register_scalar_pnpn_bc
  use time_state, only : time_state_t
  use utils, only : neko_error
  implicit none
  private

  !> A user-defined scalar Pn/Pn boundary condition.
  type, extends(bc_t) :: scalar_pnpn_bc_template_t
   contains
     !> Apply the boundary condition to a scalar field on the CPU.
     procedure :: apply_scalar => scalar_pnpn_bc_template_apply_scalar
     !> Apply the boundary condition to a vector field on the CPU.
     procedure :: apply_vector => scalar_pnpn_bc_template_apply_vector
     !> Apply the boundary condition to a scalar field on the device.
     procedure :: apply_scalar_dev => &
          scalar_pnpn_bc_template_apply_scalar_dev
     !> Apply the boundary condition to a vector field on the device.
     procedure :: apply_vector_dev => &
          scalar_pnpn_bc_template_apply_vector_dev
     !> Constructor from JSON.
     procedure :: init => scalar_pnpn_bc_template_init
     !> Destructor.
     procedure :: free => scalar_pnpn_bc_template_free
     !> Finalize the boundary condition.
     procedure :: finalize => scalar_pnpn_bc_template_finalize
  end type scalar_pnpn_bc_template_t

  public :: scalar_pnpn_bc_template_register_types

contains

  !> Constructor from JSON.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine scalar_pnpn_bc_template_init(this, coef, json)
    class(scalar_pnpn_bc_template_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%free()
    call this%init_base(coef)
    this%bc_type = -1

    ! TODO: Read parameters from json and set this%bc_type to the appropriate
    !       BC_* classification from the bc module.
  end subroutine scalar_pnpn_bc_template_init

  !> Apply the boundary condition to a scalar field on the CPU.
  !! @param[inout] x The scalar field values.
  !! @param[in] n The number of field values.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether strong or weak conditions are being applied.
  subroutine scalar_pnpn_bc_template_apply_scalar(this, x, n, time, strong)
    class(scalar_pnpn_bc_template_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    ! TODO: Apply the condition to the entries in this%msk.
    call neko_error("scalar_pnpn_bc_template: apply_scalar is not implemented")
  end subroutine scalar_pnpn_bc_template_apply_scalar

  !> Apply the boundary condition to a vector field on the CPU.
  !! @param[inout] x The first vector component.
  !! @param[inout] y The second vector component.
  !! @param[inout] z The third vector component.
  !! @param[in] n The number of values in each component.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether strong or weak conditions are being applied.
  subroutine scalar_pnpn_bc_template_apply_vector(this, x, y, z, n, time, &
       strong)
    class(scalar_pnpn_bc_template_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n), y(n), z(n)
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    ! Required by bc_t, but not used by the scalar scheme.
    call neko_error("scalar_pnpn_bc_template: apply_vector is not implemented")
  end subroutine scalar_pnpn_bc_template_apply_vector

  !> Apply the boundary condition to a scalar field on the device.
  !! @param[inout] x_d Device pointer to the scalar field values.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether strong or weak conditions are being applied.
  !! @param[inout] strm The device stream.
  subroutine scalar_pnpn_bc_template_apply_scalar_dev(this, x_d, time, &
       strong, strm)
    class(scalar_pnpn_bc_template_t), target, intent(inout) :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    ! TODO: Launch the device implementation of the scalar condition.
    call neko_error( &
         "scalar_pnpn_bc_template: apply_scalar_dev is not implemented")
  end subroutine scalar_pnpn_bc_template_apply_scalar_dev

  !> Apply the boundary condition to a vector field on the device.
  !! @param[inout] x_d Device pointer to the first vector component.
  !! @param[inout] y_d Device pointer to the second vector component.
  !! @param[inout] z_d Device pointer to the third vector component.
  !! @param[in] time The current time state.
  !! @param[in] strong Whether strong or weak conditions are being applied.
  !! @param[inout] strm The device stream.
  subroutine scalar_pnpn_bc_template_apply_vector_dev(this, x_d, y_d, z_d, &
       time, strong, strm)
    class(scalar_pnpn_bc_template_t), target, intent(inout) :: this
    type(c_ptr), intent(inout) :: x_d, y_d, z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    ! Required by bc_t, but not used by the scalar scheme.
    call neko_error( &
         "scalar_pnpn_bc_template: apply_vector_dev is not implemented")
  end subroutine scalar_pnpn_bc_template_apply_vector_dev

  !> Destructor.
  subroutine scalar_pnpn_bc_template_free(this)
    class(scalar_pnpn_bc_template_t), target, intent(inout) :: this

    call this%free_base()
  end subroutine scalar_pnpn_bc_template_free

  !> Finalize the boundary condition.
  subroutine scalar_pnpn_bc_template_finalize(this)
    class(scalar_pnpn_bc_template_t), target, intent(inout) :: this

    call this%finalize_base()
  end subroutine scalar_pnpn_bc_template_finalize

  !> Register the custom boundary condition type with scalar Pn/Pn.
  subroutine scalar_pnpn_bc_template_register_types()
    procedure(scalar_pnpn_bc_allocate), pointer :: allocator

    allocator => scalar_pnpn_bc_template_allocate
    call register_scalar_pnpn_bc("scalar_pnpn_bc_template", allocator)
  end subroutine scalar_pnpn_bc_template_register_types

  !> Allocate the custom boundary condition.
  !! @param[inout] obj The polymorphic boundary condition to allocate.
  subroutine scalar_pnpn_bc_template_allocate(obj)
    class(bc_t), pointer, intent(inout) :: obj

    allocate(scalar_pnpn_bc_template_t::obj)
  end subroutine scalar_pnpn_bc_template_allocate
end module scalar_pnpn_bc_template
