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
!
!> Implements `artificial_viscosity_t`.
module artificial_viscosity
  use num_types, only : rp
  use viscous_regularization, only : viscous_regularization_t
  use registry, only : neko_registry
  use field, only : field_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use coefs, only : coef_t
  use field_math, only : field_copy, field_add2
  use dofmap, only : dofmap_t
  use utils, only : neko_error
  implicit none
  private

  !> The handler for the artificial viscosity.
  type, public, extends(viscous_regularization_t) :: artificial_viscosity_t
   contains
     procedure, pass(this) :: init => artificial_viscosity_init_from_json
     procedure, pass(this) :: free => artificial_viscosity_free
     procedure, pass(this) :: update => artificial_viscosity_update
  end type artificial_viscosity_t

contains
  !> Constructor
  subroutine artificial_viscosity_init_from_json(this, json, coef, dof)
    class(artificial_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof

    call this%init_base(json, coef, dof)

    call json_get(json, "reg_coeff_name", &
         this%reg_coeff_name)

    call neko_registry%add_field(dof, this%reg_coeff_name)
    this%reg_coeff => neko_registry%get_field(this%reg_coeff_name)

  end subroutine artificial_viscosity_init_from_json

  !> Destructor
  subroutine artificial_viscosity_free(this)
    class(artificial_viscosity_t), intent(inout) :: this

    nullify(this%reg_coeff)
    call this%free_base()

  end subroutine artificial_viscosity_free

  !> Update of the artificial viscosity
  subroutine artificial_viscosity_update(this, effective_visc, mu)
    class(artificial_viscosity_t), intent(inout) :: this
    type(field_t), intent(inout) :: effective_visc
    type(field_t), intent(in), optional :: mu

    call field_copy(effective_visc, this%reg_coeff, this%coef%dof%size())
    if (present(mu)) then
       call field_add2(effective_visc, mu, this%coef%dof%size())
    end if

  end subroutine artificial_viscosity_update

end module artificial_viscosity
