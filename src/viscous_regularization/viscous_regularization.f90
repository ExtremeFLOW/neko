! Copyright (c) 2025-2026, The Neko Authors
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
!> Defines the base interface for viscous regularization methods.
module viscous_regularization
  use num_types, only : rp
  use json_module, only : json_file
  use field, only : field_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Base abstract type for viscous regularization methods.
  type, abstract, public :: viscous_regularization_t
     !> Spatially varying regularization coefficient.
     type(field_t), pointer :: reg_coeff => null()
     !> Name of the regularization coefficient field.
     character(len=:), allocatable :: reg_coeff_name
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> SEM map of degrees of freedom.
     type(dofmap_t), pointer :: dof => null()
   contains
     !> Initialize state shared by all regularization methods.
     procedure, pass(this) :: init_base => viscous_regularization_init_base
     !> Free state shared by all regularization methods.
     procedure, pass(this) :: free_base => viscous_regularization_free_base
     !> Initialize a regularization method.
     procedure(reg_init), pass(this), deferred :: init
     !> Free a regularization method.
     procedure(reg_free), pass(this), deferred :: free
     !> Update the effective viscosity field.
     procedure(reg_update), pass(this), deferred :: update
  end type viscous_regularization_t

  abstract interface
     !> Initialize a viscous regularization method.
     !! @param this The viscous regularization method.
     !! @param json The regularization configuration.
     !! @param coef The SEM coefficients.
     !! @param dof The SEM map of degrees of freedom.
     subroutine reg_init(this, json, coef, dof)
       import viscous_regularization_t, json_file, coef_t, dofmap_t, field_t
       class(viscous_regularization_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in), target :: coef
       type(dofmap_t), intent(in), target :: dof
     end subroutine reg_init
  end interface

  abstract interface
     !> Free a viscous regularization method.
     !! @param this The viscous regularization method.
     subroutine reg_free(this)
       import viscous_regularization_t
       class(viscous_regularization_t), intent(inout) :: this
     end subroutine reg_free
  end interface

  abstract interface
     !> Update an effective viscosity field.
     !! @param this The viscous regularization method.
     !! @param effective_visc The effective viscosity field to update.
     !! @param mu The optional physical viscosity field.
     subroutine reg_update(this, effective_visc, mu)
       import viscous_regularization_t, field_t
       class(viscous_regularization_t), intent(inout) :: this
       type(field_t), intent(inout) :: effective_visc
       type(field_t), intent(in), optional :: mu
     end subroutine reg_update
  end interface

  interface
     !> Allocate and initialize a viscous regularization method.
     !! @param object The regularization object to allocate.
     !! @param type_name The regularization type name.
     !! @param json The regularization configuration.
     !! @param coef The SEM coefficients.
     !! @param dof The SEM map of degrees of freedom.
     module subroutine viscous_regularization_factory(object, type_name, json, &
          coef, dof)
       class(viscous_regularization_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in), target :: coef
       type(dofmap_t), intent(in), target :: dof
     end subroutine viscous_regularization_factory
  end interface

  public :: viscous_regularization_factory

contains

  !> Initialize state shared by all viscous regularization methods.
  !! @param this The viscous regularization method.
  !! @param json The regularization configuration.
  !! @param coef The SEM coefficients.
  !! @param dof The SEM map of degrees of freedom.
  subroutine viscous_regularization_init_base(this, json, coef, dof)
    class(viscous_regularization_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof

    this%coef => coef
    this%dof => dof

  end subroutine viscous_regularization_init_base

  !> Free state shared by all viscous regularization methods.
  !! @param this The viscous regularization method.
  subroutine viscous_regularization_free_base(this)
    class(viscous_regularization_t), intent(inout) :: this

    nullify(this%reg_coeff)
    if (allocated(this%reg_coeff_name)) deallocate(this%reg_coeff_name)
    nullify(this%coef)
    nullify(this%dof)

  end subroutine viscous_regularization_free_base

end module viscous_regularization
