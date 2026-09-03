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
!> Implements `avm_t`.
module artificial_viscosity_model
  use num_types, only : rp
  use case, only : case_t
  use field, only : field_t
  use json_module, only : json_file
  use registry, only : neko_registry
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Base abstract type for artificial viscosity models.
  type, abstract, public :: avm_t
     !> Artificial viscosity.
     type(field_t), pointer :: reg_coeff => null()
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> SEM dofmap.
     type(dofmap_t), pointer :: dof => null()
     !> The size of the `les_model_registry`
     integer :: les_model_registry_size = 0
   contains
     !> Constructor for the avm_t (base) class.
     procedure, pass(this) :: init_base => avm_init_base
     !> Destructor for the avm_t (base) class.
     procedure, pass(this) :: free_base => avm_free_base
     !> The common constructor.
     procedure(avm_init), pass(this), deferred :: init
     !> Destructor.
     procedure(avm_free), pass(this), deferred :: free
     !> Perform artificial viscosity related computations before the time step.
     procedure(avm_preprocess), pass(this), deferred :: preprocess
     !> Perform artificial viscosity related computations after the time step.
     procedure(avm_compute), pass(this), deferred :: compute
     !> Restore model state after loading a checkpoint.
     procedure, pass(this) :: restart => avm_restart
  end type avm_t

  abstract interface
     !> Perform artificial viscosity related computations before the time step.
     !! @param time The time state.
     subroutine avm_preprocess(this, time)
       import avm_t, time_state_t
       class(avm_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine avm_preprocess
  end interface

  abstract interface
     !> Perform artificial viscosity related computations after the time step.
     !! @param time The time state.
     subroutine avm_compute(this, time)
       import avm_t, time_state_t
       class(avm_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine avm_compute
  end interface

  abstract interface
     !> Common constructor.
     !! @param case The case_t object.
     !! @param json A dictionary with parameters.
     subroutine avm_init(this, case, json)
       import avm_t, json_file, case_t
       class(avm_t), intent(inout) :: this
       class(case_t), intent(inout), target :: case
       type(json_file), intent(inout) :: json
     end subroutine avm_init
  end interface

  abstract interface
     !> Destructor.
     subroutine avm_free(this)
       import avm_t
       class(avm_t), intent(inout) :: this
     end subroutine avm_free
  end interface

  interface
     !> Artificial viscosity model allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the artificial viscosity model.
     module subroutine avm_allocator(object, type_name)
       class(avm_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
     end subroutine avm_allocator
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for an object allocator.
  !! Implemented in the user modules, should allocate the `obj` to the custom
  !! user type.
  abstract interface
     subroutine avm_allocate(obj)
       import avm_t
       class(avm_t), allocatable, intent(inout) :: obj
     end subroutine avm_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_avm(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(avm_allocate), pointer, intent(in) :: allocator
     end subroutine register_avm
  end interface

  ! A name-allocator pair for user-defined types. A helper type to define a
  ! registry of custom allocators.
  type allocator_entry
     character(len=20) :: type_name
     procedure(avm_allocate), pointer, nopass :: allocator
  end type allocator_entry

  !> Registry of artificial viscosity model allocators for user-defined types
  type(allocator_entry), allocatable :: avm_registry(:)

  !> The size of the `avm_registry`
  integer :: avm_registry_size = 0



  interface
     !> Artificial viscosity model factory. Both allocates and initializes
     !! the object.
     !! @param object The object to be allocated.
     !! @param type_name The name of the artificial viscosity model.
     !! @param case The case_t object.
     !! @param dofmap SEM map of degrees of freedom.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     module subroutine avm_factory(object, type_name, case, json)
       class(avm_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
       class(case_t), intent(inout), target :: case
       type(json_file), intent(inout) :: json
     end subroutine avm_factory
  end interface

  public :: avm_factory, avm_allocator, register_avm, &
       avm_allocate


contains
  !> Constructor for the avm_t (base) class.
  !! @param dof Map of degrees of freedom.
  !! @param coef The SEM coefficients.
  !! @param reg_coeff_name The name of the artificial viscosity field.
  subroutine avm_init_base(this, dof, coef, reg_coeff_name)
    class(avm_t), intent(inout) :: this
    class(dofmap_t), intent(in), target :: dof
    class(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: reg_coeff_name

    call this%free_base()

    call neko_registry%add_field(dof, trim(reg_coeff_name), .true.)

    this%reg_coeff => neko_registry%get_field(trim(reg_coeff_name))
    this%coef => coef
    this%dof => dof

  end subroutine avm_init_base

  !> Restore model state after loading a checkpoint.
  !! Models with additional history should override this method.
  !! @param time The restored time state.
  subroutine avm_restart(this, time)
    class(avm_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%compute(time)

  end subroutine avm_restart

  !> Destructor for the avm_t (base) class.
  subroutine avm_free_base(this)
    class(avm_t), intent(inout) :: this

    nullify(this%reg_coeff)
    nullify(this%coef)
    nullify(this%dof)

  end subroutine avm_free_base

end module artificial_viscosity_model
