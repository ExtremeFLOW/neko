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
!> Implements `avm_model_t`.
module avm_model
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

  !> Base abstract type for avm models based on the Boussinesq approximation.
  type, abstract, public :: avm_model_t
     !> Artificial viscosity.
     type(field_t), pointer :: reg_coeff => null()
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> SEM dofmap.
     type(dofmap_t), pointer :: dof => null()
     !> The size of the `les_model_registry`
     integer :: les_model_registry_size = 0
   contains
     !> Constructor for the avm_model_t (base) class.
     procedure, pass(this) :: init_base => avm_model_init_base
     !> Destructor for the avm_model_t (base) class.
     procedure, pass(this) :: free_base => avm_model_free_base
     !> The common constructor.
     procedure(avm_model_init), pass(this), deferred :: init
     !> Destructor.
     procedure(avm_model_free), pass(this), deferred :: free
     !> Perform artificial viscosity related computations before the time step.
     procedure(avm_model_preprocess), pass(this), deferred :: preprocess
     !> Perform artificial viscosity related computations after the time step.
     procedure(avm_model_compute), pass(this), deferred :: compute
  end type avm_model_t

  abstract interface
     !> Perform artificial viscosity related computations before the time step.
     !! @param time The time state.
     subroutine avm_model_preprocess(this, time)
       import avm_model_t, time_state_t
       class(avm_model_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine avm_model_preprocess
  end interface

  abstract interface
     !> Perform artificial viscosity related computations after the time step.
     !! @param time The time state.
     subroutine avm_model_compute(this, time)
       import avm_model_t, time_state_t
       class(avm_model_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine avm_model_compute
  end interface

  abstract interface
     !> Common constructor.
     !! @param case The case_t object.
     !! @param json A dictionary with parameters.
     subroutine avm_model_init(this, case, json)
       import avm_model_t, json_file, case_t
       class(avm_model_t), intent(inout) :: this
       class(case_t), intent(inout), target :: case
       type(json_file), intent(inout) :: json
     end subroutine avm_model_init
  end interface

  abstract interface
     !> Destructor.
     subroutine avm_model_free(this)
       import avm_model_t
       class(avm_model_t), intent(inout) :: this
     end subroutine avm_model_free
  end interface

  interface
     !> AVM model allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the AVM model.
     module subroutine avm_model_allocator(object, type_name)
       class(avm_model_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
     end subroutine avm_model_allocator
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for an object allocator.
  !! Implemented in the user modules, should allocate the `obj` to the custom
  !! user type.
  abstract interface
     subroutine avm_model_allocate(obj)
       import avm_model_t
       class(avm_model_t), allocatable, intent(inout) :: obj
     end subroutine avm_model_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_avm_model(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(avm_model_allocate), pointer, intent(in) :: allocator
     end subroutine register_avm_model
  end interface

  ! A name-allocator pair for user-defined types. A helper type to define a
  ! registry of custom allocators.
  type allocator_entry
     character(len=20) :: type_name
     procedure(avm_model_allocate), pointer, nopass :: allocator
  end type allocator_entry

  !> Registry of artificial viscosity model allocators for user-defined types
  type(allocator_entry), allocatable :: avm_model_registry(:)

  !> The size of the `avm_model_registry`
  integer :: avm_model_registry_size = 0



  interface
     !> avm model factory. Both allocates and initializes the object.
     !! @param object The object to be allocated.
     !! @param type_name The name of the avm model.
     !! @param case The case_t object.
     !! @param dofmap SEM map of degrees of freedom.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     module subroutine avm_model_factory(object, type_name, case, json)
       class(avm_model_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type_name
       class(case_t), intent(inout), target :: case
       type(json_file), intent(inout) :: json
     end subroutine avm_model_factory
  end interface

  public :: avm_model_factory, avm_model_allocator, register_avm_model, &
       avm_model_allocate


contains
  !> Constructor for the avm_model_t (base) class.
  !! @param case The case_t object.
  !! @param reg_coeff_name The name of the artificial viscosity field.
  subroutine avm_model_init_base(this, case, reg_coeff_name)
    class(avm_model_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    character(len=*), intent(in) :: reg_coeff_name

    associate(dofmap => case%fluid%dm_Xh, &
         coef => case%fluid%c_Xh)

      call neko_registry%add_field(dofmap, trim(reg_coeff_name), .true.)
      call neko_registry%add_field(dofmap, "avm_delta", .true.)
      this%reg_coeff => neko_registry%get_field(trim(reg_coeff_name))
      this%coef => case%fluid%c_Xh
      this%dof => case%fluid%dm_Xh

    end associate
  end subroutine avm_model_init_base

  !> Destructor for the avm_model_t (base) class.
  subroutine avm_model_free_base(this)
    class(avm_model_t), intent(inout) :: this

    nullify(this%reg_coeff)
    nullify(this%coef)

  end subroutine avm_model_free_base

end module avm_model
