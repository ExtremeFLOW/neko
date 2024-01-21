! Copyright (c) 2023, The Neko Authors
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
!> Implements `les_model_t`.
module les_model
  use num_types, only : rp
  use field, only : field_t, field_ptr_t
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  implicit none
  private

  !> Base abstract type for LES models based on the Boussinesq approximation.
  type, abstract, public :: les_model_t
     !> Subgrid kinematic viscosity.
     type(field_t), pointer :: nut => null()

     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
   contains
     !> Constructor for the les_model_t (base) class.
     procedure, pass(this) :: init_base => les_model_init_base
     !> Destructor for the les_model_t (base) class.
     procedure, pass(this) :: free_base => les_model_free_base
     !> The common constructor.
     procedure(les_model_init), pass(this), deferred :: init
     !> Destructor.
     procedure(les_model_free), pass(this), deferred :: free
     !> Compute eddy viscosity.
     procedure(les_model_compute), pass(this), deferred :: compute
  end type les_model_t

  abstract interface
     !> Compute eddy viscosity.
     !! @param t The time value.
     !! @param tstep The current time-step.
     subroutine les_model_compute(this, t, tstep)
       import les_model_t, rp
       class(les_model_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
     end subroutine les_model_compute
  end interface

  abstract interface
     !> Common constructor.
     !! @param dofmap SEM map of degrees of freedom.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     subroutine les_model_init(this, dofmap, coef, json)
       import les_model_t, json_file, dofmap_t, coef_t
       class(les_model_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       type(dofmap_t), intent(in) :: dofmap
       type(json_file), intent(inout) :: json
     end subroutine les_model_init
  end interface

  abstract interface
     !> Destructor.
     subroutine les_model_free(this)
       import les_model_t
       class(les_model_t), intent(inout) :: this
     end subroutine les_model_free
  end interface

contains
  !> Constructor for the les_model_t (base) class.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param nu_name The name of the turbulent viscosity field.
  subroutine les_model_init_base(this, dofmap, coef, nut_name)
    class(les_model_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: nut_name

    if (.not. neko_field_registry%field_exists(trim(nut_name))) then
       call neko_field_registry%add_field(dofmap, trim(nut_name))
    end if
    this%nut => neko_field_registry%get_field(trim(nut_name))
    this%coef => coef
  end subroutine les_model_init_base

  !> Destructor for the les_model_t (base) class.
  subroutine les_model_free_base(this)
    class(les_model_t), intent(inout) :: this

    nullify(this%nut)
    nullify(this%coef)
  end subroutine les_model_free_base

end module les_model