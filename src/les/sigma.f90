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
!> Implements `sigma_t`.
module sigma
  use num_types, only : rp
  use field, only : field_t
  use les_model, only : les_model_t
  use dofmap , only : dofmap_t
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use sigma_cpu, only : sigma_compute_cpu
  use coefs, only : coef_t
  implicit none
  private

  !> Implements the Sigma LES model.
  !! @note Reference DOI: 10.1063/1.3623274
  type, public, extends(les_model_t) :: sigma_t
     !> Model constant, default to 1.35.
     real(kind=rp) :: c
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => sigma_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => sigma_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => sigma_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => sigma_compute
  end type sigma_t

contains
  !> Constructor.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine sigma_init(this, dofmap, coef, json)
    class(sigma_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name
    real(kind=rp) :: c
    character(len=:), allocatable :: delta_type

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    ! Based on  C = 1.35 as default values
    call json_get_or_default(json, "c", c, 1.35_rp)

    call sigma_init_from_components(this, dofmap, coef, c, nut_name, delta_type)
  end subroutine sigma_init

  !> Constructor from components.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param c The model constant.
  !! @param nut_name The name of the SGS viscosity field.
  subroutine sigma_init_from_components(this, dofmap, coef, c, nut_name, &
       delta_type)
    class(sigma_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: c
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: delta_type

    call this%free()

    call this%init_base(dofmap, coef, nut_name, delta_type)

    this%c = c

  end subroutine sigma_init_from_components

  !> Destructor for the les_model_t (base) class.
  subroutine sigma_free(this)
    class(sigma_t), intent(inout) :: this

    call this%free_base()
  end subroutine sigma_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine sigma_compute(this, t, tstep)
    class(sigma_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    if (NEKO_BCKND_DEVICE .eq. 1) then
        call neko_error("Sigma model not implemented on accelarators.")
    else
        call sigma_compute_cpu(t, tstep, this%coef, this%nut, this%delta,&
                                this%c)
    end if

  end subroutine sigma_compute

end module sigma
