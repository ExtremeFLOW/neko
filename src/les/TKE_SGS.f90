! Copyright (c) 2025, The Neko Authors
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
!> Implements `TKE_SGS_t`.
module TKE_SGS
  use utils, only : neko_error
  use num_types, only : rp
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use les_model, only : les_model_t
  use json_utils, only : json_get_or_default, json_get
  use json_module, only : json_file
  use neko_config, only : NEKO_BCKND_DEVICE
  use TKE_SGS_cpu, only : TKE_SGS_compute_cpu
  ! use TKE_SGS_device, only : TKE_SGS_compute_device
  use field_registry, only : neko_field_registry
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Implements the TKE_SGS LES model.
  !! @note Reference DOI: 10.1175/1520-0493(1963)091<0099:GCEWTP>2.3.CO;2
  type, public, extends(les_model_t) :: TKE_SGS_t
     !> Model constant, defaults to 0.10.
     real(kind=rp) :: c_k
     !> The reference temperature
     real(kind=rp) :: T0
     !> The gravitational acceleration
     real(kind=rp) :: g
     !> Vertical direction
     character(len=:), allocatable :: vertical_dir
     !> Eddy diffusivity for temperature and TKE
     type(field_t), pointer :: nutheta => null()
     type(field_t), pointer :: nue => null()
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => TKE_SGS_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          TKE_SGS_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => TKE_SGS_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => TKE_SGS_compute
  end type TKE_SGS_t

contains
  !> Constructor.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param json A dictionary with parameters.
  subroutine TKE_SGS_init(this, fluid, json)
    class(TKE_SGS_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name, nutheta_name, nue_name
    character(len=:), allocatable :: vertical_dir
    character(len=:), allocatable :: delta_type
    logical :: if_ext
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "nutheta_field", nutheta_name, "nutheta")
    call json_get_or_default(json, "nue_field", nue_name, "nue")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "c_k", this%c_k, 0.10_rp)
    call json_get(json, "T0", this%T0)
    call json_get(json, "g", this%g)
    call json_get_or_default(json, "vertical_direction", vertical_dir, "z")
    this%vertical_dir = trim(vertical_dir)
    call json_get_or_default(json, "extrapolation", if_ext, .false.)

    if (if_ext) then
       call neko_error("Velocity extrapolation is not supported by TKE SGS model.")
    end if

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : TKE_SGS'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Delta evaluation : ', delta_type
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'c_k : ', this%c_k
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'extrapolation : ', if_ext
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call TKE_SGS_init_from_components(this, fluid, nut_name, &
         nutheta_name, nue_name, delta_type, if_ext)

  end subroutine TKE_SGS_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param nut_name The name of the Eddy viscosity field.
  !! @param nutheta_name The name of the Eddy diffusivity field for temperature.
  !! @param nue_name The name of the Eddy diffusivity field for TKE.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether trapolate the velocity.
  subroutine TKE_SGS_init_from_components(this, fluid, &
       nut_name, nutheta_name, nue_name, delta_type, if_ext)
    class(TKE_SGS_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    real(kind=rp) :: c_k
    character(len=*), intent(in) :: nut_name, nutheta_name, nue_name
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext

    call this%free()

    call this%init_base(fluid, nut_name, delta_type, if_ext)
    this%nutheta => neko_field_registry%get_field(trim(nutheta_name))
    this%nue => neko_field_registry%get_field(trim(nue_name))

  end subroutine TKE_SGS_init_from_components

  !> Destructor for the les_model_t (base) class.
  subroutine TKE_SGS_free(this)
    class(TKE_SGS_t), intent(inout) :: this

    call this%free_base()
  end subroutine TKE_SGS_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine TKE_SGS_compute(this, t, tstep)
    class(TKE_SGS_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Compute the eddy viscosity field
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("TKE SGS model is not implemented on device yet.")
    else
       call TKE_SGS_compute_cpu(t, tstep, this%coef, &
            this%nut, this%nutheta, this%nue, &
            this%delta, this%c_k, this%T0, this%g, &
            this%vertical_dir)
    end if

  end subroutine TKE_SGS_compute

end module TKE_SGS
