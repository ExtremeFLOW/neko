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
!
!> Implements `deardorff_t`.
module deardorff
  use utils, only : neko_error
  use num_types, only : rp
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use les_model, only : les_model_t
  use json_utils, only : json_get_or_default, json_get, json_get_or_lookup
  use json_module, only : json_file
  use neko_config, only : NEKO_BCKND_DEVICE
  use deardorff_cpu, only : deardorff_compute_cpu
  use deardorff_device, only : deardorff_compute_device
  use registry, only : neko_registry
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Implements the deardorff LES model.
  !! @note Reference DOI: 10.1007/BF00119502
  type, public, extends(les_model_t) :: deardorff_t
     !> Model constant, defaults to 0.10.
     real(kind=rp) :: c_k
     !> The reference temperature
     real(kind=rp) :: T0
     !> The gravitational acceleration
     real(kind=rp) :: g(3)
     !> Temperature field name
     character(len=:), allocatable :: temperature_field_name
     !> TKE field name
     character(len=:), allocatable :: TKE_field_name
     !> Eddy diffusivity for temperature and TKE
     type(field_t), pointer :: temperature_alphat => null()
     type(field_t), pointer :: TKE_alphat => null()
     !> Source term for TKE equation
     type(field_t), pointer :: TKE_source => null()
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => deardorff_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          deardorff_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => deardorff_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => deardorff_compute
  end type deardorff_t

contains
  !> Constructor.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param json A dictionary with parameters.
  subroutine deardorff_init(this, fluid, json)
    class(deardorff_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: temperature_field_name
    character(len=:), allocatable :: TKE_field_name
    character(len=:), allocatable :: nut_name
    character(len=:), allocatable :: temperature_alphat_name, TKE_alphat_name
    character(len=:), allocatable :: TKE_source_name
    character(len=:), allocatable :: delta_type
    real(kind=rp) :: c_k, T0
    real(kind=rp), allocatable :: g(:)
    logical :: if_ext
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "temperature_field", &
         temperature_field_name, "temperature")
    call json_get_or_default(json, "TKE_field", TKE_field_name, "TKE")
    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "temperature_alphat_field", &
         temperature_alphat_name, "temperature_alphat")
    call json_get_or_default(json, "TKE_alphat_field", TKE_alphat_name, &
         "TKE_alphat")
    call json_get_or_default(json, "TKE_source_field", TKE_source_name, &
         "TKE_source")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "c_k", c_k, 0.10_rp)
    call json_get(json, "T0", T0)
    call json_get_or_lookup(json, "g", g)
    if (.not. size(g) == 3) then
       call neko_error("The gravity vector should have 3 components")
    end if
    call json_get_or_default(json, "extrapolation", if_ext, .false.)

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : deardorff'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Delta evaluation : ', delta_type
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'c_k : ', c_k
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'extrapolation : ', if_ext
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call deardorff_init_from_components(this, fluid, c_k, T0, &
         temperature_field_name, TKE_field_name, nut_name, &
         temperature_alphat_name, TKE_alphat_name, TKE_source_name, g, &
         delta_type, if_ext)

    deallocate(temperature_field_name)
    deallocate(TKE_field_name)
    deallocate(nut_name)
    deallocate(temperature_alphat_name)
    deallocate(TKE_alphat_name)
    deallocate(TKE_source_name)
    deallocate(delta_type)
    deallocate(g)

  end subroutine deardorff_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param c_k The deardorff model constant.
  !! @param T0 The reference temperature.
  !! @param temperature_field_name The name of the temperature field.
  !! @param TKE_field_name The name of the TKE field.
  !! @param nut_name The name of the eddy viscosity field.
  !! @param temperature_alphat_name The name of the eddy diffusivity field for
  !! temperature.
  !! @param TKE_alphat_name The name of the eddy diffusivity field for TKE.
  !! @param TKE_source_name The name of the source term in the TKE equation
  !! @param g The gravitational acceleration vector.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether to extrapolate the velocity.
  subroutine deardorff_init_from_components(this, fluid, &
       c_k, T0, temperature_field_name, TKE_field_name, nut_name, &
       temperature_alphat_name, TKE_alphat_name, TKE_source_name, g, &
       delta_type, if_ext)
    class(deardorff_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    real(kind=rp), intent(in) :: c_k, T0
    character(len=*), intent(in) :: temperature_field_name
    character(len=*), intent(in) :: TKE_field_name
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: temperature_alphat_name, TKE_alphat_name
    character(len=*), intent(in) :: TKE_source_name
    real(kind=rp), intent(in) :: g(3)
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext

    call this%free()
    call this%init_base(fluid, nut_name, delta_type, if_ext)

    call neko_registry%add_field(fluid%dm_Xh, &
         trim(temperature_alphat_name), .true.)
    call neko_registry%add_field(fluid%dm_Xh, trim(TKE_alphat_name), .true.)
    call neko_registry%add_field(fluid%dm_Xh, trim(TKE_source_name), .true.)

    this%temperature_alphat => &
         neko_registry%get_field(trim(temperature_alphat_name))
    this%TKE_alphat => neko_registry%get_field(trim(TKE_alphat_name))
    this%TKE_source => neko_registry%get_field(trim(TKE_source_name))
    this%c_k = c_k
    this%T0 = T0
    this%temperature_field_name = temperature_field_name
    this%TKE_field_name = TKE_field_name
    this%g = -g

  end subroutine deardorff_init_from_components

  !> Destructor for the deardorff_t class.
  subroutine deardorff_free(this)
    class(deardorff_t), intent(inout) :: this

    nullify(this%temperature_alphat)
    nullify(this%TKE_alphat)
    nullify(this%TKE_source)

    if (allocated(this%temperature_field_name)) then
       deallocate(this%temperature_field_name)
    end if

    if (allocated(this%TKE_field_name)) then
       deallocate(this%TKE_field_name)
    end if

    call this%free_base()
  end subroutine deardorff_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine deardorff_compute(this, t, tstep)
    class(deardorff_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    type(field_t), pointer :: u, v, w, u_e, v_e, w_e

    if (this%if_ext .eqv. .true.) then
       ! Extrapolate the velocity fields
       associate(ulag => this%ulag, vlag => this%vlag, &
            wlag => this%wlag, ext_bdf => this%ext_bdf)

         u => neko_registry%get_field_by_name("u")
         v => neko_registry%get_field_by_name("v")
         w => neko_registry%get_field_by_name("w")
         u_e => neko_registry%get_field_by_name("u_e")
         v_e => neko_registry%get_field_by_name("v_e")
         w_e => neko_registry%get_field_by_name("w_e")

         call this%sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
              ulag, vlag, wlag, ext_bdf%advection_coeffs, ext_bdf%nadv)

       end associate
    end if

    ! Compute the eddy viscosity field
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call deardorff_compute_device(this%if_ext, t, tstep, this%coef, &
            this%temperature_field_name, this%TKE_field_name, &
            this%nut, this%temperature_alphat, this%TKE_alphat, &
            this%TKE_source, this%delta, this%c_k, this%T0, this%g)
    else
       call deardorff_compute_cpu(this%if_ext, t, tstep, this%coef, &
            this%temperature_field_name, this%TKE_field_name, &
            this%nut, this%temperature_alphat, this%TKE_alphat, &
            this%TKE_source, this%delta, this%c_k, this%T0, this%g)
    end if

  end subroutine deardorff_compute

end module deardorff
