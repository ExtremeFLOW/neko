! Copyright (c) 2024, The Neko Authors
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
!> Implements `wale_t`.
module wale
  use num_types, only : rp
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use les_model, only : les_model_t
  use json_utils, only : json_get_or_default
  use json_module, only : json_file
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use wale_cpu, only : wale_compute_cpu
  use wale_device, only : wale_compute_device
  use field_registry, only : neko_field_registry
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Implements the Wale LES model.
  !! @note Reference DOI: 10.1023/A:1009995426001
  type, public, extends(les_model_t) :: wale_t
     !> Model constant, defaults to 0.55.
     real(kind=rp) :: c_w
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => wale_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          wale_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => wale_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => wale_compute
  end type wale_t

contains
  !> Constructor.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param json A dictionary with parameters.
  subroutine wale_init(this, fluid, json)
    class(wale_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name
    real(kind=rp) :: c_w
    character(len=:), allocatable :: delta_type
    logical :: if_ext
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "c_w", c_w, 0.55_rp)
    call json_get_or_default(json, "extrapolation", if_ext, .false.)

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : Wale'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Delta evaluation : ', delta_type
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'c_w : ', c_w
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'extrapolation : ', if_ext
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call wale_init_from_components(this, fluid, c_w, nut_name, &
         delta_type, if_ext)

  end subroutine wale_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param c_w The model constant.
  !! @param nut_name The name of the SGS viscosity field.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether trapolate the velocity.
  subroutine wale_init_from_components(this, fluid, c_w, &
       nut_name, delta_type, if_ext)
    class(wale_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    real(kind=rp) :: c_w
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext

    call this%free()

    call this%init_base(fluid, nut_name, delta_type, if_ext)
    this%c_w = c_w

  end subroutine wale_init_from_components

  !> Destructor for the les_model_t (base) class.
  subroutine wale_free(this)
    class(wale_t), intent(inout) :: this

    call this%free_base()
  end subroutine wale_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine wale_compute(this, t, tstep)
    class(wale_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    type(field_t), pointer :: u, v, w, u_e, v_e, w_e

    if (this%if_ext .eqv. .true.) then
       ! Extrapolate the velocity fields
       associate(ulag => this%ulag, vlag => this%vlag, &
            wlag => this%wlag, ext_bdf => this%ext_bdf)

         u => neko_field_registry%get_field_by_name("u")
         v => neko_field_registry%get_field_by_name("v")
         w => neko_field_registry%get_field_by_name("w")
         u_e => neko_field_registry%get_field_by_name("u_e")
         v_e => neko_field_registry%get_field_by_name("v_e")
         w_e => neko_field_registry%get_field_by_name("w_e")

         call this%sumab%compute_fluid(u_e, v_e, w_e, u, v, w, &
              ulag, vlag, wlag, ext_bdf%advection_coeffs, ext_bdf%nadv)

       end associate
    end if

    ! Compute the eddy viscosity field
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call wale_compute_device(this%if_ext, t, tstep, this%coef, &
            this%nut, this%delta, this%c_w)
    else
       call wale_compute_cpu(this%if_ext, t, tstep, this%coef, &
            this%nut, this%delta, this%c_w)
    end if

  end subroutine wale_compute

end module wale
