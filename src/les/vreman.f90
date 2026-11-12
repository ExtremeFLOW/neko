! Copyright (c) 2023-2026, The Neko Authors
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
!> Implements `vreman_t`.
module vreman
  use num_types, only : rp
  use les_model, only : les_model_t
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use json_utils, only : json_get_or_default, json_get, json_get_or_lookup, &
       json_get_or_lookup_or_default
  use json_module, only : json_file
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use vreman_cpu, only : vreman_compute_cpu
  use vreman_device, only : vreman_compute_device
  use registry, only : neko_registry
  use logger, only : LOG_SIZE, neko_log
  implicit none
  private

  !> Implements the Vreman LES model.
  !! @note Reference DOI: 10.1063/1.1785131
  type, public, extends(les_model_t) :: vreman_t
     !> Model constant, defaults to 0.07.
     real(kind=rp) :: c
     !> Critical Richardson number and reference temperature.
     real(kind=rp) :: ri_c, ref_temp
     !> Gravity vector
     real(kind=rp) :: g(3)
     !> Whether to apply the Richardson-number correction.
     logical :: if_corr = .false.
     !> Whether to use anisotropic filter widths.
     logical :: if_aniso = .false.
     !> Name of the scalar field used for the stability treatment.
     character(len=:), allocatable :: scalar_name
     !> Filter-width tensor (only used if if_aniso).
     type(field_t) :: g11, g12, g13, g22, g23, g33
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => vreman_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          vreman_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => vreman_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => vreman_compute
  end type vreman_t

contains
  !> Constructor.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param json A dictionary with parameters.
  subroutine vreman_init(this, fluid, json)
    class(vreman_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name

    real(kind=rp) :: c, ri_c, ref_temp
    real(kind=rp), allocatable :: g(:)
    character(len=:), allocatable :: scalar_name
    character(len=:), allocatable :: delta_type
    logical :: if_ext, if_corr, if_aniso
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    ! Based on the Smagorinsky Cs = 0.17.
    call json_get_or_lookup_or_default(json, "c", c, 0.07_rp)
    call json_get_or_default(json, "extrapolation", if_ext, .false.)
    call json_get_or_default(json, "anisotropic", if_aniso, .false.)
    call json_get_or_default(json, "buoyancy_correction", if_corr, .false.)

    ri_c = 0.25_rp
    if (if_corr) then
       call json_get_or_lookup(json, "Ri_c", ri_c)
       if (ri_c .le. 0.0_rp) then
          call neko_error("The critical Richardson number should be positive.")
       end if
    end if

    if (if_corr) then
       call json_get(json, "scalar_field", scalar_name)
       call json_get_or_lookup(json, "reference_temperature", ref_temp)
       call json_get_or_lookup(json, "g", g)
       if (.not. size(g) == 3) then
          call neko_error("VREMAN: The gravity vector should have 3 components")
       end if
       if (ref_temp .le. 0.0_rp) then
          call neko_error("The reference temperature should be positive.")
       end if
    else
       ! These values are not used
       g = [0.0_rp, 0.0_rp, 0.0_rp]
       ref_temp = 0.0_rp
       scalar_name = ""
    end if

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : Vreman'
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Delta evaluation : ', delta_type
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'c : ', c
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'extrapolation : ', if_ext
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'anisotropic : ', if_aniso
    call neko_log%message(log_buf)
    write(log_buf, '(A, L1)') 'buoyancy correction : ', if_corr
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call vreman_init_from_components(this, fluid, c, nut_name, &
         delta_type, if_ext, if_corr, scalar_name, ri_c, ref_temp, g, &
         if_aniso)

  end subroutine vreman_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param c The model constant.
  !! @param nut_name The name of the SGS viscosity field.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether to extrapolate the velocity.
  !! @param if_corr Whether to apply buoyancy correction.
  !! @param scalar_name The name of the scalar field for the stability
  !! treatment.
  !! @param ri_c Critical Richardson number.
  !! @param ref_temp Reference temperature for Richardson number.
  !! @param g The gravity vector, pointing in the direction of gravity.
  !! @param if_aniso Whether to use anisotropic filter widths.
  subroutine vreman_init_from_components(this, fluid, c, nut_name, &
       delta_type, if_ext, if_corr, scalar_name, ri_c, ref_temp, g, &
       if_aniso)
    class(vreman_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    real(kind=rp) :: c, ri_c, ref_temp
    real(kind=rp) :: g(3)
    character(len=*), intent(in) :: scalar_name
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext, if_corr
    logical, intent(in), optional :: if_aniso

    call this%free()

    call this%init_base(fluid, nut_name, delta_type, if_ext)
    this%c = c
    this%if_corr = if_corr
    this%ri_c = ri_c
    this%ref_temp = ref_temp
    ! Reversed so that g points upwards inside the model, as in the
    ! Deardorff model; the case file gives g in the direction of gravity.
    this%g = -g
    this%scalar_name = scalar_name

    this%if_aniso = .false.
    if (present(if_aniso)) this%if_aniso = if_aniso

    if (this%if_aniso) then
       call this%g11%init(this%coef%dof, "vreman_g11")
       call this%g12%init(this%coef%dof, "vreman_g12")
       call this%g13%init(this%coef%dof, "vreman_g13")
       call this%g22%init(this%coef%dof, "vreman_g22")
       call this%g23%init(this%coef%dof, "vreman_g23")
       call this%g33%init(this%coef%dof, "vreman_g33")
       call this%compute_delta_tensor(this%g11, this%g12, this%g13, &
            this%g22, this%g23, this%g33)
    end if

    if (this%if_corr .and. &
         sqrt(sum(this%g**2)) .le. 0.0_rp) then
       call neko_error("VREMAN: The gravity vector must have at least " // &
            "one nonzero component")
    end if

  end subroutine vreman_init_from_components

  !> Destructor for the les_model_t (base) class.
  subroutine vreman_free(this)
    class(vreman_t), intent(inout) :: this

    call this%g11%free()
    call this%g12%free()
    call this%g13%free()
    call this%g22%free()
    call this%g23%free()
    call this%g33%free()
    call this%free_base()
  end subroutine vreman_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine vreman_compute(this, t, tstep)
    class(vreman_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    type(field_t), pointer :: u, v, w, u_e, v_e, w_e

    if (this%if_ext) then
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
              ulag, vlag, wlag, ext_bdf%advection_coeffs%x, ext_bdf%nadv)

       end associate
    end if

    ! Compute the eddy viscosity field
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call vreman_compute_device(this%if_ext, t, tstep, this%coef, &
            this%nut, this%delta, this%c, this%if_corr, &
            this%scalar_name, this%ri_c, this%ref_temp, this%g, &
            this%if_aniso, this%g11, this%g12, this%g13, &
            this%g22, this%g23, this%g33)
    else
       call vreman_compute_cpu(this%if_ext, t, tstep, this%coef, &
            this%nut, this%delta, this%c, this%if_corr, &
            this%scalar_name, this%ri_c, this%ref_temp, this%g, &
            this%if_aniso, this%g11, this%g12, this%g13, &
            this%g22, this%g23, this%g33)
    end if

  end subroutine vreman_compute

end module vreman
