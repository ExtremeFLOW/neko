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
!!> Implements the `subsidence_source_term_t` type.
module subsidence_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use registry, only : neko_registry
  use time_state, only : time_state_t
  use import_field_utils, only : import_fields
  use scratch_registry, only : neko_scratch_registry
  use subsidence_source_term_cpu, only : subsidence_source_term_compute_cpu
  use subsidence_source_term_device, only : &
       subsidence_source_term_compute_device
  implicit none
  private

  !> Subsidence source term: applies vertical velocity profile w_sub 
  !! to the fields via w_sub * d(phi)/dz with phi = u,v and w_sub 
  !! prescribed in the case.f90 file.
  type, public, extends(source_term_t) :: subsidence_source_term_t
     !> The subsidence vertical velocity profile w_sub(z)
     type(field_t), pointer :: w_sub => null()
     !> Velocity fields
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     !> Vertical direction indication
     real(kind=rp) :: vertical_dir(3)
      !> Registry name where the w_sub profile is stored
      character(len=1024) :: profile_registry_name
      !> Flag to check and retrieve registry entries at first compute
      logical :: check = .true.
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => subsidence_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          subsidence_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => subsidence_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => subsidence_source_term_compute
  end type subsidence_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine subsidence_source_term_init_from_json(this, json, fields, coef, &
      variable_name)
    class(subsidence_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp) :: start_time, end_time
    type(json_file) :: interp_subdict
    real(kind=rp) :: vertical_dir(3)
    character(len=1024) :: profile_registry_name

   !!! add various checks?

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))
    call json_get_or_default(json, "vertical_dir", vertical_dir, [0.0_rp, 0.0_rp, 1.0_rp])
    call json_get_or_default(json, "profile_registry_name", profile_registry_name, "w_sub")

    if (.not. size(vertical_dir) == 3) then
       call neko_error("Vertical direction must be specified with a three-component &
        vector (e.g. [0.0, 0.0, 1.0] to set z as vertical).")
    end if

        call subsidence_source_term_init_from_components(this, fields, coef, &
          vertical_dir, start_time, end_time, profile_registry_name)

  end subroutine subsidence_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param vertical_dir The vertical direction vector (e.g. [0,0,1] for z).
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine subsidence_source_term_init_from_components(this, fields, coef, &
   vertical_dir, start_time, end_time, profile_registry_name)
    class(subsidence_source_term_t), intent(inout) :: this
    class(field_list_t), intent(in), target :: fields
    type(coef_t) :: coef
    real(kind=rp), intent(in) :: vertical_dir(3)
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time
    character(len=*), intent(in) :: profile_registry_name

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

     ! Defer retrieval of the `w_sub` profile from the registry until compute
     this%profile_registry_name = trim(profile_registry_name)
     this%check = .true.

     this%u => neko_registry%get_field("u")
     this%v => neko_registry%get_field("v")
     this%vertical_dir = vertical_dir
  end subroutine subsidence_source_term_init_from_components

  !> Destructor.
  subroutine subsidence_source_term_free(this)
    class(subsidence_source_term_t), intent(inout) :: this

    call this%free_base()
    nullify(this%w_sub)
    nullify(this%u)
    nullify(this%v)
  end subroutine subsidence_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time value.
  subroutine subsidence_source_term_compute(this, time)
    class(subsidence_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
     ! On first compute, retrieve the user-provided `w_sub` profile from the
     ! field registry using the configured `profile_registry_name`.
     if (this%check) then
       if (.not. neko_registry%field_exists(trim(this%profile_registry_name))) then
         call neko_error("SUBSIDENCE: No w_sub profile set (searching for " // &
            trim(this%profile_registry_name) // ") - add it to the registry in your case file")
       end if

       this%w_sub => neko_registry%get_field(trim(this%profile_registry_name))
       this%check = .false.
     end if

     if (NEKO_BCKND_DEVICE .eq. 1) then
       call subsidence_source_term_compute_device(this%u, this%v, this%fields, &
          this%coef, this%w_sub, this%vertical_dir)
     else
       call subsidence_source_term_compute_cpu(this%u, this%v, this%fields, &
          this%coef, this%w_sub, this%vertical_dir)
     end if

  end subroutine subsidence_source_term_compute

end module subsidence_source_term