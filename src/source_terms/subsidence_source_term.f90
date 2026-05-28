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
  use device, only : HOST_TO_DEVICE
  implicit none
  private

  !> Subsidence source term: applies vertical velocity profile w_sub
  !! to the fields via w_sub * d(phi)/dz with phi = u,v or scalar T, and
  !! w_sub prescribed in the case.f90 file.
  type, public, extends(source_term_t) :: subsidence_source_term_t
     !> The subsidence vertical velocity profile w_sub(z)
     type(field_t), pointer :: w_sub => null()
     !> Velocity fields for fluid subsidence.
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     !> Scalar field for temperature subsidence.
     type(field_t), pointer :: scalar => null()
     !> Registry name and w_sub values
     character(len=:), allocatable :: profile_registry_name
     character(len=:), allocatable :: method_str
     real(kind=rp) :: div_rate
     real(kind=rp) :: w_sub_max
     real(kind=rp) :: w_sub_max_height
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
    real(kind=rp) :: div_rate = 0.0_rp, w_sub_max = 0.0_rp, w_sub_max_height = 0.0_rp
    character(len=:), allocatable :: profile_registry_name
    character(len=:), allocatable :: method_str

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call json_get_or_default(json, "profile_registry_name", profile_registry_name, "w_sub")
    this%profile_registry_name = trim(profile_registry_name)

    call json_get_or_default(json, "method", method_str, "user")
    this%method_str = trim(method_str)

    select case (trim(method_str))
    case ("linear")
       call json_get_or_default(json, "div_rate", div_rate, 1.0e-5_rp)
    case ("linear_constant")
       call json_get_or_default(json, "w_sub_max", w_sub_max, 1.0e-5_rp)
       call json_get_or_default(json, "w_sub_max_height", w_sub_max_height, 500.0_rp)
    case ("user")
       ! No extra properties to collect
    case default
       call neko_error("SUBSIDENCE: Unknown method: " // trim(method_str))
    end select

    call subsidence_source_term_init_from_components(this, fields, coef, &
         start_time, end_time, method_str, div_rate, w_sub_max, w_sub_max_height, &
         profile_registry_name, variable_name)

  end subroutine subsidence_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine subsidence_source_term_init_from_components(this, fields, coef, &
       start_time, end_time, method, div_rate, w_sub_max, w_sub_max_height, &
       profile_registry_name, variable_name)
    class(subsidence_source_term_t), intent(inout) :: this
    class(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time
    character(len=*), intent(in) :: method
    real(kind=rp), intent(in) :: div_rate
    real(kind=rp), intent(in) :: w_sub_max
    real(kind=rp), intent(in) :: w_sub_max_height
    character(len=*), intent(in) :: profile_registry_name
    character(len=*), intent(in) :: variable_name

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    this%method_str = trim(method)
    this%div_rate = div_rate
    this%w_sub_max = w_sub_max
    this%w_sub_max_height = w_sub_max_height
    this%profile_registry_name = trim(profile_registry_name)
    this%check = .true.

    if (fields%size() == 1) then
       ! Scalar subsidence source term for temperature or another scalar.
       if (.not. neko_registry%field_exists(trim(variable_name))) then
          call neko_error("SUBSIDENCE: No scalar field found for " // &
               trim(variable_name))
       end if
       this%scalar => neko_registry%get_field(trim(variable_name))
    else if (fields%size() == 3) then
       ! Fluid subsidence source term for momentum.
       this%u => neko_registry%get_field("u")
       this%v => neko_registry%get_field("v")
    else
       call neko_error("SUBSIDENCE: Unexpected number of fields for source term.")
    end if
  end subroutine subsidence_source_term_init_from_components

  !> Destructor.
  subroutine subsidence_source_term_free(this)
    class(subsidence_source_term_t), intent(inout) :: this

    call this%free_base()
    nullify(this%w_sub)
    nullify(this%u)
    nullify(this%v)
    nullify(this%scalar)
  end subroutine subsidence_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time value.
  subroutine subsidence_source_term_compute(this, time)
    class(subsidence_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i
    real(kind=rp) :: z

    if (this%check) then
       if (trim(this%method_str) .eq. "user") then
          if (.not. neko_registry%field_exists(trim(this%profile_registry_name))) then
             call neko_error("SUBSIDENCE: No w_sub profile set (searching for " // &
                  trim(this%profile_registry_name) // "). Add it to the registry in your case file")
          end if
          this%w_sub => neko_registry%get_field(trim(this%profile_registry_name))
       else
          if (associated(this%u)) then
             call neko_registry%add_field(this%u%dof, trim(this%profile_registry_name))
          else
             call neko_registry%add_field(this%scalar%dof, trim(this%profile_registry_name))
          end if

          this%w_sub => neko_registry%get_field(trim(this%profile_registry_name))

          select case (trim(this%method_str))
          case ("linear")
             do i = 1, this%w_sub%size()
                z = this%w_sub%dof%z(i,1,1,1)
                this%w_sub%x(i,1,1,1) = this%div_rate * z
             end do

          case ("linear_constant")
             do i = 1, this%w_sub%size()
                z = this%w_sub%dof%z(i,1,1,1)
                if (z .le. this%w_sub_max_height) then
                   this%w_sub%x(i,1,1,1) = (this%w_sub_max / this%w_sub_max_height) * z
                else
                   this%w_sub%x(i,1,1,1) = this%w_sub_max
                end if
             end do
          end select

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call this%w_sub%copy_from(HOST_TO_DEVICE, .true.)
          end if
       end if
       this%check = .false.
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call subsidence_source_term_compute_device(this%u, this%v, this%scalar, &
            this%fields, this%coef, this%w_sub)
    else
       call subsidence_source_term_compute_cpu(this%u, this%v, this%scalar, &
            this%fields, this%coef, this%w_sub)
    end if

  end subroutine subsidence_source_term_compute

end module subsidence_source_term
