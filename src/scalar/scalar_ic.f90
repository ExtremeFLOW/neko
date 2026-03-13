! Copyright (c) 2021-2025, The Neko Authors
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
!> Scalar initial condition
module scalar_ic_m
  use gather_scatter_m, only : gs_t, GS_OP_ADD
  use neko_config_m, only : NEKO_BCKND_DEVICE
  use num_types_m, only : rp
  use device_math_m, only : device_col2
  use device_m, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use field_m, only : field_t
  use utils_m, only : neko_error, filename_chsuffix, filename_suffix, &
       neko_warning, NEKO_FNAME_LEN, extract_fld_file_index
  use coefs_m, only : coef_t
  use math_m, only : col2, cfill, cfill_mask
  use user_intf_m, only : user_initial_conditions_intf
  use json_module, only : json_file
  use json_utils_m, only : json_get, json_get_or_default, &
       json_get_or_lookup, json_get_or_lookup_or_default
  use point_zone_m, only : point_zone_t
  use point_zone_registry_m, only : neko_point_zone_registry
  use logger_m, only : neko_log, LOG_SIZE
  use fld_file_data_m, only : fld_file_data_t
  use fld_file_m, only : fld_file_t
  use checkpoint_m, only : chkp_t
  use file_m, only : file_t
  use global_interpolation_m, only : global_interpolation_t
  use interpolation_m, only : interpolator_t
  use space_m, only : space_t, GLL
  use field_list_m, only : field_list_t
  use import_field_utils_m, only : import_fields
  implicit none
  private

  interface set_scalar_ic
     module procedure set_scalar_ic_int, set_scalar_ic_usr
  end interface set_scalar_ic

  public :: set_scalar_ic

contains

  !> Set scalar initial condition (builtin)
  !! @details Set scalar initial condition using one of the builtin types
  !! currently supported:
  !! - uniform
  !! - point zone
  !! - field
  !! @param s Scalar field.
  !! @param coef Coefficient.
  !! @param gs Gather-Scatter object.
  !! @param type Type of initial condition.
  !! @param params JSON parameters.
  !! @param i Index of the scalar field.
  subroutine set_scalar_ic_int(s, coef, gs, type, params, i)
    type(field_t), intent(inout) :: s
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    character(len=*) :: type
    type(json_file), intent(inout) :: params
    integer, intent(in) :: i

    ! Variables for retrieving JSON parameters
    real(kind=rp) :: ic_value
    character(len=:), allocatable :: read_str
    character(len=NEKO_FNAME_LEN) :: fname, mesh_fname
    real(kind=rp) :: zone_value, tol
    logical :: interpolate
    integer :: tgt_scal_idx

    if (trim(type) .eq. 'uniform') then

       call json_get_or_lookup(params, 'value', ic_value)
       call set_scalar_ic_uniform(s, ic_value)

    else if (trim(type) .eq. 'point_zone') then

       call json_get_or_lookup(params, 'base_value', ic_value)
       call json_get(params, 'zone_name', read_str)
       call json_get_or_lookup(params, 'zone_value', zone_value)

       call set_scalar_ic_point_zone(s, ic_value, read_str, zone_value)

    else if (trim(type) .eq. 'field') then

       call json_get(params, 'file_name', read_str)
       fname = trim(read_str)
       call json_get_or_default(params, 'interpolate', interpolate, &
            .false.)
       call json_get_or_lookup_or_default(params, 'tolerance', tol, 0.000001_rp)
       call json_get_or_default(params, 'mesh_file_name', read_str, &
            "none")
       mesh_fname = trim(read_str)

       ! Give the user the option to select which scalar they want to import
       ! the values from, in the fld file. 0 corresponds to temperature.
       call json_get_or_default(params, 'target_index', tgt_scal_idx, i)

       call set_scalar_ic_fld(s, fname, interpolate, tol, mesh_fname, i, &
            tgt_scal_idx)

    else
       call neko_error('Invalid initial condition')
    end if

    call set_scalar_ic_common(s, coef, gs)

  end subroutine set_scalar_ic_int

  !> Set scalar intial condition (user defined)
  !! @details Set scalar initial condition using a user defined function.
  !! @param scheme_name Name of the scheme calling the user routine.
  !! @param s Scalar field.
  !! @param coef Coefficient.
  !! @param gs Gather-Scatter object.
  !! @param user_proc User defined initial condition function.
  subroutine set_scalar_ic_usr(scheme_name, s, coef, gs, user_proc)
    character(len=*), intent(in) :: scheme_name
    type(field_t), target, intent(inout) :: s
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    procedure(user_initial_conditions_intf) :: user_proc
    type(field_list_t) :: fields

    call neko_log%message("Type: user")

    call fields%init(1)
    call fields%assign_to_field(1, s)

    call user_proc(scheme_name, fields)
    call set_scalar_ic_common(s, coef, gs)

  end subroutine set_scalar_ic_usr

  !> Set scalar initial condition (common)
  !! @details Finalize scalar initial condition by distributing the initial
  !! condition across elements and multiplying by the coefficient (if any).
  !! @param s Scalar field.
  !! @param coef Coefficient.
  !! @param gs Gather-Scatter object.
  subroutine set_scalar_ic_common(s, coef, gs)
    type(field_t), intent(inout) :: s
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n

    n = s%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, n, HOST_TO_DEVICE, sync = .false.)
    end if

    ! Ensure continuity across elements for initial conditions
    call gs%op(s%x, n, GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(s%x_d, coef%mult_d, n)
    else
       call col2(s%x, coef%mult, n)
    end if

  end subroutine set_scalar_ic_common

  !> Uniform initial condition
  !! @details Set scalar initial condition to a uniform value across the domain.
  !! @param s Scalar field.
  !! @param ic_value Desired value of the scalar field.
  subroutine set_scalar_ic_uniform(s, ic_value)
    type(field_t), intent(inout) :: s
    real(kind=rp), intent(in) :: ic_value
    integer :: n
    character(len=LOG_SIZE) :: log_buf

    call neko_log%message("Type : uniform")
    write (log_buf, '(A,ES12.6)') "Value: ", ic_value
    call neko_log%message(log_buf)

    s = ic_value
    n = s%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cfill(s%x, ic_value, n)
    end if

  end subroutine set_scalar_ic_uniform

  !> Point zone initial condition
  !! @details Set scalar initial condition to a uniform value across a point
  !! zone.
  !! @param s Scalar field.
  !! @param base_value Base value of the scalar field.
  !! @param zone_name Name of the point zone.
  !! @param zone_value Desired value of the scalar field in the point zone.
  subroutine set_scalar_ic_point_zone(s, base_value, zone_name, zone_value)
    type(field_t), intent(inout) :: s
    real(kind=rp), intent(in) :: base_value
    character(len=*), intent(in) :: zone_name
    real(kind=rp), intent(in) :: zone_value

    ! Internal variables
    character(len=LOG_SIZE) :: log_buf
    class(point_zone_t), pointer :: zone
    integer :: size

    call neko_log%message("Type       : point_zone")
    write (log_buf, '(A,ES12.6)') "Base value: ", base_value
    call neko_log%message(log_buf)
    call neko_log%message("Zone name : " // trim(zone_name))
    write (log_buf, '(A,ES12.6)') "Zone value: ", zone_value
    call neko_log%message(log_buf)

    size = s%dof%size()
    zone => neko_point_zone_registry%get_point_zone(trim(zone_name))

    call set_scalar_ic_uniform(s, base_value)
    call cfill_mask(s%x, zone_value, size, zone%mask%get(), zone%size)

  end subroutine set_scalar_ic_point_zone

  !> Set the initial condition of the scalar based on a field.
  !! @detail The field is read from an `fld` file. If enabled, interpolation
  !! is also possible. In that case, the mesh coordinates can be read from
  !! another file in the `fld` field series.
  !! @param s The scalar field.
  !! @param file_name The name of the "fld" file series.
  !! @param sample_idx index of the field file .f000* to read, default is
  !! -1.
  !! @param interpolate Flag to indicate wether or not to interpolate the
  !! values onto the current mesh.
  !! @param tolerance If interpolation is enabled, tolerance for finding the
  !! points in the mesh.
  !! @param mesh_file_name If interpolation is enabled, name of the field
  !! file series where the mesh coordinates are located.
  !! @param i Index of the scalar field. 0 corresponds to temperature.
  !! @param target_idx The index of the scalar field to import from the
  !! fld file. 0 corresponds to temperature.
  subroutine set_scalar_ic_fld(s, file_name, &
       interpolate, tolerance, mesh_file_name, i, target_idx)
    type(field_t), target, intent(inout) :: s
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: interpolate
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(inout) :: mesh_file_name
    integer, intent(in) :: i
    integer, intent(in) :: target_idx

    character(len=LOG_SIZE) :: log_buf
    type(field_t), pointer :: ss
    type(field_list_t) :: s_tgt_list

    if (i .ne. target_idx) then
       write (log_buf, '(A,I0,A,I0)') "Loading scalar #", target_idx, &
            " into scalar #", i
       call neko_log%message(log_buf)
    end if

    ! use a pointer since import_fields needs a pointer as input
    ss => s

    ! Put ss in a field list of 1 element
    call s_tgt_list%init(1)
    call s_tgt_list%assign(1, ss)

    call import_fields(file_name, mesh_file_name, &
         s_target_list = s_tgt_list, & ! The target field
         s_index_list = [target_idx], & ! Take values from target scalar
         interpolate = interpolate, tolerance = tolerance)

    call s_tgt_list%free()

    nullify(ss)

    ! If we are on GPU we need to move s back to the host
    ! since set_scalar_ic_common copies it again to the device.
    call s%copy_from(device_to_host, .true.)

  end subroutine set_scalar_ic_fld

end module scalar_ic_m
