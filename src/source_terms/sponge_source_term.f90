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
!> Implements the `sponge_t` type.
module sponge_source_term
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use json_utils, only : json_get, json_get_or_default, json_get
  use utils, only: neko_error
  use device, only: device_memcpy, HOST_TO_DEVICE
  use device_math, only: device_sub3, device_col2, device_add2s2
  use time_state, only: time_state_t
  use math, only: sub3, col2, add2s2
  use logger, only: neko_log, NEKO_LOG_DEBUG
  use neko_config, only: NEKO_BCKND_DEVICE
  use source_term, only: source_term_t
  use case, only: case_t
  use simcomp_executor, only: neko_simcomps
  use flow_ic, only: set_flow_ic_fld, set_flow_ic
  use field_list, only: field_list_t
  use coefs, only: coef_t
  use utils, only: NEKO_FNAME_LEN
  use file, only: file_t
  use comm, only: pe_rank
  use fld_file_output, only: fld_file_output_t
  implicit none
  private

  type, public, extends(source_term_t) :: sponge_source_term_t

     !> Velocity components at current time.
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     !> Base flow components.
     type(field_t) :: u_bf
     type(field_t) :: v_bf
     type(field_t) :: w_bf
     !> Work array.
     type(field_t) :: wk
     !> Flag indicating wether the baseflow has been set by the user.
     logical :: baseflow_set = .false.
     !> Indicates wether the baseflow should be applied from the initial
     !! condition.
     logical :: baseflow_is_ic = .false.
     !> Fringe field. This field need to be added to the registry and filled
     !! by the user.
     type(field_t), pointer :: fringe
     !> Fringe amplitude in the 3 cartesian directions.
     real(kind=rp) :: amplitudes(3)
     !> Name of the fringe field in the neko field registry.
     character(len=1024) :: fringe_registry_name
     !> Wether or not to write the fringe field to disk.
     logical :: dump_fields = .false.
     !> Name of the fld file in which to dump the fields
     character(NEKO_FNAME_LEN) :: dump_fname

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => sponge_init_from_json
     !> Initialize a sponge with a constant baseflow.
     procedure, pass(this) :: init_constant => &
          sponge_init_constant
     !> Initialize a sponge with a baseflow imported from the initial condition.
     procedure, pass(this) :: init_ic => &
          sponge_init_ic
     !> Initialize a sponge with a baseflow imported from a field file.
     procedure, pass(this) :: init_field => &
          sponge_init_field
     !> Common constructor.
     procedure, pass(this) :: init_common => &
          sponge_init_common
     !> Destructor.
     procedure, pass(this) :: free => sponge_free
     !> Compute the sponge field.
     procedure, pass(this) :: compute_ => sponge_compute
  end type sponge_source_term_t

contains

  !> Constructor from json.
  subroutine sponge_init_from_json(this, json, fields, coef, variable_name)
    class(sponge_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    real(kind=rp), allocatable :: amplitudes(:)
    real(kind=rp), allocatable :: ramp(:), pad(:)

    character(len=:), allocatable :: baseflow_file, baseflow_method
    character(len=:), allocatable :: read_str, fringe_registry_name, dump_fname
    character(len=NEKO_FNAME_LEN) :: fname, mesh_fname
    logical :: interpolate, dump_fields
    real(kind=rp), allocatable :: constant_value(:)
    real(kind=rp) :: start_time, end_time, tolerance
    integer :: nzones

    type(json_file) :: baseflow_subdict
    integer :: i,izone

    call neko_log%section("SPONGE SOURCE TERM", LVL=NEKO_LOG_DEBUG)

    call json_get_or_default(json, "dump_fields", dump_fields, .false.)
    call json_get_or_default(json, "dump_file_name", dump_fname, &
         "spng_fields")
    call json_get_or_default(json, "fringe_registry_name", &
         fringe_registry_name, "sponge_fringe")

    call json_get(json, "amplitudes", amplitudes)
    if (size(amplitudes) .ne. 3) &
         call neko_error("(SPONGE) Expected 3 elements for 'amplitudes'")

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))


    call json_get(json, 'baseflow', &
         baseflow_subdict)
    call json_get(baseflow_subdict, "method", baseflow_method)

    select case (trim(baseflow_method))

       ! Import the base flow from an fld file. The same parameters as field
       ! initial condition apply since we use the same subroutine to read adn
       ! potentially interpolate from a field file.
    case ("field")

       ! The lines below are just copy pasted from set_flow_ic_int, because we
       ! want to make sure to separate the JSON constructor
       ! sponge_init_from_json and the constructor from attributes
       ! sponge_init_field.
       call json_get(baseflow_subdict, 'file_name', read_str)
       fname = trim(read_str)
       call json_get_or_default(baseflow_subdict, 'interpolate', interpolate, &
            .false.)
       call json_get_or_default(baseflow_subdict, 'tolerance', tolerance, &
            0.000001_rp)
       call json_get_or_default(baseflow_subdict, 'mesh_file_name', read_str, &
            "none")
       mesh_fname = trim(read_str)

       call this%init_field(fields, coef, start_time, end_time, amplitudes, &
            fringe_registry_name, dump_fields, dump_fname, fname, interpolate, &
            tolerance, mesh_fname)

       ! Constant base flow
    case ("constant")

       call json_get(baseflow_subdict, "value", constant_value)
       if (size(constant_value) .lt. 3) &
            call neko_error("(SPONGE) Expected 3 elements for 'value'")

       call this%init_constant(fields, coef, start_time, end_time, &
            amplitudes, fringe_registry_name, dump_fields, dump_fname, &
            constant_value)

       ! Apply the initial condition as base flow.
    case ("initial_condition")

       call this%init_ic(fields, coef, start_time, end_time, amplitudes, &
            fringe_registry_name, dump_fields, dump_fname)

    case default
       call neko_error("(SPONGE)" // trim(baseflow_method) // &
            " is not a valid method")
    end select

    call neko_log%end_section(lvl=NEKO_LOG_DEBUG)

  end subroutine sponge_init_from_json

  !> Initialize a sponge with a constant baseflow.
  subroutine sponge_init_constant(this, fields, coef, start_time, end_time, &
       amplitudes, fringe_registry_name, dump_fields, dump_fname, &
       constant_values)
    class(sponge_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time, end_time
    real(kind=rp), intent(in) :: amplitudes(:)
    character(len=*), intent(in) :: fringe_registry_name, dump_fname
    logical, intent(in) :: dump_fields
    real(kind=rp), intent(in) :: constant_values(:)

    !
    ! Common constructor
    !
    call sponge_init_common(this, fields, coef, start_time, end_time, &
         amplitudes, fringe_registry_name, dump_fields, dump_fname)

    !
    ! Assign constant values
    !
    this%u_bf = constant_values(1)
    this%v_bf = constant_values(2)
    this%w_bf = constant_values(3)

    this%baseflow_set = .true.

  end subroutine sponge_init_constant

  !> Initialize a sponge with a baseflow imported from a field file.
  subroutine sponge_init_field(this, fields, coef, start_time, end_time, &
       amplitudes, fringe_registry_name, dump_fields, dump_fname, file_name, &
       interpolate, tolerance, mesh_file_name)
    class(sponge_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time, end_time
    real(kind=rp), intent(in) :: amplitudes(:)
    character(len=*), intent(in) :: fringe_registry_name, dump_fname
    logical, intent(in) :: dump_fields
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: interpolate
    real(kind=rp), intent(in) :: tolerance
    character(len=*), intent(inout) :: mesh_file_name

    !
    ! Common constructor
    !
    call sponge_init_common(this, fields, coef, start_time, end_time, &
         amplitudes, fringe_registry_name, dump_fields, dump_fname)

    !
    ! Use the initial condition field subroutine to set a field as baseflow
    !
    call set_flow_ic_fld(this%u_bf, this%v_bf, this%w_bf, this%wk, &
         file_name, interpolate, tolerance, mesh_file_name)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%u_bf%x, this%u_bf%x_d, this%u_bf%size(), &
            HOST_TO_DEVICE, .false.)
       call device_memcpy(this%v_bf%x, this%v_bf%x_d, this%v_bf%size(), &
            HOST_TO_DEVICE, .false.)
       call device_memcpy(this%w_bf%x, this%w_bf%x_d, this%w_bf%size(), &
            HOST_TO_DEVICE, .false.)
    end if

    this%baseflow_set = .true.

  end subroutine sponge_init_field

  !> Initialize a sponge with a baseflow imported from the initial condition.
  !! This method of initializing a sponge is JSON dependent and requires
  !! the neko case object to have a valid JSON object loaded.
  subroutine sponge_init_ic(this, fields, coef, start_time, end_time, &
       amplitudes, fringe_registry_name, dump_fields, dump_fname)
    class(sponge_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time, end_time
    real(kind=rp), intent(in) :: amplitudes(:)
    character(len=*), intent(in) :: fringe_registry_name, dump_fname
    logical, intent(in) :: dump_fields

    !
    ! Common constructor
    !
    call sponge_init_common(this, fields, coef, start_time, end_time, &
         amplitudes, fringe_registry_name, dump_fields, dump_fname)

    ! NOTE: The initial condition can only be obtained through
    ! the case file, which can only be accessed through neko_simcomps.
    ! However at this point neko_simcomps is not yet initialized, therefore
    ! we will do this at the very beginning of compute_ based on the flag below.
    this%baseflow_is_ic = .true.

  end subroutine sponge_init_ic

  !> Common constructor.
  subroutine sponge_init_common(this, fields, coef, start_time, end_time, &
       amplitudes, fringe_registry_name, dump_fields, dump_fname)
    class(sponge_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time, end_time
    real(kind=rp), intent(in) :: amplitudes(:)
    character(len=*), intent(in) :: fringe_registry_name, dump_fname
    logical, intent(in) :: dump_fields

    integer :: i

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    this%amplitudes(1) = amplitudes(1)
    this%amplitudes(2) = amplitudes(2)
    this%amplitudes(3) = amplitudes(3)

    this%fringe_registry_name = trim(fringe_registry_name)
    this%dump_fields = dump_fields
    this%dump_fname = trim(dump_fname)

    call neko_log%message("Initializing sponge", lvl = NEKO_LOG_DEBUG)

    call neko_log%message("Pointing at fields u,v,w", &
         lvl = NEKO_LOG_DEBUG)
    this%u => neko_field_registry%get_field_by_name("u")
    this%v => neko_field_registry%get_field_by_name("v")
    this%w => neko_field_registry%get_field_by_name("w")

    call neko_log%message("Initializing bf fields", &
         lvl = NEKO_LOG_DEBUG)

    call this%u_bf%init(this%u%dof, "sponge_u_bf")
    call this%v_bf%init(this%u%dof, "sponge_v_bf")
    call this%w_bf%init(this%u%dof, "sponge_w_bf")
    call this%wk%init(this%u%dof, "sponge_wk")

  end subroutine sponge_init_common

  !> Destructor.
  subroutine sponge_free(this)
    class(sponge_source_term_t), intent(inout) :: this

    call this%free_base()

    call this%u_bf%free()
    call this%v_bf%free()
    call this%w_bf%free()
    call this%wk%free

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%fringe)

    this%fringe_registry_name = ""
    this%baseflow_set = .false.

  end subroutine sponge_free

  !> Compute the sponge field.
  !! The sponge is applied according to the following definition:
  !!                    f_x = amp_x * fringe(x,y,z) * (u_bf - u)
  !!                    f_y = amp_y * fringe(x,y,z) * (v_bf - v)
  !!                    f_z = amp_z * fringe(x,y,z) * (w_bf - w)
  !!
  !! The fringe field must be added to the regisry and filled by the user
  !! according to where they want to apply the sponge.
  !! NOTE: Unfortunately this function is JSON dependent due to the initial
  !! condition being JSON dependent quite deep, there should be a nice way
  !! to get rid of that but it should take some effort.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine sponge_compute(this, time)
    class(sponge_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    integer :: i, nfields
    type(field_t), pointer :: fu, fv, fw

    type(json_file) :: json_subdict
    character(len=:), allocatable :: string_val

    type(fld_file_output_t) :: fout
    class(case_t), pointer :: case
    case => null()

    !
    ! If the baseflow comes from initial condition, call the initial condition
    ! function as they are called in simulation.f90
    !
    if (this%baseflow_is_ic .and. .not. this%baseflow_set) then
       case => neko_simcomps%case

       ! Here we recreate what is in the initial condition part of case.f90
       call json_get(case%params, 'case.fluid.initial_condition.type', &
            string_val)
       call json_get(case%params, 'case.fluid.initial_condition', &
            json_subdict)

       call neko_log%section("Fluid initial condition ")

       if (trim(string_val) .ne. 'user') then
          call set_flow_ic(this%u_bf, this%v_bf, this%w_bf, &
               this%wk, case%fluid%c_Xh, case%fluid%gs_Xh, string_val, &
               json_subdict)
       else
          call json_get(case%params, 'case.fluid.scheme', string_val)
          if (trim(string_val) .eq. 'compressible') then
             call set_flow_ic(case%fluid%rho, &
                  this%u_bf, this%v_bf, this%w_bf, this%wk, &
                  case%fluid%c_Xh, case%fluid%gs_Xh, &
                  case%user%initial_conditions, case%fluid%name)
          else
             !! NOTE: this%fringe is used instead of pressure as a work array
             call set_flow_ic(this%u_bf, this%v_bf, this%w_bf, &
                  this%wk, case%fluid%c_Xh, case%fluid%gs_Xh, &
                  case%user%initial_conditions, case%fluid%name)
          end if
       end if
       this%baseflow_set = .true.
       nullify(case)
    end if

    ! Do some checks at the first timestep
    if (time%tstep .eq. 1) then
       if (.not. this%baseflow_set) call neko_error("SPONGE: No baseflow set")

       ! Check if the user has added the fringe field in the registry
       if (.not. &
            neko_field_registry%field_exists(trim(this%fringe_registry_name))) &
            call neko_error("SPONGE: No fringe field set (" // &
                this%fringe_registry_name // " not found)")

       ! This will throw an error if the user hasn't added 'sponge_fringe'
       ! to the registry.
       this%fringe => &
            neko_field_registry%get_field(trim(this%fringe_registry_name))

       !
       ! Dump the fringe and/or baseflow fields for visualization
       !
       if (this%dump_fields) then
          call fout%init(sp, trim(this%dump_fname), 4)
          call fout%fields%assign_to_ptr(1, this%fringe)
          call fout%fields%assign_to_field(2, this%u_bf)
          call fout%fields%assign_to_field(3, this%v_bf)
          call fout%fields%assign_to_field(4, this%w_bf)
          call fout%sample(time%t)
       end if
    end if

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! wk = u_bf - u
       call device_sub3(this%wk%x_d, this%u_bf%x_d, this%u%x_d, this%u%size())
       ! wk = fringe * wk = fringe * (u_bf - u)
       call device_col2(this%wk%x_d, this%fringe%x_d, this%fringe%size())
       ! fu = fu + amplitude(1)*wk = fu + amplitude(1)*fringe*(u_bf - u)
       call device_add2s2(fu%x_d, this%wk%x_d, this%amplitudes(1), &
            fu%dof%size())

       call device_sub3(this%wk%x_d, this%v_bf%x_d, this%v%x_d, this%v%size())
       call device_col2(this%wk%x_d, this%fringe%x_d, this%fringe%size())
       call device_add2s2(fv%x_d, this%wk%x_d, this%amplitudes(2), &
            fv%dof%size())

       call device_sub3(this%wk%x_d, this%w_bf%x_d, this%w%x_d, this%w%size())
       call device_col2(this%wk%x_d, this%fringe%x_d, this%fringe%size())
       call device_add2s2(fw%x_d, this%wk%x_d, this%amplitudes(3), &
            fw%dof%size())
    else
       call sub3(this%wk%x, this%u_bf%x, this%u%x, this%u%size())
       call col2(this%wk%x, this%fringe%x, this%fringe%size())
       call add2s2(fu%x, this%wk%x, this%amplitudes(1), fu%dof%size())

       call sub3(this%wk%x, this%v_bf%x, this%v%x, this%v%size())
       call col2(this%wk%x, this%fringe%x, this%fringe%size())
       call add2s2(fv%x, this%wk%x, this%amplitudes(2), fv%dof%size())

       call sub3(this%wk%x, this%w_bf%x, this%w%x, this%w%size())
       call col2(this%wk%x, this%fringe%x, this%fringe%size())
       call add2s2(fw%x, this%wk%x, this%amplitudes(3), fw%dof%size())
    end if

  end subroutine sponge_compute

end module sponge_source_term
