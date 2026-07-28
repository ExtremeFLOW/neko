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
!
!> Implements the `wall_shear_stress_t` type.
module wall_shear_stress_simcomp
  use num_types, only : rp
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use simulation_component, only : simulation_component_t
  use time_based_controller, only : time_based_controller_t
  use time_state, only : time_state_t
  use case, only : case_t
  use field, only : field_t
  use field_math, only : field_copy
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use coefs, only : coef_t
  use vector, only : vector_t
  use boundary_data, only : boundary_data_t
  use operators, only : strain_rate
  use drag_torque, only : calc_force_array, device_calc_force_array
  use math, only : vdot3, sqrt_inplace
  use device_math, only : device_vdot3, device_sqrt_inplace
  use neko_config, only : NEKO_BCKND_DEVICE
  use ale_manager, only : neko_ale
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error, NEKO_VARNAME_LEN
  implicit none
  private

  !> A simulation component that computes the wall shear stress on one or more
  !! labelled boundary zones and registers it in the `neko_registry`.
  !! The result is written into the registry fields `<computed_field>_x`,
  !! `_y`, `_z` and `_mag`, which are zero away from the marked zones.
  type, public, extends(simulation_component_t) :: wall_shear_stress_t
     !> X velocity component.
     type(field_t), pointer :: u => null()
     !> Y velocity component.
     type(field_t), pointer :: v => null()
     !> Z velocity component.
     type(field_t), pointer :: w => null()
     !> Pressure.
     type(field_t), pointer :: p => null()
     !> Dynamic viscosity.
     type(field_t), pointer :: mu => null()
     !> X component of the wall shear stress.
     type(field_t), pointer :: tau_x => null()
     !> Y component of the wall shear stress.
     type(field_t), pointer :: tau_y => null()
     !> Z component of the wall shear stress.
     type(field_t), pointer :: tau_z => null()
     !> Magnitude of the wall shear stress.
     type(field_t), pointer :: tau_mag => null()
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Boundary points of the chosen zones, with their geometry. Built with
     !! the `coef` normal convention, pointing out of the fluid.
     type(boundary_data_t) :: bdata
     !> Labelled zones to include.
     integer, allocatable :: zone_indices(:)
     !> Masked dynamic viscosity.
     type(vector_t) :: mu_msk
     !> Masked pressure.
     type(vector_t) :: p_msk
     !> Masked strain rate components.
     type(vector_t) :: s11, s22, s33, s12, s13, s23
     !> Viscous traction components at the masked points.
     type(vector_t) :: t1, t2, t3
     !> Pressure traction.
     type(vector_t) :: pt1, pt2, pt3
     !> Which traction components are registered as fields.
     logical :: want_x = .true.
     logical :: want_y = .true.
     logical :: want_z = .true.
     logical :: want_mag = .true.
     !> Whether ALE or not.
     logical :: ale_enabled = .false.
     !> Base name of the registered fields.
     character(len=:), allocatable :: computed_field
   contains
     !> Constructors.
     procedure, pass(this) :: init => wall_shear_stress_init_from_json
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     procedure, pass(this) :: init_from_controllers => &
          wall_shear_stress_init_from_controllers
     procedure, pass(this) :: init_from_controllers_properties => &
          wall_shear_stress_init_from_controllers_properties
     procedure, private, pass(this) :: init_common => &
          wall_shear_stress_init_common
     !> Destructor.
     procedure, pass(this) :: free => wall_shear_stress_free
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute_ => wall_shear_stress_compute
  end type wall_shear_stress_t

contains

  !> Constructor from json.
  !! @param json JSON object with the parameters.
  !! @param case The case object.
  subroutine wall_shear_stress_init_from_json(this, json, case)
    class(wall_shear_stress_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name, computed_field, fluid_name
    character(len=:), allocatable :: viscosity_field
    integer, allocatable :: zone_indices(:)
    character(len=NEKO_VARNAME_LEN), allocatable :: fields(:)
    character(len=NEKO_VARNAME_LEN), allocatable :: components(:)
    logical :: want_x, want_y, want_z, want_mag
    integer :: i

    call json_get_or_default(json, "name", name, "wall_shear_stress")
    call json_get_or_default(json, "computed_field", computed_field, "tau")
    call json_get_or_default(json, "fluid_name", fluid_name, "fluid")
    call json_get_or_default(json, "viscosity_field", viscosity_field, &
         trim(fluid_name) // "_mu_tot")
    call json_get(json, "zone_indices", zone_indices)

    ! Parse the component selection. Default is all four.
    if (json%valid_path("components")) then
       call json_get(json, "components", components)
    else
       allocate(components(1))
       components(1) = "all"
    end if
    call wall_shear_stress_parse_components(components, &
         want_x, want_y, want_z, want_mag)

    call wall_shear_stress_build_field_list(computed_field, want_x, want_y, &
         want_z, want_mag, fields)

    call this%init_base(json, case)

    ! Register the selected fields so they are available in the registry.
    do i = 1, size(fields)
       call neko_registry%add_field(case%fluid%c_Xh%dof, trim(fields(i)), &
            ignore_existing = .true.)
    end do

    call this%init_common(name, computed_field, viscosity_field, &
         zone_indices, want_x, want_y, want_z, want_mag, &
         case%fluid%c_Xh)

  end subroutine wall_shear_stress_init_from_json

  !> Translate a `components` list into the four selection flags.
  !! @param components The requested component names. Accepted values are
  !! "x", "y", "z", "mag" and "all".
  !! @param want_x Set true if the x component is requested.
  !! @param want_y Set true if the y component is requested.
  !! @param want_z Set true if the z component is requested.
  !! @param want_mag Set true if the magnitude is requested.
  subroutine wall_shear_stress_parse_components(components, &
       want_x, want_y, want_z, want_mag)
    character(len=*), intent(in) :: components(:)
    logical, intent(out) :: want_x, want_y, want_z, want_mag
    integer :: i
    character(len=:), allocatable :: c

    want_x = .false.
    want_y = .false.
    want_z = .false.
    want_mag = .false.

    if (size(components) .eq. 0) then
       call neko_error("wall_shear_stress: 'components' must not be empty")
    end if

    do i = 1, size(components)
       c = trim(components(i))
       select case (c)
       case ("all")
          want_x = .true.
          want_y = .true.
          want_z = .true.
          want_mag = .true.
       case ("x")
          want_x = .true.
       case ("y")
          want_y = .true.
       case ("z")
          want_z = .true.
       case ("mag")
          want_mag = .true.
       case default
          call neko_error("wall_shear_stress: unknown component '" // c // &
               "'. Use x, y, z, mag, or all.")
       end select
    end do

    if (.not. (want_x .or. want_y .or. want_z .or. want_mag)) then
       call neko_error("wall_shear_stress: no components selected")
    end if

  end subroutine wall_shear_stress_parse_components

  !> Build the array of selected field names in x, y, z, mag order.
  !! @param computed_field The base name of the fields.
  !! @param want_x Whether the x component is selected.
  !! @param want_y Whether the y component is selected.
  !! @param want_z Whether the z component is selected.
  !! @param want_mag Whether the magnitude is selected.
  !! @param fields Allocated with the selected names.
  subroutine wall_shear_stress_build_field_list(computed_field, want_x, &
       want_y, want_z, want_mag, fields)
    character(len=*), intent(in) :: computed_field
    logical, intent(in) :: want_x, want_y, want_z, want_mag
    character(len=NEKO_VARNAME_LEN), allocatable, intent(out) :: fields(:)
    integer :: n, k

    n = 0
    if (want_x) n = n + 1
    if (want_y) n = n + 1
    if (want_z) n = n + 1
    if (want_mag) n = n + 1

    allocate(fields(n))
    k = 0
    if (want_x) then
       k = k + 1
       fields(k) = trim(computed_field) // "_x"
    end if
    if (want_y) then
       k = k + 1
       fields(k) = trim(computed_field) // "_y"
    end if
    if (want_z) then
       k = k + 1
       fields(k) = trim(computed_field) // "_z"
    end if
    if (want_mag) then
       k = k + 1
       fields(k) = trim(computed_field) // "_mag"
    end if

  end subroutine wall_shear_stress_build_field_list

  !> Constructor from components, passing controllers.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution order priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param computed_field The base name of the registered fields.
  !! @param viscosity_field The name of the registered viscosity field.
  !! @param zone_indices Labelled zones to include.
  !! @param coef The SEM coefficients.
  subroutine wall_shear_stress_init_from_controllers(this, name, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       computed_field, viscosity_field, zone_indices, coef, &
       want_x, want_y, want_z, want_mag)
    class(wall_shear_stress_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*), intent(in) :: computed_field
    character(len=*), intent(in) :: viscosity_field
    integer, intent(in) :: zone_indices(:)
    type(coef_t), intent(inout), target :: coef
    logical, intent(in), optional :: want_x, want_y, want_z, want_mag
    character(len=NEKO_VARNAME_LEN), allocatable :: fields(:)
    logical :: wx, wy, wz, wm
    integer :: i

    ! Default to all four components when not specified.
    wx = .true.
    wy = .true.
    wz = .true.
    wm = .true.
    if (present(want_x)) wx = want_x
    if (present(want_y)) wy = want_y
    if (present(want_z)) wz = want_z
    if (present(want_mag)) wm = want_mag
    if (.not. (wx .or. wy .or. wz .or. wm)) then
       call neko_error("wall_shear_stress: no components selected")
    end if

    call wall_shear_stress_build_field_list(computed_field, wx, wy, wz, wm, &
         fields)

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)

    do i = 1, size(fields)
       call neko_registry%add_field(coef%dof, trim(fields(i)), &
            ignore_existing = .true.)
    end do

    call this%init_common(name, computed_field, viscosity_field, &
         zone_indices, wx, wy, wz, wm, coef)

  end subroutine wall_shear_stress_init_from_controllers

  !> Constructor from components, passing properties to the
  !! `time_based_controller` components in the base type.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution order priority of the simcomp.
  !! @param preprocess_control Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_control Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_control Control mode for output.
  !! @param output_value Value parameter for output.
  !! @param computed_field The base name of the registered fields.
  !! @param viscosity_field The name of the registered viscosity field.
  !! @param zone_indices Labelled zones to include.
  !! @param coef The SEM coefficients.
  subroutine wall_shear_stress_init_from_controllers_properties(this, name, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, computed_field, &
       viscosity_field, zone_indices, coef, &
       want_x, want_y, want_z, want_mag)
    class(wall_shear_stress_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=*), intent(in) :: computed_field
    character(len=*), intent(in) :: viscosity_field
    integer, intent(in) :: zone_indices(:)
    type(coef_t), intent(inout), target :: coef
    logical, intent(in), optional :: want_x, want_y, want_z, want_mag
    character(len=NEKO_VARNAME_LEN), allocatable :: fields(:)
    logical :: wx, wy, wz, wm
    integer :: i

    wx = .true.
    wy = .true.
    wz = .true.
    wm = .true.
    if (present(want_x)) wx = want_x
    if (present(want_y)) wy = want_y
    if (present(want_z)) wz = want_z
    if (present(want_mag)) wm = want_mag
    if (.not. (wx .or. wy .or. wz .or. wm)) then
       call neko_error("wall_shear_stress: no components selected")
    end if

    call wall_shear_stress_build_field_list(computed_field, wx, wy, wz, wm, &
         fields)

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)

    do i = 1, size(fields)
       call neko_registry%add_field(coef%dof, trim(fields(i)), &
            ignore_existing = .true.)
    end do

    call this%init_common(name, computed_field, viscosity_field, &
         zone_indices, wx, wy, wz, wm, coef)

  end subroutine wall_shear_stress_init_from_controllers_properties

  !> Common part of all constructors.
  !! @param name The unique name of the simcomp.
  !! @param computed_field The base name of the registered fields.
  !! @param viscosity_field The name of the registered viscosity field.
  !! @param zone_indices Labelled zones to include.
  !! @param coef The SEM coefficients.
  subroutine wall_shear_stress_init_common(this, name, computed_field, &
       viscosity_field, zone_indices, want_x, want_y, &
       want_z, want_mag, coef)
    class(wall_shear_stress_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: computed_field
    character(len=*), intent(in) :: viscosity_field
    integer, intent(in) :: zone_indices(:)
    logical, intent(in) :: want_x, want_y, want_z, want_mag
    type(coef_t), intent(inout), target :: coef
    character(len=LOG_SIZE) :: log_buf
    integer :: n_pts, glb_n_pts

    this%name = name
    this%coef => coef
    this%computed_field = computed_field
    this%want_x = want_x
    this%want_y = want_y
    this%want_z = want_z
    this%want_mag = want_mag

    allocate(this%zone_indices(size(zone_indices)))
    this%zone_indices = zone_indices

    this%ale_enabled = .false.
    if (associated(neko_ale)) then
       if (neko_ale%active) this%ale_enabled = .true.
    end if

    this%u => neko_registry%get_field_by_name("u")
    this%v => neko_registry%get_field_by_name("v")
    this%w => neko_registry%get_field_by_name("w")
    this%p => neko_registry%get_field_by_name("p")

    if (.not. neko_registry%field_exists(trim(viscosity_field))) then
       call neko_error("wall_shear_stress: the viscosity field '" // &
            trim(viscosity_field) // "' is not in the registry")
    end if
    this%mu => neko_registry%get_field_by_name(trim(viscosity_field))

    if (want_x) this%tau_x => &
         neko_registry%get_field_by_name(trim(computed_field) // "_x")
    if (want_y) this%tau_y => &
         neko_registry%get_field_by_name(trim(computed_field) // "_y")
    if (want_z) this%tau_z => &
         neko_registry%get_field_by_name(trim(computed_field) // "_z")
    if (want_mag) this%tau_mag => &
         neko_registry%get_field_by_name(trim(computed_field) // "_mag")

    ! `only_facets = .true.`
    ! The normals are requested in the `coef` convention, pointing out of the
    ! fluid domain, because that is what `calc_force_array` expects. The
    ! outward convention would flip the sign of the computed traction.
    call this%bdata%init(this%coef, this%zone_indices, &
         outward_normals = .false.)

    n_pts = this%bdata%n_local

    if (n_pts .gt. 0) then
       call this%mu_msk%init(n_pts)
       call this%p_msk%init(n_pts)
       call this%s11%init(n_pts)
       call this%s22%init(n_pts)
       call this%s33%init(n_pts)
       call this%s12%init(n_pts)
       call this%s13%init(n_pts)
       call this%s23%init(n_pts)
       call this%t1%init(n_pts)
       call this%t2%init(n_pts)
       call this%t3%init(n_pts)
       call this%pt1%init(n_pts)
       call this%pt2%init(n_pts)
       call this%pt3%init(n_pts)
    end if

    glb_n_pts = this%bdata%n_global

    call neko_log%section("Wall shear stress")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,*(I0,:,", "))') "Zone indices: ", this%zone_indices
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Global number of masked points: ", glb_n_pts
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Viscosity field: ", trim(viscosity_field)
    call neko_log%message(log_buf)
    log_buf = "Registered fields: "
    if (this%want_x) log_buf = trim(log_buf) // " " // &
         trim(computed_field) // "_x"
    if (this%want_y) log_buf = trim(log_buf) // " " // &
         trim(computed_field) // "_y"
    if (this%want_z) log_buf = trim(log_buf) // " " // &
         trim(computed_field) // "_z"
    if (this%want_mag) log_buf = trim(log_buf) // " " // &
         trim(computed_field) // "_mag"
    call neko_log%message(log_buf)
    write(log_buf, '(A,L1)') "Moving mesh (ALE): ", this%ale_enabled
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine wall_shear_stress_init_common

  !> Destructor.
  subroutine wall_shear_stress_free(this)
    class(wall_shear_stress_t), intent(inout) :: this

    call this%bdata%free()

    call this%mu_msk%free()
    call this%p_msk%free()
    call this%s11%free()
    call this%s22%free()
    call this%s33%free()
    call this%s12%free()
    call this%s13%free()
    call this%s23%free()
    call this%t1%free()
    call this%t2%free()
    call this%t3%free()
    call this%pt1%free()
    call this%pt2%free()
    call this%pt3%free()

    if (allocated(this%zone_indices)) deallocate(this%zone_indices)
    if (allocated(this%computed_field)) deallocate(this%computed_field)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    nullify(this%mu)
    nullify(this%tau_x)
    nullify(this%tau_y)
    nullify(this%tau_z)
    nullify(this%tau_mag)
    nullify(this%coef)

    call this%free_base()

  end subroutine wall_shear_stress_free

  !> Compute the wall shear stress.
  !! @param time The current time state.
  subroutine wall_shear_stress_compute(this, time)
    class(wall_shear_stress_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    type(field_t), pointer :: fx, fy, fz
    integer :: temp_indices(6)
    integer :: work_indices(3)
    integer :: n, n_pts
    character(len=LOG_SIZE) :: log_buf

    n = this%coef%dof%size()
    n_pts = this%bdata%n_local

    ! A no-op on a static mesh, so no ALE test is needed here.
    call this%bdata%update_geometry(to_host = .false.)

    ! Strain rate over the whole field, then gather it at the mask.
    call neko_scratch_registry%request_field(s11, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(s22, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(s33, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(s12, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(s13, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(s23, temp_indices(6), .false.)

    call strain_rate(s11, s22, s33, s12, s13, s23, this%u, this%v, &
         this%w, this%coef)

    if (n_pts .gt. 0) then
       call this%bdata%get(this%mu, this%mu_msk, to_host = .false.)
       call this%bdata%get(this%p, this%p_msk, to_host = .false.)
       call this%bdata%get(s11, this%s11, to_host = .false.)
       call this%bdata%get(s22, this%s22, to_host = .false.)
       call this%bdata%get(s33, this%s33, to_host = .false.)
       call this%bdata%get(s12, this%s12, to_host = .false.)
       call this%bdata%get(s13, this%s13, to_host = .false.)
       call this%bdata%get(s23, this%s23, to_host = .false.)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_calc_force_array(this%pt1%x_d, this%pt2%x_d, &
               this%pt3%x_d, this%t1%x_d, this%t2%x_d, this%t3%x_d, &
               this%s11%x_d, this%s22%x_d, this%s33%x_d, &
               this%s12%x_d, this%s13%x_d, this%s23%x_d, &
               this%p_msk%x_d, this%bdata%n_x%x_d, &
               this%bdata%n_y%x_d, this%bdata%n_z%x_d, &
               this%mu_msk%x_d, n_pts)
       else
          call calc_force_array(this%pt1%x, this%pt2%x, this%pt3%x, &
               this%t1%x, this%t2%x, this%t3%x, &
               this%s11%x, this%s22%x, this%s33%x, &
               this%s12%x, this%s13%x, this%s23%x, &
               this%p_msk%x, this%bdata%n_x%x, &
               this%bdata%n_y%x, this%bdata%n_z%x, &
               this%mu_msk%x, n_pts)
       end if

       ! Remove the wall-normal part, leaving the tangential traction.
       ! The projection is invariant to the sign of the normal.
       call this%bdata%tangential(this%t1, this%t2, this%t3)
    end if

    call neko_scratch_registry%relinquish_field(temp_indices)


    call neko_scratch_registry%request_field(fx, work_indices(1), .false.)
    call neko_scratch_registry%request_field(fy, work_indices(2), .false.)
    call neko_scratch_registry%request_field(fz, work_indices(3), .false.)

    ! Scatter the vectors into fields.
    call this%bdata%scatter(this%t1, fx)
    call this%bdata%scatter(this%t2, fy)
    call this%bdata%scatter(this%t3, fz)

    ! Copy the work fields into the selected registered fields.
    if (this%want_x) call field_copy(this%tau_x, fx)
    if (this%want_y) call field_copy(this%tau_y, fy)
    if (this%want_z) call field_copy(this%tau_z, fz)

    ! The magnitude
    if (this%want_mag) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_vdot3(this%tau_mag%x_d, fx%x_d, fy%x_d, fz%x_d, &
               fx%x_d, fy%x_d, fz%x_d, n)
          call device_sqrt_inplace(this%tau_mag%x_d, n)
       else
          call vdot3(this%tau_mag%x, fx%x, fy%x, fz%x, fx%x, fy%x, fz%x, n)
          call sqrt_inplace(this%tau_mag%x, n)
       end if
    end if

    call neko_scratch_registry%relinquish_field(work_indices)

    write(log_buf, '(A,A,A,E15.7,A,*(I0,:,", "))') &
         "WSS: '", trim(this%name), "' computed at t = ", &
         time%t, ", zones: ", this%zone_indices
    call neko_log%message(log_buf)

  end subroutine wall_shear_stress_compute

end module wall_shear_stress_simcomp
