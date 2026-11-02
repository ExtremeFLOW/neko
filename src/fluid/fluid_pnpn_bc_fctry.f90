
! Copyright (c) 2024-2026, The Neko Authors
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
!> Defines the boundary condition factories, allocators and user type
!! registries for `fluid_pnpn_t`.
submodule(fluid_pnpn) fluid_pnpn_bc_fctry
  ! We explicitly import all necessary modules to work around an ifx 2026.1
  ! internal compiler error even when they are only available through the
  ! parent. Generic interfaces (neko_error, json_get, json_file, ...) are the
  ! exception, since gfortran rejects importing them a second time, so they
  ! are only accessed through the parent.
  use vector_bc_projector, only : segregated_vector_bc_projector_t, &
       coupled_vector_bc_projector_t
  use user_intf, only : user_t
  use utils, only : neko_type_error, neko_type_registration_error
  use json_module, only : json_array
  use neumann, only : neumann_t
  use field_dirichlet, only : field_dirichlet_t
  use inflow, only : inflow_t
  use blasius, only : blasius_t
  use dirichlet, only : dirichlet_t
  use dong_outflow, only : dong_outflow_t
  use symmetry_aligned, only : symmetry_aligned_t
  use symmetry, only : symmetry_t
  use non_normal_aligned, only : non_normal_aligned_t
  use non_normal, only : non_normal_t
  use no_slip, only : no_slip_t
  use zero_dirichlet, only : zero_dirichlet_t
  use shear_stress, only : shear_stress_t
  use wall_model_bc, only : wall_model_bc_t
  use field_dirichlet_vector, only : field_dirichlet_vector_t
  use expression_dirichlet, only : expression_dirichlet_t
  use expression_dirichlet_vector, only : expression_dirichlet_vector_t
  use overset_interface, only : overset_interface_t
  use overset_interface_vector, only : overset_interface_vector_t
  implicit none

  !> A joint boundary condition type, i.e. a value of the `type` keyword in
  !! a boundary condition entry of the case file, and the velocity and pressure
  !! condition types it expands to. An empty string means that no condition is
  !! created for that field.
  type :: fluid_pnpn_bc_preset_t
     character(len=25) :: name
     character(len=25) :: velocity
     character(len=25) :: pressure
  end type fluid_pnpn_bc_preset_t

  !> The joint boundary condition types accepted by the `type` keyword.
  !! To add a new one, append a row with the type name, the velocity condition
  !! type and the pressure condition type. When Neumann is zero valued
  !! it amounts to do nothing, and is encoded as "".
  type(fluid_pnpn_bc_preset_t), parameter :: FLUID_PNPN_BC_PRESETS(17) = [ &
       fluid_pnpn_bc_preset_t("symmetry", "symmetry", ""), &
       fluid_pnpn_bc_preset_t("velocity_value", "dirichlet", ""), &
       fluid_pnpn_bc_preset_t("expression_velocity", "expression", ""), &
       fluid_pnpn_bc_preset_t("expression_pressure", "", "expression"), &
       fluid_pnpn_bc_preset_t("no_slip", "no_slip", ""), &
       fluid_pnpn_bc_preset_t("outflow", "", "zero_dirichlet"), &
       fluid_pnpn_bc_preset_t("normal_outflow", "non_normal", &
       "zero_dirichlet"), &
       fluid_pnpn_bc_preset_t("outflow+user", "", "user_dirichlet"), &
       fluid_pnpn_bc_preset_t("normal_outflow+user", "non_normal", &
       "user_dirichlet"), &
       fluid_pnpn_bc_preset_t("outflow+dong", "", "dong"), &
       fluid_pnpn_bc_preset_t("normal_outflow+dong", "non_normal", "dong"), &
       fluid_pnpn_bc_preset_t("shear_stress", "shear_stress", ""), &
       fluid_pnpn_bc_preset_t("user_velocity", "user_dirichlet", ""), &
       fluid_pnpn_bc_preset_t("user_pressure", "", "user_dirichlet"), &
       fluid_pnpn_bc_preset_t("blasius_profile", "blasius_profile", ""), &
       fluid_pnpn_bc_preset_t("wall_model", "wall_model", ""), &
       fluid_pnpn_bc_preset_t("overset_interface", "overset_interface", &
       "overset_interface")]

  !> The built-in velocity boundary condition types.
  character(len=25), parameter :: FLUID_PNPN_VELOCITY_BCS(11) = [ &
       character(len=25) :: &
       "dirichlet", &
       "expression", &
       "no_slip", &
       "symmetry", &
       "non_normal", &
       "neumann", &
       "shear_stress", &
       "wall_model", &
       "blasius_profile", &
       "user_dirichlet", &
       "overset_interface"]

  !> The built-in pressure boundary condition types.
  character(len=25), parameter :: FLUID_PNPN_PRESSURE_BCS(7) = [ &
       character(len=25) :: &
       "dirichlet", &
       "zero_dirichlet", &
       "expression", &
       "dong", &
       "neumann", &
       "user_dirichlet", &
       "overset_interface"]

contains

  !> Split a boundary condition entry of the case file into a velocity part
  !! and a pressure part.
  !! @details An entry either has a joint `type` keyword, which is expanded
  !! using `FLUID_PNPN_BC_PRESETS`, or both a `velocity` and a `pressure`
  !! object, each with its own `type`. Mixing the two forms is an error. All
  !! other keywords of the entry, e.g. `zone_indices` and `name`, are copied
  !! to both parts.
  !! @param[inout] json The boundary condition entry.
  !! @param[inout] velocity_json The velocity part, valid if `has_velocity`.
  !! @param[inout] pressure_json The pressure part, valid if `has_pressure`.
  !! @param[out] has_velocity Whether the entry defines a velocity condition.
  !! @param[out] has_pressure Whether the entry defines a pressure condition.
  module subroutine fluid_pnpn_bc_split(json, velocity_json, pressure_json, &
       has_velocity, has_pressure)
    type(json_file), intent(inout) :: json
    type(json_file), intent(inout) :: velocity_json
    type(json_file), intent(inout) :: pressure_json
    logical, intent(out) :: has_velocity
    logical, intent(out) :: has_pressure
    character(len=:), allocatable :: type
    character(len=25) :: velocity_type, pressure_type
    logical :: couple_pressure
    integer :: i

    velocity_type = ""
    pressure_type = ""
    has_velocity = json%valid_path("velocity")
    has_pressure = json%valid_path("pressure")

    if (json%valid_path("type")) then
       if (has_velocity .or. has_pressure) then
          call neko_error("A fluid boundary condition cannot have both a " // &
               "'type' keyword and a 'velocity' or 'pressure' object.")
       end if

       call json_get(json, "type", type)

       do i = 1, size(FLUID_PNPN_BC_PRESETS)
          if (trim(type) .eq. trim(FLUID_PNPN_BC_PRESETS(i)%name)) exit
       end do

       if (i .gt. size(FLUID_PNPN_BC_PRESETS)) then
          call neko_type_error("fluid_pnpn boundary conditions", type, &
               FLUID_PNPN_BC_PRESETS%name)
       end if

       velocity_type = FLUID_PNPN_BC_PRESETS(i)%velocity
       pressure_type = FLUID_PNPN_BC_PRESETS(i)%pressure

       ! The overset interface only constrains the pressure on request.
       if (trim(type) .eq. "overset_interface") then
          call json_get_or_default(json, "couple_pressure", couple_pressure, &
               .false.)
          if (.not. couple_pressure) pressure_type = ""
       end if

       has_velocity = len_trim(velocity_type) .gt. 0
       has_pressure = len_trim(pressure_type) .gt. 0
    else if (.not. (has_velocity .and. has_pressure)) then
       call neko_error("A fluid boundary condition needs either a 'type' " // &
            "keyword, or both a 'velocity' and a 'pressure' object.")
    end if

    if (has_velocity) then
       call fluid_pnpn_bc_extract_part(json, "velocity", velocity_type, &
            velocity_json)
    end if

    if (has_pressure) then
       call fluid_pnpn_bc_extract_part(json, "pressure", pressure_type, &
            pressure_json)
    end if
  end subroutine fluid_pnpn_bc_split

  !> Build the JSON object for one field of a boundary condition entry.
  !! @details `part` receives all keywords of `entry` except `type`,
  !! `velocity` and `pressure`. Then either `type` is set to `preset_type`,
  !! for an entry with a joint type, or the contents of the `key` object in
  !! `entry` are copied into `part`, overriding any keyword already present.
  !! `entry` itself is not modified.
  !! @param[inout] entry The boundary condition entry in the case file.
  !! @param[in] key Either `velocity` or `pressure`.
  !! @param[in] preset_type The condition type from the joint type, or empty.
  !! @param[inout] part The JSON object to build.
  subroutine fluid_pnpn_bc_extract_part(entry, key, preset_type, part)
    type(json_file), intent(inout) :: entry
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: preset_type
    type(json_file), intent(inout) :: part
    type(json_core) :: core
    type(json_value), pointer :: root, sub, child, existing, child_clone
    character(len=:), allocatable :: buffer, child_name
    logical :: found
    integer :: i, n

    ! Start from an independent copy of the entry
    call entry%print_to_string(buffer)
    call part%destroy()
    call part%initialize(strict_type_checking = .true.)
    call part%load_from_string(buffer)

    ! Retain only the top-level stuff like the zone indices
    call part%remove("type")
    call part%remove("velocity")
    call part%remove("pressure")

    if (len_trim(preset_type) .gt. 0) then
       call part%add("type", trim(preset_type))
    end if

    ! The sub-object is read from the unmodified entry, not from the copy.
    call entry%get(key, sub, found)
    if (found) then
       call part%get_core(core)
       call part%get(root)
       n = core%count(sub)
       do i = 1, n
          call core%get_child(sub, i, child, found)
          call core%info(child, name = child_name)

          call core%get_child(root, child_name, existing, found)
          if (found) call core%remove(existing, destroy = .true.)

          nullify(child_clone)
          call core%clone(child, child_clone)
          call core%add(root, child_clone)
       end do
    end if
  end subroutine fluid_pnpn_bc_extract_part

  !> Allocator for velocity boundary conditions.
  !! @param[inout] object The object to be allocated.
  !! @param[in] type_name The name of the boundary condition type.
  !! @param[in] full_stress_formulation Whether the scheme uses the full
  !! viscous stress formulation, which selects the coupled implementations of
  !! the mixed conditions.
  module subroutine fluid_pnpn_velocity_bc_allocator(object, type_name, &
       full_stress_formulation)
    class(bc_t), pointer, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    logical, intent(in) :: full_stress_formulation
    integer :: i

    if (associated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ("dirichlet")
       allocate(inflow_t::object)
    case ("expression")
       allocate(expression_dirichlet_vector_t::object)
    case ("neumann")
       allocate(neumann_t::object)
    case ("no_slip")
       allocate(no_slip_t::object)
    case ("symmetry")
       if (full_stress_formulation) then
          allocate(symmetry_t::object)
       else
          allocate(symmetry_aligned_t::object)
       end if
    case ("non_normal")
       if (full_stress_formulation) then
          allocate(non_normal_t::object)
       else
          allocate(non_normal_aligned_t::object)
       end if
    case ("shear_stress")
       allocate(shear_stress_t::object)
    case ("wall_model")
       allocate(wall_model_bc_t::object)
    case ("blasius_profile")
       allocate(blasius_t::object)
    case ("user_dirichlet")
       allocate(field_dirichlet_vector_t::object)
    case ("overset_interface")
       allocate(overset_interface_vector_t::object)
    case default
       do i = 1, fluid_pnpn_velocity_bc_registry_size
          if (trim(type_name) .eq. &
               trim(fluid_pnpn_velocity_bc_registry(i)%type_name)) then
             call fluid_pnpn_velocity_bc_registry(i)%allocator(object)
             return
          end if
       end do

       call neko_type_error("fluid_pnpn velocity boundary conditions", &
            type_name, FLUID_PNPN_VELOCITY_BCS)
    end select
  end subroutine fluid_pnpn_velocity_bc_allocator

  !> Allocator for pressure boundary conditions.
  !! @param[inout] object The object to be allocated.
  !! @param[in] type_name The name of the boundary condition type.
  module subroutine fluid_pnpn_pressure_bc_allocator(object, type_name)
    class(bc_t), pointer, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    integer :: i

    if (associated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ("dirichlet")
       allocate(dirichlet_t::object)
    case ("zero_dirichlet")
       allocate(zero_dirichlet_t::object)
    case ("expression")
       allocate(expression_dirichlet_t::object)
    case ("dong")
       allocate(dong_outflow_t::object)
    case ("user_dirichlet")
       allocate(field_dirichlet_t::object)
    case ("overset_interface")
       allocate(overset_interface_t::object)
    case default
       do i = 1, fluid_pnpn_pressure_bc_registry_size
          if (trim(type_name) .eq. &
               trim(fluid_pnpn_pressure_bc_registry(i)%type_name)) then
             call fluid_pnpn_pressure_bc_registry(i)%allocator(object)
             return
          end if
       end do

       call neko_type_error("fluid_pnpn pressure boundary conditions", &
            type_name, FLUID_PNPN_PRESSURE_BCS)
    end select
  end subroutine fluid_pnpn_pressure_bc_allocator

  !> Register a custom velocity boundary condition allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param[in] type_name The name of the boundary condition type.
  !! @param[in] allocator The allocator for the custom user type.
  module subroutine register_fluid_pnpn_velocity_bc(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(fluid_pnpn_bc_allocate), pointer, intent(in) :: allocator

    call fluid_pnpn_bc_register(fluid_pnpn_velocity_bc_registry, &
         fluid_pnpn_velocity_bc_registry_size, FLUID_PNPN_VELOCITY_BCS, &
         "fluid_pnpn velocity boundary condition", type_name, allocator)
  end subroutine register_fluid_pnpn_velocity_bc

  !> Register a custom pressure boundary condition allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param[in] type_name The name of the boundary condition type.
  !! @param[in] allocator The allocator for the custom user type.
  module subroutine register_fluid_pnpn_pressure_bc(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(fluid_pnpn_bc_allocate), pointer, intent(in) :: allocator

    call fluid_pnpn_bc_register(fluid_pnpn_pressure_bc_registry, &
         fluid_pnpn_pressure_bc_registry_size, FLUID_PNPN_PRESSURE_BCS, &
         "fluid_pnpn pressure boundary condition", type_name, allocator)
  end subroutine register_fluid_pnpn_pressure_bc

  !> Add an allocator to one of the boundary condition registries.
  !! @param[inout] registry The registry to add to.
  !! @param[inout] registry_size The number of entries in `registry`.
  !! @param[in] known_bcs The built-in type names of the same field, which a
  !! custom type must not clash with.
  !! @param[in] base_type Description of the registry for error messages.
  !! @param[in] type_name The name of the boundary condition type.
  !! @param[in] allocator The allocator for the custom user type.
  subroutine fluid_pnpn_bc_register(registry, registry_size, known_bcs, &
       base_type, type_name, allocator)
    type(fluid_pnpn_bc_allocator_entry), allocatable, intent(inout) :: &
         registry(:)
    integer, intent(inout) :: registry_size
    character(len=*), intent(in) :: known_bcs(:)
    character(len=*), intent(in) :: base_type
    character(len=*), intent(in) :: type_name
    procedure(fluid_pnpn_bc_allocate), pointer, intent(in) :: allocator
    type(fluid_pnpn_bc_allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(known_bcs)
       if (trim(type_name) .eq. trim(known_bcs(i))) then
          call neko_type_registration_error(base_type, type_name, .true.)
       end if
    end do

    do i = 1, registry_size
       if (trim(type_name) .eq. trim(registry(i)%type_name)) then
          call neko_type_registration_error(base_type, type_name, .false.)
       end if
    end do

    if (registry_size .eq. 0) then
       allocate(registry(1))
    else
       allocate(temp(registry_size + 1))
       temp(1:registry_size) = registry
       call move_alloc(temp, registry)
    end if

    registry_size = registry_size + 1
    registry(registry_size)%type_name = type_name
    registry(registry_size)%allocator => allocator
  end subroutine fluid_pnpn_bc_register

  !> Factory routine for pressure boundary conditions.
  !! @details Allocates the condition based on the `type` keyword, connects
  !! the user routines where needed, initializes the object from `json`, marks
  !! the zones in `zone_indices`, sets the name and finalizes. The `neumann`
  !! type is the Neumann condition of the Pn/Pn scheme, with the normal
  !! gradient given by the momentum equation. It needs no object, so `object`
  !! is left unassociated for it.
  !! @param[inout] object The boundary condition to be allocated.
  !! @param[in] scheme The `fluid_pnpn_t` scheme.
  !! @param[inout] json The pressure part of the boundary condition entry.
  !! @param[in] coef The SEM coeffcients.
  !! @param[in] user The user interface.
  module subroutine pressure_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(fluid_pnpn_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    type(user_t), target, intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i
    integer, allocatable :: zone_indices(:)
    character(len=:), allocatable :: default_name
    character(len=64) :: buf

    call json_get(json, "type", type)

    if (trim(type) .eq. "neumann") then
       ! The flux is set by the scheme and is not a user input.
       if (json%valid_path("flux")) then
          call neko_error("The 'neumann' pressure boundary condition " // &
               "does not accept a 'flux' keyword, the flux is set by the " // &
               "scheme.")
       end if

       if (associated(object)) then
          call object%free()
          deallocate(object)
       end if
       return
    end if

    call fluid_pnpn_pressure_bc_allocator(object, type)

    select case (trim(type))
    case ("user_dirichlet")
       select type (obj => object)
       type is (field_dirichlet_t)
          obj%update => user%dirichlet_conditions
          call json%add("field_name", scheme%p%name)
       end select
    case ("overset_interface")
       select type (obj => object)
       type is (overset_interface_t)
          call json%add("field_name", scheme%p%name)
          obj%morph_interface => user%morph_interface
       end select
    end select

    call json_get_or_lookup(json, "zone_indices", zone_indices)
    call object%init(coef, json)

    do i = 1, size(zone_indices)
       call object%mark_labeled_zone(zone_indices(i))
    end do

    write(buf, '("pressure_bc_", I0)') zone_indices(1)
    default_name = trim(buf)
    call json_get_or_default(json, "name", object%name, default_name)
    object%zone_indices = zone_indices
    call object%finalize()

    deallocate(type)
    deallocate(zone_indices)
    deallocate(default_name)
  end subroutine pressure_bc_factory

  !> Factory routine for velocity boundary conditions.
  !! @details Allocates the condition based on the `type` keyword, connects
  !! the user routines where needed, initializes the object from `json`, marks
  !! the zones in `zone_indices`, sets the name and finalizes. The `neumann`
  !! type without a `flux` keyword is a homogeneous Neumann condition, i.e. no
  !! constraint, which needs no object, so `object` is left unassociated for
  !! it. With a `flux` keyword, which must be an array of 3 reals, a
  !! `neumann_t` is created that adds the flux to the momentum right-hand side.
  !! @param[inout] object The boundary condition to be allocated.
  !! @param[in] scheme The `fluid_pnpn_t` scheme.
  !! @param[inout] json The velocity part of the boundary condition entry.
  !! @param[in] coef The SEM coeffcients.
  !! @param[in] user The user interface.
  module subroutine velocity_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(fluid_pnpn_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    type(user_t), target, intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i
    integer, allocatable :: zone_indices(:)
    character(len=:), allocatable :: default_name
    character(len=:), allocatable :: bc_name
    character(len=64) :: buf
    integer :: var_type, n_children

    call json_get(json, "type", type)

    if (trim(type) .eq. "neumann") then
       if (.not. json%valid_path("flux")) then
          if (associated(object)) then
             call object%free()
             deallocate(object)
          end if
          return
       end if

       ! neumann_t applies the flux component-wise to a vector, so a scalar
       ! flux is not enough.
       call json%info("flux", var_type = var_type, n_children = n_children)
       if (var_type .ne. json_array .or. n_children .ne. 3) then
          call neko_error("The 'flux' of the 'neumann' velocity boundary " // &
               "condition must be an array of 3 reals.")
       end if
    end if

    call fluid_pnpn_velocity_bc_allocator(object, type, &
         scheme%full_stress_formulation)

    select case (trim(type))
    case ("wall_model")
       ! Kind of hack, but  OK for now
       call json%add("scheme_name", scheme%name)

       select type (wall_bc => object)
       type is (wall_model_bc_t)
          wall_bc%user => user
       end select
    case ("user_dirichlet")
       select type (obj => object)
       type is (field_dirichlet_vector_t)
          obj%update => user%dirichlet_conditions
       end select
    case ("overset_interface")
       select type (obj => object)
       type is (overset_interface_vector_t)
          obj%morph_interface => user%morph_interface
       end select
    end select

    call json_get_or_lookup(json, "zone_indices", zone_indices)
    write(buf,'("velocity_bc_",I0)') zone_indices(1)
    default_name = trim(buf)
    call json_get_or_default(json, "name", bc_name, default_name)

    call object%init(coef, json)
    do i = 1, size(zone_indices)
       call object%mark_labeled_zone(zone_indices(i))
    end do

    object%name = bc_name
    object%zone_indices = zone_indices

    call object%finalize()

    deallocate(type)
    deallocate(zone_indices)
    deallocate(default_name)
    deallocate(bc_name)
  end subroutine velocity_bc_factory

end submodule fluid_pnpn_bc_fctry
