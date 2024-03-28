! Copyright (c) 2023, The Neko Authors
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
!> Implements the `brinkman_source_term_t` type.
module brinkman_source_term
  use num_types, only: rp
  use field, only: field_t
  use field_list, only: field_list_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use field_registry, only: neko_field_registry
  use source_term, only: source_term_t
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use utils, only: neko_error
  use brinkman_source_term_cpu, only: brinkman_source_term_compute_cpu
  use brinkman_source_term_device, only: brinkman_source_term_compute_device
  use brinkman_source_term_cpu, only: implicit_brinkman_source_term_compute_cpu, brinkman_source_term_compute_cpu
  use brinkman_source_term_device, only: brinkman_source_term_compute_device,  implicit_brinkman_source_term_compute_device
  implicit none
  private

  !> A Brinkman source term.
  !! The region and strength are controlled by assigning regions types and
  !! brinkman limits to the source term.
  type, public, extends(source_term_t) :: brinkman_source_term_t
     private

     !> The value of the source term.
     type(field_t), pointer :: indicator => null()
     !> Brinkman permeability field.
     type(field_t), pointer :: brinkman => null()
     ! > Flag for implicit vs explicit
     logical ::if_implicit
   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init => brinkman_source_term_init_from_json
     !> Destructor.
     procedure, public, pass(this) :: free => brinkman_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: compute_ => brinkman_source_term_compute

     ! ----------------------------------------------------------------------- !
     ! Private methods
     procedure, pass(this) :: init_boundary_mesh
     procedure, pass(this) :: init_point_zone
  end type brinkman_source_term_t

contains

  ! ========================================================================== !
  ! Public methods

  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine brinkman_source_term_init_from_json(this, json, fields, coef)
    use file, only: file_t
    use tri_mesh, only: tri_mesh_t
    use device, only: device_memcpy, HOST_TO_DEVICE
    use filters, only: smooth_step_field, step_function_field, permeability_field
    use signed_distance, only: signed_distance_field
    use profiler, only: profiler_start_region, profiler_end_region
    use json_module, only: json_core, json_value
    implicit none

    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: start_time, end_time

    character(len=:), allocatable :: filter_type
    real(kind=rp), dimension(:), allocatable :: brinkman_limits
    logical :: brinkman_implicit
    real(kind=rp) :: brinkman_penalty

    type(json_value), pointer :: json_object_list
    type(json_core) :: core
    type(json_value), pointer :: json_object_pointer
    character(len=:), allocatable :: json_buffer

    character(len=:), allocatable :: object_type
    type(json_file) :: object_settings
    integer :: n_regions
    integer :: i
    logical :: found



    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    ! Read the options for the permeability field
    call json_get(json, 'brinkman.limits', brinkman_limits)
    call json_get(json, 'brinkman.penalty', brinkman_penalty)

    ! Read implicit vs explicit formulation
    call json_get(json, 'brinkman.implicit', brinkman_implicit)
    this%if_implicit = brinkman_implicit

    if (size(brinkman_limits) .ne. 2) then
       call neko_error('brinkman_limits must be a 2 element array of reals')
    end if

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    ! ------------------------------------------------------------------------ !
    ! Allocate the permeability and indicator field

    if (neko_field_registry%field_exists('brinkman_indicator') &
        .or. neko_field_registry%field_exists('brinkman')) then
       call neko_error('Brinkman field already exists.')
    end if

    call neko_field_registry%add_field(coef%dof, 'brinkman_indicator')
    this%indicator => neko_field_registry%get_field_by_name('brinkman_indicator')

    call neko_field_registry%add_field(coef%dof, 'brinkman')
    this%brinkman => neko_field_registry%get_field_by_name('brinkman')

    ! ------------------------------------------------------------------------ !
    ! Select which constructor should be called

    call json%get('objects', json_object_list)
    n_regions = core%count(json_object_list)

    do i=1, n_regions

       ! Create a new json containing just the subdict for this object.
       call core%get_child(json_object_list, i, json_object_pointer)
       call core%print_to_string(json_object_pointer, json_buffer)
       call object_settings%load_from_string(json_buffer)


       call json_get(object_settings, 'type', object_type)
       select case (object_type)
         case ('boundary_mesh')
          call this%init_boundary_mesh(object_settings)
         case ('point_zone')
          call this%init_point_zone(object_settings)
         case default
          call neko_error('Unknown region type')
       end select

    end do

    ! Run filter on the full indicator field to smooth it out.
    call json_get_or_default(json, 'filter.type', filter_type, 'none')

    select case (filter_type)
      case ('none')
       ! Do nothing
      case default
       call neko_error('Unknown filter type')
    end select

    ! ------------------------------------------------------------------------ !
    ! Compute the permeability field

    call permeability_field(this%brinkman, this%indicator, &
      & brinkman_limits(1), brinkman_limits(2), brinkman_penalty)

    ! Copy the permeability field to the device
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%brinkman%x, this%brinkman%x_d, &
                          this%brinkman%dof%size(), HOST_TO_DEVICE, .true.)
    end if

  end subroutine brinkman_source_term_init_from_json

  !> Destructor.
  subroutine brinkman_source_term_free(this)
    class(brinkman_source_term_t), intent(inout) :: this

    this%brinkman => null()
    call this%free_base()
  end subroutine brinkman_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine brinkman_source_term_compute(this, t, tstep)
    use device, only: device_memcpy, HOST_TO_DEVICE
    implicit none

    class(brinkman_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

	if(this%if_implicit) then
   ! this is implicit
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call implicit_brinkman_source_term_compute_device(this%fields, this%brinkman)
    else
       call implicit_brinkman_source_term_compute_cpu(this%fields, this%brinkman)
    end if
   else
	! this is all explicit
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call brinkman_source_term_compute_device(this%fields, this%brinkman)
    else
       call brinkman_source_term_compute_cpu(this%fields, this%brinkman)
    end if
   endif
  end subroutine brinkman_source_term_compute

  ! ========================================================================== !
  ! Private methods

  !> Initializes the source term from a boundary mesh.
  subroutine init_boundary_mesh(this, json)
    use file, only: file_t
    use tri_mesh, only: tri_mesh_t
    use device, only: device_memcpy, HOST_TO_DEVICE
    use filters, only: smooth_step_field, step_function_field, permeability_field
    use signed_distance, only: signed_distance_field
    use profiler, only: profiler_start_region, profiler_end_region
    implicit none

    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    ! Options
    character(len=:), allocatable :: mesh_file_name
    character(len=:), allocatable :: distance_transform
    character(len=:), allocatable :: filter_type

    type(file_t) :: mesh_file
    type(tri_mesh_t) :: boundary_mesh
    real(kind=rp) :: scalar

    type(field_t) :: temp_field

    ! ------------------------------------------------------------------------ !
    ! Read the options for the boundary mesh

    call json_get(json, 'name', mesh_file_name)

    ! Settings on how to filter the design field
    call json_get(json, 'distance_transform.type', distance_transform)

    ! ------------------------------------------------------------------------ !
    ! Load the immersed boundary mesh

    mesh_file = file_t(mesh_file_name)
    call mesh_file%read(boundary_mesh)

    if (boundary_mesh%nelv .eq. 0) then
       call neko_error('No elements in the boundary mesh')
    end if

    ! ------------------------------------------------------------------------ !
    ! Compute the permeability field

    ! Assign the signed distance field to all GLL points in the permeability
    ! field. Initally we just run a brute force loop over all GLL points and
    ! compute the signed distance function. This should be replaced with a
    ! more efficient method, such as a tree search.

    call temp_field%init(this%indicator%dof)

    ! Select how to transform the distance field to a design field
    select case (distance_transform)
      case ('smooth_step')
       call json_get(json, 'distance_transform.value', scalar)

       call signed_distance_field(temp_field, boundary_mesh, scalar)
       call smooth_step_field(temp_field, scalar, 0.0_rp)

      case ('step')

       call json_get(json, 'distance_transform.value', scalar)

       call signed_distance_field(temp_field, boundary_mesh, scalar)
       call step_function_field(temp_field, scalar, 1.0_rp, 0.0_rp)

      case default
       call neko_error('Unknown distance transform')
    end select

    ! ------------------------------------------------------------------------ !
    ! Run filter on the temporary indicator field to smooth it out.
    call json_get_or_default(json, 'filter.type', filter_type, 'none')

    select case (filter_type)
      case ('none')
       ! Do nothing
      case default
       call neko_error('Unknown filter type')
    end select

    ! Update the global indicator field by max operator
    this%indicator%x = max(this%indicator%x, temp_field%x)

  end subroutine init_boundary_mesh

  !> Initializes the source term from a point zone.
  subroutine init_point_zone(this, json)
    use filters, only: smooth_step_field, step_function_field, permeability_field
    use signed_distance, only: signed_distance_field
    use profiler, only: profiler_start_region, profiler_end_region
    use point_zone, only: point_zone_t
    use device, only: device_memcpy, HOST_TO_DEVICE
    use point_zone_registry, only: neko_point_zone_registry
    implicit none

    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    ! Options
    character(len=:), allocatable :: zone_name
    character(len=:), allocatable :: filter_type

    type(field_t) :: temp_field
    class(point_zone_t), pointer :: my_point_zone
    integer :: i

    ! ------------------------------------------------------------------------ !
    ! Read the options for the point zone

    call json_get(json,'name', zone_name)
    call json_get_or_default(json, 'filter.type', filter_type, 'none')

    ! Compute the indicator field

    call temp_field%init(this%indicator%dof)

    my_point_zone => neko_point_zone_registry%get_point_zone(zone_name)

    do i = 1, my_point_zone%size
       temp_field%x(my_point_zone%mask(i), 1, 1, 1) = 1.0_rp
    end do

    ! Run filter on the temporary indicator field to smooth it out.

    select case (filter_type)
      case ('none')
       ! Do nothing
      case default
       call neko_error('Unknown filter type')
    end select

    ! Update the global indicator field by max operator
    this%indicator%x = max(this%indicator%x, temp_field%x)

  end subroutine init_point_zone

end module brinkman_source_term
