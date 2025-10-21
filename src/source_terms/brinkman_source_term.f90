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
  use aabb, only : aabb_t, get_aabb
  use coefs, only : coef_t
  use device, only : device_memcpy, HOST_TO_DEVICE
  use field, only : field_t
  use field_list, only : field_list_t
  use math, only : cfill_mask, pwmax2
  use device_math, only : device_cfill_mask, device_pwmax2
  use field_math, only : field_pwmax2, field_subcol3, field_copy
  use field_registry, only : neko_field_registry
  use mappings, only : smooth_step_field, step_function_field, &
       permeability_field
  use file, only : file_t
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use tri_mesh, only : tri_mesh_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp, dp
  use point_zone, only : point_zone_t
  use point_zone_registry, only : neko_point_zone_registry
  use profiler, only : profiler_start_region, profiler_end_region
  use signed_distance, only : signed_distance_field
  use source_term, only : source_term_t
  use utils, only : neko_error
  use filter, only : filter_t
  use PDE_filter, only : PDE_filter_t
  use fld_file_output, only : fld_file_output_t
  use fld_file_data, only : fld_file_data_t
  use num_types, only : sp, dp
  use time_state, only : time_state_t

  use global_interpolation, only: global_interpolation_t
  use interpolation, only: interpolator_t
  use space, only: space_t, GLL
  implicit none
  private

  !> A Brinkman source term.
  !! The region and strength are controlled by assigning regions types and
  !! brinkman limits to the source term.
  type, public, extends(source_term_t) :: brinkman_source_term_t
     private

     !> The unfiltered indicator field
     type(field_t), pointer :: indicator_unfiltered
     !> The value of the source term.
     type(field_t), pointer :: indicator
     !> Brinkman permeability field.
     type(field_t), pointer :: brinkman
     !> Filter
     class(filter_t), allocatable :: filter
   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init => &
          brinkman_source_term_init_from_json
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
  !! @param variable_name The name of the variable where the source term acts.
  subroutine brinkman_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp) :: start_time, end_time

    character(len=:), allocatable :: filter_type
    real(kind=rp) :: filter_radius
    real(kind=rp), dimension(:), allocatable :: brinkman_limits
    real(kind=rp) :: brinkman_penalty

    type(json_value), pointer :: json_object_list
    type(json_core) :: core

    character(len=:), allocatable :: object_type
    type(json_file) :: object_settings
    integer :: n_regions
    integer :: i
    type(fld_file_output_t) :: output


    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    ! Read the options for the permeability field
    call json_get(json, 'brinkman.limits', brinkman_limits)
    call json_get(json, 'brinkman.penalty', brinkman_penalty)

    if (size(brinkman_limits) .ne. 2) then
       call neko_error('brinkman_limits must be a 2 element array of reals')
    end if

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    ! ------------------------------------------------------------------------ !
    ! Allocate the permeability and indicator field

    call neko_field_registry%add_field(coef%dof, 'brinkman_indicator', .true.)
    call neko_field_registry%add_field(coef%dof, 'brinkman_indicator_unfiltered', &
         .true.)
    call neko_field_registry%add_field(coef%dof, 'brinkman_permeability', &
         .true.)

    this%indicator => neko_field_registry%get_field('brinkman_indicator')
    this%indicator_unfiltered => &
         neko_field_registry%get_field('brinkman_indicator_unfiltered')
    this%brinkman => neko_field_registry%get_field('brinkman_permeability')

    ! ------------------------------------------------------------------------ !
    ! Select which constructor should be called

    call json%info('objects', n_children = n_regions)

    do i = 1, n_regions
       call json_extract_item(json, "objects", i, object_settings)
       call json_get_or_default(object_settings, 'type', object_type, 'none')

       select case (object_type)
       case ('boundary_mesh')
          call this%init_boundary_mesh(object_settings)
       case ('point_zone')
          call this%init_point_zone(object_settings)

       case ('none')
          call object_settings%print()
          call neko_error('Brinkman source term objects require a region type')
       case default
          call neko_error('Brinkman source term unknown region type')
       end select

    end do

    ! ------------------------------------------------------------------------ !
    ! Filter the indicator field

    call json_get_or_default(json, 'filter.type', filter_type, 'none')
    select case (filter_type)
    case ('PDE')
       ! Initialize the unfiltered design field
       call this%indicator_unfiltered%init(coef%dof)

       ! Allocate a PDE filter
       allocate(PDE_filter_t::this%filter)

       ! Initialize the filter
       call this%filter%init(json, coef)

       ! Copy the current indicator to unfiltered (essentially a rename)
       call field_copy(this%indicator_unfiltered, this%indicator)

       ! Apply the filter
       call this%filter%apply(this%indicator, this%indicator_unfiltered)

       ! Set up sampler to include the unfiltered and filtered fields
       call output%init(sp, 'brinkman', 3)
       call output%fields%assign_to_field(1, this%indicator_unfiltered)
       call output%fields%assign_to_field(2, this%indicator)
       call output%fields%assign_to_field(3, this%brinkman)

    case ('none')
       ! Set up sampler to include the unfiltered field
       call output%init(sp, 'brinkman', 2)
       call output%fields%assign_to_field(1, this%indicator)
       call output%fields%assign_to_field(2, this%brinkman)

    case default
       call neko_error('Brinkman source term unknown filter type')
    end select

    ! ------------------------------------------------------------------------ !
    ! Compute the permeability field

    this%brinkman = this%indicator
    call permeability_field(this%brinkman, &
         brinkman_limits(1), brinkman_limits(2), brinkman_penalty)

    ! Sample the Brinkman field
    call output%sample(0.0_rp)

  end subroutine brinkman_source_term_init_from_json

  !> Destructor.
  subroutine brinkman_source_term_free(this)
    class(brinkman_source_term_t), intent(inout) :: this

    nullify(this%indicator)
    nullify(this%indicator_unfiltered)
    nullify(this%brinkman)

    if (allocated(this%filter)) then
       call this%filter%free()
       deallocate(this%filter)
    end if
    call this%free_base()

  end subroutine brinkman_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine brinkman_source_term_compute(this, time)
    class(brinkman_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, fu, fv, fw
    integer :: n

    n = this%fields%item_size(1)

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    call field_subcol3(fu, u, this%brinkman, n)
    call field_subcol3(fv, v, this%brinkman, n)
    call field_subcol3(fw, w, this%brinkman, n)

  end subroutine brinkman_source_term_compute

  ! ========================================================================== !
  ! Private methods

  !> Initializes the source term from a boundary mesh.
  subroutine init_boundary_mesh(this, json)
    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    ! Options
    character(len=:), allocatable :: mesh_file_name
    character(len=:), allocatable :: distance_transform
    character(len=:), allocatable :: filter_type
    character(len=:), allocatable :: mesh_transform

    ! Read the options for the boundary mesh
    type(file_t) :: mesh_file
    type(tri_mesh_t) :: boundary_mesh
    real(kind=rp) :: scalar_r
    real(kind=dp) :: scalar_d
    logical :: cache, cache_exist
    character(len=:), allocatable :: cache_filename
    type(file_t) :: cache_file
    type(fld_file_output_t) :: cache_output
    type(fld_file_data_t) :: cache_data
    type(global_interpolation_t) :: global_interp
    type(space_t) :: prev_Xh
    type(interpolator_t) :: space_interp

    ! Mesh transform options variables
    real(kind=dp), dimension(:), allocatable :: box_min, box_max
    logical :: keep_aspect_ratio
    real(kind=dp), dimension(3) :: scaling
    real(kind=dp), dimension(3) :: translation
    type(field_t) :: temp_field
    type(aabb_t) :: mesh_box, target_box
    integer :: idx_p
    character(len=LOG_SIZE) :: log_msg

    ! ------------------------------------------------------------------------ !
    ! Read the options for the boundary mesh

    call json_get(json, 'name', mesh_file_name)
    call json_get_or_default(json, 'cache', cache, .false.)

    ! Settings on how to filter the design field
    call json_get(json, 'distance_transform.type', distance_transform)

    ! ------------------------------------------------------------------------ !
    ! Check if we can load from cache
    if (cache) then
       call json_get(json, 'cache_file', cache_filename)

       inquire(file=trim(cache_filename) // "0.nek5000", exist=cache_exist)
       write(log_msg, '(A)') "Checking for Brinkman source term cache."
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)

       if (cache_exist) then
          write(log_msg, '(A)') "Loading Brinkman source term from cache."
          call neko_log%message(log_msg, NEKO_LOG_DEBUG)

          call cache_data%init()
          call temp_field%init(this%coef%dof)

          call cache_file%init(cache_filename // "0.fld")
          call cache_file%set_counter(0)
          call cache_file%read(cache_data)

          !
          ! Check that the data in the fld file matches the current case.
          ! Note that this is a safeguard and there are corner cases where
          ! two different meshes have the same dimension and same # of elements
          ! but this should be enough to cover obvious cases.
          !
          if (cache_data%glb_nelv .ne. temp_field%msh%glb_nelv .or. &
               cache_data%gdim .ne. temp_field%msh%gdim) then
             call neko_error("The fld file must match the current mesh! " // &
                  "Use 'interpolate': 'true' to enable interpolation.")
          end if

          ! Do the space-to-space interpolation
          call prev_Xh%init(GLL, cache_data%lx, cache_data%ly, cache_data%lz)
          call space_interp%init(temp_field%Xh, prev_Xh)
          call space_interp%map_host(temp_field%x, cache_data%p%x, &
               cache_data%nelv, temp_field%Xh)
          call space_interp%free()
          call prev_Xh%free()

          ! Synchronize to device if needed
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_memcpy(temp_field%x, temp_field%x_d, &
                  temp_field%size(), HOST_TO_DEVICE, sync = .true.)
          end if

          ! Update the global indicator field by max operator
          call field_pwmax2(this%indicator, temp_field)

          ! Clean up
          call cache_data%free()
          call temp_field%free()
          call cache_file%free()
          return
       end if
    end if

    ! ------------------------------------------------------------------------ !
    ! Load the immersed boundary mesh

    call mesh_file%init(mesh_file_name)
    call mesh_file%read(boundary_mesh)

    if (boundary_mesh%nelv .eq. 0) then
       call neko_error('No elements in the boundary mesh')
    end if

    ! ------------------------------------------------------------------------ !
    ! Transform the mesh if specified.

    call json_get_or_default(json, 'mesh_transform.type', &
         mesh_transform, 'none')

    select case (mesh_transform)
    case ('none')
       ! Do nothing
    case ('bounding_box')
       call json_get(json, 'mesh_transform.box_min', box_min)
       call json_get(json, 'mesh_transform.box_max', box_max)
       call json_get_or_default(json, 'mesh_transform.keep_aspect_ratio', &
            keep_aspect_ratio, .true.)

       if (size(box_min) .ne. 3 .or. size(box_max) .ne. 3) then
          call neko_error('Case file: mesh_transform. &
          &box_min and box_max must be 3 element arrays of reals')
       end if

       call target_box%init(box_min, box_max)

       mesh_box = get_aabb(boundary_mesh)

       scaling = target_box%get_diagonal() / mesh_box%get_diagonal()
       if (keep_aspect_ratio) then
          scaling = minval(scaling)
       end if

       translation = - scaling * mesh_box%get_min() + target_box%get_min()

       do idx_p = 1, boundary_mesh%mpts
          boundary_mesh%points(idx_p)%x = &
               scaling * boundary_mesh%points(idx_p)%x + translation
       end do

       ! Report the transformation applied
       write(log_msg, '(A)') "The following transformation was applied:"
       call neko_log%message(log_msg)
       write(log_msg, '(A, 3F12.6)') "Scaling: ", scaling
       call neko_log%message(log_msg)
       write(log_msg, '(A, 3F12.6)') "Translation: ", translation
       call neko_log%message(log_msg)

    case default
       call neko_error('Unknown mesh transform')
    end select

    ! ------------------------------------------------------------------------ !
    ! Compute the permeability field

    ! Assign the signed distance field to all GLL points in the permeability
    ! field. Initially we just run a brute force loop over all GLL points and
    ! compute the signed distance function. This should be replaced with a
    ! more efficient method, such as a tree search.

    call temp_field%init(this%coef%dof)

    ! Select how to transform the distance field to a design field
    select case (distance_transform)
    case ('smooth_step')
       call json_get(json, 'distance_transform.value', scalar_d)
       scalar_r = real(scalar_d, kind=rp)

       call signed_distance_field(temp_field, boundary_mesh, scalar_d)
       call smooth_step_field(temp_field, scalar_r, 0.0_rp)

    case ('step')

       call json_get(json, 'distance_transform.value', scalar_d)
       scalar_r = real(scalar_d, kind=rp)

       call signed_distance_field(temp_field, boundary_mesh, scalar_d)
       call step_function_field(temp_field, scalar_r, 1.0_rp, 0.0_rp)

    case default
       call neko_error('Unknown distance transform')
    end select

    ! Write the field to cache
    if (cache) then
       write(log_msg, '(A)') "Writing Brinkman source term to cache."
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
       call cache_output%init(dp, cache_filename, 1)
       call cache_output%fields%assign_to_field(1, temp_field)
       call cache_output%sample(0.0_rp)
    end if

    ! Update the global indicator field by max operator
    call field_pwmax2(this%indicator, temp_field)

    call temp_field%free()

  end subroutine init_boundary_mesh

  !> Initializes the source term from a point zone.
  subroutine init_point_zone(this, json)
    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    ! Options
    character(len=:), allocatable :: zone_name

    type(field_t) :: temp_field
    class(point_zone_t), pointer :: zone
    integer :: i

    ! ------------------------------------------------------------------------ !
    ! Read the options for the point zone

    call json_get(json, 'name', zone_name)

    ! Compute the indicator field
    call temp_field%init(this%coef%dof)

    zone => neko_point_zone_registry%get_point_zone(zone_name)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill_mask(temp_field%x_d, 1.0_rp, temp_field%size(), &
            zone%mask%get_d(), zone%size)
    else
       call cfill_mask(temp_field%x, 1.0_rp, temp_field%size(), &
            zone%mask%get(), zone%size)
    end if

    ! Update the global indicator field by max operator
    call field_pwmax2(this%indicator, temp_field)

  end subroutine init_point_zone

end module brinkman_source_term
