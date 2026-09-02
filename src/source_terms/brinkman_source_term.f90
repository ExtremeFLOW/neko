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
!> Implements the `brinkman_source_term_t` type.
module brinkman_source_term
  use aabb, only : aabb_t, get_aabb
  use coefs, only : coef_t
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use field, only : field_t
  use field_list, only : field_list_t
  use math, only : cfill_mask, copy
  use device_math, only : device_cfill_mask
  use field_math, only : field_pwmax2, field_subcol3, field_copy, field_rzero
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use mappings, only : smooth_step_field, step_function_field, &
       permeability_field
  use file, only : file_t
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item, &
       json_get_or_lookup, json_get_or_lookup_or_default
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_DEBUG
  use tri_mesh, only : tri_mesh_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp, dp, sp
  use point_zone, only : point_zone_t
  use point_zone_registry, only : neko_point_zone_registry
  use profiler, only : profiler_start_region, profiler_end_region
  use signed_distance, only : signed_distance_field
  use source_term, only : source_term_t
  use utils, only : neko_error, filename_suffix
  use filter, only : filter_t
  use PDE_filter, only : PDE_filter_t
  use field_output, only : field_output_t
  use fld_file_output, only : fld_file_output_t
  use fld_file_data, only : fld_file_data_t
  use time_state, only : time_state_t

  use global_interpolation, only: global_interpolation_t
  use interpolation, only: interpolator_t
  use space, only: space_t, GLL
  use dofmap, only : dofmap_t
  use amr_restart_component, only : amr_restart_component_t
  use amr_reconstruct, only : amr_reconstruct_t
  implicit none
  private

  !> Abstract type for Brinkman object
  type, abstract, extends(amr_restart_component_t) :: brinkman_object_t

   contains
     !> Free base
     procedure, pass(this) :: free_base => brinkman_object_free_base
     !> Object initialisation
     procedure(brinkman_object_init), pass(this), deferred :: init
     !> Object destructor.
     procedure(brinkman_object_free), pass(this), deferred :: free
     !> Getting indicator
     procedure(brinkman_object_get), pass(this), deferred :: get
  end type brinkman_object_t

  !> Brinkman object wrapper
  type :: brinkman_object_entry_t
     class(brinkman_object_t), allocatable :: obj
  end type brinkman_object_entry_t

  !> Boundary mesh object
  type, extends(brinkman_object_t) :: mesh_object_t
     !> Boundary mesh
     type(tri_mesh_t) :: boundary_mesh
     !> Distance filter function
     character(len=:), allocatable :: distance_transform
     !> Transformation distance
     real(kind=dp) :: scalar_d
   contains
     !> Object initialisation
     procedure, public, pass(this) :: init => mesh_object_init
     !> Getting indicator
     procedure, public, pass(this) :: free => mesh_object_free
     !> Getting indicator
     procedure, public, pass(this) :: get => mesh_object_get
     !> AMR restart
     procedure, pass(this) :: amr_restart => mesh_object_amr_restart
  end type mesh_object_t

  !> Point zone object
  type, extends(brinkman_object_t) :: point_object_t
     !> Point zone pointer
     class(point_zone_t), pointer :: zone => null()
   contains
     !> Object initialisation
     procedure, public, pass(this) :: init => point_object_init
     !> Getting indicator
     procedure, public, pass(this) :: free => point_object_free
     !> Getting indicator
     procedure, public, pass(this) :: get => point_object_get
     !> AMR restart
     procedure, pass(this) :: amr_restart => point_object_amr_restart
  end type point_object_t

  !> Field object
  type, extends(brinkman_object_t) :: field_object_t
     !> Source field for indicator
     type(field_t), pointer :: source_field => null()
   contains
     !> Object initialisation
     procedure, public, pass(this) :: init => field_object_init
     !> Getting indicator
     procedure, public, pass(this) :: free => field_object_free
     !> Getting indicator
     procedure, public, pass(this) :: get => field_object_get
     !> AMR restart
     procedure, pass(this) :: amr_restart => field_object_amr_restart
  end type field_object_t

  !> File object
  type, extends(brinkman_object_t) :: file_object_t
     !> Source field for indicator
     type(field_t), pointer :: source_field => null()
   contains
     !> Object initialisation
     procedure, public, pass(this) :: init => file_object_init
     !> Getting indicator
     procedure, public, pass(this) :: free => file_object_free
     !> Getting indicator
     procedure, public, pass(this) :: get => file_object_get
     !> AMR restart
     procedure, pass(this) :: amr_restart => file_object_amr_restart
  end type file_object_t

  !> A Brinkman source term.
  !! The region and strength are controlled by assigning regions types and
  !! brinkman limits to the source term.
  type, public, extends(source_term_t) :: brinkman_source_term_t
     private
     !> The unfiltered indicator field
     type(field_t), pointer :: unfiltered => null()
     !> The value of the source term.
     type(field_t), pointer :: indicator => null()
     !> Brinkman permeability field.
     type(field_t), pointer :: brinkman => null()
     !> Filter
     class(filter_t), allocatable :: filter

     ! Variables needed for forcing restart
     !> Number of objects
     integer :: n_objects
     !> Objects used to generate indicator
     type(brinkman_object_entry_t), dimension(:), allocatable :: objects
     !> Limits for permeability field
     real(kind=rp), dimension(2) :: brinkman_limits
     !> Penalty for permeability field
     real(kind=rp) :: brinkman_penalty

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init => &
          brinkman_source_term_init_from_json
     !> Destructor.
     procedure, public, pass(this) :: free => brinkman_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: compute_ => brinkman_source_term_compute
     !> AMR restart
     procedure, pass(this) :: amr_restart => brinkman_source_term_amr_restart

     ! ----------------------------------------------------------------------- !
     ! Internal subroutines for adding a given object type to the Brinkman
     ! source term

!     procedure, pass(this) :: add_boundary_mesh
!     procedure, pass(this) :: add_point_zone
!     procedure, pass(this) :: add_file
!     procedure, pass(this) :: add_field

  end type brinkman_source_term_t

  abstract interface
     !> Object initialisation
     subroutine brinkman_object_init(this, json, dof)
       import brinkman_object_t, json_file, dofmap_t
       class(brinkman_object_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       type(dofmap_t), intent(in) :: dof
     end subroutine brinkman_object_init

     !> Object destructor.
     subroutine brinkman_object_free(this)
       import brinkman_object_t
       class(brinkman_object_t), intent(inout) :: this
     end subroutine brinkman_object_free

     !> Getting indicator
     subroutine brinkman_object_get(this, field)
       import brinkman_object_t, field_t
       class(brinkman_object_t), intent(inout) :: this
       type(field_t), intent(inout) :: field
     end subroutine brinkman_object_get
  end interface

contains

  !> Free base Brinkman object
  subroutine brinkman_object_free_base(this)
    class(brinkman_object_t), intent(inout) :: this
    call this%free_amr_base()
  end subroutine brinkman_object_free_base

  !> Mesh object initialisation
  !! @param[inout]  json    runtime parameters
  !! @param[in]     dof     external dofmap for the field
  subroutine mesh_object_init(this, json, dof)
    class(mesh_object_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(in) :: dof
    ! Options
    character(len=:), allocatable :: mesh_file_name, mesh_transform
    ! Read the options for the boundary mesh
    type(file_t) :: mesh_file
    ! Mesh transform options variables
    real(kind=dp), dimension(:), allocatable :: box_min, box_max
    logical :: keep_aspect_ratio
    real(kind=dp), dimension(3) :: scaling
    real(kind=dp), dimension(3) :: translation
    type(aabb_t) :: mesh_box, target_box
    integer :: idx_p
    character(len=LOG_SIZE) :: log_msg

    call this%free()

    ! ------------------------------------------------------------------------ !
    ! Read the options for the boundary mesh

    call json_get(json, 'name', mesh_file_name)

    ! ------------------------------------------------------------------------ !
    ! Load the immersed boundary mesh

    call mesh_file%init(mesh_file_name)
    call mesh_file%read(this%boundary_mesh)

    if (this%boundary_mesh%nelv .eq. 0) then
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
       call json_get_or_lookup(json, 'mesh_transform.box_min', box_min)
       call json_get_or_lookup(json, 'mesh_transform.box_max', box_max)
       call json_get_or_default(json, 'mesh_transform.keep_aspect_ratio', &
            keep_aspect_ratio, .true.)

       if (size(box_min) .ne. 3 .or. size(box_max) .ne. 3) then
          call neko_error('Case file: mesh_transform. &
          &box_min and box_max must be 3 element arrays of reals')
       end if

       call target_box%init(box_min, box_max)

       mesh_box = get_aabb(this%boundary_mesh)

       scaling = target_box%get_diagonal() / mesh_box%get_diagonal()
       if (keep_aspect_ratio) then
          scaling = minval(scaling)
       end if

       translation = - scaling * mesh_box%get_min() + target_box%get_min()

       do idx_p = 1, this%boundary_mesh%gpts
          this%boundary_mesh%points(idx_p)%x = &
               scaling * this%boundary_mesh%points(idx_p)%x + translation
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

    ! Settings on how to filter the design field
    call json_get(json, 'distance_transform.type', this%distance_transform)

    select case (this%distance_transform)
    case ('smooth_step')
       call json_get_or_lookup(json, 'distance_transform.value', this%scalar_d)
    case ('step')
       call json_get_or_lookup(json, 'distance_transform.value', this%scalar_d)
    case default
       call neko_error('Unknown distance transform')
    end select

  end subroutine mesh_object_init

  !> Mesh object destructor.
  subroutine mesh_object_free(this)
    class(mesh_object_t), intent(inout) :: this

    call this%boundary_mesh%free()
    if (allocated(this%distance_transform)) deallocate(this%distance_transform)
    this%scalar_d = 0.0_dp

  end subroutine mesh_object_free

  !> Getting indicator
  !! @param[inout]  field    indicator field
  subroutine mesh_object_get(this, field)
    class(mesh_object_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    real(kind=rp) :: scalar_r
    type(field_t), pointer :: temp_field
    integer :: temp_idx

    ! ------------------------------------------------------------------------ !
    ! Compute the permeability field

    ! Assign the signed distance field to all GLL points in the permeability
    ! field. Initially we just run a brute force loop over all GLL points and
    ! compute the signed distance function. This should be replaced with a
    ! more efficient method, such as a tree search.

    call neko_scratch_registry%request_field(temp_field, temp_idx, .true.)

    ! Select how to transform the distance field to a design field
    select case (this%distance_transform)
    case ('smooth_step')
       scalar_r = real(this%scalar_d, kind=rp)
       call signed_distance_field(temp_field, this%boundary_mesh, this%scalar_d)
       call smooth_step_field(temp_field, scalar_r, 0.0_rp)

    case ('step')
       scalar_r = real(this%scalar_d, kind=rp)
       call signed_distance_field(temp_field, this%boundary_mesh, this%scalar_d)
       call step_function_field(temp_field, scalar_r, 1.0_rp, 0.0_rp)

    case default
       call neko_error('Unknown distance transform')
    end select

    ! Update the global indicator field by max operator
    call field_pwmax2(field, temp_field)

    call neko_scratch_registry%relinquish(temp_idx)

  end subroutine mesh_object_get

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine mesh_object_amr_restart(this, reconstruct, counter, time)
    class(mesh_object_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time
  end subroutine mesh_object_amr_restart

  !> Point object initialisation
  !! @param[inout]  json    runtime parameters
  !! @param[in]     dof     external dofmap for the field
  subroutine point_object_init(this, json, dof)
    class(point_object_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(in) :: dof
    character(len=:), allocatable :: zone_name

    call this%free()

    ! ------------------------------------------------------------------------ !
    ! Read the options for the point zone

    call json_get(json, 'name', zone_name)
    this%zone => neko_point_zone_registry%get_point_zone(zone_name)
    deallocate(zone_name)

  end subroutine point_object_init

  !> Point object destructor.
  subroutine point_object_free(this)
    class(point_object_t), intent(inout) :: this

    if (associated(this%zone)) nullify(this%zone)

  end subroutine point_object_free

  !> Getting indicator
  !! @param[inout]  field    indicator field
  subroutine point_object_get(this, field)
    class(point_object_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    type(field_t), pointer :: temp_field
    integer :: temp_idx

    ! Compute the indicator field
    call neko_scratch_registry%request_field(temp_field, temp_idx, .true.)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill_mask(temp_field%x_d, 1.0_rp, temp_field%size(), &
            this%zone%mask%get_d(), this%zone%size)
    else
       call cfill_mask(temp_field%x, 1.0_rp, temp_field%size(), &
            this%zone%mask%get(), this%zone%size)
    end if

    ! Update the global indicator field by max operator
    call field_pwmax2(field, temp_field)

    call neko_scratch_registry%relinquish(temp_idx)

  end subroutine point_object_get

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine point_object_amr_restart(this, reconstruct, counter, time)
    class(point_object_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time

    write(*, *) "Nothing done for point zone"
  end subroutine point_object_amr_restart

  !> Field object initialisation
  !! @param[inout]  json    runtime parameters
  !! @param[in]     dof     external dofmap for the field
  subroutine field_object_init(this, json, dof)
    class(field_object_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(in) :: dof
    character(len=:), allocatable :: field_name

    call this%free()

    ! Read the field name
    call json_get(json, 'name', field_name)
    this%source_field => neko_registry%get_field(field_name)
    deallocate(field_name)

  end subroutine field_object_init

  !> Field object destructor.
  subroutine field_object_free(this)
    class(field_object_t), intent(inout) :: this

    if (associated(this%source_field)) nullify(this%source_field)

  end subroutine field_object_free

  !> Getting indicator
  !! @param[inout]  field    indicator field
  subroutine field_object_get(this, field)
    class(field_object_t), intent(inout) :: this
    type(field_t), intent(inout) :: field

    ! Update the global indicator field by max operator
    call field_pwmax2(field, this%source_field)

  end subroutine field_object_get

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine field_object_amr_restart(this, reconstruct, counter, time)
    class(field_object_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time
  end subroutine field_object_amr_restart

  !> File object initialisation
  !! @param[inout]  json    runtime parameters
  !! @param[in]     dof     external dofmap for the field
  subroutine file_object_init(this, json, dof)
    class(file_object_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(dofmap_t), intent(in) :: dof
    type(file_t) :: file
    character(len=:), allocatable :: file_name, field_name
    character(len=80) :: suffix
    integer :: file_idx

    call this%free()

    call json_get(json, 'file_name', file_name)
    call json_get_or_default(json, 'field_name', field_name, &
         'brinkman_indicator')
    call json_get_or_default(json, 'file_index', file_idx, 0)


    call neko_registry%add_field(dof, trim(field_name))
    this%source_field => neko_registry%get_field(field_name)

    call file%init(file_name)
    call file%set_counter(file_idx)

    call filename_suffix(file_name, suffix)
    select case (trim(suffix))
    case ('fld')
       block
         type(fld_file_data_t) :: fld_data
         type(field_list_t) :: fld_fields
         integer :: idx(1)

         call fld_data%init()
         call file%read(fld_data)
         select case (field_name(1:1))
         case ('p')
            call fld_data%import_fields(p = this%source_field)
         case ('u')
            call fld_data%import_fields(u = this%source_field)
         case ('v')
            call fld_data%import_fields(v = this%source_field)
         case ('w')
            call fld_data%import_fields(w = this%source_field)
         case ('t')
            call fld_data%import_fields(t = this%source_field)
         case ('s')

            if (len_trim(field_name) .eq. 3) then
               read(field_name(2:3), '(I2)') idx(1)
            else if (len_trim(field_name) .eq. 2) then
               read(field_name(2:2), '(I1)') idx(1)
            else
               call neko_error('For fields with prefix s, the field name ' // &
                    'must be in the format sXX, where XX is the index of ' // &
                    'the field in the fld file')
            end if

            call fld_fields%init(1)
            call fld_fields%assign(1, this%source_field)

            call fld_data%import_fields(s_target_list = fld_fields, &
                 s_index_list = idx)

            call fld_fields%free()
         case default
            call neko_error('Unknown field prefix in field name: ' // &
                 trim(field_name))
         end select

         call fld_data%free()
       end block

    case ('vtkhdf')

       ! VTKHDF will read the name of the field object.
       call file%read(this%source_field)

    case default
       call neko_error("Brinkman cannot read file: " // trim(file_name))
    end select

    call file%free()

  end subroutine file_object_init

  !> File object destructor.
  subroutine file_object_free(this)
    class(file_object_t), intent(inout) :: this

    if (associated(this%source_field)) nullify(this%source_field)

  end subroutine file_object_free

  !> Getting indicator
  !! @param[inout]  field    indicator field
  subroutine file_object_get(this, field)
    class(file_object_t), intent(inout) :: this
    type(field_t), intent(inout) :: field

    ! Update the global indicator field by max operator
    call field_pwmax2(field, this%source_field)

  end subroutine file_object_get

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine file_object_amr_restart(this, reconstruct, counter, time)
    class(file_object_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time
  end subroutine file_object_amr_restart

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
    character(len=LOG_SIZE) :: error_msg

    character(len=:), allocatable :: filter_type
    real(kind=rp), dimension(:), allocatable :: brinkman_limits
    real(kind=rp) :: brinkman_penalty

    character(len=:), allocatable :: object_type
    type(json_file) :: object_settings, filter_settings
    integer :: n_objects, i

    type(file_t) :: output
    type(field_list_t) :: output_fields
    logical :: output_enable
    character(len=:), allocatable :: output_path, output_format, &
         output_precision, fname, suffix
    integer :: precision

    ! ------------------------------------------------------------------------ !
    ! Read the general options for the Brinkman source term

    call neko_log%section('Brinkman source term')

    ! Read the options for the permeability field
    call json_get_or_lookup(json, 'limits', brinkman_limits)
    call json_get_or_lookup(json, 'penalty', brinkman_penalty)

    ! Mandatory fields for the general source term
    call json_get_or_lookup_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_lookup_or_default(json, "end_time", end_time, huge(0.0_rp))

    ! Output settings
    call json_get_or_default(json, 'output.enable', output_enable, .false.)
    call json_get_or_default(json, 'output.precision', output_precision, 'sp')
    call json_get_or_default(json, 'output.path', output_path, './')
    call json_get_or_default(json, 'output.format', output_format, 'fld')

    if (size(brinkman_limits) .ne. 2) then
       call neko_error('brinkman_limits must be a 2 element array of reals')
    end if

    if (output_path(len_trim(output_path):len_trim(output_path)) .ne. '/') then
       output_path = trim(output_path) // '/'
    end if

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    ! Store limits and penalty
    this%brinkman_limits(:) = brinkman_limits(:)
    this%brinkman_penalty = brinkman_penalty

    ! ------------------------------------------------------------------------ !
    ! Allocate the permeability and indicator field

    call neko_registry%add_field(coef%dof, 'brinkman_indicator', .true.)
    call neko_registry%add_field(coef%dof, 'brinkman_permeability', .true.)

    this%indicator => neko_registry%get_field('brinkman_indicator')
    this%brinkman => neko_registry%get_field('brinkman_permeability')

    ! ------------------------------------------------------------------------ !
    ! Initialise objects; select which constructor should be called

    call json%info('objects', n_children = n_objects)
    this%n_objects = n_objects
    allocate(this%objects(this%n_objects))
    do i = 1, n_objects
       call json_extract_item(json, "objects", i, object_settings)
       call json_get(object_settings, 'type', object_type)

       select case (object_type)
       case ('boundary_mesh')
          allocate(mesh_object_t :: this%objects(i)%obj)
       case ('point_zone')
          allocate(point_object_t :: this%objects(i)%obj)
       case ('field')
          allocate(field_object_t :: this%objects(i)%obj)
       case ('file')
          allocate(file_object_t :: this%objects(i)%obj)
       case default
          write(error_msg, '(A,I0,A,A)') 'Brinkman object: ', i, &
               ' unknown type: ', trim(object_type)
          call neko_error(trim(error_msg))
       end select

       call this%objects(i)%obj%init(object_settings, this%coef%dof)
    end do

    ! ------------------------------------------------------------------------ !
    ! Fill the indicator field

    call field_rzero(this%indicator)
    do i = 1, n_objects
       if (allocated(this%objects(i)%obj)) &
            call this%objects(i)%obj%get(this%indicator)
    end do

    do i = 1, n_objects
       call json_extract_item(json, "objects", i, object_settings)
       call json_get(object_settings, 'type', object_type)

       select case (object_type)
       case ('boundary_mesh')
!          call this%add_boundary_mesh(object_settings)
       case ('point_zone')
!          call this%add_point_zone(object_settings)
       case ('field')
!          call this%add_field(object_settings)
       case ('file')
!          call this%add_file(object_settings)

       case default
          write(error_msg, '(A,I0,A,A)') 'Brinkman object: ', i, &
               ' unknown type: ', trim(object_type)
          call neko_error(trim(error_msg))
       end select

    end do

    ! ------------------------------------------------------------------------ !
    ! Filter the indicator field

    call json_get_or_default(json, 'filter.type', filter_type, 'none')
    select case (filter_type)
    case ('PDE')

       ! Copy the current indicator to unfiltered (essentially a rename)
       call neko_registry%add_field(coef%dof, 'brinkman_unfiltered', .true.)
       this%unfiltered => neko_registry%get_field('brinkman_unfiltered')
       call field_copy(this%unfiltered, this%indicator)

       ! Allocate a PDE filter
       allocate(PDE_filter_t::this%filter)
       call json_get(json, 'filter', filter_settings)
       call this%filter%init(filter_settings, coef)
       call this%filter%apply(this%indicator, this%unfiltered)

    case ('none')
       ! Do nothing
    case default
       call neko_error('Brinkman source term unknown filter type')
    end select

    ! ------------------------------------------------------------------------ !
    ! Compute the permeability field

    this%brinkman = this%indicator
    call permeability_field(this%brinkman, &
         brinkman_limits(1), brinkman_limits(2), brinkman_penalty)

    ! ------------------------------------------------------------------------ !
    ! Set up output the brinkman fields if enabled

    if (output_enable) then
       fname = trim(output_path) // 'brinkman'
       select case (trim(output_format))
       case ('nek5000')
          suffix = '.fld'
       case ('adios2')
          suffix = '.bp'
       case ('hdf5')
          suffix = '.h5'
       case default
          suffix = '.' // trim(output_format)
       end select
       fname = trim(fname) // trim(suffix)

       select case (trim(output_precision))
       case ('sp', 'single')
          precision = sp
       case ('dp', 'double')
          precision = dp
       case default
          call neko_error('Unknown output precision')
       end select

       call output%init(fname, precision = precision)

       call output_fields%init(2)
       call output_fields%assign_to_field(1, this%indicator)
       call output_fields%assign_to_field(2, this%brinkman)

       call neko_log%message('Brinkman output')
       call neko_log%message('  1: Indicator')
       call neko_log%message('  2: Permeability')
       if (associated(this%unfiltered)) then
          call output_fields%append(this%unfiltered)
          call neko_log%message('  3: Unfiltered Indicator')
       end if

       call output%write(output_fields)

       call output%free()
       call output_fields%free()
    end if

    call neko_log%end_section()
  end subroutine brinkman_source_term_init_from_json

  !> Destructor.
  subroutine brinkman_source_term_free(this)
    class(brinkman_source_term_t), intent(inout) :: this
    integer :: il

    if (allocated(this%objects)) then
       do il = 1, this%n_objects
          if (allocated(this%objects(il)%obj)) then
             call this%objects(il)%obj%free()
             deallocate(this%objects(il)%obj)
          end if
       end do
       deallocate(this%objects)
    end if
    this%n_objects = 0
    this%brinkman_limits(:) = 0.0_rp
    this%brinkman_penalty = 0.0_rp

    if (associated(this%indicator)) nullify(this%indicator)
    if (associated(this%unfiltered)) nullify(this%unfiltered)
    if (associated(this%brinkman)) nullify(this%brinkman)

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

    u => neko_registry%get_field('u')
    v => neko_registry%get_field('v')
    w => neko_registry%get_field('w')

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
  subroutine add_boundary_mesh(this, json)
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
    type(field_t), pointer :: temp_field
    type(aabb_t) :: mesh_box, target_box
    integer :: idx_p, temp_idx
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

          call neko_scratch_registry%request_field(temp_field, temp_idx, .true.)

          call cache_data%init()
          call cache_file%init(cache_filename // "0.fld")
          call cache_file%set_counter(0)
          call cache_file%read(cache_data)
          call cache_data%import_fields(p = temp_field)

          ! Update the global indicator field by max operator
          call field_pwmax2(this%indicator, temp_field)

          ! Clean up
          call neko_scratch_registry%relinquish(temp_idx)
          call cache_data%free()
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
       call json_get_or_lookup(json, 'mesh_transform.box_min', box_min)
       call json_get_or_lookup(json, 'mesh_transform.box_max', box_max)
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

       do idx_p = 1, boundary_mesh%gpts
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

    call neko_scratch_registry%request_field(temp_field, temp_idx, .true.)

    ! Select how to transform the distance field to a design field
    select case (distance_transform)
    case ('smooth_step')
       call json_get_or_lookup(json, 'distance_transform.value', scalar_d)
       scalar_r = real(scalar_d, kind=rp)

       call signed_distance_field(temp_field, boundary_mesh, scalar_d)
       call smooth_step_field(temp_field, scalar_r, 0.0_rp)

    case ('step')

       call json_get_or_lookup(json, 'distance_transform.value', scalar_d)
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

    call neko_scratch_registry%relinquish(temp_idx)

  end subroutine add_boundary_mesh

  !> Initializes the source term from a point zone.
  subroutine add_point_zone(this, json)
    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    ! Options
    character(len=:), allocatable :: zone_name

    type(field_t), pointer :: temp_field
    class(point_zone_t), pointer :: zone
    integer :: i, temp_idx

    ! ------------------------------------------------------------------------ !
    ! Read the options for the point zone

    call json_get(json, 'name', zone_name)

    ! Compute the indicator field
    call neko_scratch_registry%request_field(temp_field, temp_idx, .true.)

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

    call neko_scratch_registry%relinquish(temp_idx)
  end subroutine add_point_zone

  !> Initializes the source term from a file.
  subroutine add_file(this, json)
    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(file_t) :: file
    type(field_t), pointer :: temp_field
    character(len=:), allocatable :: file_name, field_name, tmp_str
    character(len=80) :: suffix
    integer :: file_idx, temp_idx

    call json_get(json, 'file_name', file_name)
    call json_get_or_default(json, 'field_name', field_name, &
         'brinkman_indicator')
    call json_get_or_default(json, 'file_index', file_idx, 0)

    call neko_scratch_registry%request_field(temp_field, temp_idx, .true.)

    call file%init(file_name)
    call file%set_counter(file_idx)

    call filename_suffix(file_name, suffix)
    select case (trim(suffix))
    case ('fld')
       block
         type(fld_file_data_t) :: fld_data
         type(field_list_t) :: fld_fields
         integer :: idx(1)

         call fld_data%init()
         call file%read(fld_data)
         select case (field_name(1:1))
         case ('p')
            call fld_data%import_fields(p = temp_field)
         case ('u')
            call fld_data%import_fields(u = temp_field)
         case ('v')
            call fld_data%import_fields(v = temp_field)
         case ('w')
            call fld_data%import_fields(w = temp_field)
         case ('t')
            call fld_data%import_fields(t = temp_field)
         case ('s')

            if (len_trim(field_name) .eq. 3) then
               read(field_name(2:3), '(I2)') idx(1)
            else if (len_trim(field_name) .eq. 2) then
               read(field_name(2:2), '(I1)') idx(1)
            else
               call neko_error('For fields with prefix s, the field name ' // &
                    'must be in the format sXX, where XX is the index of ' // &
                    'the field in the fld file')
            end if

            call fld_fields%init(1)
            call fld_fields%assign(1, temp_field)

            call fld_data%import_fields(s_target_list = fld_fields, &
                 s_index_list = idx)

            call fld_fields%free()
         case default
            call neko_error('Unknown field prefix in field name: ' // &
                 trim(field_name))
         end select

         call fld_data%free()
       end block

    case ('vtkhdf')

       ! VTKHDF will read the name of the field object.
       tmp_str = trim(temp_field%name)
       temp_field%name = trim(field_name)

       call file%read(temp_field)

       temp_field%name = trim(tmp_str)

    case default
       call neko_error("Brinkman cannot read file: " // trim(file_name))
    end select

    ! Update the global indicator field by max operator
    call field_pwmax2(this%indicator, temp_field)

    call neko_scratch_registry%relinquish(temp_idx)
    call file%free()

  end subroutine add_file

  !> Initializes the source term from a field.
  subroutine add_field(this, json)
    class(brinkman_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: field_name
    type(field_t), pointer :: temp_field

    ! Read the field name
    call json_get(json, 'name', field_name)
    temp_field => neko_registry%get_field(field_name)

    ! Update the global indicator field by max operator
    call field_pwmax2(this%indicator, temp_field)

    if (associated(temp_field)) nullify(temp_field)

  end subroutine add_field

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     time          time state
  subroutine brinkman_source_term_amr_restart(this, reconstruct, counter, time)
    class(brinkman_source_term_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter
    type(time_state_t), intent(in) :: time
    integer :: il

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    call this%amr_restart_base(reconstruct, counter, time)

    ! Reallocate Brinkman filtered and unfiltered indicator, permeability
    if (reconstruct%nold .ne. reconstruct%nnew) then
       if (associated(this%indicator)) &
            call this%indicator%amr_reallocate(reconstruct, counter, time)
       if (associated(this%unfiltered)) &
            call this%unfiltered%amr_reallocate(reconstruct, counter, time)
       if (associated(this%brinkman)) &
            call this%brinkman%amr_reallocate(reconstruct, counter, time)
    end if

    ! Regenerate indicator
    call field_rzero(this%indicator)
    do il = 1, this%n_objects
       if (allocated(this%objects(il)%obj)) then
          call this%objects(il)%obj%amr_restart(reconstruct, counter, time)
!          call this%objects(il)%obj%get(this%indicator)
       end if
    end do
    
    ! testing
    block
      integer :: lx, il, jl, kl, el
      real(kind=rp) :: x, y, z
      lx = this%coef%dof%Xh%lx
      do el = 1, reconstruct%nnew
         do kl = 1, lx
            do jl = 1, lx
               do il = 1, lx
                  x = this%coef%dof%x(il, jl, kl, el)
                  y = this%coef%dof%y(il, jl, kl, el)
                  z = this%coef%dof%z(il, jl, kl, el)
                  if (sqrt((x - 1.0_rp)**2 + (y - 0.5_rp)**2 + &
                       (z - 0.5_rp)**2) .lt. 0.1) &
                       this%indicator%x(il, jl, kl, el) = 1.0_rp
               end do
            end do
         end do
      end do
    end block
    

    ! Perform filtering
    if (allocated(this%filter)) then
       select type (filter => this%filter)
       type is (PDE_filter_t)
          if (.not. associated(this%unfiltered)) &
               call neko_error('Brinkman source term missing unfiltered field')
          call field_copy(this%unfiltered, this%indicator)

          ! Restart PDE filter
          call this%filter%amr_restart(reconstruct, counter, time)
          ! apply filter
          call this%filter%apply(this%indicator, this%unfiltered)
       class default
          call neko_error('Brinkman source term unknown filter type')
       end select
    end if

    ! ------------------------------------------------------------------------ !
    ! Compute the permeability field

    this%brinkman = this%indicator
    call permeability_field(this%brinkman, this%brinkman_limits(1), &
         this%brinkman_limits(2), this%brinkman_penalty)

    
    block
      type(file_t) :: output
      type(field_list_t) :: output_fields
      character(len=:), allocatable :: fname
      integer :: precision
      fname = './brinkman_ref.fld'
      precision = sp
      call output%init(fname, precision = precision)
      call output_fields%init(2)
      call output_fields%assign_to_field(1, this%indicator)
      call output_fields%assign_to_field(2, this%brinkman)
      call output%write(output_fields)
      call output%free()
      call output_fields%free()

      write(*,*) "TESTfilterPDE", allocated(this%filter), &
           associated(this%unfiltered)
      call neko_log%message('Working on brinkman')
    end block
    

  end subroutine brinkman_source_term_amr_restart

end module brinkman_source_term
