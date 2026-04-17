! Copyright (c) 2022, The Neko Authors
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
!> Implements type field_subsampler_t.
module field_subsampler
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use time_state, only : time_state_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use logger, only : neko_log, NEKO_LOG_DEBUG
  use device, only : DEVICE_TO_HOST, HOST_TO_DEVICE, device_memcpy
  use math, only : masked_gather_copy
  use device_math, only : device_masked_gather_copy_aligned
  use comm, only : pe_rank, NEKO_COMM
  use mpi_f08, only : MPI_ALLREDUCE, MPI_IN_PLACE, MPI_SUM, MPI_INTEGER
  use utils, only : neko_warning, neko_error
  use simulation_component, only : simulation_component_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use case, only : case_t
  use point_zone, only : point_zone_t
  use time_based_controller, only : time_based_controller_t
  use point_zone_registry, only : neko_point_zone_registry
  use registry, only : neko_registry, registry_t
  use dofmap, only : dofmap_t
  use space, only : space_t, GLL
  use mesh, only : mesh_t
  use field_list, only : field_list_t
  use interpolation, only : interpolator_t
  use scratch_registry, only : neko_scratch_registry
  use field_writer, only : field_writer_t
  use math, only : glsum

  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Implements the field_subsampler_t simulation components, which allows
  !! for masking regions of the domain with a point_zone and to re-sample
  !! fields with a different polynomial order.
  type, public, extends(simulation_component_t) :: field_subsampler_t

     !> Abstract interface for the compute() subroutine, which will be assigned
     !! at runtime depending on the subsampling method used.
     procedure(compute_intrf), pass(this), pointer :: compute_impl => &
          dummy_compute

     !> Fields to subsample.
     character(len=20), allocatable :: field_names(:)
     !> Number of fields to subsample (size of field_names).
     integer :: n_fields = 0
     !> Pointers to the subsampled fields in the registry.
     type(field_list_t) :: fields
     !> Pointers to the source fields in the registry.
     type(field_list_t) :: source_fields
     !> Internal field writer object for output
     type(field_writer_t) :: writer
     !> The dofmap to use for the newly created masked/re-sampled fields.
     type(dofmap_t) :: dof

     ! =====================================================================
     ! Variables for point zone masking.

     !> Point zone to use for subsampling.
     class(point_zone_t), pointer :: point_zone => null()
     !> Points to the new mesh for element masking. If no point zone, is
     !! used, this%msh will point to the same mesh as the current sim.
     ! Otherwise, it will point to its own internal masked mesh_t.
     type(mesh_t), pointer :: msh => null()
     !> Indicates if the instance of msh is internally managed
     logical :: internal_mesh = .false.

     ! =====================================================================
     ! Variables for space-to-space interpolation

     !> New polynomial order to use for subsampling.
     integer :: lx = -1
     !> Space-to-space interpolator object.
     type(interpolator_t) :: interpolator
     !> Points to the new space for resampling. If no point zone is provided,
     !! will point to the same space as the current sim. Otherwise, will
     !! point to its own masked one.
     type(space_t), pointer :: Xh => null()
     !> Indicates if the instance of msh is internally managed. Also serves
     !! as a flag for determining whether or not to perform space-to-space
     !! interpolation
     logical :: internal_space = .false.

     !> Flag to indicate whether this%check() has been called and no errors
     !! were found.
     logical :: checked = .false.

   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => field_subsampler_init_json
     !> Common constructor.
     procedure, pass(this) :: init_common => &
          field_subsampler_init_common
     !> Destructor.
     procedure, pass(this) :: free => field_subsampler_free
     !> Compute the indicator (called according to the simcomp controller).
     procedure, pass(this) :: compute_ => compute_wrapper
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          field_subsampler_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          field_subsampler_init_from_controllers_properties
     !> Check valid initialization and assignment of variables
     procedure, pass(this) :: check => field_subsampler_check
     !> Dummy compute routine
     procedure, pass(this) :: dummy_compute
  end type field_subsampler_t

  !> Abstract interface for the compute() subroutine, which will be assigned
  !! at runtime depending on the subsampling method used.
  abstract interface
     subroutine compute_intrf(this, time)
       import :: field_subsampler_t, time_state_t
       class(field_subsampler_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine compute_intrf
  end interface

contains

  !> Checks the validity of the setup before calling compute()
  subroutine field_subsampler_check(this)
    class(field_subsampler_t), intent(inout) :: this

    integer :: n, ierr

    logical :: is_valid
    is_valid = .false.

    ! Internal mesh means we are doing point zone masking
    if (this%internal_mesh) then

       ! Check if point zone exists and mask is set up
       is_valid = (associated(this%point_zone) .and. &
            this%point_zone%mask%is_set())
       if (.not. is_valid) call neko_error("The point zone is missing " // &
            " or has not been set up properly")

       n = this%point_zone%size
       call MPI_Allreduce(MPI_IN_PLACE, n, 1, MPI_INTEGER, MPI_SUM, &
            NEKO_COMM, ierr)

       is_valid = n .gt. 0
       if (.not. is_valid) call neko_warning("Point zone is empty")

    end if

    ! Internal space means we are doing space-to-space interpolation
    if (this%internal_space) then
       is_valid = this%lx .gt. 0
       if (.not. is_valid) call neko_error("lx has not been set up properly")

       is_valid = allocated(this%interpolator%Xh_to_Yh)
       if (.not. is_valid) call neko_error("The interpolator has not been " // &
            "initialized properly")
    end if

    ! Check the dofmap
    is_valid = (this%dof%global_size() .gt. 0 .and. allocated(this%dof%x))
    if (.not. is_valid) call neko_error("Dofmap not initialized or empty")

    ! Check the internal field list, pointing to the fields in the registry
    ! Technically the size of the fields should be checked but since they are
    ! all based on this%dofmap we assume if we get here we are fine.
    is_valid = (this%fields%size() .gt. 0 .and. allocated(this%fields%items))
    if (.not. is_valid) call neko_error("Internal field_list not initialized")

    ! Lastly, check that we have properly associated compute_impl
    is_valid = associated(this%compute_impl)
    if (.not. is_valid) call neko_error("compute_impl not associated")

    this%checked = .true.

  end subroutine field_subsampler_check

  !> Wrapper for the run-time-assigned subroutine
  !! compute_impl, which will be assigned at runtime to the correct compute().
  subroutine compute_wrapper(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (.not. this%checked) call this%check()
    call this%compute_impl(time)

  end subroutine compute_wrapper

  !> Dummy subroutine assigned by default to compute_impl.
  subroutine dummy_compute(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call neko_error("field_subsampler must be initialized first!")

  end subroutine dummy_compute

  !> Constructor.
  subroutine field_subsampler_init_json(this, json, case)
    class(field_subsampler_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    character(len=:), allocatable :: name
    character(len=20), allocatable :: which_fields(:)
    integer :: lx

    logical :: do_space_interp, do_point_zone_masking

    character(len=:), allocatable :: pz_name
    class(point_zone_t), pointer :: point_zone
    point_zone => null()

    call this%init_base(json, case)

    call json_get_or_default(json, "name", name, "field_subsampler")
    call json_get(json, "source_fields", which_fields)

    ! ========================================================================
    ! Check if the user has provided a polynomial order for space-to-space
    ! interpolation and initialize accordingly
    do_space_interp = json%valid_path('polynomial_order')
    if (do_space_interp) then
       call json_get_or_default(json, "polynomial_order", lx, &
            case%fluid%Xh%lx - 1)
       lx = lx + 1
    end if

    ! ========================================================================
    ! Check if the user has provided a point zone for masking and initialize
    ! accordingly
    do_point_zone_masking = json%valid_path('point_zone')
    if (do_point_zone_masking) then
       call json_get(json, "point_zone", pz_name)
       point_zone => neko_point_zone_registry%get_point_zone(trim(pz_name))
    end if

    ! ========================================================================
    ! Call the common constructor with the appropriate optional arguments
    if (do_space_interp .and. do_point_zone_masking) then
       call field_subsampler_init_common(this, name, which_fields, &
            lx = lx, point_zone = point_zone)
    else if (do_space_interp .and. .not. do_point_zone_masking) then
       call field_subsampler_init_common(this, name, which_fields, &
            lx = lx)
    else if (.not. do_space_interp .and. do_point_zone_masking) then
       call field_subsampler_init_common(this, name, which_fields, &
            point_zone = point_zone)
    else
       call neko_error("Invalid configuration: please pass either " // &
            "a point zone or a valid polynomial order.")
    end if

    ! =======================================================================
    ! Initialize the field writer based on the new subsampled fields!
    ! We cannot use the JSON because "fields" is already taken by the field

    block
      character(len=1024) :: new_field_names(this%n_fields)
      integer :: i

      do i = 1, this%n_fields
         new_field_names(i) = trim(this%fields%name(i))
      end do

      ! Will be picked up by the field_writer.
      call json%add("fields", new_field_names)

      ! This is needed so the field_writer doesn't pick up the point zone!
      call json%remove("point_zone")

      ! Force a different file name, as we cannot add these subsampled
      ! fields to the regular fluid output (different # of elements and
      ! polynomial order not supported)
      if (.not. json%valid_path("output_filename")) then
         call json%add("output_filename", this%name)
      end if

      call this%writer%init(json, case)

      ! Put the point zone back to have the JSON back as it was :)
      call json%add("point_zone", trim(this%point_zone%name))

    end block

  end subroutine field_subsampler_init_json

  !> Actual constructor.
  !! @param name The unique name of the simcomp.
  !! @param which_fields The names of the fields to be subsampled.
  !! @param lx The number of GLL points to use for subsampling.
  !! @param point_zone The point zone to use for subsampling.
  subroutine field_subsampler_init_common(this, name, which_fields, &
       lx, point_zone)
    class(field_subsampler_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=20), intent(in) :: which_fields(:)
    integer, intent(in), optional :: lx
    class(point_zone_t), pointer, intent(in), optional :: point_zone

    this%name = name
    this%field_names = which_fields

    this%n_fields = size(which_fields)

    ! ========================================================================
    ! Polynomial order / space interpolation

    this%internal_space = .false.
    if (present(lx)) then

       this%lx = lx

       if (this%lx .eq. this%case%fluid%Xh%lx) then
          if (pe_rank .eq. 0) then
             call neko_warning("No change in polynomial order")
             call neko_warning("Space-to-space interpolation disabled.")
          end if
       else

          ! Create the subsampled space, with the (strong) assumption that
          ! lx = ly = lz.
          allocate(this%Xh)
          call this%Xh%init(GLL, this%lx, this%lx, this%lx)
          this%internal_space = .true.

          call this%interpolator%init(this%Xh, this%case%fluid%Xh)
       end if

    end if

    ! If we haven't initialized our own space above, retrieve it from the
    ! fluid.
    if (.not. this%internal_space) then
       this%lx = this%case%fluid%Xh%lx
       this%Xh => this%case%fluid%Xh
    end if

    ! ========================================================================
    ! Point zone masking

    if (present(point_zone)) then
       this%point_zone => point_zone

       if (.not. point_zone%full_elements) call neko_error("full_elements" // &
            "must be enabled when sampling a point zone")

       ! Create the subsampled mesh
       allocate(this%msh)
       call this%case%fluid%msh%subset_by_mask(this%msh, &
            point_zone%mask, &
            this%case%fluid%Xh%lx, & ! Yes, lx, ly and lz should come from the
            this%case%fluid%Xh%ly, & ! fluid object, since the point_zone
            this%case%fluid%Xh%lz) ! mask is built from those values.
       this%internal_mesh = .true.

    else
       this%msh => this%case%fluid%msh
       this%internal_mesh = .false.
    end if

    ! Create the subsampled dofmap
    call this%dof%init(this%msh, this%Xh)

    ! =======================================================================
    ! Register the new fields in the registry and initialize the field_lists

    block
      integer :: i
      character(len=2048) :: field_name

      call this%fields%init(this%n_fields)
      call this%source_fields%init(this%n_fields)

      do i = 1, this%n_fields

         ! Point the new, subsampled fields to the registry
         field_name = this%name // "_" // this%field_names(i)
         call neko_registry%add_field(this%dof, field_name)
         this%fields%items(i)%ptr => neko_registry%get_field(trim(field_name))

         ! Point the source fields from the registry to avoid searches
         this%source_fields%items(i)%ptr => &
              neko_registry%get_field(this%field_names(i))
      end do
    end block

    !========================================================================
    ! Assign the correct compute() depending on whether interpolation is
    ! enabled or not and if a point zone is provided or not.

    if (this%internal_space .and. .not. this%internal_mesh) then
       this%compute_impl => field_subsampler_compute_Xh
    else if (.not. this%internal_space .and. this%internal_mesh) then
       this%compute_impl => field_subsampler_compute_pz
    else if (this%internal_space .and. this%internal_mesh) then
       this%compute_impl => field_subsampler_compute_pz_Xh
    else
       call neko_error("Please pass either a point zone or a " // &
            "valid polynomial order.")
    end if

  end subroutine field_subsampler_init_common

  !> Destructor
  subroutine field_subsampler_free(this)
    class(field_subsampler_t), intent(inout) :: this

    this%lx = -1

    if (allocated(this%field_names)) deallocate(this%field_names)
    nullify(this%point_zone)

    if (this%internal_space) then
       call this%Xh%free()
       deallocate(this%Xh)
       this%internal_space = .false.
    end if

    if (this%internal_mesh) then
       call this%msh%free()
       deallocate(this%msh)
       this%internal_mesh = .false.
    end if

    call this%dof%free()
    call this%interpolator%free()

    nullify(this%Xh)
    nullify(this%msh)

    call this%source_fields%free()
    call this%fields%free()
    call this%writer%free()

    call this%free_base()

    this%compute_impl => dummy_compute

  end subroutine field_subsampler_free

  !> Subsample the fields based only on a point zone, no space-to-space
  !! interpolation.
  subroutine field_subsampler_compute_pz(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    integer :: i, n, n_mask

    n = this%dof%size()
    n_mask = this%point_zone%mask%size()

    do i = 1, size(this%field_names)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_gather_copy_aligned(this%fields%x_d(i), &
               this%source_fields%x_d(i), &
               this%point_zone%mask%get_d(), n, n_mask)
       else
          call masked_gather_copy(this%fields%x(i), &
               this%source_fields%x(i), &
               this%point_zone%mask%get(), n, n_mask)
       end if

    end do

  end subroutine field_subsampler_compute_pz

  !> Subsample fields wihout any point zone, only space-to-space
  !! interpolation is enabled.
  subroutine field_subsampler_compute_Xh(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    integer :: i

    do i = 1, size(this%field_names)
       call this%interpolator%map(this%fields%x(i), &
            this%source_fields%x(i), &
            this%msh%nelv, this%Xh)
    end do

  end subroutine field_subsampler_compute_Xh

  !> Subsample the fields based a point zone, and with space-to-space
  !! interpolation enabled.
  !! @details This routine uses the work array wk as a buffer to first
  !! do the masked copy.
  subroutine field_subsampler_compute_pz_Xh(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: wk
    integer :: i, n, n_mask, tmp_index

    n = this%dof%size()
    n_mask = this%point_zone%mask%size()

    call neko_scratch_registry%request_field(wk, tmp_index, .false.)

    do i = 1, size(this%field_names)

       ! First, do a masked copy onto the submesh but that has the same poly.
       ! order
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_gather_copy_aligned(wk%x_d, &
               this%source_fields%x_d(i), &
               this%point_zone%mask%get_d(), n, n_mask)
       else
          call masked_gather_copy(wk%x, &
               this%source_fields%x(i), &
               this%point_zone%mask%get(), n, n_mask)
       end if

       ! Then, map the two different spaces
       call this%interpolator%map(this%fields%x(i), wk%x, &
            this%msh%nelv, this%Xh)

    end do

    call neko_scratch_registry%relinquish_field(tmp_index)

  end subroutine field_subsampler_compute_pz_Xh

  !> Constructor from components, passing controllers.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param which_fields List of field names to subsample.
  !! @param lx Number of GLL points to use for the space-to-space
  !! interpolation.
  !! @param point_zone Point zone to use for the subsampling.
  subroutine field_subsampler_init_from_controllers(this, name, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       which_fields, lx, point_zone)
    class(field_subsampler_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=20), intent(in) :: which_fields(:)
    integer, intent(in) :: lx
    class(point_zone_t), intent(in), pointer, optional :: point_zone

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)

    call this%init_common(name, which_fields, lx = lx, &
         point_zone = point_zone)

  end subroutine field_subsampler_init_from_controllers

  !> Constructor from components, passing properties to the
  !! time_based_controller` components in the base type.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_controller Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_controller Control mode for output.
  !! @param output_value Value parameter for output.
  !! @param which_fields List of field names to subsample.
  !! @param lx Number of GLL points to use for the space-to-space
  !! interpolation.
  !! @param point_zone Point zone to use for the subsampling.
  subroutine field_subsampler_init_from_controllers_properties(this, name, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, which_fields, lx, &
       point_zone)
    class(field_subsampler_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=20), intent(in) :: which_fields(:)
    integer, intent(in), optional :: lx
    class(point_zone_t), pointer, intent(in), optional :: point_zone

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)

    call this%init_common(name, which_fields, lx = lx, &
         point_zone = point_zone)

  end subroutine field_subsampler_init_from_controllers_properties

end module field_subsampler
