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
  use comm, only : pe_rank
  use utils, only : neko_warning, neko_error
  use simulation_component, only : simulation_component_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use case, only : case_t
  use point_zone, only : point_zone_t
  use point_zone_registry, only : neko_point_zone_registry
  use registry, only : neko_registry, registry_t
  use dofmap, only : dofmap_t
  use space, only : space_t, GLL
  use mesh, only : mesh_t
  use field_list, only : field_list_t
  use interpolation, only : interpolator_t

  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Provides tools to calculate the spectral error indicator
  !! @details
  !! This is a posteriori error measure, based on the local properties of
  !! the spectral solution, which was developed by Mavriplis. This method
  !! formally only gives an indication of the error.
  type, public, extends(simulation_component_t) :: field_subsampler_t

     procedure(compute_intrf), pass(this), pointer :: compute_impl => &
          dummy_compute

     !> Fields to subsample.
     character(len=20), allocatable :: field_names(:)
     !> Pointers to the subsampled fields in the registry
     type(field_list_t) :: fields
     !> Point zone to use for subsampling.
     class(point_zone_t), pointer :: point_zone => null()

     !> Points to the mesh our fields are based on. If no point zone, is
     !! used, this%msh will point to the same mesh as the current sim.
     ! Otherwise, it will point to its own internal masked mesh_t.
     type(mesh_t), pointer :: msh => null()
     !> Indicates if the instance of msh is internally managed
     logical :: internal_mesh = .false.

     !> The dofmap to use for our fields. We know for sure this dofmap is
     !! going to be different from case%fluid%dm_Xh because we force the
     !! user to provide either a point_zone or a different polynomial order.
     type(dofmap_t) :: dof

     !> Points to the space our fields are based on. If no point zone, will
     !! point to the same space as the current sim. Otherwise, will point to its
     !! own masked one.
     type(space_t), pointer :: Xh => null()
     !> Indicates if the instance of msh is internally managed. Also serves
     !! as a flag for determining whether or not to perform space-to-space
     !! interpolation
     logical :: internal_space = .false.

     !> Space-to-space interpolation
     type(interpolator_t) :: interpolator
     !> Temporary work field for space-to-space interpolation
     type(field_t) :: wk
     !> New polynomial order to use for subsampling.
     integer :: lx = -1

   contains
     !> Constructor.
     procedure, pass(this) :: init => field_subsampler_init
     !> Destructor.
     procedure, pass(this) :: free => field_subsampler_free
     !> Compute the indicator (called according to the simcomp controller).
     procedure, pass(this) :: compute_ => compute_wrapper
  end type field_subsampler_t

  abstract interface
     subroutine compute_intrf(this, time)
       import :: field_subsampler_t, time_state_t
       class(field_subsampler_t), intent(inout) :: this
       type(time_state_t), intent(in) :: time
     end subroutine compute_intrf
  end interface


contains

  subroutine compute_wrapper(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%compute_impl(time)

  end subroutine compute_wrapper

  subroutine dummy_compute(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call neko_error("field_subsampler must be initialized first!")

  end subroutine dummy_compute

  !> Constructor.
  subroutine field_subsampler_init(this, json, case)
    class(field_subsampler_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    character(len=:), allocatable :: name
    character(len=20), allocatable :: which_fields(:)
    integer :: p

    character(len=:), allocatable :: pz_name
    class(point_zone_t), pointer :: point_zone
    point_zone => null()

    call json_get_or_default(json, "name", name, "field_subsampler")
    call json_get(json, "fields", which_fields)
    call json_get_or_default(json, "polynomial_order", p, case%fluid%Xh%lx - 1)

    call this%init_base(json, case)

    if (json%valid_path('point_zone')) then

       call json_get(json, "point_zone", pz_name)
       point_zone => neko_point_zone_registry%get_point_zone(trim(pz_name))

       call field_subsampler_init_from_components(this, name, which_fields, &
            p, point_zone = point_zone)

    else

       ! Throw an error if the user doesn't subsample by point zone or by
       ! different polynomial order. (in that case there is no subsampling
       ! done at all)
       if (p .eq. case%fluid%Xh%lx - 1) then
          call neko_error("No subsampling strategy defined. Please " // &
               "specify either a different polynomial order or a " // &
               "point zone for subsampling.")
       end if

       call field_subsampler_init_from_components(this, name, which_fields, p)

    end if

  end subroutine field_subsampler_init

  !> Actual constructor.
  !! @param name The unique name of the simcomp.
  !! @param which_fields The names of the fields to be subsampled.
  !! @param p The new polynomial order to use for subsampling. NOTE: this is
  !! NOT `lx`! `lx = p + 1`.
  !! @param point_zone The point zone to use for subsampling.
  subroutine field_subsampler_init_from_components(this, name, which_fields, &
       p, point_zone)
    class(field_subsampler_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=20), intent(in) :: which_fields(:)
    integer, intent(in), optional :: p
    class(point_zone_t), pointer, intent(in), optional :: point_zone

    this%name = name
    this%field_names = which_fields

    ! ========================================================================
    ! Polynomial order / space interpolation

    this%internal_space = .false.
    if (present(p)) then

       this%lx = p + 1

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

       if (.not. point_zone%full_elements) &
            call neko_error("full_elements must be enabled for subsampling")

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
    ! Register fields in the registry and initialize the field_list

    block
      integer :: i, n_fields
      character(len=2048) :: field_name
      n_fields = size(this%field_names)

      call this%fields%init(n_fields)

      do i = 1, n_fields
         field_name = this%name // "_" // this%field_names(i)
         call neko_registry%add_field(this%dof, field_name)

         this%fields%items(i)%ptr => neko_registry%get_field(trim(field_name))
      end do
    end block

    !========================================================================
    ! Assign the correct compute() depending on whether interpolation is
    ! enabled or not

    ! Only space-to-space interpolation
    if (this%internal_space .and. .not. this%internal_mesh) then
       this%compute_impl => field_subsampler_compute_Xh

       ! Only masked copying based on the point_zone
    else if (.not. this%internal_space .and. this%internal_mesh) then
       this%compute_impl => field_subsampler_compute_pz

       ! Both masked copying and space-to-space interpolation
    else

       ! initialize a work array with masked mesh but same poly. order as
       ! the simulation
       call this%wk%init(this%msh, this%case%fluid%Xh)

       this%compute_impl => field_subsampler_compute_pz_Xh
    end if

  end subroutine field_subsampler_init_from_components

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

    call this%wk%free()
    call this%fields%free()
    call this%free_base()

  end subroutine field_subsampler_free

  !> Subsample the fields based only on a point zone, no space-to-space
  !! interpolation.
  subroutine field_subsampler_compute_pz(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: f
    integer :: i, n, n_mask


    n = this%dof%size()
    n_mask = this%point_zone%mask%size()

    do i = 1, size(this%field_names)

       f => neko_registry%get_field(this%field_names(i))

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_gather_copy_aligned(this%fields%x_d(i), f%x_d, &
               this%point_zone%mask%get_d(), n, n_mask)
       else
          call masked_gather_copy(this%fields%x(i), f%x, &
               this%point_zone%mask%get(), n, n_mask)
       end if

    end do

  end subroutine field_subsampler_compute_pz

  !> Subsample fields wihout any point zone, only space-to-space
  !! interpolation is enabled.
  subroutine field_subsampler_compute_Xh(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: f
    integer :: i

    do i = 1, size(this%field_names)
       f => neko_registry%get_field(this%field_names(i))
       call this%interpolator%map(this%fields%x(i), f%x, &
            this%msh%nelv, this%Xh)
    end do

  end subroutine field_subsampler_compute_Xh

  !> Subsample the fields based a point zone, and with space-to-space
  !! interpolation enabled.
  subroutine field_subsampler_compute_pz_Xh(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: f
    integer :: i, n, n_mask

    n = this%dof%size()
    n_mask = this%point_zone%mask%size()

    do i = 1, size(this%field_names)

       f => neko_registry%get_field(this%field_names(i))

       ! =====================================================================
       ! First, do a masked copy onto the submesh but that has the same poly.
       ! order

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_gather_copy_aligned(this%wk%x_d, f%x_d, &
               this%point_zone%mask%get_d(), n, n_mask)
       else
          call masked_gather_copy(this%wk%x, f%x, &
               this%point_zone%mask%get(), n, n_mask)
       end if

       ! =====================================================================
       ! Then, map the two different spaces
       call this%interpolator%map(this%fields%x(i), this%wk%x, &
            this%msh%nelv, this%Xh)

    end do

  end subroutine field_subsampler_compute_pz_Xh

end module field_subsampler
