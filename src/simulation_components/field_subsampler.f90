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
  use logger, only : neko_log
  use device, only : DEVICE_TO_HOST, HOST_TO_DEVICE, device_memcpy
  use comm, only : pe_rank
  use utils, only : neko_warning, neko_error
  use field_writer, only : field_writer_t
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

  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Provides tools to calculate the spectral error indicator
  !! @details
  !! This is a posteriori error measure, based on the local properties of
  !! the spectral solution, which was developed by Mavriplis. This method
  !! formally only gives an indication of the error.
  type, public, extends(simulation_component_t) :: field_subsampler_t

     !> Fields to subsample.
     character(len=20), allocatable :: field_names(:)

     !> Point zone to use for subsampling.
     class(point_zone_t), pointer :: point_zone => null()
     !> Points to the mesh our fields are based on. If no point zone, will
     !! point to the same mesh as the current sim. Otherwise, will point to its
     !! own masked one.
     type(mesh_t), pointer :: msh => null()
     !> Indicates if the instance of msh is internally managed
     logical :: internal_mesh

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
     logical :: internal_space

     !> New polynomial order to use for subsampling.
     integer :: lx = -1

     !> Field writer controller for the output
     type(field_writer_t) :: writer

   contains
     !> Constructor.
     procedure, pass(this) :: init => field_subsampler_init
     !> Destructor.
     procedure, pass(this) :: free => field_subsampler_free
     !> Compute the indicator (called according to the simcomp controller).
     procedure, pass(this) :: compute_ => field_subsampler_compute
     !> Update the masked fields.
     !procedure, pass(this) :: update => field_subsampler_update

  end type field_subsampler_t

contains

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

       ! Warn the user in case they don't subsample by point zone or by
       ! different polynomial order, in that case we duplicate entire fields
       if (p .eq. case%fluid%Xh%lx - 1) then
          call neko_error("No subsampling strategy defined. Please " // &
               "specify either a lower polynomial order or a point zone " // &
               "for subsampling.")
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
    if (present(p)) then
       ! If we detect the same polynomial order, skip interpolation
       if (p .ne. this%case%fluid%Xh%lx - 1) then
          this%lx = p + 1

          ! Create the subsampled space
          allocate(this%Xh)
          call this%Xh%init(GLL, this%lx, this%lx, this%lx)
          this%internal_space = .true.
       else
          call neko_warning("Same polynomial order, no interpolation performed.")
       end if

    else

       ! Point to the same space as the fluid simulation
       this%Xh => this%case%fluid%Xh
       this%internal_space = .false.

    end if

    ! ========================================================================
    ! Point zone masking
    if (present(point_zone)) then
       this%point_zone => point_zone

       ! Create the subsampled mesh
       allocate(this%msh)
       call this%case%fluid%dm_Xh%msh%subset_by_mask(this%msh, &
            point_zone%mask, &
            this%lx, this%lx, this%lx)
    else
       this%msh => this%case%fluid%msh
    end if

    ! Create the subsampled dofmap
    call this%dof%init(this%msh, this%Xh)

    ! ========================================================================
    ! Field writer
    
    ! - Initialize the field_writer with the masked fields AND lower poly order
    !   => this will also add them to the registry!

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

    nullify(this%Xh)
    nullify(this%msh)

    call this%writer%free()
    call this%free_base()

  end subroutine field_subsampler_free

  !> Compute the spectral error indicator.
  subroutine field_subsampler_compute(this, time)
    class(field_subsampler_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

  end subroutine field_subsampler_compute


end module field_subsampler
