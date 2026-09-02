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
!
!> Implements `wall_model_t`.
module wall_model
  use num_types, only : rp
  use field, only : field_t
  use json_module, only : json_file
  use registry, only : neko_registry
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use vector, only : vector_t
  use utils, only : neko_error, nonlinear_index
  use wall_sampler, only : wall_sampler_t
  use wall_sampler_fctry, only : wall_sampler_factory
  use user_intf, only : user_t
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  use device, only : device_get_ptr
  use wall_model_device, only : wall_model_compute_mag_field_device
  use json_utils, only : json_get, json_get_or_lookup
  implicit none
  private

  !> Base abstract type for wall-stress models for wall-modelled LES.
  type, abstract, public :: wall_model_t
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Map of degrees of freedom.
     type(dofmap_t), pointer :: dof => null()
     !> The dynamic viscosity of the fluid
     type(field_t), pointer :: mu => null()
     !> The density of the fluid
     type(field_t), pointer :: rho => null()
     !> The name of the scheme for which the wall model is used.
     !! This is used to identify the name of the mu an rho fields in the
     !! registry.
     character(len=:), allocatable :: scheme_name
     !> The boundary condition mask. Stores the array size at index zero!
     integer, pointer :: msk(:) => null()
     type(c_ptr) :: msk_d = C_NULL_PTR
     !> The boundary condition facet ids. Stores the array size at index zero!
     integer, pointer :: facet(:) => null()
     !> The x component of the shear stress.
     type(vector_t) :: tau_x
     !> The y component of the shear stress.
     type(vector_t) :: tau_y
     !> The z component of the shear stress.
     type(vector_t) :: tau_z
     !> The x component of the normal.
     type(vector_t) :: n_x
     !> The y component of the normal.
     type(vector_t) :: n_y
     !> The z component of the normal.
     type(vector_t) :: n_z
     !> Strategy used to sample solution fields away from the wall.
     class(wall_sampler_t), allocatable :: sampler
     !> Number of nodes in the boundary
     integer :: n_nodes = 0
     !> The 3D field with the computed stress magnitude at the boundary.
     type(field_t), pointer :: tau_field => null()
   contains
     !> Constructor for the wall_model_t (base) class.
     procedure, pass(this) :: init_base => wall_model_init_base
     !> Base type implementation of the deferred `partial_init`.
     procedure, pass(this) :: partial_init_base => wall_model_partial_init_base
     !> Base type implementation of the deferred `finilize`.
     procedure, pass(this) :: finalize_base => wall_model_finalize_base
     !> Destructor for the wall_model_t (base) class.
     procedure, pass(this) :: free_base => wall_model_free_base
     !> Compute the wall shear stress's magnitude.
     procedure, pass(this) :: compute_mag_field => wall_model_compute_mag_field
     !> The common constructor.
     procedure(wall_model_init), pass(this), deferred :: init
     !> A part of the constructor that parses the JSON without initializing the
     !! base `wall_model_t` type. Used in `wall_model_bc_t` during the init
     !! stage of the bc construction. So, in this routine you cannot use the
     !! bc mask and facets. The construction can be finished by calling
     !! `finalize_base`. This generally follows the same pattern as in `bc_t`.
     procedure(wall_model_partial_init), pass(this), deferred :: partial_init
     !> Finalization of the partial construction, similar to `bc_t`.
     procedure(wall_model_finalize), pass(this), deferred :: finalize
     !> Destructor.
     procedure(wall_model_free), pass(this), deferred :: free
     !> Compute the wall shear stress.
     procedure(wall_model_compute), pass(this), deferred :: compute
     !> Build geometry shared by wall models and samplers.
     procedure, pass(this) :: setup_geometry => wall_model_setup_geometry
     !> Check that n_samples is one.
     procedure, pass(this) :: validate_single_sample => &
          wall_model_validate_single_sample
  end type wall_model_t

  abstract interface
     !> Compute wall shear stress.
     !! @param t The time value.
     !! @param tstep The current time-step.
     subroutine wall_model_compute(this, t, tstep)
       import wall_model_t, rp
       class(wall_model_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
     end subroutine wall_model_compute
  end interface

  abstract interface
     !> Common constructor.
     !! @param scheme_name The name of the scheme for which the wall model is
     !! used.
     !! @param coef SEM coefficients.
     !! @param msk The boundary mask.
     !! @param facet The boundary facets.
     !! @param json A dictionary with parameters.
     subroutine wall_model_init(this, scheme_name, coef, msk, facet, &
          json)
       import wall_model_t, json_file, dofmap_t, coef_t, rp
       class(wall_model_t), intent(inout) :: this
       character(len=*), intent(in) :: scheme_name
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
       type(json_file), intent(inout) :: json
     end subroutine wall_model_init
  end interface

  abstract interface
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     !! @param coef SEM coefficients.
     !! @param scheme_name Name of the owning scheme.
     !! @param json A dictionary with parameters.
     subroutine wall_model_partial_init(this, coef, scheme_name, json)
       import wall_model_t, json_file, dofmap_t, coef_t, rp, user_t
       class(wall_model_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       character(len=*), intent(in) :: scheme_name
       type(json_file), intent(inout) :: json
     end subroutine wall_model_partial_init
  end interface

  abstract interface
     !> Finilzation of partial construction, similar to `bc_t`
     !! @param msk The boundary mask.
     !! @param facet The boundary facets.
     subroutine wall_model_finalize(this, msk, facet, bc_name, user)
       import wall_model_t, user_t
       class(wall_model_t), intent(inout) :: this
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
       character(len=*), optional, intent(in) :: bc_name
       type(user_t), target, optional, intent(in) :: user
     end subroutine wall_model_finalize
  end interface

  abstract interface
     !> Destructor.
     subroutine wall_model_free(this)
       import wall_model_t
       class(wall_model_t), intent(inout) :: this
     end subroutine wall_model_free
  end interface

  interface
     !> Wall model factory. Both constructs and initializes the object.
     !! @param object The object to be allocated.
     !! @param scheme_name The name of the scheme for which the wall model is
     !! used.
     !! @param coef SEM coefficients.
     !! @param msk The boundary mask.
     !! @param facet The boundary facets.
     !! @param json A dictionary with parameters.
     module subroutine wall_model_factory(object, scheme_name, coef, msk, &
          facet, json)
       class(wall_model_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: scheme_name
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
       type(json_file), intent(inout) :: json
     end subroutine wall_model_factory
  end interface

  interface
     !> Wall model allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the type to allocate.
     module subroutine wall_model_allocator(object, type_name)
       class(wall_model_t), allocatable, intent(inout) :: object
       character(len=:), allocatable, intent(in) :: type_name
     end subroutine wall_model_allocator
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for an object allocator.
  !! Implemented in the user modules, should allocate the `obj` to the custom
  !! user type.
  abstract interface
     subroutine wall_model_allocate(obj)
       import wall_model_t
       class(wall_model_t), allocatable, intent(inout) :: obj
     end subroutine wall_model_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_wall_model(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(wall_model_allocate), pointer, intent(in) :: allocator
     end subroutine register_wall_model
  end interface

  ! A name-allocator pair for user-defined types. A helper type to define a
  ! registry of custom allocators.
  type allocator_entry
     character(len=20) :: type_name
     procedure(wall_model_allocate), pointer, nopass :: allocator
  end type allocator_entry

  !> Registry of wall model allocators for user-defined types
  type(allocator_entry), allocatable :: wall_model_registry(:)

  !> The size of the `wall_model_registry`
  integer :: wall_model_registry_size = 0

  public :: wall_model_factory, wall_model_allocator, register_wall_model, &
       wall_model_allocate

contains
  !> Constructor for the wall_model_t (base) class.
  !! @param coef SEM coefficients.
  !! @param msk The underlying mask of the boundary condition.
  !! @param facet, The underlying facet index list of the boundary condition.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param sampler Configured sampling strategy; ownership is transferred.
  subroutine wall_model_init_base(this, scheme_name, coef, msk, facet, sampler)
    class(wall_model_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    integer, target, intent(in) :: msk(0:)
    integer, target, intent(in) :: facet(0:)
    character(len=*) :: scheme_name
    class(wall_sampler_t), allocatable, intent(inout) :: sampler

    call this%free_base

    this%coef => coef
    this%dof => coef%dof
    if (.not. allocated(sampler)) then
       call neko_error('Wall model sampler has not been initialized')
    end if
    call move_alloc(sampler, this%sampler)
    this%scheme_name = trim(scheme_name)
    this%mu => neko_registry%get_field_by_name(this%scheme_name // "_mu")
    this%rho => neko_registry%get_field_by_name(this%scheme_name // &
         "_rho")

    call neko_registry%add_field(this%dof, "tau", &
         ignore_existing = .true.)
    this%tau_field => neko_registry%get_field("tau")

    call this%finalize_base(msk, facet)
  end subroutine wall_model_init_base

  !> Partial initialization based on JSON, prior to knowing the mask and facets.
  !! @param coef SEM coefficients.
  !! @param The name of the scheme for which the wall model is used.
  !! @param json A dictionary with parameters.
  subroutine wall_model_partial_init_base(this, coef, scheme_name, json)
    class(wall_model_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: scheme_name
    type(json_file), intent(inout) :: json

    call this%free_base()

    this%coef => coef
    this%dof => coef%dof
    call wall_sampler_factory(this%sampler, json)
    this%scheme_name = trim(scheme_name)

    this%mu => neko_registry%get_field_by_name(this%scheme_name // "_mu")
    this%rho => neko_registry%get_field_by_name(this%scheme_name // &
         "_rho")

    call neko_registry%add_field(this%dof, "tau", &
         ignore_existing = .true.)
    this%tau_field => neko_registry%get_field("tau")
  end subroutine wall_model_partial_init_base

  subroutine wall_model_finalize_base(this, msk, facet, bc_name, user)
    class(wall_model_t), intent(inout) :: this
    integer, target, intent(in) :: msk(0:)
    integer, target, intent(in) :: facet(:)
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user

    this%msk(0:msk(0)) => msk
    if (NEKO_BCKND_DEVICE .eq. 1) this%msk_d = device_get_ptr(msk)
    this%facet(0:msk(0)) => facet

    call this%tau_x%init(this%msk(0))
    call this%tau_y%init(this%msk(0))
    call this%tau_z%init(this%msk(0))

    call this%n_x%init(this%msk(0))
    call this%n_y%init(this%msk(0))
    call this%n_z%init(this%msk(0))

    call this%setup_geometry()
    if (.not. allocated(this%sampler)) then
       call neko_error('Wall model sampler has not been initialized')
    end if
    if (present(bc_name)) then
       if (present(user)) then
          call this%sampler%finalize(this%coef, this%msk, this%facet, &
               this%n_x, this%n_y, this%n_z, bc_name, user)
       else
          call this%sampler%finalize(this%coef, this%msk, this%facet, &
               this%n_x, this%n_y, this%n_z, bc_name)
       end if
    else
       if (this%sampler%user_values) then
          call neko_error('Wall model user sampler has no boundary context')
       end if
       call this%sampler%finalize(this%coef, this%msk, this%facet, &
            this%n_x, this%n_y, this%n_z)
    end if

  end subroutine wall_model_finalize_base

  !> Destructor for the wall_model_t (base) class.
  subroutine wall_model_free_base(this)
    class(wall_model_t), intent(inout) :: this

    nullify(this%coef)
    nullify(this%msk)
    nullify(this%facet)
    nullify(this%tau_field)
    nullify(this%mu)
    nullify(this%rho)
    this%msk_d = C_NULL_PTR

    call this%tau_x%free()
    call this%tau_y%free()
    call this%tau_z%free()

    if (allocated(this%sampler)) then
       call this%sampler%free()
       deallocate(this%sampler)
    end if

    if (allocated(this%scheme_name)) then
       deallocate(this%scheme_name)
    end if

    call this%n_x%free()
    call this%n_y%free()
    call this%n_z%free()

    nullify(this%dof)
  end subroutine wall_model_free_base

  !> Build wall-node normals independently of the sampling strategy.
  subroutine wall_model_setup_geometry(this)
    class(wall_model_t), intent(inout) :: this
    integer :: n_nodes, fid, idx(4), i, linear
    real(kind=rp) :: normal(3)

    call this%coef%require_facets('wall_model')

    n_nodes = this%msk(0)
    this%n_nodes = n_nodes

    do i = 1, n_nodes
       linear = this%msk(i)
       fid = this%facet(i)
       idx = nonlinear_index(linear, this%coef%Xh%lx, this%coef%Xh%ly,&
            this%coef%Xh%lz)
       normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), fid)

       this%n_x%x(i) = normal(1)
       this%n_y%x(i) = normal(2)
       this%n_z%x(i) = normal(3)

    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%n_x%x, this%n_x%x_d, n_nodes, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%n_y%x, this%n_y%x_d, n_nodes, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%n_z%x, this%n_z%x_d, n_nodes, HOST_TO_DEVICE, &
            sync = .true.)
    end if

  end subroutine wall_model_setup_geometry

  subroutine wall_model_validate_single_sample(this)
    class(wall_model_t), intent(inout) :: this

    if (this%sampler%n_samples /= 1) then
       call neko_error('This wall model requires exactly one sample per node')
    end if
  end subroutine wall_model_validate_single_sample

  subroutine wall_model_compute_mag_field(this)
    class(wall_model_t), intent(inout) :: this
    integer :: i, m
    real(kind=rp) :: magtau

    m = this%msk(0)
    if (m > 0) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call wall_model_compute_mag_field_device(this%tau_x%x_d, &
               this%tau_y%x_d, &
               this%tau_z%x_d, &
               this%tau_field%x_d, &
               this%msk_d, m)
       else
          do i = 1, m
             magtau = sqrt(this%tau_x%x(i)**2 + &
                  this%tau_y%x(i)**2 + &
                  this%tau_z%x(i)**2)
             this%tau_field%x(this%msk(i),1,1,1) = magtau
          end do
       end if
    end if

  end subroutine wall_model_compute_mag_field

end module wall_model
