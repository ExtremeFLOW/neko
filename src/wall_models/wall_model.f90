! Copyright (c) 2024, The Neko Authors
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
  use field_registry, only : neko_field_registry
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use vector, only : vector_t
  use utils, only : neko_error, nonlinear_index
  use math, only : glmin, glmax
  use comm, only : pe_rank
  use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  use file, only : file_t
  use field_registry, only : neko_field_registry
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  use device, only : device_map, device_free, device_get_ptr
  use wall_model_device, only : wall_model_compute_mag_field_device
  use json_utils, only : json_get
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
     !> The r indices of the sampling points
     integer, allocatable :: ind_r(:)
     type(c_ptr) :: ind_r_d = C_NULL_PTR
     !> The s indices of the sampling points
     integer, allocatable :: ind_s(:)
     type(c_ptr) :: ind_s_d = C_NULL_PTR
     !> The t indices of the sampling points
     integer, allocatable :: ind_t(:)
     type(c_ptr) :: ind_t_d = C_NULL_PTR
     !> The element indices of the sampling points
     integer, allocatable :: ind_e(:)
     type(c_ptr) :: ind_e_d = C_NULL_PTR
     !> The sampling height
     type(vector_t) :: h
     !> Sampling index
     integer :: h_index = 0
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
     !! This constructor gets the scheme_name from the JSON. This is currently
     !! hacked in by the `pnpn_bc_factory`.
     procedure(wall_model_partial_init), pass(this), deferred :: partial_init
     !> Finalization of the partial construction, similar to `bc_t`.
     procedure(wall_model_finalize), pass(this), deferred :: finalize
     !> Destructor.
     procedure(wall_model_free), pass(this), deferred :: free
     !> Compute the wall shear stress.
     procedure(wall_model_compute), pass(this), deferred :: compute
     !> Find the sampling points based on the value of `h_index`.
     procedure, pass(this) :: find_points => wall_model_find_points
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
     !! @param h_index The off-wall index of the sampling cell.
     !! @param json A dictionary with parameters.
     subroutine wall_model_init(this, scheme_name, coef, msk, facet, &
          h_index, json)
       import wall_model_t, json_file, dofmap_t, coef_t, rp
       class(wall_model_t), intent(inout) :: this
       character(len=*), intent(in) :: scheme_name
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
       integer, intent(in) :: h_index
       type(json_file), intent(inout) :: json
     end subroutine wall_model_init
  end interface

  abstract interface
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     subroutine wall_model_partial_init(this, coef, json)
       import wall_model_t, json_file, dofmap_t, coef_t, rp
       class(wall_model_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       type(json_file), intent(inout) :: json
     end subroutine wall_model_partial_init
  end interface

  abstract interface
     !> Finilzation of partial construction, similar to `bc_t`
     !! @param msk The boundary mask.
     !! @param facet The boundary facets.
     subroutine wall_model_finalize(this, msk, facet)
       import wall_model_t
       class(wall_model_t), intent(inout) :: this
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
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
     !! @param h_index The off-wall index of the sampling cell.
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
  !! @param index The off-wall index of the sampling point.
  subroutine wall_model_init_base(this, scheme_name, coef, msk, facet, index)
    class(wall_model_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    integer, target, intent(in) :: msk(0:)
    integer, target, intent(in) :: facet(0:)
    character(len=*) :: scheme_name
    integer, intent(in) :: index

    call this%free_base

    this%coef => coef
    this%dof => coef%dof
    this%h_index = index
    this%scheme_name = trim(scheme_name)
    this%mu => neko_field_registry%get_field_by_name(this%scheme_name // "_mu")
    this%rho => neko_field_registry%get_field_by_name(this%scheme_name // &
         "_rho")

    call neko_field_registry%add_field(this%dof, "tau", &
         ignore_existing = .true.)
    this%tau_field => neko_field_registry%get_field("tau")

    call this%finalize_base(msk, facet)
  end subroutine wall_model_init_base

  !> Partial initialization based on JSON, prior to knowing the mask and facets.
  !! @param coef SEM coefficients.
  !! @param The name of the scheme for which the wall model is used.
  !! @param json A dictionary with parameters.
  subroutine wall_model_partial_init_base(this, coef, json)
    class(wall_model_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%free_base()

    this%coef => coef
    this%dof => coef%dof
    call json_get(json, "h_index", this%h_index)
    call json_get(json, "scheme_name", this%scheme_name)

    this%mu => neko_field_registry%get_field_by_name(this%scheme_name // "_mu")
    this%rho => neko_field_registry%get_field_by_name(this%scheme_name // &
         "_rho")

    call neko_field_registry%add_field(this%dof, "tau", &
         ignore_existing = .true.)
    this%tau_field => neko_field_registry%get_field("tau")
  end subroutine wall_model_partial_init_base

  subroutine wall_model_finalize_base(this, msk, facet)
    class(wall_model_t), intent(inout) :: this
    integer, target, intent(in) :: msk(0:)
    integer, target, intent(in) :: facet(:)

    this%msk(0:msk(0)) => msk
    if (NEKO_BCKND_DEVICE .eq. 1) this%msk_d = device_get_ptr(msk)
    this%facet(0:msk(0)) => facet

    call this%tau_x%init(this%msk(0))
    call this%tau_y%init(this%msk(0))
    call this%tau_z%init(this%msk(0))

    allocate(this%ind_r(this%msk(0)))
    allocate(this%ind_s(this%msk(0)))
    allocate(this%ind_t(this%msk(0)))
    allocate(this%ind_e(this%msk(0)))

    call this%h%init(this%msk(0))
    call this%n_x%init(this%msk(0))
    call this%n_y%init(this%msk(0))
    call this%n_z%init(this%msk(0))

    call this%find_points()

    ! Initialize pointers for device
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%ind_r, this%ind_r_d, this%n_nodes)
       call device_map(this%ind_s, this%ind_s_d, this%n_nodes)
       call device_map(this%ind_t, this%ind_t_d, this%n_nodes)
       call device_map(this%ind_e, this%ind_e_d, this%n_nodes)
       call device_memcpy(this%ind_r, this%ind_r_d, this%n_nodes, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%ind_s, this%ind_s_d, this%n_nodes, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%ind_t, this%ind_t_d, this%n_nodes, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%ind_e, this%ind_e_d, this%n_nodes, &
            HOST_TO_DEVICE, sync = .false.)
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

    call this%tau_x%free()
    call this%tau_y%free()
    call this%tau_z%free()

    if (allocated(this%ind_r)) then
       deallocate(this%ind_r)
    end if
    if (allocated(this%ind_s)) then
       deallocate(this%ind_s)
    end if
    if (allocated(this%ind_t)) then
       deallocate(this%ind_t)
    end if
    if (allocated(this%ind_e)) then
       deallocate(this%ind_e)
    end if

    if (c_associated(this%msk_d)) then
       call device_free(this%msk_d)
    end if
    if (c_associated(this%ind_r_d)) then
       call device_free(this%ind_r_d)
    end if
    if (c_associated(this%ind_s_d)) then
       call device_free(this%ind_s_d)
    end if
    if (c_associated(this%ind_t_d)) then
       call device_free(this%ind_t_d)
    end if
    if (c_associated(this%ind_e_d)) then
       call device_free(this%ind_e_d)
    end if

    if (allocated(this%scheme_name)) then
       deallocate(this%scheme_name)
    end if


    call this%h%free()
    call this%n_x%free()
    call this%n_y%free()
    call this%n_z%free()

    nullify(this%dof)
  end subroutine wall_model_free_base

  !> Find sampling points based on the requested index.
  subroutine wall_model_find_points(this)
    class(wall_model_t), intent(inout) :: this
    integer :: n_nodes, fid, idx(4), i, linear
    real(kind=rp) :: normal(3), p(3), x, y, z, xw, yw, zw, magp
    real(kind=rp) :: hmin, hmax
    type(field_t), pointer :: h_field
    type(file_t) :: h_file
    character(len=LOG_SIZE), allocatable :: log_msg

    n_nodes = this%msk(0)
    this%n_nodes = n_nodes

    call neko_field_registry%add_field(this%coef%dof, "sampling_height", &
         ignore_existing=.true.)

    h_field => neko_field_registry%get_field_by_name("sampling_height")

    do i = 1, n_nodes
       linear = this%msk(i)
       fid = this%facet(i)
       idx = nonlinear_index(linear, this%coef%Xh%lx, this%coef%Xh%ly,&
            this%coef%Xh%lz)
       normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), fid)

       this%n_x%x(i) = normal(1)
       this%n_y%x(i) = normal(2)
       this%n_z%x(i) = normal(3)

       ! inward normal
       normal = -normal

       select case (fid)
       case (1)
          this%ind_r(i) = idx(1) + this%h_index
          this%ind_s(i) = idx(2)
          this%ind_t(i) = idx(3)
       case (2)
          this%ind_r(i) = idx(1) - this%h_index
          this%ind_s(i) = idx(2)
          this%ind_t(i) = idx(3)
       case (3)
          this%ind_r(i) = idx(1)
          this%ind_s(i) = idx(2) + this%h_index
          this%ind_t(i) = idx(3)
       case (4)
          this%ind_r(i) = idx(1)
          this%ind_s(i) = idx(2) - this%h_index
          this%ind_t(i) = idx(3)
       case (5)
          this%ind_r(i) = idx(1)
          this%ind_s(i) = idx(2)
          this%ind_t(i) = idx(3) + this%h_index
       case (6)
          this%ind_r(i) = idx(1)
          this%ind_s(i) = idx(2)
          this%ind_t(i) = idx(3) - this%h_index
       case default
          call neko_error("The face index is not correct ")
       end select
       this%ind_e(i) = idx(4)

       ! Location of the wall node
       xw = this%dof%x(idx(1), idx(2), idx(3), idx(4))
       yw = this%dof%y(idx(1), idx(2), idx(3), idx(4))
       zw = this%dof%z(idx(1), idx(2), idx(3), idx(4))

       ! Location of the sampling point
       x = this%dof%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
            this%ind_e(i))
       y = this%dof%y(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
            this%ind_e(i))
       z = this%dof%z(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
            this%ind_e(i))

       ! Vector from the sampling point to the wall
       p(1) = x - xw
       p(2) = y - yw
       p(3) = z - zw

       ! Total distance to the sampling point
       magp = sqrt(p(1)**2 + p(2)**2 + p(3)**2)

       ! Project on the normal direction to get h
       this%h%x(i) = p(1)*normal(1) + p(2)*normal(2) + p(3)*normal(3)

       h_field%x(linear,1,1,1) = this%h%x(i)

       ! Look at how much the total distance distance from the normal and warn
       ! if significant
       if ((this%h%x(i) - magp) / magp > 0.1) then
          write(log_msg,*) "Significant misalignment between wall normal and"
          call neko_log%message(log_msg, NEKO_LOG_DEBUG)
          write(log_msg,*) "sampling point direction at wall node", xw, yw, zw
          call neko_log%message(log_msg, NEKO_LOG_DEBUG)
       end if
    end do

!    hmin = glmin(this%h%x, n_nodes)
!    hmax = glmax(this%h%x, n_nodes)
!    if (pe_rank .eq. 0) then
!       write(*, "(A, F10.4, F10.4)") "   h min / max:", hmin, hmax
!    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%h%x, this%h%x_d, n_nodes, HOST_TO_DEVICE,&
            sync = .false.)
       call device_memcpy(this%n_x%x, this%n_x%x_d, n_nodes, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%n_y%x, this%n_y%x_d, n_nodes, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%n_z%x, this%n_z%x_d, n_nodes, HOST_TO_DEVICE, &
            sync = .true.)
    end if

    ! Each wall_model bc will do a write unfortunately... But very helpful
    ! for setup debugging.
    call h_file%init("sampling_height.fld")
    call h_file%write(h_field)
  end subroutine wall_model_find_points

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
