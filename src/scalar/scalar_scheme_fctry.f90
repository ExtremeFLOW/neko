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
!> Defines a factory subroutine for scalar schemes.
submodule (scalar_scheme_m) scalar_scheme_fctry
  use scalar_pnpn_m, only : scalar_pnpn_t
  use utils_m, only : neko_type_error, neko_type_registration_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: SCALAR_KNOWN_TYPES(1) = [character(len=20) :: &
       "pnpn"]

contains

  !> Scalar scheme factory. Both constructs and initializes the object.
  !! @param object The object to be created and initialized.
  !! @param msh The mesh.
  !! @param coef The coefficients.
  !! @param gs The gather-scatter.
  !! @param params The parameter dictionary in json.
  !! @param numerics_params The numerical parameter dictionary in json.
  !! @param user Type with user-defined procedures.
  !! @param chkp Checkpoint for restarts.
  !! @param ulag, vlag, wlag The lagged velocity fields.
  !! @param time_scheme The time scheme controller.
  !! @param rho The density field.
  module subroutine scalar_scheme_factory(object, msh, coef, gs, params, &
       numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
    class(scalar_scheme_t), allocatable, intent(inout) :: object
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    type(field_t), target, intent(in) :: rho
    character(len=:), allocatable :: type_name

    call json_get_or_default(params, "scheme", type_name, "pnpn")

    call scalar_scheme_allocator(object, type_name)

    call object%init(msh, coef, gs, params, numerics_params, user, chkp, &
         ulag, vlag, wlag, time_scheme, rho)

  end subroutine scalar_scheme_factory

  !> Scalar scheme allocator.
  !! @param object The object to be allocated.
  !! @param type_name The name of the scalar scheme type.
  module subroutine scalar_scheme_allocator(object, type_name)
    class(scalar_scheme_t), allocatable, intent(inout) :: object
    character(len=*), intent(in):: type_name
    integer :: i

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ("pnpn")
       allocate(scalar_pnpn_t::object)
    case default
       do i = 1, scalar_scheme_registry_size
          if (trim(type_name) == &
               trim(scalar_scheme_registry(i)%type_name)) then
             call scalar_scheme_registry(i)%allocator(object)
             return
          end if
       end do
       call neko_type_error("scalar scheme", trim(type_name), &
            SCALAR_KNOWN_TYPES)
    end select

  end subroutine scalar_scheme_allocator

  !> Register a custom scalar scheme allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_scalar_scheme(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(scalar_scheme_allocate), pointer, intent(in) :: allocator
    type(scalar_scheme_allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(SCALAR_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(SCALAR_KNOWN_TYPES(i))) then
          call neko_type_registration_error("scalar scheme", type_name, &
               .true.)
       end if
    end do

    do i = 1, scalar_scheme_registry_size
       if (trim(type_name) .eq. &
            trim(scalar_scheme_registry(i)%type_name)) then
          call neko_type_registration_error("scalar scheme", type_name, &
               .false.)
       end if
    end do

    ! Expand registry
    if (scalar_scheme_registry_size == 0) then
       allocate(scalar_scheme_registry(1))
    else
       allocate(temp(scalar_scheme_registry_size + 1))
       temp(1:scalar_scheme_registry_size) = scalar_scheme_registry
       call move_alloc(temp, scalar_scheme_registry)
    end if

    scalar_scheme_registry_size = scalar_scheme_registry_size + 1
    scalar_scheme_registry(scalar_scheme_registry_size)%type_name = type_name
    scalar_scheme_registry(scalar_scheme_registry_size)%allocator => allocator
  end subroutine register_scalar_scheme

end submodule scalar_scheme_fctry
