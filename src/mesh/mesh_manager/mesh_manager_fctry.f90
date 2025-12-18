! Copyright (c) 2019-2025, The Neko Authors
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
!> Defines a factory subroutine for mesh manager.
submodule (mesh_manager) mesh_manager_fctry
  use utils, only : neko_type_error
  use json_utils, only : json_get
  use mesh_manager_p4est, only : mesh_manager_p4est_t
  use mesh_manager_transfer_p4est, only : mesh_manager_transfer_p4est_t

  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: MESHMNG_KNOWN_TYPES(1) = [character(len=20) :: &
       "p4est"]
contains

  !> Mesh manager factory. Both constructs and initializes the object.
  !! @param[inout]  object          The object to be initialised.
  !! @param[inout]  json            JSON object initialising the mesh manager.
  !! @param[in]     ifpartition     partitioning flag
  module subroutine mesh_manager_factory(object, json, ifpartition)
    class(mesh_manager_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    logical, intent(in) :: ifpartition
    character(len=:), allocatable :: type_name

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if
    call json_get(json, "type", type_name)
    ! Allocate
    call mesh_manager_allocator(object, type_name)

    ! Initialise base types
    call object%init_base(type_name)
    deallocate(type_name)

    ! data redistribution
    call object%transfer%init_base(ifpartition)

  end subroutine mesh_manager_factory

  !> Mesh manager allocator.
  !! @param[inout]  object      The object to be allocated.
  !! @param[in]     type_name   The name of the type to allocate.
  module subroutine mesh_manager_allocator(object, type_name)
    class(mesh_manager_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name

    select case (trim(type_name))
    case ("p4est")
       ! allocate type itself
       allocate(mesh_manager_p4est_t::object)
       ! allocate redistribution type
       allocate(mesh_manager_transfer_p4est_t::object%transfer)
    case default
       call neko_type_error("mesh manager", type_name, MESHMNG_KNOWN_TYPES)
    end select
  end subroutine mesh_manager_allocator

end submodule mesh_manager_fctry
