! Copyright (c) 2021-2024, The Neko Authors
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
submodule (precon) precon_fctry
  use identity, only : ident_t
  use device_identity, only : device_ident_t
  use jacobi, only : jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use device_jacobi, only : device_jacobi_t
  use hsmg, only : hsmg_t
  use phmg, only : phmg_t
  use utils, only : neko_type_error, neko_type_registration_error
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: PC_KNOWN_TYPES(4) = [character(len=20) :: &
       "jacobi", &
       "hsmg", &
       "phmg", &
       "ident"]

contains

  !> Allocate a preconditioner
  module subroutine precon_allocator(pc, type_name)
    class(pc_t), allocatable, intent(inout) :: pc
    character(len=*), intent(in) :: type_name
    integer :: i

    if (allocated(pc)) then
       call precon_destroy(pc)
       deallocate(pc)
    end if

    select case (trim(type_name))
    case ('jacobi')
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_jacobi_t::pc)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(device_jacobi_t::pc)
       else
          allocate(jacobi_t::pc)
       end if
    case ('hsmg')
       allocate(hsmg_t::pc)
    case ('phmg')
       allocate(phmg_t::pc)
    case('ident')
       if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(device_ident_t::pc)
       else
          allocate(ident_t::pc)
       end if
    case default
       do i = 1, precon_registry_size
          if (trim(type_name) .eq. trim(precon_registry(i)%type_name)) then
             call precon_registry(i)%allocator(pc)
             return
          end if
       end do

       call neko_type_error("preconditioner", type_name, PC_KNOWN_TYPES)
    end select

  end subroutine precon_allocator

  !> Destroy a preconditioner
  module subroutine precon_destroy(pc)
    class(pc_t), allocatable, intent(inout) :: pc

    if (allocated(pc)) then
       select type (pcp => pc)
       type is (jacobi_t)
          call pcp%free()
       type is (sx_jacobi_t)
          call pcp%free()
       type is (device_jacobi_t)
          call pcp%free()
       type is (hsmg_t)
          call pcp%free()
       type is (phmg_t)
          call pcp%free()
       end select
    end if

  end subroutine precon_destroy

  !> Register a custom preconditioner allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_precon(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(precon_allocate), pointer, intent(in) :: allocator
    type(precon_allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(PC_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(PC_KNOWN_TYPES(i))) then
          call neko_type_registration_error("preconditioner", type_name, &
               .true.)
       end if
    end do

    do i = 1, precon_registry_size
       if (trim(type_name) .eq. trim(precon_registry(i)%type_name)) then
          call neko_type_registration_error("preconditioner", type_name, &
               .false.)
       end if
    end do

    if (precon_registry_size .eq. 0) then
       allocate(precon_registry(1))
    else
       allocate(temp(precon_registry_size + 1))
       temp(1:precon_registry_size) = precon_registry
       call move_alloc(temp, precon_registry)
    end if

    precon_registry_size = precon_registry_size + 1
    precon_registry(precon_registry_size)%type_name = type_name
    precon_registry(precon_registry_size)%allocator => allocator
  end subroutine register_precon

end submodule precon_fctry
