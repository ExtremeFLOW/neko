! Copyright (c) 2021-2026, The Neko Authors
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
submodule (ax_product) ax_helm_fctry
  use neko_config, only : NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
       NEKO_BCKND_DEVICE
  use ax_helm_device, only : ax_helm_device_t
  use ax_helm_xsmm, only : ax_helm_xsmm_t
  use ax_helm_sx, only : ax_helm_sx_t
  use ax_helm, only : ax_helm_t
  use ax_helm_cpu, only : ax_helm_cpu_t
  use ax_helm_full_cpu, only : ax_helm_full_cpu_t
  use ax_helm_full_device, only : ax_helm_full_device_t
  use utils, only : neko_error, neko_type_error, &
       neko_type_registration_error
  implicit none

  ! List of all possible types created by the allocator routine
  character(len=20) :: AX_HELM_KNOWN_TYPES(2) = [character(len=20) :: &
       "standard", &
       "full"]

contains

  !> Allocate a Helmholtz problem matrix-vector product.
  !! The implementation is selected by name and compute backend.
  !! @param object The matrix-vector product type to be allocated.
  !! @param type_name The name of the matrix-vector product type.
  module subroutine ax_helm_allocator(object, type_name)
    class(ax_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    integer :: i

    if (allocated(object)) then
       deallocate(object)
    end if

    select case (trim(type_name))
    case ("standard")
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(ax_helm_sx_t::object)
       else if (NEKO_BCKND_XSMM .eq. 1) then
          allocate(ax_helm_xsmm_t::object)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(ax_helm_device_t::object)
       else
          allocate(ax_helm_cpu_t::object)
       end if
    case ("full")
       if (NEKO_BCKND_XSMM .eq. 1) then
          call neko_error("Full stress formulation is only available " // &
               "on the CPU and device")
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(ax_helm_full_device_t::object)
       else
          allocate(ax_helm_full_cpu_t::object)
       end if
    case default
       do i = 1, ax_helm_registry_size
          if (trim(type_name) .eq. &
               trim(ax_helm_registry(i)%type_name)) then
             call ax_helm_registry(i)%allocator(object)
             return
          end if
       end do

       call neko_type_error("matrix-vector product", type_name, &
            AX_HELM_KNOWN_TYPES)
    end select

  end subroutine ax_helm_allocator

  !> Register a custom matrix-vector product allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_ax_helm(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(ax_helm_allocate), pointer, intent(in) :: allocator
    type(ax_helm_allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(AX_HELM_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(AX_HELM_KNOWN_TYPES(i))) then
          call neko_type_registration_error("matrix-vector product", &
               type_name, .true.)
       end if
    end do

    do i = 1, ax_helm_registry_size
       if (trim(type_name) .eq. trim(ax_helm_registry(i)%type_name)) then
          call neko_type_registration_error("matrix-vector product", &
               type_name, .false.)
       end if
    end do

    if (ax_helm_registry_size .eq. 0) then
       allocate(ax_helm_registry(1))
    else
       allocate(temp(ax_helm_registry_size + 1))
       temp(1:ax_helm_registry_size) = ax_helm_registry
       call move_alloc(temp, ax_helm_registry)
    end if

    ax_helm_registry_size = ax_helm_registry_size + 1
    ax_helm_registry(ax_helm_registry_size)%type_name = type_name
    ax_helm_registry(ax_helm_registry_size)%allocator => allocator
  end subroutine register_ax_helm

end submodule ax_helm_fctry
