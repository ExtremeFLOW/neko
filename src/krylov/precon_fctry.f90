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
  use utils, only : concat_string_array, neko_error
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: PC_KNOWN_TYPES(3) = [character(len=20) :: &
     "jacobi", &
     "hsmg", &
     "ident"]

contains

  !> Create a preconditioner
  module subroutine precon_factory(pc, type_name)
    class(pc_t), target, allocatable, intent(inout) :: pc
    character(len=*), intent(in) :: type_name
    character(len=:), allocatable :: type_string

    if (allocated(pc)) then
       call precon_destroy(pc)
       deallocate(pc)
    end if

    if (trim(type_name) .eq. 'jacobi') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_jacobi_t::pc)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(device_jacobi_t::pc)
       else
          allocate(jacobi_t::pc)
       end if
    else if (type_name(1:4) .eq. 'hsmg') then
       allocate(hsmg_t::pc)
    else if(trim(type_name) .eq. 'ident') then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(device_ident_t::pc)
       else
          allocate(ident_t::pc)
       end if
    else
       type_string =  concat_string_array(PC_KNOWN_TYPES, &
            NEW_LINE('A') // "-  ",  .true.)
       call neko_error("Unknown preconditioner type: " &
                       // trim(type_name) // ".  Known types are: " &
                       // type_string)
    end if

  end subroutine precon_factory

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
       end select
    end if

  end subroutine precon_destroy

end submodule precon_fctry
