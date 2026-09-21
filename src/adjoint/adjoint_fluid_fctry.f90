!> @file adjoint_fluid_fctry.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Factory for all adjoint fluid schemes
module adjoint_fluid_fctry
  use adjoint_fluid_scheme, only: adjoint_fluid_scheme_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use utils, only: neko_type_error
  implicit none
  private

  public :: adjoint_fluid_scheme_factory

  ! List of all possible types created by the factory routine
  character(len=20) :: KNOWN_TYPES(1) = [character(len=20) :: &
       "pnpn"]

contains

  !> Initialise a adjoint fluid scheme
  subroutine adjoint_fluid_scheme_factory(object, type_name)
    class(adjoint_fluid_scheme_t), intent(inout), allocatable :: object
    character(len=*), intent(in) :: type_name

    select case (trim(type_name))
    case ('pnpn')
       allocate(adjoint_fluid_pnpn_t::object)
    case default
       call neko_type_error("adjoint fluid scheme", type_name, KNOWN_TYPES)
    end select

  end subroutine adjoint_fluid_scheme_factory

end module adjoint_fluid_fctry
