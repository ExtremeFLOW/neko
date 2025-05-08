! Copyright (c) 2025, The Neko Authors
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
!> Factory for all fluid schemes
submodule (fluid_scheme_base) fluid_base_fctry
  use fluid_pnpn, only : fluid_pnpn_t
  use fluid_scheme_compressible_euler, only : fluid_scheme_compressible_euler_t
  use utils, only : neko_type_error

  ! List of all possible types created by the factory routine
  character(len=20) :: FLUID_KNOWN_TYPES(2) = [character(len=20) :: &
       "pnpn", "compressible"]

contains

  !> Initialise a fluid scheme
  module subroutine fluid_scheme_base_factory(object, type_name)
    class(fluid_scheme_base_t), intent(inout), allocatable :: object
    character(len=*) :: type_name
    character(len=:), allocatable :: type_string

    select case (trim(type_name))
    case ('pnpn')
       allocate(fluid_pnpn_t::object)
    case ('compressible')
       allocate(fluid_scheme_compressible_euler_t::object)
    case default
       call neko_type_error("fluid scheme base", type_name, FLUID_KNOWN_TYPES)
    end select

  end subroutine fluid_scheme_base_factory

end submodule fluid_base_fctry
