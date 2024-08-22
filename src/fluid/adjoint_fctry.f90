! Copyright (c) 2021, The Neko Authors
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
!> Factory for all adjoint schemes
module adjoint_fctry
  use adjoint_scheme, only : adjoint_scheme_t
  use adjoint_pnpn, only: adjoint_pnpn_t
  use utils, only : concat_string_array, neko_error
  implicit none
  private

  public :: adjoint_scheme_factory

  ! List of all possible types created by the factory routine
  character(len=20) :: KNOWN_TYPES(1) = [character(len=20) :: &
       "pnpn"]

contains

  !> Initialise a adjoint scheme
  subroutine adjoint_scheme_factory(object, type_name)
    class(adjoint_scheme_t), intent(inout), allocatable :: object
    character(len=*) :: type_name
    character(len=:), allocatable :: type_string

    if (trim(type_name) .eq. 'pnpn') then
       allocate(adjoint_pnpn_t::object)
    else
       type_string = concat_string_array(KNOWN_TYPES, NEW_LINE('A') // "-  ", &
            .true.)
       call neko_error("Unknown adjoint scheme type: " &
            // trim(type_name) // ". Known types are: " &
            // type_string)
    end if

  end subroutine adjoint_scheme_factory

end module adjoint_fctry
