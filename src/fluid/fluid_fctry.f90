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
!> Factory for all fluid schemes
module fluid_fctry
  use fluid_method, only : fluid_scheme_t
  use fluid_plan1, only : fluid_plan1_t
  use fluid_pnpn, only : fluid_pnpn_t    
  use utils, only : neko_error
  use neko_config
  implicit none

contains

  !> Initialise a fluid scheme
  subroutine fluid_scheme_factory(fluid, fluid_scheme)
    class(fluid_scheme_t), intent(inout), allocatable :: fluid
    character(len=*) :: fluid_scheme

    if (trim(fluid_scheme) .eq. 'plan1') then
       allocate(fluid_plan1_t::fluid)
    else if (trim(fluid_scheme) .eq. 'pnpn') then
       allocate(fluid_pnpn_t::fluid)
    else
       call neko_error('Invalid fluid scheme')
    end if
    
  end subroutine fluid_scheme_factory

end module fluid_fctry
