! Copyright (c) 2021-2025, The Neko Authors
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
submodule (les_model) les_model_fctry
  use vreman, only : vreman_t
  use smagorinsky, only : smagorinsky_t
  use dynamic_smagorinsky, only : dynamic_smagorinsky_t
  use sigma, only : sigma_t
  use wale, only : wale_t
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: LES_KNOWN_TYPES(5) = [character(len=20) :: &
       "vreman", &
       "smagorinsky", &
       "dymamic_smagorinsky", &
       "sigma", &
       "wale"]

contains
  !> LES model factory. Both constructs and initializes the object.
  !! @param object The object to be allocated.
  !! @param type_name The name of the LES model.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  module subroutine les_model_factory(object, type_name, dofmap, coef, json)
    class(les_model_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: type_string

    if (allocated(object)) deallocate(object)

    select case (trim(type_name))
    case ('vreman')
       allocate(vreman_t::object)
    case ('smagorinsky')
       allocate(smagorinsky_t::object)
    case ('dynamic_smagorinsky')
       allocate(dynamic_smagorinsky_t::object)
    case ('sigma')
       allocate(sigma_t::object)
    case ('wale')
       allocate(wale_t::object)
    case default
       call neko_type_error("LES model", type_name, LES_KNOWN_TYPES)
    end select

  end subroutine les_model_factory

end submodule les_model_fctry
