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
submodule (les_model) les_model_fctry
  use vreman, only : vreman_t
  use smagorinsky, only : smagorinsky_t
  use dynamic_smagorinsky, only : dynamic_smagorinsky_t
  use sigma, only : sigma_t
  use utils, only : concat_string_array, neko_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: LES_KNOWN_TYPES(4) = [character(len=20) :: &
     "vreman", &
     "smagorinsky", &
     "dymamic_smagorinsky", &
     "sigma"]

contains
  !> LES model factory. Both constructs and initializes the object.
  !! @param object The object to be allocated.
  !! @param type_name The name of the LES model.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  module subroutine les_model_factory(object, type_name, dofmap, coef, json)
    class(les_model_t), allocatable, target, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: type_string

    if (allocated(object)) then
       deallocate(object)
    else if (trim(type_name) .eq. 'vreman') then
       allocate(vreman_t::object)
    else if (trim(type_name) .eq. 'smagorinsky') then
       allocate(smagorinsky_t::object)
    else if (trim(type_name) .eq. 'dynamic_smagorinsky') then
       allocate(dynamic_smagorinsky_t::object)
    else if (trim(type_name) .eq. 'sigma') then
       allocate(sigma_t::object)
    else
       type_string =  concat_string_array(LES_KNOWN_TYPES, &
            NEW_LINE('A') // "-  ", .true.)
       call neko_error("Unknown LES model type: " &
                       // trim(type_name) // ".  Known types are: " &
                       // type_string)
       stop

    end if
    call object%init(dofmap, coef, json)
  end subroutine les_model_factory

end submodule les_model_fctry
