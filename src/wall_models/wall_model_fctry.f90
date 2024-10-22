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
submodule (wall_model) wall_model_fctry
  use vreman, only : vreman_t
  use spalding, only : spalding_t
  use rough_log_law, only : rough_log_law_t
  use utils, only : concat_string_array
  use json_utils, only : json_get
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: WALLM_KNOWN_TYPES(2) = [character(len=20) :: &
     "spalding", &
     "rough_log_law"]

contains

  !> Wall model factory. Both constructs and initializes the object.
  !! @param object The object to be allocated.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param nu The molecular kinematic viscosity.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  module subroutine wall_model_factory(object, coef, msk, facet, nu, &
       h_index, json)
    class(wall_model_t), allocatable, target, intent(inout) :: object
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    real(kind=rp), intent(in) :: nu
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string

    type_string =  concat_string_array(WALLM_KNOWN_TYPES, &
         NEW_LINE('A') // "-  ", prepend=.true.)

    call json_get(json, "model", type_name)

    if (trim(type_name) .eq. "spalding") then
       allocate(spalding_t::object)
    else if (trim(type_name) .eq. "rough_log_law") then
       allocate(rough_log_law_t::object)
    else
       call neko_error("Unknown wall model type: " // trim(type_name) // &
          trim(type_name) // ".  Known types are: "  // type_string)
       stop
    end if

    ! Initialize
    call object%init(coef, msk, facet, nu, h_index, json)

  end subroutine wall_model_factory

end submodule wall_model_fctry
