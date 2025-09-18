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
submodule (filter) filter_fctry
  use elementwise_filter, only : elementwise_filter_t
  use PDE_filter, only : PDE_filter_t
  use utils, only : concat_string_array, neko_type_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: FILTER_KNOWN_TYPES(2) = [character(len=20) :: &
       "elementwise", &
       "PDE"]

contains
  !> Filter factory. Both constructs and initializes the object.
  !! @param object The object to be allocated.
  !! @param type_name The name of the filter.
  !! @param json A dictionary with parameters.
  !! @param coef SEM coefficients.
  module subroutine filter_factory(object, type_name, json, coef)
    class(filter_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: type_string

    if (allocated(object)) then
       deallocate(object)
    else if (trim(type_name) .eq. 'elementwise') then
       allocate(elementwise_filter_t::object)
    else if (trim(type_name) .eq. 'PDE') then
       allocate(pde_filter_t::object)
    else
       call neko_type_error("filter type: ", type_name, &
            FILTER_KNOWN_TYPES)
    end if

    call object%init(json, coef)
  end subroutine filter_factory

end submodule filter_fctry
