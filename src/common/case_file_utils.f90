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
!> Utilities for handling the case file, which are not pure JSON. Therefore,
!! split from the JSON utilities parent module.
submodule (json_utils) case_file_utils
  use registry, only : neko_registry
  use vector, only : vector_t
  implicit none

contains

!> Retrieves a real either from the json or from the corresponding scalar
!! in the `neko_registry`.
!! @details First tries to retrieve the array from the JSON, looking for a real
!! entry under the given name. If not found, the name is looked up in the
!! `neko_registry` as a scalar, and the data is copied from there.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_real_from_registry_or_entry(json, name, val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=rp), intent(out) :: val
    real(kind=rp), pointer :: scalar_ptr

    logical :: found

    call json%get(name, val, found)
    if (found) return

    scalar_ptr => neko_registry%get_scalar(name)
    val = scalar_ptr
  end subroutine json_get_real_from_registry_or_entry

!> Retrieves a real either from the json or from the corresponding scalar
!! in the `neko_registry`, otherwise sets the default value.
!! @details First tries to retrieve the array from the JSON, looking for a real
!! entry unde the given name. If not found, the name is looked up in the
!! `neko_registry` as a scalar, and the data is copied from there.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
!! @param[in] defalt The default value to be used if the parameter is not found
  module subroutine json_get_real_from_registry_or_entry_or_default(json, name,&
       val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=rp), intent(out) :: val
    real(kind=rp), intent(in) :: default

    real(kind=rp), pointer :: scalar_ptr
    logical :: found

    call json%get(name, val, found)
    if (found) return

    found = neko_registry%scalar_exists(name)
    if (found) then
       scalar_ptr => neko_registry%get_scalar(name)
       val = scalar_ptr
    else
       val = default
    end if
  end subroutine json_get_real_from_registry_or_entry_or_default

!> Retrieves a real array either from the json or from the corresponding vector
!! in the `neko_registry`.
!! @details First tries to retrieve the array from the JSON, looking for a real
!! vector entry under the given name. If not found, the name is looked up in the
!! `neko_registry` as a vector, and the data is copied from there.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_real_array_from_registry_or_entry(json, name, val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=rp), allocatable, intent(out) :: val(:)
    type(vector_t), pointer :: vec_ptr

    logical :: found

    call json%get(name, val, found)
    if (found) return

    vec_ptr => neko_registry%get_vector(name)
    val = vec_ptr%x
  end subroutine json_get_real_array_from_registry_or_entry

end submodule case_file_utils
