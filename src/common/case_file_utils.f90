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
  use registry, only : neko_const_registry
  use vector, only : vector_t
  implicit none

contains

!> Retrieves a real either from the json or from the corresponding scalar
!! in the `neko_const_registry`.
!! @details First tries to retrieve the value from the JSON, looking for a real
!! entry under the given name. If not found, it looks for the same entry, but
!! with a string value. The retrieved string is the name looked up in the
!! `neko_const_registry` as a scalar, and the data is copied from there.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_or_lookup_real(json, name, val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=sp), intent(out) :: val

    character(len=:), allocatable :: reg_name
    logical :: found
    integer :: var_type

    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 6) .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither a real nor a string.")
    end if

    ! Try to find a real
    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Try to find a string. It must exist
    call json_get(json, name, reg_name)

    ! Retrieve the value from the registry
    val = real(neko_const_registry%get_real_scalar(reg_name), kind=sp)
  end subroutine json_get_or_lookup_real

!> Same as `json_get_or_lookup_real`, but for double precision.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_or_lookup_double(json, name, val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: val

    character(len=:), allocatable :: reg_name
    logical :: found
    integer :: var_type

    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 6) .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither a real nor a string.")
    end if

    ! Try to find a real
    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Try to find a string. It must exist
    call json_get(json, name, reg_name)

    call json%print()
    ! Retrieve the value from the registry
    val = real(neko_const_registry%get_real_scalar(reg_name), kind=dp)
  end subroutine json_get_or_lookup_double

!> Retrieves an integer either from the json or from the corresponding scalar
!! in the `neko_const_registry`.
!! @details First tries to retrieve the value from the JSON, looking for an int
!! entry under the given name. If not found, it looks for the same entry, but
!! with a string value. The retrieved string is the name looked up in the
!! `neko_const_registry` as a scalar, and the data is copied from there,
!! using explicit conversion to an int.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_or_lookup_integer(json, name, val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    integer, intent(out) :: val

    character(len=:), allocatable :: reg_name
    logical :: found
    integer :: var_type

    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 5) .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither an integer nor a string.")
    end if

    ! Try to find an int
    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Try to find a string. It must exist
    call json_get(json, name, reg_name)

    ! Retrieve the value from the registry
    val = neko_const_registry%get_integer_scalar(reg_name)
  end subroutine json_get_or_lookup_integer

!> Retrieves a real either from the json or from the corresponding scalar
!! in the `neko_registry`, otherwise sets the default value.
!! @details First tries to retrieve the value from the JSON, looking for a real
!! entry under the given name. If not found, it looks for the same entry, but
!! with a string value. If this fails, the default value is assigned. Otherwise,
!! the retrieved string is the name looked up in the `neko_const_registry` as a
!! scalar, and the data is copied from there.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
!! @param[in] defalt The default value to be used if the parameter is not found
  module subroutine json_get_or_lookup_or_default_real(json, name,&
       val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=sp), intent(out) :: val
    real(kind=sp), intent(in) :: default

    character(len=:), allocatable :: reg_name
    logical :: found
    integer :: var_type

    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 6) .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither a real nor a string.")
    end if

    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Try to find a string.
    call json%get(name, reg_name, found)
    if (.not. found) then
       ! We use another call here instead of just assigning the defaut value
       ! to handle whether defaults are allowed as per `json_no_defaults`.
       call json_get_or_default(json, name, val, default)
    else
       ! Retrieve the array from the registry
       val = real(neko_const_registry%get_real_scalar(reg_name), kind=sp)
    end if
  end subroutine json_get_or_lookup_or_default_real

!> Same as `json_get_or_lookup_or_default_real`, but for doubles.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
!! @param[in] defalt The default value to be used if the parameter is not found
  module subroutine json_get_or_lookup_or_default_double(json, name, &
       val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: val
    real(kind=dp), intent(in) :: default

    character(len=:), allocatable :: reg_name
    logical :: found
    integer :: var_type

    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 6) .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither a real nor a string.")
    end if

    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Try to find a string.
    call json%get(name, reg_name, found)
    if (.not. found) then
       ! We use another call here instead of just assigning the defaut value
       ! to handle whether defaults are allowed as per `json_no_defaults`.
       call json_get_or_default(json, name, val, default)
    else
       ! Retrieve the array from the registry
       val = real(neko_const_registry%get_real_scalar(reg_name), kind=dp)
    end if
  end subroutine json_get_or_lookup_or_default_double

!> Retrieves an integer either from the json or from the corresponding scalar
!! in the `neko_registry`, otherwise sets the default value.
!! @details First tries to retrieve the value from the JSON, looking for an int
!! entry under the given name. If not found, it looks for the same entry, but
!! with a string value. If this fails, the default value is assigned. Otherwise,
!! the retrieved string is the name looked up in the `neko_const_registry` as a
!! scalar, and the data is copied from there, using explicit convresion to an
!! int.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
!! @param[in] defalt The default value to be used if the parameter is not found
  module subroutine json_get_or_lookup_or_default_integer(json, &
       name, val, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer, intent(in) :: default

    character(len=:), allocatable :: reg_name
    logical :: found
    integer :: var_type

    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 5) .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither an integer nor a string.")
    end if

    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Try to find a string.
    call json%get(name, reg_name, found)
    if (.not. found) then
       ! We use another call here instead of just assigning the defaut value
       ! to handle whether defaults are allowed as per `json_no_defaults`.
       call json_get_or_default(json, name, val, default)
    else
       ! Retrieve the array from the registry
       val = neko_const_registry%get_integer_scalar(reg_name)
    end if
  end subroutine json_get_or_lookup_or_default_integer

!> Retrieves a real array either from the json or from the corresponding vector
!! in the `neko_registry`.
!! @details First tries to retrieve the values from the JSON, looking for a real
!! array entry under the given name. If not found, it looks for the same entry,
!! but with a string value. The retrieved string is the name looked up in the
!! `neko_const_registry` as a scalar, and the data is copied from there.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_or_lookup_real_array(json, name, val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=sp), allocatable, intent(inout) :: val(:)

    type(vector_t), pointer :: vec_ptr
    logical :: found
    character(len=:), allocatable :: reg_name
    integer :: var_type

    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 3) .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither an array nor a string.")
    end if

    ! Try to find a real array
    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Finding an array failed. Check if it is a string, otherwise error.
    ! If found is false here, json_get will emit the correct error.
    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither a real array nor a string.")
    end if

    ! Try to find a string. It must exist
    call json_get(json, name, reg_name)

    ! Retrieve the array from the registry
    vec_ptr => neko_const_registry%get_vector(reg_name)
    val = real(vec_ptr%x, kind=sp)
  end subroutine json_get_or_lookup_real_array

!> Sampe as `json_get_or_lookup_real_array`, but for double precision.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_or_lookup_double_array(json, name, val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=dp), allocatable, intent(inout) :: val(:)

    type(vector_t), pointer :: vec_ptr
    logical :: found
    character(len=:), allocatable :: reg_name
    integer :: var_type

    ! Try to find a real array
    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Finding an array failed. Check if it is a string, otherwise error.
    ! If found is false here, json_get will emit the correct error.
    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither a real array nor a string.")
    end if

    ! Try to find a string. It must exist
    call json_get(json, name, reg_name)

    ! Retrieve the array from the registry
    vec_ptr => neko_const_registry%get_vector(reg_name)
    val = real(vec_ptr%x, kind=dp)
  end subroutine json_get_or_lookup_double_array

!> Retrieves an int array either from the json or from the corresponding vector
!! in the `neko_registry`.
!! @details First tries to retrieve the values from the JSON, looking for an int
!! array entry under the given name. If not found, it looks for the same entry,
!! but with a string value. The retrieved string is the name looked up in the
!! `neko_const_registry` as a scalar, and the data is copied from there.
!! @param[inout] json The json to retrieve the parameter from.
!! @param[in] name The full path to the parameter.
!! @param[out] value The variable to be populated with the retrieved parameter
  module subroutine json_get_or_lookup_integer_array(json, name, &
       val)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    integer, allocatable, intent(inout) :: val(:)

    type(vector_t), pointer :: vec_ptr
    logical :: found
    character(len=:), allocatable :: reg_name
    integer :: i
    integer :: var_type

    ! Try to find an integer array
    call json%get(name, val, found)
    ! The value is retrieved from the JSON keyword, we are done
    if (found) return

    ! Finding an array failed. Check if it is a string, otherwise error.
    ! If found is false here, json_get will emit the correct error.
    call json%info(name, found = found, var_type = var_type)
    if (found .and. (var_type .ne. 7)) then
       call neko_error("Parameter " // name // &
            " is neither an integer array nor a string.")
    end if

    ! Try to find a string. It must exist
    call json_get(json, name, reg_name)

    ! Retrieve the array from the registry
    vec_ptr => neko_const_registry%get_vector(reg_name)
    val = int(vec_ptr%x)

    do i = 1, size(val)
       if (.not. abscmp(real(val(i), kind=rp), vec_ptr%x(i))) then
          call neko_error("json_get_or_lookup_integer_array: " &
               // "Value retrieved from registry entry '" // reg_name &
               // "' is not an integer array.")
       end if
    end do
  end subroutine json_get_or_lookup_integer_array

end submodule case_file_utils
