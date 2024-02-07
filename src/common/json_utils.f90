! Copyright (c) 2019-2021, The Neko Authors
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
!> Utilities for retrieving parameters from the case files.
module json_utils
  use num_types, only : rp
  use json_module, only : json_file, json_value, json_core
  use utils, only : neko_error
  implicit none
  private

  public :: json_get, json_get_or_default, json_extract_subdict

  !> Retrieves a parameter by name or throws an error
  interface json_get
     module procedure json_get_real, json_get_integer, json_get_logical, &
                     json_get_string, json_get_real_array, &
                     json_get_integer_array, json_get_logical_array, &
                     json_get_string_array
  end interface json_get

  !> Retrieves a parameter by name or assigns a provided default value.
  !! In the latter case also adds the missing paramter to the json
  interface json_get_or_default
     module procedure json_get_or_default_real, json_get_or_default_integer,&
                     json_get_or_default_string, json_get_or_default_logical
  end interface json_get_or_default

contains

  !> Retrieves a real parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_real(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=rp), intent(out) :: value
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
  end subroutine json_get_real

  !> Retrieves an integer parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_integer(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
  end subroutine json_get_integer

  !> Retrieves a logical parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_logical(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    logical, intent(out) :: value
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
  end subroutine json_get_logical

  !> Retrieves a string parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_string(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: value
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
  end subroutine json_get_string

  !> Retrieves a real array parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_real_array(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=rp), allocatable, intent(out) :: value(:)
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
  end subroutine json_get_real_array

  !> Retrieves a integer array parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_integer_array(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    integer, allocatable, intent(out) :: value(:)
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
  end subroutine json_get_integer_array

  !> Retrieves a logical array parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_logical_array(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    logical, allocatable, intent(out) :: value(:)
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
  end subroutine json_get_logical_array

  !> Retrieves a string array parameter by name or throws an error
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_string_array(json, name, value)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    character(len=*), allocatable, intent(out) :: value(:)
    logical :: found
    type(json_value), pointer :: json_val, val_ptr
    type(json_core) :: core
    character(len=:), allocatable :: string_value
    integer :: i, n_children

    if (.not. json%valid_path(name)) then
       call neko_error("Parameter "//name//" missing from the case file")
    end if
    call json%info(name, n_children=n_children)

    if (.not. allocated(value)) then
       allocate(value(n_children))
    else if (len(value) .lt. n_children) then
       deallocate(value)
       allocate(value(n_children))
    end if

    call json%get(name, json_val, found)
    call json%get_core(core)

    do i = 1, n_children
       call core%get_child(json_val, i, val_ptr, found)
       call core%get(val_ptr, string_value)

       if (len(string_value) .gt. 0) then
          value(i) = string_value
       end if
    end do

  end subroutine json_get_string_array

  !> Retrieves a real parameter by name or assigns a provided default value.
  !! In the latter case also adds the missing paramter to the json.
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_or_default_real(json, name, value, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    real(kind=rp), intent(out) :: value
    real(kind=rp), intent(in) :: default
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       value = default
       call json%add(name, value)
    end if
  end subroutine json_get_or_default_real

  !> Retrieves an integer parameter by name or assigns a provided default value.
  !! In the latter case also adds the missing paramter to the json.
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_or_default_integer(json, name, value, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    integer, intent(in) :: default
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       value = default
       call json%add(name, value)
    end if
  end subroutine json_get_or_default_integer

  !> Retrieves a logical parameter by name or assigns a provided default value.
  !! In the latter case also adds the missing paramter to the json.
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_or_default_logical(json, name, value, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    logical, intent(out) :: value
    logical, intent(in) :: default
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       value = default
       call json%add(name, value)
    end if
  end subroutine json_get_or_default_logical

  !> Retrieves a string parameter by name or assigns a provided default value.
  !! In the latter case also adds the missing paramter to the json.
  !! @param json The json to retrieve the parameter from.
  !! @param name The full path to the parameter.
  !! @value value The variable to be populated with the retrieved parameter.
  subroutine json_get_or_default_string(json, name, value, default)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: value
    character(len=*), intent(in) :: default
    logical :: found

    call json%get(name, value, found)

    if (.not. found) then
       value = default
       call json%add(name, value)
    end if
  end subroutine json_get_or_default_string

  !> Extract `i`th json object from a json array.
  !! @param core JSON core object.
  !! @param parent The parent JSON object with the array.
  !! @param i The index of the object to extract.
  !! @param subdict JSON object object to be filled with the subdict.
  subroutine json_extract_subdict(core, parent, i, subdict)
    type(json_core), intent(inout) :: core
    type(json_value), pointer, intent(in) :: parent
    integer, intent(in) :: i
    type(json_file), intent(inout) :: subdict
    type(json_value), pointer :: ptr
    logical :: found
    character(len=:), allocatable :: buffer

    call core%get_child(parent, i, ptr, found)
    call core%print_to_string(ptr, buffer)
    call subdict%load_from_string(buffer)

  end subroutine json_extract_subdict

end module json_utils
