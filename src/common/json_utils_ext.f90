!> @file json_utils_ext.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
module json_utils_ext
  use mpi_f08, only: MPI_Comm_rank, MPI_Initialized, MPI_Bcast, &
       MPI_INTEGER, MPI_CHARACTER
  use comm, only: NEKO_COMM
  use json_file_module, only: json_file
  use json_value_module, only: json_value
  use utils, only: neko_error, filename_suffix

  implicit none
  private

  public :: json_key_fallback, json_read_file

  interface json_key_fallback
     module procedure json_key_fallback_string
     module procedure json_key_fallback_json
  end interface json_key_fallback

contains

  !> Create a json_string based on fallback logic.
  !! If the lookup key is present in the json object, return it.
  !! If the fallback key is present in the json object, return it.
  !! Otherwise, return the lookup key.
  function json_key_fallback_string(json, lookup, fallback) result(string)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: lookup
    character(len=*), intent(in) :: fallback
    character(len=:), allocatable :: string

    if ((lookup .in. json)) then
       string = lookup
    else if (fallback .in. json) then
       string = fallback
    else
       string = lookup
    end if

  end function json_key_fallback_string

  !> Create a json object based on fallback logic.
  !! If the key is present in the lookup json object, point to this object.
  !! If the key is present in the fallback json object, point to this object.
  !! Otherwise, point to the lookup object.
  function json_key_fallback_json(lookup_json, fallback_json, key) &
       result(json_pointer)
    type(json_file), target, intent(inout) :: lookup_json
    type(json_file), target, intent(inout) :: fallback_json
    character(len=*), intent(in) :: key
    type(json_file), pointer :: json_pointer

    if ((key .in. lookup_json)) then
       json_pointer => lookup_json
    else if (key .in. fallback_json) then
       json_pointer => fallback_json
    else
       json_pointer => lookup_json
    end if

  end function json_key_fallback_json

  !> Read a json file taking mpi into account.
  !! This function reads a json file and broadcasts it to all ranks.
  !! @param[in] filename The name of the file to read.
  !! @return The json object.
  function json_read_file(filename) result(json)
    character(len=*), intent(in) :: filename
    type(json_file) :: json

    logical :: mpi_is_initialized
    integer :: rank, ierr, length
    character(len=:), allocatable :: json_buffer
    character(len=4) :: suffix

    ! Check if the file is a json or Neko case file.
    call filename_suffix(filename, suffix)

    if (trim(suffix) .ne. 'json' .and. trim(suffix) .ne. 'case' ) then
       call neko_error('Invalid case file')
    end if

    ! Initialize the mpi variables.
    rank = 0
    mpi_is_initialized = .false.

    ! Check if MPI is initialized and get the rank if it is.
    call MPI_Initialized(mpi_is_initialized, ierr)
    if (mpi_is_initialized) call MPI_Comm_rank(NEKO_COMM, rank, ierr)

    ! Read the json file and broadcast it to all ranks.
    if (rank .eq. 0) then
       call json%load_file(filename = trim(filename))
    end if

    ! Serialize the json object to a string so it can be broadcast.
    if (mpi_is_initialized) then
       if (rank .eq. 0) call json%print_to_string(json_buffer)

       length = len(json_buffer)
       call MPI_Bcast(length, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

       if (rank .ne. 0) allocate(character(len=length) :: json_buffer)
       call MPI_Bcast(json_buffer, length, MPI_CHARACTER, 0, NEKO_COMM, ierr)

       if (rank .ne. 0) then
          call json%load_from_string(json_buffer)
          deallocate(json_buffer)
       end if
    end if

  end function json_read_file

end module json_utils_ext
