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
!> Utilities
!! @details Various utility functions
module utils
  use, intrinsic :: iso_fortran_env, only : error_unit, output_unit
  use iso_c_binding
  use num_types, only: rp, dp
  implicit none
  private

  integer, parameter :: NEKO_FNAME_LEN = 1024
  integer, parameter :: NEKO_VARNAME_LEN = 256

  interface neko_error
     module procedure neko_error_plain, neko_error_msg
  end interface neko_error

  interface read_duration
     module procedure read_duration_scalar
     module procedure read_duration_components
  end interface read_duration

  public :: neko_error, neko_warning, nonlinear_index, filename_chsuffix, &
       filename_path, filename_name, filename_suffix, &
       filename_suffix_pos, filename_tslash_pos, filename_split, &
       linear_index, split_string, NEKO_FNAME_LEN, index_is_on_facet, &
       concat_string_array, extract_fld_file_index, neko_type_error, &
       neko_type_registration_error, NEKO_VARNAME_LEN, mkdir, read_duration

  interface
     function c_mkdir(path, mode) bind(C, name="mkdir")
       import :: c_char, c_int
       character(kind=c_char), dimension(*) :: path
       integer(c_int), value :: mode
       integer(c_int) :: c_mkdir
     end function c_mkdir
  end interface

contains

  !> Find position (in the string) of a filename's suffix
  pure function filename_suffix_pos(fname) result(suffix_pos)
    character(len=*), intent(in) :: fname
    integer :: suffix_pos
    suffix_pos = scan(trim(fname), '.', back = .true.)
  end function filename_suffix_pos

  !> Find position (in the string) of a filename's trailing slash
  pure function filename_tslash_pos(fname) result(tslash_pos)
    character(len=*), intent(in) :: fname
    integer :: tslash_pos
    tslash_pos = scan(trim(fname), '/', back = .true.)
  end function filename_tslash_pos

  !> Extract the path to a file
  subroutine filename_path(fname, path)
    character(len=*), intent(in) :: fname
    character(len=*), intent(out) :: path
    integer :: tslash_pos

    tslash_pos = filename_tslash_pos(fname)
    if (tslash_pos .gt. 0) then
       path = trim(fname(1:tslash_pos))
    else
       path = './'
    end if

  end subroutine filename_path

  !> Extract the base name of a file (without path and suffix)
  subroutine filename_name(fname, name)
    character(len=*), intent(in) :: fname
    character(len=*), intent(out) :: name
    integer :: tslash_pos, suffix_pos, start, end

    tslash_pos = filename_tslash_pos(fname)
    suffix_pos = filename_suffix_pos(fname)
    if (tslash_pos .eq. 0) then
       start = 1
    else
       start = tslash_pos + 1
    end if
    if (suffix_pos .eq. 0) then
       end = len_trim(fname)
    else
       end = suffix_pos - 1
    end if
    name = trim(fname(start:end))
  end subroutine filename_name

  !> Extract a filename's suffix
  subroutine filename_suffix(fname, suffix)
    character(len=*) :: fname
    character(len=*) :: suffix
    suffix = trim(fname(filename_suffix_pos(fname) + 1:len_trim(fname)))
  end subroutine filename_suffix

  !> Extract file name components
  subroutine filename_split(fname, path, name, suffix)
    character(len=*), intent(in) :: fname
    character(len=*), intent(out) :: path, name, suffix
    integer :: tslash_pos, suffix_pos

    tslash_pos = filename_tslash_pos(fname)
    suffix_pos = filename_suffix_pos(fname)

    if (tslash_pos .gt. 0) then
       path = trim(fname(1:tslash_pos))
    else
       path = './'
    end if

    if (suffix_pos .gt. 0) then
       name = trim(fname(tslash_pos + 1:suffix_pos - 1))
       suffix = trim(fname(suffix_pos:len_trim(fname)))
    else
       name = trim(fname(tslash_pos + 1:len_trim(fname)))
       suffix = ''
    end if

  end subroutine filename_split

  !> Change a filename's suffix
  subroutine filename_chsuffix(fname, new_fname, new_suffix)
    character(len=*) :: fname
    character(len=*) :: new_fname
    character(len=*) :: new_suffix
    integer :: suffix_pos

    suffix_pos = filename_suffix_pos(fname)
    new_fname = trim(fname(1:suffix_pos))//new_suffix

  end subroutine filename_chsuffix

  !> Recursively create a directory and all parent directories if they do not
  !! exist. This should be safer than the
  !! `execute_command_line('mkdir -p ' // path)` approach, which can be used to
  !! execute arbitrary code if `path` is not properly sanitized.
  !! @param path The path of the directory to create.
  !! @param mode The octal mode to create the directory with, will adjust by
  !! umask.
  recursive subroutine mkdir(path, mode)
    character(len=*), intent(in) :: path
    integer, intent(in), optional :: mode
    integer :: slash_pos, i, path_len
    character(kind=c_char), allocatable :: c_path(:)
    integer(c_int) :: dir_mode, ierr

    if (present(mode)) then
       dir_mode = int(mode, kind=c_int)
    else
       dir_mode = int(o'777', kind=c_int)
    end if

    slash_pos = scan(path, '/', back = .true.)
    if (slash_pos .gt. 0) then
       call mkdir(trim(path(1:slash_pos-1)), dir_mode)
    end if

    path_len = len_trim(path)
    allocate(c_path(path_len + 1))
    do i = 1, path_len
       c_path(i) = path(i:i)
    end do
    c_path(path_len + 1) = c_null_char

    ierr = c_mkdir(c_path, dir_mode)
    deallocate(c_path)
  end subroutine mkdir

  !> Extracts the index of a field file. For example, "myfield.f00045"
  !! will return `45`. If the suffix of the file name is invalid, returns
  !! a default index value.
  !! @param fld_filename Name of the fld file, e.g. `myfield0.f00035`.
  !! @param default_index The index to return in case the suffix of
  !! `fld_filename` is invalid.
  function extract_fld_file_index(fld_filename, default_index) result(index)
    character(len=*), intent(in) :: fld_filename
    integer, intent(in) :: default_index

    character(len=80) :: suffix
    integer :: index, fpos, i
    logical :: valid

    call filename_suffix(fld_filename, suffix)

    valid = .true.

    ! This value will be modified when reading the file name extension
    ! e.g. "field0.f00035" will set sample_idx = 35
    index = default_index

    !
    ! Try to extract the index of the field file from the suffix "fxxxxx"
    !
    fpos = scan(trim(suffix), 'f')
    if (fpos .eq. 1) then
       ! Make sure that the suffix only contains integers from 0 to 9
       do i = 2, len(trim(suffix))
          if (.not. (iachar(suffix(i:i)) >= iachar('0') &
               .and. iachar(suffix(i:i)) <= iachar('9'))) then
             valid = .false.
          end if
       end do
    else
       valid = .false.
    end if

    ! Must be exactly 6 characters long, i.e. an 'f' with 5 integers after
    if (len(trim(suffix)) .ne. 6) valid = .false.

    if (valid) read (suffix(2:), "(I5.5)") index

  end function extract_fld_file_index

  !> Split a string based on delimiter (tokenizer)
  !! OBS: very hacky, this should really be improved, it is rather embarrasing
  !! code.
  function split_string(string, delimiter) result(split_str)
    character(len=*) :: string
    character(len=*) :: delimiter
    character(len=100), allocatable :: split_str(:)
    integer :: length, i, i2, offset, j
    i = 0
    offset = 1
    length = 1
    if (len(trim(string)) .eq. 0) then
       allocate(split_str(1))
       split_str(1) = trim(string)
       return
    end if
    do while (.true.)
       i = scan(string(offset:), delimiter, back = .false.)
       if (i .eq. 0) exit
       length = length + 1
       offset = offset + i
    end do

    allocate(split_str(length))
    i = 0
    j = 1
    offset = 1
    do while (.true.)
       i2 = scan(trim(string(offset:)), delimiter, back = .false.)
       if (i2 .eq. 0) then
          split_str(j) = trim(string(offset:))
          exit
       end if
       split_str(j) = trim(string(offset:offset+i2-2))
       offset = offset+i2
       j = j + 1
    end do
  end function split_string

  !> Compute the address of a (i,j,k,l) array
  !! with sizes (1:lx, 1:ly, 1:lz, :)
  pure function linear_index(i, j, k, l, lx, ly, lz) result(index)
    integer, intent(in) :: i, j, k, l, lx, ly, lz
    integer :: index

    index = (i + lx * ((j - 1) + ly * ((k - 1) + lz * ((l - 1)))))
  end function linear_index

  !> Compute (i,j,k,l) array given linear index
  !! with sizes (1:lx, 1:ly, 1:lz, :)
  pure function nonlinear_index(linear_index, lx, ly, lz) result(index)
    integer, intent(in) :: linear_index, lx, ly, lz
    integer :: index(4)
    integer :: lin_idx
    lin_idx = linear_index -1
    index(4) = lin_idx/(lx*ly*lz)
    index(3) = (lin_idx-(lx*ly*lz)*index(4))/(lx*ly)
    index(2) = (lin_idx-(lx*ly*lz)*index(4)-(lx*ly)*index(3))/lx
    index(1) = (lin_idx-(lx*ly*lz)*index(4)-(lx*ly)*index(3)-lx*index(2))
    index(1) = index(1) + 1
    index(2) = index(2) + 1
    index(3) = index(3) + 1
    index(4) = index(4) + 1

  end function nonlinear_index

  pure function index_is_on_facet(i, j, k, lx, ly, lz, facet) result(is_on)
    integer, intent(in) :: i, j, k, lx, ly, lz, facet
    logical :: is_on

    is_on = .false.
    select case (facet)
    case (1)
       if (i .eq. 1) is_on = .true.
    case (2)
       if (i .eq. lx) is_on = .true.
    case (3)
       if (j .eq. 1) is_on = .true.
    case (4)
       if (j .eq. ly) is_on = .true.
    case (5)
       if (k .eq. 1) is_on = .true.
    case (6)
       if (k .eq. lz) is_on = .true.
    end select

  end function index_is_on_facet

  !> Reports an error and stops execution
  !! @param[optional] error_code The error code to report.
  subroutine neko_error_plain(error_code)
    integer, optional :: error_code

    if (present(error_code)) then
       write(error_unit, *) '*** ERROR ***', error_code
       error stop
    else
       write(error_unit, *) '*** ERROR ***'
       error stop
    end if

  end subroutine neko_error_plain

  !> Reports an error and stops execution
  !! @param error_msg The error message to report.
  subroutine neko_error_msg(error_msg)
    character(len=*) :: error_msg
    write(error_unit, *) '*** ERROR: ', trim(error_msg), ' ***'
    error stop
  end subroutine neko_error_msg

  !> Reports an error allocating a type for a particular base pointer class.
  !! @details Should be used in factories.
  !! @param base_type The base type of the object, which the factory tried to
  !! construct.
  !! @param wrong_type The type that was attempted to construct.
  !! @param known_types A list of the types that are known.
  subroutine neko_type_error(base_type, wrong_type, known_types)
    character(len=*), intent(in) :: base_type
    character(len=*), intent(in) :: wrong_type
    character(len=*), intent(in) :: known_types(:)
    integer :: i

    write(error_unit, *) '*** ERROR WHEN SELECTING TYPE ***'
    write(error_unit, *) 'Type ', wrong_type, ' does not exist for ', base_type
    write(error_unit, *) 'Valid types are:'
    do i = 1, size(known_types)
       write(error_unit, *) "    ", known_types(i)
    end do
    error stop
  end subroutine neko_type_error

  subroutine neko_type_registration_error(base_type, wrong_type, known)
    character(len=*), intent(in) :: base_type
    character(len=*), intent(in) :: wrong_type
    logical, intent(in) :: known

    write(error_unit, *) '*** ERROR WHEN REGISTERING TYPE ***'
    write(error_unit, *) 'Type name ', wrong_type, &
         ' conflicts with and already existing ', base_type, " type"
    if (known) then
       write(error_unit, *) 'Please rename your custom type.'
    else
       write(error_unit, *) 'The already existing type is also custom.' // &
            ' Make all custom type names unique!'
    end if
    error stop
  end subroutine neko_type_registration_error

  !> Reports a warning to standard output
  subroutine neko_warning(warning_msg)
    character(len=*) :: warning_msg
    write(output_unit, *) '*** WARNING: ', trim(warning_msg), ' ***'
  end subroutine neko_warning

  !> Concatenate an array of strings into one string with array items
  !! separated by spaces.
  !! @param array The array of strings.
  !! @param sep The separator put between the strings in the array.
  !! @param prepend Whether to also prepend the string with the separator.
  function concat_string_array(array, sep, prepend) result(result)
    character(len=*), intent(in) :: array(:)
    character(len=*), intent(in) :: sep
    logical, intent(in) :: prepend
    character(:), allocatable :: result
    integer :: i

    result = trim(array(1))
    do i = 2, size(array)
       result = result // sep // trim(array(i))
    end do

    if (prepend) then
       result = sep // result
    end if

  end function concat_string_array

  !> Parse runtime string to total seconds.
  !! Supported formats are:
  !! - seconds as an integer or real value (for example "5400" or "5400.0")
  !! - mm:ss
  !! - hh:mm:ss
  !! - dd-hh:mm:ss
  subroutine read_duration_scalar(runtime_string, runtime_seconds, &
       ierr)
    character(len=*), intent(in) :: runtime_string
    real(kind=rp), intent(inout) :: runtime_seconds
    integer, optional, intent(out) :: ierr
    real(kind=dp) :: parsed_seconds

    if (present(ierr)) ierr = 0
    if (len_trim(runtime_string) .eq. 0) return

    parsed_seconds = read_duration_internal(runtime_string, ierr)

    if (present(ierr)) then
       if (ierr .ne. 0) return
    end if

    runtime_seconds = parsed_seconds
  end subroutine read_duration_scalar

  !> Parse runtime string to components.
  !! The output array maps to the largest unit available based on size:
  !! - size 4: [dd, hh, mm, ss]
  !! - size 3: [hh, mm, ss]
  !! - size 2: [mm, ss]
  !! - size 1: [ss]
  subroutine read_duration_components(runtime_string, runtime_values, ierr)
    character(len=*), intent(in) :: runtime_string
    real(kind=rp), intent(inout) :: runtime_values(:)
    integer, optional, intent(out) :: ierr
    integer, parameter :: i64 = selected_int_kind(18)
    integer :: n_values
    integer(kind=i64) :: total_whole, days, hours, minutes, seconds_whole
    real(kind=dp) :: parsed_seconds, second_value, frac_seconds

    if (present(ierr)) ierr = 0
    if (len_trim(runtime_string) .eq. 0) return

    n_values = size(runtime_values)
    if (n_values .lt. 1 .or. n_values .gt. 4) then
       call set_error_or_throw( &
            'Error parsing duration: output array size must be 1 to 4', ierr)
       return
    end if

    parsed_seconds = read_duration_internal(runtime_string, ierr)

    if (present(ierr)) then
       if (ierr .ne. 0) return
    end if

    total_whole = int(parsed_seconds, kind=i64)
    frac_seconds = parsed_seconds - real(total_whole, kind=dp)
    if (frac_seconds .lt. 0.0_dp) frac_seconds = 0.0_dp

    days = total_whole / 86400_i64
    total_whole = total_whole - 86400_i64 * days
    hours = total_whole / 3600_i64
    total_whole = total_whole - 3600_i64 * hours
    minutes = total_whole / 60_i64
    seconds_whole = total_whole - 60_i64 * minutes

    second_value = real(seconds_whole, kind=dp) + frac_seconds
    if (second_value .ge. 60.0_dp) then
       second_value = second_value - 60.0_dp
       minutes = minutes + 1_i64
       if (minutes .ge. 60_i64) then
          minutes = minutes - 60_i64
          hours = hours + 1_i64
          if (hours .ge. 24_i64) then
             hours = hours - 24_i64
             days = days + 1_i64
          end if
       end if
    end if

    select case (n_values)
    case (4)
       runtime_values(1) = real(days, kind=rp)
       runtime_values(2) = real(hours, kind=rp)
       runtime_values(3) = real(minutes, kind=rp)
       runtime_values(4) = real(second_value, kind=rp)
    case (3)
       runtime_values(1) = real(days * 24_i64 + hours, kind=rp)
       runtime_values(2) = real(minutes, kind=rp)
       runtime_values(3) = real(second_value, kind=rp)
    case (2)
       runtime_values(1) = real((days * 24_i64 + hours) * 60_i64 + &
            minutes, kind=rp)
       runtime_values(2) = real(second_value, kind=rp)
    case (1)
       runtime_values(1) = real(parsed_seconds, kind=rp)
    end select
  end subroutine read_duration_components

  !> Parse runtime string to total seconds.
  real(kind=dp) function read_duration_internal(runtime_string, ierr) &
       result(runtime_seconds)
    character(len=*), intent(in) :: runtime_string
    integer, optional, intent(out) :: ierr
    character(len=:), allocatable :: time_string
    real(kind=dp) :: parsed_seconds, read_real
    logical :: has_minutes, has_hours, has_days
    integer :: read_int, ios, sep

    if (present(ierr)) ierr = 0

    runtime_seconds = 0.0_dp
    time_string = trim(adjustl(runtime_string))

    sep = index(time_string, ':')
    has_minutes = sep .gt. 0
    has_hours = .false.
    if (has_minutes .and. sep .lt. len_trim(time_string)) then
       has_hours = index(time_string(sep + 1:len_trim(time_string)), ':') &
            .gt. 0
    end if
    has_days = index(time_string, '-') .gt. 0

    if ((has_days .and. .not. has_hours) .or. &
         (has_hours .and. .not. has_minutes)) then
       call set_error_or_throw('Error parsing duration: Bad format', ierr)
       return
    end if

    parsed_seconds = 0.0_dp

    ! Read the days field.
    if (has_days) then
       sep = index(time_string, '-')
       read(time_string(1:sep - 1), *, iostat=ios) read_int
       if (ios .ne. 0 .or. read_int .lt. 0) then
          call set_error_or_throw( &
               'Error parsing duration: Invalid days value', ierr)
          return
       end if

       time_string = time_string(sep + 1:)
       parsed_seconds = parsed_seconds + real(86400 * read_int, kind=dp)
    end if

    ! Read the hours.
    if (has_hours) then
       sep = index(time_string, ':')
       read(time_string(1:sep - 1), *, iostat=ios) read_int
       if (ios .ne. 0 .or. read_int .lt. 0 .or. &
            (has_days .and. read_int .gt. 23)) then
          call set_error_or_throw( &
               'Error parsing duration: Invalid hours value', ierr)
          return
       end if

       time_string = time_string(sep + 1:)
       parsed_seconds = parsed_seconds + real(3600 * read_int, kind=dp)
    end if

    ! Read the minutes.
    if (has_minutes) then
       sep = index(time_string, ':')
       read(time_string(1:sep - 1), *, iostat=ios) read_int
       if (ios .ne. 0 .or. read_int .lt. 0 .or. &
            (has_hours .and. read_int .gt. 59)) then
          call set_error_or_throw( &
               'Error parsing duration: Invalid minutes value', ierr)
          return
       end if

       time_string = time_string(sep + 1:)
       parsed_seconds = parsed_seconds + real(60 * read_int, kind=dp)
    end if

    ! Read the seconds.
    read(time_string, *, iostat=ios) read_real
    if (ios .ne. 0 .or. read_real .lt. 0.0_dp .or. &
         (has_minutes .and. read_real .ge. 60.0_dp)) then
       call set_error_or_throw( &
            'Error parsing duration: Invalid seconds value', ierr)
       return
    end if

    parsed_seconds = parsed_seconds + read_real
    runtime_seconds = parsed_seconds

    if (allocated(time_string)) deallocate(time_string)
  end function read_duration_internal

  !> Raise parser error or set ierr, depending on call mode.
  subroutine set_error_or_throw(message, ierr)
    character(len=*), intent(in) :: message
    integer, optional, intent(out) :: ierr

    if (present(ierr)) then
       ierr = 1
    else
       call neko_error(trim(message))
    end if
  end subroutine set_error_or_throw
end module utils
