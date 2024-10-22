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
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  implicit none
  private

  integer, parameter :: NEKO_FNAME_LEN = 1024

  interface neko_error
     module procedure neko_error_plain, neko_error_msg
  end interface neko_error

  public :: neko_error, neko_warning, nonlinear_index, filename_chsuffix, &
       filename_suffix, filename_suffix_pos, filename_tslash_pos, &
       linear_index, split_string, NEKO_FNAME_LEN, index_is_on_facet, &
       concat_string_array, extract_fld_file_index


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

  !> Extract a filename's suffix
  subroutine filename_suffix(fname, suffix)
    character(len=*) :: fname
    character(len=*) :: suffix
    suffix = trim(fname(filename_suffix_pos(fname) + 1:len_trim(fname)))
  end subroutine filename_suffix

  !> Change a filename's suffix
  subroutine filename_chsuffix(fname, new_fname, new_suffix)
    character(len=*) :: fname
    character(len=*) :: new_fname
    character(len=*) :: new_suffix
    integer :: suffix_pos

    suffix_pos = filename_suffix_pos(fname)
    new_fname = trim(fname(1:suffix_pos))//new_suffix

  end subroutine filename_chsuffix

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
    write(error_unit, *) '*** ERROR: ', error_msg, ' ***'
    error stop
  end subroutine neko_error_msg

  !> Reports a warning to standard output
  subroutine neko_warning(warning_msg)
    character(len=*) :: warning_msg
    write(output_unit, *) '*** WARNING: ', warning_msg, ' ***'
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

    if (prepend .eqv. .true.) then
       result = sep // result
    end if

  end function concat_string_array

end module utils
