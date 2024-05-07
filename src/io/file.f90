! Copyright (c) 2019-2024, The Neko Authors
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
!> Module for file I/O operations.
module file
  use utils, only : neko_error, neko_warning, filename_suffix
  use num_types, only : rp
  use generic_file, only : generic_file_t
  use nmsh_file, only : nmsh_file_t
  use chkp_file, only : chkp_file_t
  use map_file, only : map_file_t
  use rea_file, only : rea_file_t
  use re2_file, only : re2_file_t
  use fld_file, only : fld_file_t
  use fld_file_data, only : fld_file_data_t
  use vtk_file, only : vtk_file_t
  use stl_file, only : stl_file_t
  use csv_file, only : csv_file_t
  use hdf5_file, only : hdf5_file_t
  implicit none

  !> A wrapper around a polymorphic `generic_file_t` that handles its init.
  !! This is essentially a factory for `generic_file_t` descendants additionally
  !! handling special CSV file parameters (header and precision).
  type file_t
     class(generic_file_t), allocatable :: file_type
   contains
     !> Writes data to a file.
     procedure :: write => file_write
     !> Read @a data from a file.
     procedure :: read => file_read
     !> Set a file's counter.
     procedure :: set_counter => file_set_counter
     !> Set a file's start counter.
     procedure :: set_start_counter => file_set_start_counter
     !> Set a file's header.
     procedure :: set_header => file_set_header
     !> Set a file's output precision.
     procedure :: set_precision => file_set_precision
     !> File operation destructor.
     final :: file_free
  end type file_t

  interface file_t
     module procedure file_init
  end interface file_t

contains

  !> File reader/writer constructor.
  !! @param fname Filename.
  function file_init(fname, header, precision) result(this)
    character(len=*) :: fname
    character(len=*), optional :: header
    integer, optional :: precision
    type(file_t), target :: this
    character(len=80) :: suffix
    class(generic_file_t), pointer :: q

    call filename_suffix(fname, suffix)

    if (allocated(this%file_type)) then
       deallocate(this%file_type)
    end if

    if (suffix .eq. "rea") then
       allocate(rea_file_t::this%file_type)
    else if (suffix .eq. "re2") then
       allocate(re2_file_t::this%file_type)
    else if (suffix .eq. "map") then
       allocate(map_file_t::this%file_type)
    else if (suffix .eq. "vtk") then
       allocate(vtk_file_t::this%file_type)
    else if (suffix .eq. "nmsh") then
       allocate(nmsh_file_t::this%file_type)
    else if (suffix .eq. "fld") then
       allocate(fld_file_t::this%file_type)
    else if (suffix .eq. "chkp") then
       allocate(chkp_file_t::this%file_type)
    else if (suffix .eq. "stl") then
       allocate(stl_file_t::this%file_type)
    else if (suffix .eq. "csv") then
       allocate(csv_file_t::this%file_type)
       this%file_type%serial = .true.
    else if ((suffix .eq. "hdf5") .or. (suffix .eq. "h5")) then
       allocate(hdf5_file_t::this%file_type)
    else
       call neko_error('Unknown file format')
    end if

    call this%file_type%init(fname)

    if (present(header)) then
       call this%set_header(header)
    end if

    if (present(precision)) then
       call this%set_precision(precision)
    end if

  end function file_init

  !> File operation destructor.
  subroutine file_free(this)
    type(file_t), intent(inout) :: this

    if (allocated(this%file_type)) then
       deallocate(this%file_type)
    end if

  end subroutine file_free

  !> Writes data to a file.
  !! @param data Data to be written.
  subroutine file_write(this, data, t)
    class(file_t), intent(inout) :: this
    class(*), intent(inout) :: data
    real(kind=rp), intent(in), optional :: t

    if (present(t)) then
       call this%file_type%write(data, t)
    else
       call this%file_type%write(data)
    end if

  end subroutine file_write

  !> Read @a data from a file.
  !! @param data Read data.
  subroutine file_read(this, data)
    class(file_t), intent(in) :: this
    class(*), intent(inout) :: data

    call this%file_type%read(data)

  end subroutine file_read

  !> Set a file's counter.
  subroutine file_set_counter(this, n)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: n

    select type(ft => this%file_type)
    class is (generic_file_t)
       call ft%set_counter(n)
    end select

  end subroutine file_set_counter

  !> Set a file's start counter.
  subroutine file_set_start_counter(this, n)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: n

    select type(ft => this%file_type)
    class is (generic_file_t)
       call ft%set_start_counter(n)
    end select

  end subroutine file_set_start_counter

  !> Set a file's header.
  subroutine file_set_header(this, hd)
    class(file_t), intent(inout) :: this
    character(len=*), intent(in) :: hd

    character(len=80) :: suffix

    select type(ft => this%file_type)
    class is (csv_file_t)
       call ft%set_header(hd)
    class default
       call filename_suffix(this%file_type%fname, suffix)
       call neko_warning("No set_header defined for " // trim(suffix) // " yet!")
    end select

  end subroutine file_set_header

  !> Set a file's output precision.
  !! @param precision Precision as defined in `num_types`.
  subroutine file_set_precision(this, precision)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: precision

    character(len=80) :: suffix

    select type(ft => this%file_type)
    type is (fld_file_t)
       call ft%set_precision(precision)
    class default
       call filename_suffix(this%file_type%fname, suffix)
       call neko_warning("No precision strategy defined for " // trim(suffix) //&
            " files!")
    end select

  end subroutine file_set_precision

end module file
