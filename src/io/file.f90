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
  use bp_file, only : bp_file_t
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
     !> Constructor
     procedure, pass (this) :: init => file_init
     !> Writes data to a file.
     procedure :: write => file_write
     !> Read @a data from a file.
     procedure :: read => file_read
     !> Get a file's name.
     procedure :: get_fname => file_get_fname
     !> Get a file's base name.
     procedure :: get_base_fname => file_get_base_fname
     !> Get a file's counter.
     procedure :: get_counter => file_get_counter
     !> Set a file's counter.
     procedure :: set_counter => file_set_counter
     !> Set a file's start counter.
     procedure :: set_start_counter => file_set_start_counter
     !> Set a file's header.
     procedure :: set_header => file_set_header
     !> Set a file's output precision.
     procedure :: set_precision => file_set_precision
     !> Set a file's output layout.
     procedure :: set_layout => file_set_layout
     !> Sets the file's overwrite flag.
     procedure, pass (this) :: set_overwrite => file_set_overwrite
     !> File operation destructor.
     procedure, pass(this) :: free => file_free
  end type file_t

contains

  !> Constructor.
  !! @param fname Filename.
  subroutine file_init(this, fname, header, precision, layout, overwrite)
    class(file_t), intent(inout) :: this
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: header
    integer, intent(in), optional :: precision
    integer, intent(in), optional :: layout
    logical, intent(in), optional :: overwrite
    character(len=80) :: suffix
    class(generic_file_t), pointer :: q

    call filename_suffix(fname, suffix)

    if (allocated(this%file_type)) then
       deallocate(this%file_type)
    end if

    select case (suffix)
    case ("rea")
       allocate(rea_file_t::this%file_type)
    case ("re2")
       allocate(re2_file_t::this%file_type)
    case ("map")
       allocate(map_file_t::this%file_type)
    case ("vtk")
       allocate(vtk_file_t::this%file_type)
    case ("nmsh")
       allocate(nmsh_file_t::this%file_type)
    case ("bp")
       allocate(bp_file_t::this%file_type)
    case ("fld")
       allocate(fld_file_t::this%file_type)
    case ("chkp")
       allocate(chkp_file_t::this%file_type)
    case ("stl")
       allocate(stl_file_t::this%file_type)
    case ("csv")
       allocate(csv_file_t::this%file_type)
       this%file_type%serial = .true.
    case ("hdf5", "h5")
       allocate(hdf5_file_t::this%file_type)
    case default
       call neko_error('Unknown file format')
    end select

    call this%file_type%init(fname)

    if (present(header)) then
       call this%set_header(header)
    end if

    if (present(precision)) then
       call this%set_precision(precision)
    end if

    if (present(layout).and. (suffix .eq. "bp")) then
       call this%set_layout(layout)
    end if

    if (present(overwrite)) then
       call this%set_overwrite(overwrite)
    end if

  end subroutine file_init

  !> File operation destructor.
  subroutine file_free(this)
    class(file_t), intent(inout) :: this

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

    call this%file_type%write(data, t = t)

  end subroutine file_write

  !> Read @a data from a file.
  !! @param data Read data.
  subroutine file_read(this, data)
    class(file_t), intent(in) :: this
    class(*), intent(inout) :: data

    call this%file_type%read(data)

  end subroutine file_read

  !> Get a file's name.
  function file_get_fname(this) result(fname)
    class(file_t), intent(in) :: this
    character(len=1024) :: fname

    fname = ""

    select type (ft => this%file_type)
    class is (generic_file_t)
       fname = ft%get_fname()
    end select

  end function file_get_fname

  !> Get a file's base name.
  function file_get_base_fname(this) result(fname)
    class(file_t), intent(in) :: this
    character(len=1024) :: fname

    fname = ""

    select type (ft => this%file_type)
    class is (generic_file_t)
       fname = ft%get_base_fname()
    end select

  end function file_get_base_fname

  !> Get a file's counter.
  function file_get_counter(this) result(n)
    class(file_t), intent(inout) :: this
    integer :: n
    n = 0

    select type (ft => this%file_type)
    class is (generic_file_t)
       n = ft%get_counter()
    end select

  end function file_get_counter

  !> Set a file's counter.
  subroutine file_set_counter(this, n)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: n

    select type (ft => this%file_type)
    class is (generic_file_t)
       call ft%set_counter(n)
    end select

  end subroutine file_set_counter

  !> Set a file's start counter.
  subroutine file_set_start_counter(this, n)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: n

    select type (ft => this%file_type)
    class is (generic_file_t)
       call ft%set_start_counter(n)
    end select

  end subroutine file_set_start_counter

  !> Set a file's header.
  subroutine file_set_header(this, hd)
    class(file_t), intent(inout) :: this
    character(len=*), intent(in) :: hd
    character(len=80) :: suffix

    select type (ft => this%file_type)
    class is (csv_file_t)
       call ft%set_header(hd)
    class default
       call filename_suffix(this%file_type%get_fname(), suffix)
       call neko_warning("No set_header defined for " // trim(suffix) // " yet")
    end select

  end subroutine file_set_header

  !> Set a file's output precision.
  !! @param precision Precision as defined in `num_types`.
  subroutine file_set_precision(this, precision)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: precision
    character(len=80) :: suffix

    select type (ft => this%file_type)
    type is (fld_file_t)
       call ft%set_precision(precision)
    type is (bp_file_t)
       call ft%set_precision(precision)
    class default
       call filename_suffix(this%file_type%get_fname(), suffix)
       call neko_warning("No precision strategy defined for " // trim(suffix) &
            // " files")
    end select

  end subroutine file_set_precision

  !> Set a file's output layout.
  !! @param layout The data layout as defined in bp_file.f90 and src/io/buffer/.
  subroutine file_set_layout(this, layout)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: layout
    character(len=80) :: suffix

    select type (ft => this%file_type)
    type is (bp_file_t)
       call ft%set_layout(layout)
    class default
       call filename_suffix(this%file_type%get_fname(), suffix)
       call neko_warning("No set_layout defined for " // trim(suffix) // " yet")
    end select

  end subroutine file_set_layout

  !> Sets the file's overwrite flag.
  subroutine file_set_overwrite(this, overwrite)
    class(file_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    character(len=80) :: suffix

    select type (ft => this%file_type)
    class is (generic_file_t)
       call ft%set_overwrite(overwrite)
    end select
  end subroutine file_set_overwrite

end module file
