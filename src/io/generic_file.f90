! Copyright (c) 2019-2023, The Neko Authors
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
module generic_file
  use num_types, only: rp
  use utils, only: neko_warning, neko_error, &
       filename_suffix_pos, filename_suffix
  use comm, only: pe_rank, NEKO_COMM
  use mpi_f08, only: MPI_Bcast, MPI_LOGICAL
  implicit none

  !> A generic file handler.
  type, abstract :: generic_file_t
     character(len=1024), private :: fname
     integer, private :: counter = -1
     integer, private :: start_counter = 0
     !> File format is serial
     logical :: serial = .false.
     logical :: overwrite = .false.
   contains
     !> Generic file constructor.
     procedure :: init => generic_file_init
     !< Write method.
     procedure(generic_file_write), deferred :: write
     !< Read method.
     procedure(generic_file_read), deferred :: read
     !> Get a file's name.
     procedure :: get_fname => generic_file_get_fname
     !> Get base name without counter.
     procedure :: get_base_fname => generic_file_get_base_fname
     !> Set the file counter to @a n.
     procedure :: set_counter => generic_file_set_counter
     !> Get the file counter.
     procedure :: get_counter => generic_file_get_counter
     !> Set the file start counter to @a n.
     procedure :: set_start_counter => generic_file_set_start_counter
     !> Get the file start counter.
     procedure :: get_start_counter => generic_file_get_start_counter
     !> Increment the file counter by 1.
     procedure :: increment_counter => generic_file_increment_counter
     !> Ensure the file exists
     procedure :: check_exists => generic_file_check_exists
     !> Set overwrite mode
     procedure :: set_overwrite => generic_file_set_overwrite
  end type generic_file_t

  abstract interface
     subroutine generic_file_write(this, data, t)
       import :: generic_file_t
       import :: rp
       class(generic_file_t), intent(inout) :: this
       class(*), target, intent(in) :: data
       real(kind=rp), intent(in), optional :: t
     end subroutine generic_file_write
  end interface

  abstract interface
     subroutine generic_file_read(this, data)
       import :: generic_file_t
       class(generic_file_t) :: this
       class(*), target, intent(inout) :: data
     end subroutine generic_file_read
  end interface

contains

  !> Generic file constructor.
  !! @param fname Filename.
  subroutine generic_file_init(this, fname)
    class(generic_file_t) :: this
    character(len=*) :: fname

    this%fname = fname
    this%counter = -1

  end subroutine generic_file_init

  !> Get a file's name.
  function generic_file_get_fname(this) result(fname)
    class(generic_file_t), intent(in) :: this
    character(len=1024) :: fname
    integer :: suffix_pos
    character(len=5) :: id_str

    if (this%counter .ne. -1) then
       suffix_pos = filename_suffix_pos(this%fname)
       write(id_str, '(i5.5)') this%counter
       fname = trim(this%fname(1:suffix_pos-1)) // id_str // &
            this%fname(suffix_pos:)
    else
       fname = this%fname
    end if

  end function generic_file_get_fname

  !> Get base name without counter.
  function generic_file_get_base_fname(this) result(fname)
    class(generic_file_t), intent(in) :: this
    character(len=1024) :: fname
    integer :: suffix_pos

    fname = trim(this%fname)

  end function generic_file_get_base_fname

  !> Set the file counter to @a n.
  subroutine generic_file_set_counter(this, n)
    class(generic_file_t), intent(inout) :: this
    integer, intent(in) :: n
    this%counter = n
  end subroutine generic_file_set_counter

  !> Increment the file counter by 1.
  subroutine generic_file_increment_counter(this)
    class(generic_file_t), intent(inout) :: this

    if (this%counter .eq. -1) then
       this%counter = this%start_counter
    else
       this%counter = this%counter + 1
    end if
  end subroutine generic_file_increment_counter

  !> Set the file start counter to @a n.
  subroutine generic_file_set_start_counter(this, n)
    class(generic_file_t), intent(inout) :: this
    integer, intent(in) :: n
    this%start_counter = n
  end subroutine generic_file_set_start_counter

  !> Get the file counter.
  pure function generic_file_get_counter(this) result(n)
    class(generic_file_t), intent(in) :: this
    integer :: n
    n = this%counter
  end function generic_file_get_counter

  !> Get the file start counter.
  pure function generic_file_get_start_counter(this) result(n)
    class(generic_file_t), intent(in) :: this
    integer :: n
    n = this%start_counter
  end function generic_file_get_start_counter

  !> check if the file exists
  subroutine generic_file_check_exists(this)
    class(generic_file_t), intent(in) :: this
    character(len=1024) :: fname
    logical :: file_exists
    integer :: neko_mpi_ierr

    fname = this%get_fname()
    file_exists = .false.

    if (pe_rank .eq. 0 .or. this%serial) then
       ! Stop if the file does not exist
       inquire(file = fname, exist = file_exists)
    end if

    if (.not. this%serial) then
       call MPI_Bcast(file_exists, 1, MPI_LOGICAL, 0, NEKO_COMM, neko_mpi_ierr)
    end if

    if (.not. file_exists) then
       call neko_error('File does not exist: ' // trim(fname))
    end if

  end subroutine generic_file_check_exists

  !> Set overwrite mode
  subroutine generic_file_set_overwrite(this, overwrite)
    class(generic_file_t), intent(inout) :: this
    logical, intent(in) :: overwrite
    character(len=80) :: suffix

    call filename_suffix(this%get_fname(), suffix)
    call neko_warning("No set_overwrite defined for " // trim(suffix) // " yet")
  end subroutine generic_file_set_overwrite
end module generic_file
