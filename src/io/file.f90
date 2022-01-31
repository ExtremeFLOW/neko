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
module file
  use utils
  use generic_file
  use nmsh_file
  use chkp_file        
  use map_file
  use rea_file
  use re2_file
  use fld_file
  use vtk_file
  implicit none
  
  type file_t
     class(generic_file_t), allocatable :: file_type
   contains
     procedure :: write => file_write
     procedure :: read => file_read
     procedure :: set_counter => file_set_counter
     final :: file_free
  end type file_t

  interface file_t
     module procedure file_init
  end interface file_t

contains

  !> File reader/writer constructor
  !! @param fname Filename
  function file_init(fname) result(this)
    character(len=*) :: fname
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
    else
       call neko_error('Unknown file format')
    end if

    call this%file_type%init(fname)

  end function file_init

  !> File operation destructor
  subroutine file_free(this)
    type(file_t), intent(inout) :: this

    if (allocated(this%file_type)) then
       deallocate(this%file_type)
    end if

  end subroutine file_free

  !> Write @a data to a file
  !! @param data Data to be written
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
   
  !> Read @a data from a file
  !! @param data Read data
  subroutine file_read(this, data)
    class(file_t), intent(in) :: this
    class(*), intent(inout) :: data

    call this%file_type%read(data)
    
  end subroutine file_read

  !> Set a file's counter
  subroutine file_set_counter(this, n)
    class(file_t), intent(inout) :: this
    integer, intent(in) :: n
    call this%file_type%set_counter(n)
  end subroutine file_set_counter

end module file
