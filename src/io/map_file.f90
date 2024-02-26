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
!> NEKTON map file
!! @details This module is used to read/write NEKTON vertex mapping data
module map_file
  use generic_File
  use utils
  use comm
  use map
  implicit none
  private

  !> Interface for NEKTON map files
  type, public, extends(generic_file_t) :: map_file_t
   contains
     procedure :: read => map_file_read
     procedure :: write => map_file_write
  end type map_file_t

contains

  !> Load NEKTON map file
  subroutine map_file_read(this, data)
    class(map_file_t) :: this
    class(*), target, intent(inout) :: data
    type(map_t), pointer :: nm
    integer :: j, k, neli, nnzi, ierr

    call this%check_exists()

    select type(data)
    type is (map_t)
       nm => data
    class default
       call neko_error("Invalid output data")
    end select

    open(unit=10, file=trim(this%fname), status='old', iostat=ierr)
    if (pe_rank .eq. 0) then
       write(*, '(A,A)') " Reading NEKTON map file ", this%fname
    end if

    read(10, *) neli, nnzi

    !> @todo Check if neli matches map%nel

    do j = 1, nm%nel
       read(10, *) nm%imap(j),(nm%vertex(k, j), k=1,nm%nlv)
    end do

    close(unit=10)

  end subroutine map_file_read

  subroutine map_file_write(this, data, t)
    class(map_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    call neko_error("Not implemented yet!")
  end subroutine map_file_write

end module map_file
