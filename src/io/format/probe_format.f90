! Copyright (c) 2020-2022, The Neko Authors
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
!> File format for probes input
!! @details this module defines the probe data type
module probe_format
  use num_types
  use point
  implicit none

  type, public :: probe_t
     integer :: npts                             !< Number of probe points
     type(point_t), allocatable :: xyz_coords(:) !< Array of point coordinates
   contains
     procedure, pass(this) :: init => probe_init
     procedure, pass(this) :: free => probe_free
  end type probe_t

contains

  !> Initialize a probe object with an array of size npts
  subroutine probe_init(this, npts)
    class(probe_t), intent(inout) :: this
    integer, intent(in) :: npts

    this%npts = npts
    allocate(this%xyz_coords(npts))

  end subroutine probe_init

  !> Deallocate a probe object
  subroutine probe_free(this)
    class(probe_t), intent(inout) :: this

    if (allocated(this%xyz_coords)) then
       deallocate(this%xyz_coords)
    end if

    this%npts = -1

  end subroutine probe_free

end module probe_format
