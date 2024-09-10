! Copyright (c) 2022-2024, The Neko Authors
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
!> Interface to CrayPat F77 API
module craypat
  implicit none
  private
  
#ifdef CRAYPAT
  include 'pat_apif.h'

  public :: craypat_record_start, craypat_record_stop, &
       craypat_region_begin, craypat_region_end
  
contains

  !> Turn on CrayPat recording
  subroutine craypat_record_start
    integer :: ierr
    call PAT_record(PAT_STATE_ON, ierr)
  end subroutine craypat_record_start

  !> Turn off CrayPat recording
  subroutine craypat_record_stop
    integer :: ierr
    call PAT_record(PAT_STATE_OFF, ierr)
  end subroutine craypat_record_stop

  !> Start a CrayPat region
  subroutine craypat_region_begin(name, region_id)
    character(len=*) :: name
    integer, intent(in) :: region_id
    integer :: ierr

    call PAT_record(PAT_STATE_QUERY, ierr)
    if (ierr .eq. PAT_STATE_ON) then
       call PAT_region_begin(region_id, trim(name), ierr)
    end if

  end subroutine craypat_region_begin

  !> End a CrayPat region
  subroutine craypat_region_end(region_id)
    integer, intent(in) :: region_id
    integer :: ierr

    call PAT_record(PAT_STATE_QUERY, ierr)
    if (ierr .eq. PAT_STATE_ON) then
       call PAT_region_end(region_id, ierr)
    end if

  end subroutine craypat_region_end

#endif

end module craypat
