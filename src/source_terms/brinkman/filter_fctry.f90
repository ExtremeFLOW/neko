! Copyright (c) 2024, The Neko Authors
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
!
!> Defines a factory subroutine for filters.
module filter_fctry
  use filter, only : filter_t
  use PDE_filter, only : PDE_filter_t
  use json_module, only : json_file
  use coefs, only : coef_t
  use json_utils, only : json_get
  use logger, only : neko_log
  implicit none
  private

  public :: filter_factory

contains

  !> Filter factory. Both constructs and initializes the object.
  !! @param json JSON object initializing the filter.
  subroutine filter_factory(filter, json, coef)
    class(filter_t), allocatable, intent(inout) :: filter
    type(json_file), intent(inout) :: json
    class(coef_t), intent(inout) :: coef
    character(len=:), allocatable :: filter_type

! fuck the json for now
!    call json_get(json, "type", filter_type)
!
!    if (trim(filter_type) .eq. "PDE") then
!       allocate(PDE_filter_t::filter)
!    else
!       call neko_log%error("Unknown filter type: " &
!                           // trim(filter_type))
!       stop
!    end if

    allocate(PDE_filter_t::filter)
    ! Initialize
    call filter%init(json, coef)

  end subroutine filter_factory

end module filter_fctry
