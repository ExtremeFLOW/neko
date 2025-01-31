! Copyright (c) 2021, The Neko Authors
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
!> Defines a mean square field
!
module mean_sqr_field
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use mean_field, only : mean_field_t
  use device_math, only : device_cmult, device_addsqr2s2
  use math, only : addsqr2s2
  implicit none
  private

  type, public, extends(mean_field_t) :: mean_sqr_field_t
   contains
     procedure, pass(this) :: update => mean_sqr_field_update
  end type mean_sqr_field_t

contains

  !> Update a mean sqr field
  subroutine mean_sqr_field_update(this, k)
    class(mean_sqr_field_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k !< Time since last sample

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cmult(this%mf%x_d, this%time, size(this%mf%x))
       call device_addsqr2s2(this%mf%x_d, this%f%x_d, k, size(this%mf%x))
       this%time = this%time + k
       call device_cmult(this%mf%x_d, 1.0_rp / this%time, size(this%mf%x))
    else
       this%mf%x = this%mf%x * this%time
       call addsqr2s2(this%mf%x, this%f%x, k, this%mf%dof%size())
       this%time = this%time + k
       this%mf%x = this%mf%x / this%time
    end if

  end subroutine mean_sqr_field_update

end module mean_sqr_field
