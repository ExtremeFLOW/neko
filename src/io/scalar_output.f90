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
!> Defines an output for a scalar
module scalar_output
  use scalar
  use neko_config
  use device
  use output
  implicit none
  private

  !> Fluid output
  type, public, extends(output_t) :: scalar_output_t
     class(scalar_scheme_t), pointer :: scalar
   contains
     procedure, pass(this) :: sample => scalar_output_sample
  end type scalar_output_t

  interface scalar_output_t
     module procedure scalar_output_init
  end interface scalar_output_t

contains

  function scalar_output_init(scalar, name, path) result(this)
    class(scalar_scheme_t), intent(in), target :: scalar
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(scalar_output_t) :: this
    character(len=1024) :: fname

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '_scalar.fld'
    else if (present(name)) then
       fname = trim(name) // '_scalar.fld'
    else if (present(path)) then
       fname = trim(path) // 'scalar.fld'
    else       
       fname = '.fld'
    end if
    
    call output_init(this, fname)    
    this%scalar => scalar
  end function scalar_output_init

  !> Sample a scalar solution at time @a t
  subroutine scalar_output_sample(this, t)
    class(scalar_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then

       associate(s=> this%scalar%s, dm_Xh => this%scalar%dm_Xh)
       
         call device_memcpy(s%x, s%x_d, dm_Xh%size(), DEVICE_TO_HOST)
         
       end associate
       
    end if
       
    call this%file_%write(this%scalar, t)

  end subroutine scalar_output_sample
  
end module scalar_output
