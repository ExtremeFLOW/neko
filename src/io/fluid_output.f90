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
!> Defines an output for a fluid
module fluid_output
  use fluid_scheme
  use neko_config
  use device
  use output
  implicit none
  private

  !> Fluid output
  type, public, extends(output_t) :: fluid_output_t
     class(fluid_scheme_t), pointer :: fluid
   contains
     procedure, pass(this) :: sample => fluid_output_sample
  end type fluid_output_t

  interface fluid_output_t
     module procedure fluid_output_init
  end interface fluid_output_t

contains

  function fluid_output_init(fluid, name, path) result(this)
    class(fluid_scheme_t), intent(in), target :: fluid
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(fluid_output_t) :: this
    character(len=1024) :: fname

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '.fld'
    else if (present(name)) then
       fname = trim(name) // '.fld'
    else if (present(path)) then
       fname = trim(path) // 'field.fld'
    else       
       fname = 'field.fld'
    end if
    
    call output_init(this, fname)    
    this%fluid => fluid
  end function fluid_output_init

  !> Sample a fluid solution at time @a t
  subroutine fluid_output_sample(this, t)
    class(fluid_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    if (NEKO_BCKND_DEVICE .eq. 1) then

       associate(p => this%fluid%p, u =>this%fluid%u, v => this%fluid%v, &
            w => this%fluid%w, dm_Xh => this%fluid%dm_Xh)
       
         call device_memcpy(p%x, p%x_d, dm_Xh%size(), DEVICE_TO_HOST)
         call device_memcpy(u%x, u%x_d, dm_Xh%size(), DEVICE_TO_HOST)
         call device_memcpy(v%x, v%x_d, dm_Xh%size(), DEVICE_TO_HOST)
         call device_memcpy(w%x, w%x_d, dm_Xh%size(), DEVICE_TO_HOST)
         
       end associate
       
    end if
       
    call this%file_%write(this%fluid, t)

  end subroutine fluid_output_sample
  
end module fluid_output
