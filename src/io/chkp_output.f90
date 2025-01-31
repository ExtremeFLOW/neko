! Copyright (c) 2021-2023, The Neko Authors
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
!> Defines an output for a checkpoint
module chkp_output
  use checkpoint, only : chkp_t
  use output
  use num_types, only : rp
  implicit none
  private

  type, public, extends(output_t) :: chkp_output_t
     type(chkp_t), pointer :: chkp
   contains
     procedure, pass(this) :: sample => chkp_output_sample
  end type chkp_output_t

  interface chkp_output_t
     module procedure chkp_output_init
  end interface chkp_output_t

contains

  function chkp_output_init(chkp, name, path, fmt) result(this)
    type(chkp_t), intent(in), target :: chkp
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    character(len=*), intent(in), optional :: fmt
    type(chkp_output_t) :: this
    character(len=1024) :: fname
    character(len=10) :: suffix

    suffix = '.chkp'
    if (present(fmt)) then
       if (fmt .eq. 'hdf5') then
          suffix = '.h5'
       end if
    end if

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // trim(suffix)
    else if (present(name)) then
       fname = trim(name) // trim(suffix)
    else if (present(path)) then
       fname = trim(path) // 'fluid' // trim(suffix)
    else
       fname= 'fluid' // trim(suffix)
    end if

    call this%init_base(fname)
    this%chkp => chkp
  end function chkp_output_init

  !> Sample a checkpoint at time @a t
  subroutine chkp_output_sample(this, t)
    class(chkp_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%chkp%sync_host()
    call this%file_%write(this%chkp, t)

  end subroutine chkp_output_sample

end module chkp_output
