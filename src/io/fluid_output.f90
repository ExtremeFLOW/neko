! Copyright (c) 2020-2023, The Neko Authors
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
  use fluid_scheme, only : fluid_scheme_t
  use scalar_scheme, only : scalar_scheme_t
  use field_list, only : field_list_t
  use neko_config
  use device
  use output
  implicit none

  !> Fluid output
  type, public, extends(output_t) :: fluid_output_t
     type(field_list_t) :: fluid
   contains
     procedure, pass(this) :: sample => fluid_output_sample
  end type fluid_output_t

  interface fluid_output_t
     module procedure fluid_output_init
  end interface fluid_output_t

contains

  function fluid_output_init(fluid, scalar, name, path) result(this)
    class(fluid_scheme_t), intent(in), target :: fluid
    class(scalar_scheme_t), intent(in), optional, target :: scalar
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

    if (allocated(this%fluid%fields)) then
       deallocate(this%fluid%fields)
    end if

    if (present(scalar)) then
       allocate(this%fluid%fields(5))
    else
       allocate(this%fluid%fields(4))
    end if

    this%fluid%fields(1)%f => fluid%p
    this%fluid%fields(2)%f => fluid%u
    this%fluid%fields(3)%f => fluid%v
    this%fluid%fields(4)%f => fluid%w

    if (present(scalar)) then
       this%fluid%fields(5)%f => scalar%s
    end if
    
  end function fluid_output_init

  !> Sample a fluid solution at time @a t
  subroutine fluid_output_sample(this, t)
    class(fluid_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then

       associate(fields => this%fluid%fields)
         do i = 1, size(fields)
            call device_memcpy(fields(i)%f%x, fields(i)%f%x_d, &
                 fields(i)%f%dof%size(), DEVICE_TO_HOST)
         end do
       end associate

    end if
       
    call this%file_%write(this%fluid, t)

  end subroutine fluid_output_sample
  
end module fluid_output
