! Copyright (c) 2026, The Neko Authors
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
!> Implements `field_output_t`.
module field_output
  use num_types, only : rp
  use field_list, only : field_list_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, DEVICE_TO_HOST
  use output, only : output_t
  implicit none
  private

  !> A simple output saving a list of fields to a file.
  type, public, extends(output_t) :: field_output_t
     !> The list of fields to save.
     type(field_list_t) :: fields
   contains
     !> Constructor.
     procedure, pass(this) :: init => field_output_init
     !> Destructor
     procedure, pass(this) :: free => field_output_free
     !> Writes the data.
     procedure, pass(this) :: sample => field_output_sample
  end type field_output_t

contains

  !> Constructor.
  !! @param name The base name of the files.
  !! @param nfields The number of field pointers to preallocate in the field
  !!        list.
  !! @param precision the precision of the reals in the file.
  !! @param path Optional path to the write folder.
  !! @param format Optional suffix for the file name.
  subroutine field_output_init(this, name, nfields, precision, path, format)
    class(field_output_t), intent (inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: nfields
    integer, intent(in), optional :: precision
    character(len=*), intent(in), optional :: path
    character(len=*), intent(in), optional :: format
    character(len=1024) :: fname, suffix

    call this%free()

    suffix = '.fld'
    if (present(format)) then
       select case (trim(format))
       case ('nek5000', 'fld')
          suffix = '.fld'
       case ('vtkhdf')
          suffix = '.vtkhdf'
       case ('adios2')
          suffix = '.bp'
       case default
          suffix = '.' // trim(format)
       end select
    end if

    if (present(path)) then
       fname = trim(path) // trim(name) // trim(suffix)
    else
       fname = trim(name) // trim(suffix)
    end if

    call this%init_base(fname, precision)

    call this%fields%init(nfields)

  end subroutine field_output_init

  !> Destructor
  subroutine field_output_free(this)
    class(field_output_t), intent(inout) :: this

    call this%free_base()
    call this%fields%free()

  end subroutine field_output_free

  !> Writes the data.
  !! @param t The time value.
  subroutine field_output_sample(this, t)
    class(field_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%fields%size()
          associate(field => this%fields%items(i)%ptr)
            call field%copy_from(DEVICE_TO_HOST, &
                 sync = i .eq. this%fields%size())
          end associate
       end do
    end if

    call this%file_%write(this%fields, t)

  end subroutine field_output_sample

end module field_output
