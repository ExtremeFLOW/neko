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
!> Implements `fld_file_output_t`.
module fld_file_output
  use num_types, only : rp
  use field_list, only : field_list_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, DEVICE_TO_HOST
  use output, only : output_t
  implicit none
  private

  !> A simple output saving a list of fields to a .fld file.
  type, public, extends(output_t) :: fld_file_output_t
     ! The list of fields to save.
     type(field_list_t) :: fields
   contains
     ! Constructor.
     procedure, pass(this) :: init => fld_file_output_init
     ! Writes the data.
     procedure, pass(this) :: sample => fld_file_output_sample
  end type fld_file_output_t

contains

  !> Constructor.
  !! @param precision the precison of the reals in the file.
  !! @param name The base name of the files.
  !! @param name The number of field pointers to preallocate in the field list.
  !! @param path Optional path to the write folder.
  subroutine fld_file_output_init(this, precision, name, nfields, path)
    class(fld_file_output_t), intent (inout) :: this
    integer, intent(in) :: precision
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: path
    integer, intent(in) :: nfields
    character(len=1024) :: fname

    if (present(path)) then
       fname = trim(path) // trim(name) // '.fld'
    else
       fname = trim(name) // '.fld'
    end if

    call this%init_base(fname, precision)

    call this%fields%init(nfields)

   end subroutine fld_file_output_init

  !> Writes the data.
  !! @param t The time value.
  subroutine fld_file_output_sample(this, t)
    class(fld_file_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(fields => this%fields%items)
         do i = 1, size(fields)
            call device_memcpy(fields(i)%ptr%x, fields(i)%ptr%x_d, &
                 fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
                 sync=(i .eq. size(fields))) ! Sync on the last field
         end do
       end associate

    end if

    call this%file_%write(this%fields, t)

  end subroutine fld_file_output_sample

end module fld_file_output
