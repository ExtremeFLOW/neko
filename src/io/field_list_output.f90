! Copyright (c) 2021-2022, The Neko Authors
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
!> Defines an output for a mean flow field
module field_list_output
  use field
  use field_list
  use num_types
  use output
  implicit none
  private

  type, public, extends(output_t) :: field_list_output_t
     type(field_list_t), pointer :: fields
     real(kind=rp) :: T_begin
   contains
     procedure, pass(this) :: sample => field_list_output_sample
  end type field_list_output_t

  interface field_list_output_t
     module procedure field_list_output_init
  end interface field_list_output_t

contains
  
  function field_list_output_init(fields, T_begin, name, path) result(this)
    type(field_list_t), intent(in), target :: fields
    real(kind=rp), intent(in) :: T_begin
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    type(field_list_output_t) :: this
    character(len=1024) :: fname

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '.fld'
    else if (present(name)) then
       fname = trim(name) // '.fld'
    else if (present(path)) then
       fname = trim(path) // 'field_list.fld'
    else
       fname = 'filed_list.fld'
    end if

    call output_init(this, fname)
    this%fields => fields
    this%T_begin = T_begin
  end function field_list_output_init

  !> Sample a mean flow field at time @a t
  subroutine field_list_output_sample(this, t)
    class(field_list_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i
    associate( out_fields => this%fields%fields)
    if (t .ge. this%T_begin) then
       if ( NEKO_BCKND_DEVICE .eq. 1) then
          do i = 1, size(out_fields)
          call device_memcpy(out_fields(i)%f%x, out_fields(i)%f%x_d,&
               out_fields(i)%f%dof%size(), DEVICE_TO_HOST)
          end do
       end if
       call this%file_%write(this%fields, t)
    end if
end associate

  end subroutine field_list_output_sample
  
end module field_list_output


