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
module mean_field_output
  use num_types, only : rp
  use fluid_scheme, only : fluid_scheme_t
  use scalar_scheme, only : scalar_scheme_t
  use field_list, only : field_list_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device
  use mean_field, only : mean_field_t
  use output, only : output_t
  implicit none
  private

  !> Output for a list of mean fields
  type, public, extends(output_t) :: mean_field_output_t
     type(mean_field_t), pointer :: mean_fields(:) 
     type(field_list_t) :: fields
     real(kind=rp) :: start_time
     integer :: n_fields
   contains
     procedure, pass(this) :: sample => mean_field_output_sample
     procedure, pass(this) :: init => mean_field_output_init
  end type mean_field_output_t

contains

  subroutine mean_field_output_init(this,mean_fields, n_fields, start_time, precision, name, path) 
    class(mean_field_output_t), intent(inout):: this
    integer, intent(in) :: precision
    integer, intent(in) :: n_fields
    class(mean_field_t), intent(inout), target :: mean_fields(n_fields)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    real(kind=rp), intent(in) :: start_time
    character(len=1024) :: fname
    integer :: i

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // '.fld'
    else if (present(name)) then
       fname = trim(name) // '.fld'
    else if (present(path)) then
       fname = trim(path) // 'mean_fields.fld'
    else
       fname = 'mean_fields.fld'
    end if
 

    call this%init_base(fname, precision)

    call this%fields%init(n_fields)
    this%n_fields = n_fields
    this%mean_fields => mean_fields
    do i = 1, n_fields
       this%fields%items(i)%ptr => this%mean_fields(i)%mf
    end do
    
  end subroutine mean_field_output_init

  !> Sample the mean solution at time @a t and reset
  subroutine mean_field_output_sample(this, t)
    class(mean_field_output_t), intent(inout) :: this
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

    do i = 1, this%n_fields
       call this%mean_fields(i)%reset()
    end do

  end subroutine mean_field_output_sample

end module mean_field_output
