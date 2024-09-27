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
  use field_list, only : field_list_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use fld_file_data, only : fld_file_data_t
  use map_2d, only : map_2d_t
  use map_1d, only : map_1d_t
  use coefs, only : coef_t
  use device
  use mean_field, only : mean_field_t
  use output, only : output_t
  use matrix, only : matrix_t
  implicit none
  private

  !> Output for a list of mean fields
  type, public, extends(output_t) :: mean_field_output_t
     !> list of mean fields
     type(mean_field_t), pointer :: mean_fields(:) 
     type(field_list_t) :: fields
     real(kind=rp) :: start_time !> Time to start output
     integer :: n_fields !> Number of fields
     !> Space averaging object for 2 homogeneous directions.
     type(map_1d_t) :: map_1d
     !> Space averaging object for 1 homogeneous direction.
     type(map_2d_t) :: map_2d
     !> The dimension of the output fields. Either 1, 2, or 3.
     integer :: output_dim
   contains
     procedure, pass(this) :: sample => mean_field_output_write
     procedure, pass(this) :: init => mean_field_output_init
  end type mean_field_output_t

contains

  subroutine mean_field_output_init(this,mean_fields, n_fields, start_time, coef, avg_dir, name, path) 
    class(mean_field_output_t), intent(inout):: this
    integer, intent(in) :: n_fields
    class(mean_field_t), intent(inout), target :: mean_fields(n_fields)
    type(coef_t), intent(inout) :: coef
    character(len=*), intent(in) :: avg_dir
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    real(kind=rp), intent(in) :: start_time
    character(len=1024) :: fname
    integer :: i
 
    if (trim(avg_dir) .eq. 'none' .or. &
        trim(avg_dir) .eq. 'x' .or.&
        trim(avg_dir) .eq. 'y' .or.&
        trim(avg_dir) .eq. 'z'&
       ) then
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.fld'
       else if (present(name)) then
          fname = trim(name) // '.fld'
       else if (present(path)) then
          fname = trim(path) // 'user_stats.fld'
       else
          fname = 'user_stats.fld'
       end if

       this%output_dim = 3

       if (trim(avg_dir) .eq. 'x' .or.&
           trim(avg_dir) .eq. 'y' .or.&
           trim(avg_dir) .eq. 'z' ) then
          call this%map_2d%init_char(coef, avg_dir, 1e-7_rp)
          this%output_dim = 2
       end if
    else
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.csv'
       else if (present(name)) then
          fname = trim(name) // '.csv'
       else if (present(path)) then
          fname = trim(path) // 'user_stats.csv'
       else
          fname = 'user_stats.csv'
       end if
       call this%map_1d%init_char(coef, avg_dir, 1e-7_rp)
       this%output_dim = 1
    end if



    call this%init_base(fname)

    call this%fields%init(n_fields)
    this%n_fields = n_fields
    this%mean_fields => mean_fields
    do i = 1, n_fields
       this%fields%items(i)%ptr => this%mean_fields(i)%mf
    end do
    
  end subroutine mean_field_output_init

  !> Sample the mean solution at time @a t and reset
  subroutine mean_field_output_write(this, t)
    class(mean_field_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i
    type(fld_file_data_t) :: output_2d
    type(matrix_t) :: avg_output_1d
    real(kind=rp) :: u, v, w, p

    !if (NEKO_BCKND_DEVICE .eq. 1) then

    !   associate(fields => this%fields%items)
    !     do i = 1, size(fields)
    !        call device_memcpy(fields(i)%ptr%x, fields(i)%ptr%x_d, &
    !             fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
    !             sync=(i .eq. size(fields))) ! Sync on the last field
    !     end do
    !   end associate
    !end if

    !call this%file_%write(this%fields, t)
    associate (out_fields => this%fields%items)
      if (t .ge. this%start_time) then
         if ( NEKO_BCKND_DEVICE .eq. 1) then
            do i = 1, size(out_fields)
               call device_memcpy(out_fields(i)%ptr%x, out_fields(i)%ptr%x_d,&
                    out_fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
                    sync = (i .eq. size(out_fields))) ! Sync on last field
            end do
         end if
         if (this%output_dim .eq. 1) then
            call this%map_1d%average_planes(avg_output_1d, &
                 this%fields)
            call this%file_%write(avg_output_1d, t)
         else if (this%output_dim .eq. 2) then
            call this%map_2d%average(output_2d, this%fields)
            !!Switch around fields to get correct orders
            !!Put average direction mean_vel in scalar45
            !do i = 1, this%map_2d%n_2d
            !   u = output_2d%v%x(i)
            !   v = output_2d%w%x(i)
            !   w = output_2d%p%x(i)
            !   p = output_2d%u%x(i)
            !   output_2d%p%x(i) = p
            !   output_2d%u%x(i) = u
            !   output_2d%v%x(i) = v
            !   output_2d%w%x(i) = w
            !end do
            
            call this%file_%write(output_2d, t)
         else
            call this%file_%write(this%fields, t)
         end if
    do i = 1, this%n_fields
       call this%mean_fields(i)%reset()
    end do
      end if
    end associate

  end subroutine mean_field_output_write

end module mean_field_output
