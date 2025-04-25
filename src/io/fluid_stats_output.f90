! Copyright (c) 2024, The Neko Authors
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
!> Implements `fluid_stats_ouput_t`.
module fluid_stats_output
  use fluid_stats, only : fluid_stats_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use map_1d, only : map_1d_t
  use map_2d, only : map_2d_t
  use fld_file_data, only : fld_file_data_t
  use device
  use output, only : output_t
  use matrix, only : matrix_t
  implicit none
  private

  !> Defines an output for the fluid statistics computed using the 
  !! `fluid_stats_t` object.
  type, public, extends(output_t) :: fluid_stats_output_t
     !> Pointer to the object computing the statistics.
     type(fluid_stats_t), pointer :: stats
     !> Space averaging object for 2 homogeneous directions.
     type(map_1d_t) :: map_1d
     !> Space averaging object for 1 homogeneous direction.
     type(map_2d_t) :: map_2d
     real(kind=rp) :: T_begin
     !> The dimension of the output fields. Either 1, 2, or 3.
     integer :: output_dim
   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_stats_output_init
     !> Samples the fields computed by the `stats` component.
     procedure, pass(this) :: sample => fluid_stats_output_sample
  end type fluid_stats_output_t


contains

  !> Constructor.
  subroutine fluid_stats_output_init(this, stats, T_begin, hom_dir, name, path)
    type(fluid_stats_t), intent(inout), target :: stats
    real(kind=rp), intent(in) :: T_begin
    character(len=*), intent(in) :: hom_dir
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    class(fluid_stats_output_t), intent(inout) :: this
    character(len=1024) :: fname

    if (trim(hom_dir) .eq. 'none' .or. &
        trim(hom_dir) .eq. 'x' .or.&
        trim(hom_dir) .eq. 'y' .or.&
        trim(hom_dir) .eq. 'z'&
       ) then
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.fld'
       else if (present(name)) then
          fname = trim(name) // '.fld'
       else if (present(path)) then
          fname = trim(path) // 'fluid_stats.fld'
       else
          fname = 'fluid_stats.fld'
       end if

       this%output_dim = 3

       if (trim(hom_dir) .eq. 'x' .or.&
           trim(hom_dir) .eq. 'y' .or.&
           trim(hom_dir) .eq. 'z' ) then
          call this%map_2d%init_char(stats%coef, hom_dir, 1e-7_rp)
          this%output_dim = 2
       end if
    else
       if (present(name) .and. present(path)) then
          fname = trim(path) // trim(name) // '.csv'
       else if (present(name)) then
          fname = trim(name) // '.csv'
       else if (present(path)) then
          fname = trim(path) // 'fluid_stats.csv'
       else
          fname = 'fluid_stats.csv'
       end if
       call this%map_1d%init_char(stats%coef, hom_dir, 1e-7_rp)
       this%output_dim = 1
    end if

    call this%init_base(fname)
    this%stats => stats
    this%T_begin = T_begin
  end subroutine fluid_stats_output_init

  !> Sample fluid_stats at time @a t
  subroutine fluid_stats_output_sample(this, t)
    class(fluid_stats_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i
    type(matrix_t) :: avg_output_1d
    type(fld_file_data_t) :: output_2d
    real(kind=rp) :: u, v, w, p
    associate (out_fields => this%stats%stat_fields%items)
      if (t .ge. this%T_begin) then
         call this%stats%make_strong_grad()
         if ( NEKO_BCKND_DEVICE .eq. 1) then
            do i = 1, size(out_fields)
               call device_memcpy(out_fields(i)%ptr%x, out_fields(i)%ptr%x_d,&
                    out_fields(i)%ptr%dof%size(), DEVICE_TO_HOST, &
                    sync = (i .eq. size(out_fields))) ! Sync on last field
            end do
         end if
         if (this%output_dim .eq. 1) then
            call this%map_1d%average_planes(avg_output_1d, &
                 this%stats%stat_fields)
            call this%file_%write(avg_output_1d, t)
         else if (this%output_dim .eq. 2) then
            call this%map_2d%average(output_2d, this%stats%stat_fields)
            !Switch around fields to get correct orders
            !Put average direction mean_vel in scalar45
            do i = 1, this%map_2d%n_2d
               u = output_2d%v%x(i)
               v = output_2d%w%x(i)
               w = output_2d%p%x(i)
               p = output_2d%u%x(i)
               output_2d%p%x(i) = p
               output_2d%u%x(i) = u
               output_2d%v%x(i) = v
               output_2d%w%x(i) = w
            end do
            
            call this%file_%write(output_2d, t)
         else
            call this%file_%write(this%stats%stat_fields, t)
         end if
         call this%stats%reset()
      end if
    end associate
  end subroutine fluid_stats_output_sample

end module fluid_stats_output


