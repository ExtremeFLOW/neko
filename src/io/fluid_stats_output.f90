! Copyright (c) 2024-2026, The Neko Authors
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
  use num_types, only : rp, dp
  use fld_file_data, only : fld_file_data_t
  use device, only : device_memcpy, DEVICE_TO_HOST
  use output, only : output_t
  use matrix, only : matrix_t
  use fld_file, only : fld_file_t
  use utils, only : neko_error
  implicit none
  private

  !> Defines an output for the fluid statistics computed using the
  !! `fluid_stats_t` object. Statistics accumulated in the averaged space
  !! are written from their accumulators, 3D mean fields are written as
  !! they are or averaged over the direction(s) of the statistics first.
  type, public, extends(output_t) :: fluid_stats_output_t
     !> Pointer to the object computing the statistics.
     type(fluid_stats_t), pointer :: stats => null()
     real(kind=dp) :: T_begin
     !> The dimension of the output fields. Either 1, 2, or 3.
     integer :: output_dim
   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_stats_output_init
     !> Destructor.
     procedure, pass(this) :: free => fluid_stats_output_free
     !> Samples the fields computed by the `stats` component.
     procedure, pass(this) :: sample => fluid_stats_output_sample
  end type fluid_stats_output_t


contains

  !> Constructor.
  !! @param stats The statistics to write. Their averaging direction, set
  !! at their initialisation, decides the dimension and format of the output.
  !! @param T_begin Time from which the statistics are written.
  !! @param hom_dir Averaging direction of the statistics, must agree with
  !! `stats`.
  !! @param name Name of the output file. Optional.
  !! @param path Path of the output file. Optional.
  subroutine fluid_stats_output_init(this, stats, T_begin, hom_dir, name, path)
    class(fluid_stats_output_t), intent(inout) :: this
    type(fluid_stats_t), intent(inout), target :: stats
    real(kind=dp), intent(in) :: T_begin
    character(len=*), intent(in) :: hom_dir
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: path
    character(len=1024) :: fname
    character(len=4) :: suffix
    integer :: expected_dim

    this%output_dim = stats%output_dim
    if (trim(hom_dir) .eq. 'none' .or. len_trim(hom_dir) .eq. 0) then
       expected_dim = 3
    else
       expected_dim = 3 - len_trim(hom_dir)
    end if
    if (expected_dim .ne. this%output_dim) then
       call neko_error('fluid_stats_output: the averaging direction does' // &
            ' not match the statistics')
    end if

    if (this%output_dim .eq. 1) then
       suffix = '.csv'
    else
       suffix = '.fld'
    end if

    if (present(name) .and. present(path)) then
       fname = trim(path) // trim(name) // suffix
    else if (present(name)) then
       fname = trim(name) // suffix
    else if (present(path)) then
       fname = trim(path) // 'fluid_stats' // suffix
    else
       fname = 'fluid_stats' // suffix
    end if

    call this%init_base(fname)

    select type (ft => this%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .false.
       ft%skip_velocity = .false.
       ft%skip_temperature = .false.
    end select

    this%stats => stats
    this%T_begin = T_begin
  end subroutine fluid_stats_output_init

  !> Destructor.
  subroutine fluid_stats_output_free(this)
    class(fluid_stats_output_t), intent(inout) :: this

    call this%free_base()

    nullify(this%stats)

  end subroutine fluid_stats_output_free

  !> Sample fluid_stats at time @a t
  subroutine fluid_stats_output_sample(this, t)
    class(fluid_stats_output_t), intent(inout) :: this
    real(kind=dp), intent(in) :: t
    integer :: i
    type(matrix_t) :: avg_output_1d
    type(fld_file_data_t) :: output_2d

    if (t .lt. this%T_begin) return

    associate (stats => this%stats, out_fields => this%stats%stat_fields)
      if (stats%avg_dim .eq. 3) then
         ! 3D mean fields, written as they are or averaged on the host.
         call stats%make_strong_grad()
         if (NEKO_BCKND_DEVICE .eq. 1) then
            do i = 1, out_fields%size()
               call device_memcpy(out_fields%items(i)%ptr%x, &
                    out_fields%items(i)%ptr%x_d, out_fields%item_size(i), &
                    DEVICE_TO_HOST, sync = (i .eq. out_fields%size()))
            end do
         end if
         select case (this%output_dim)
         case (1)
            call stats%map_1d%average_planes(avg_output_1d, out_fields)
            call this%file_%write(avg_output_1d, t)
            call avg_output_1d%free()
         case (2)
            call stats%map_2d%average(output_2d, out_fields)
            call reorder_2d(output_2d, stats%map_2d%n_2d)
            call this%file_%write(output_2d, t)
            call output_2d%free()
         case default
            call this%file_%write(out_fields, t)
         end select
      else if (this%output_dim .eq. 1) then
         call stats%map_1d%accumulated_average(avg_output_1d)
         call this%file_%write(avg_output_1d, t)
         call avg_output_1d%free()
      else
         call stats%map_2d%accumulated_average(output_2d)
         call reorder_2d(output_2d, stats%map_2d%n_2d)
         call this%file_%write(output_2d, t)
         call output_2d%free()
      end if
      call stats%reset()
    end associate

  end subroutine fluid_stats_output_sample

  !> Moves the 2D statistics into the slots of the fld file: the pressure
  !! and velocity averages into their own slots, with the mean velocity in
  !! the averaging direction in the last scalar.
  subroutine reorder_2d(output_2d, n)
    type(fld_file_data_t), intent(inout) :: output_2d
    integer, intent(in) :: n
    real(kind=rp) :: u, v, w, p
    integer :: i

    do i = 1, n
       u = output_2d%v%x(i)
       v = output_2d%w%x(i)
       w = output_2d%p%x(i)
       p = output_2d%u%x(i)
       output_2d%p%x(i) = p
       output_2d%u%x(i) = u
       output_2d%v%x(i) = v
       output_2d%w%x(i) = w
    end do

  end subroutine reorder_2d

end module fluid_stats_output
