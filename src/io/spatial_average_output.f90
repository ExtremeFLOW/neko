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
!> Output for spatially averaged fields.
module spatial_average_output
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use fld_file_data, only : fld_file_data_t
  use map_1d, only : map_1d_t
  use map_2d, only : map_2d_t
  use coefs, only : coef_t
  use device, only : DEVICE_TO_HOST
  use field_math, only : field_copy
  use matrix, only : matrix_t
  use output, only : output_t
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use utils, only : NEKO_FNAME_LEN, neko_error
  implicit none
  private

  !> Output for spatially averaged fields.
  type, public, extends(output_t) :: spatial_average_output_t
     !> Pointers to the fields to average and write.
     type(field_list_t) :: fields
     !> Space averaging object for 2 homogeneous directions.
     type(map_1d_t) :: map_1d
     !> Space averaging object for 1 homogeneous direction.
     type(map_2d_t) :: map_2d
     !> The dimension of the output fields. Either 1 or 2.
     integer :: output_dim = 0
   contains
     !> Constructor.
     procedure, pass(this) :: init => spatial_average_output_init
     !> Destructor.
     procedure, pass(this) :: free => spatial_average_output_free
     !> Sample the current fields, spatially average, and write.
     procedure, pass(this) :: sample => spatial_average_output_sample
  end type spatial_average_output_t

contains

  !> Constructor.
  !! @param fields Names of the fields to average.
  !! @param coef SEM coefficients.
  !! @param avg_direction Direction(s) to average in. Either 'x', 'y', 'z',
  !! 'xy', 'xz', or 'yz', matching the statistics output semantics.
  !! @param filename Base filename for the output file.
  subroutine spatial_average_output_init(this, fields, coef, &
       avg_direction, filename)
    class(spatial_average_output_t), intent(inout) :: this
    character(len=*), intent(in) :: fields(:)
    type(coef_t), target, intent(inout) :: coef
    character(len=*), intent(in) :: avg_direction
    character(len=*), intent(in) :: filename
    character(len=NEKO_FNAME_LEN) :: output_filename
    integer :: i

    call this%free()
    output_filename = trim(filename)

    call this%fields%init(size(fields))
    do i = 1, size(fields)
       call this%fields%assign(i, neko_registry%get_field(trim(fields(i))))
    end do

    select case (trim(avg_direction))
    case ('x', 'y', 'z')
       this%output_dim = 2
       output_filename = trim(output_filename) // '.fld'
       call this%map_2d%init_char(coef, avg_direction, 1e-7_rp)
    case ('xy', 'xz', 'yz')
       this%output_dim = 1
       output_filename = trim(output_filename) // '.csv'
       call this%map_1d%init_char(coef, avg_direction, 1e-7_rp)
    case default
       call neko_error('Unsupported avg_direction for spatial_average: ' // &
            trim(avg_direction))
    end select

    call this%init_base(output_filename)
  end subroutine spatial_average_output_init

  !> Destructor.
  subroutine spatial_average_output_free(this)
    class(spatial_average_output_t), intent(inout) :: this

    call this%free_base()
    call this%fields%free()
    call this%map_1d%free()
    call this%map_2d%free()

    this%output_dim = 0
  end subroutine spatial_average_output_free

  !> Sample the current fields, spatially average, and write.
  !! @param t The current time value.
  subroutine spatial_average_output_sample(this, t)
    class(spatial_average_output_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    type(fld_file_data_t) :: output_2d
    type(matrix_t) :: output_1d
    type(field_list_t) :: temp_fields
    type(field_t), pointer :: temp_field
    integer, allocatable :: temp_indices(:)
    integer :: i

    select case (this%output_dim)
    case (1)
       call this%fields%copy_from(DEVICE_TO_HOST, .true.)
       call this%map_1d%average_planes(output_1d, this%fields)
       call this%file_%write(output_1d, t)
       call output_1d%free()
    case (2)
       ! map_2d mutates input, so we need temporaries from the scratch registry
       allocate(temp_indices(this%fields%size()))
       call temp_fields%init(this%fields%size())
       do i = 1, this%fields%size()
          call neko_scratch_registry%request_field(temp_field, &
               temp_indices(i), .false.)
          call field_copy(temp_field, this%fields%items(i)%ptr)
          call temp_fields%assign(i, temp_field)
       end do
       call temp_fields%copy_from(DEVICE_TO_HOST, .true.)
       call this%map_2d%average(output_2d, temp_fields)
       call this%file_%write(output_2d, t)
       call output_2d%free()
       call temp_fields%free()
       call neko_scratch_registry%relinquish_field(temp_indices)
       deallocate(temp_indices)
    case default
       call neko_error('Invalid spatial_average output dimension')
    end select
  end subroutine spatial_average_output_sample

end module spatial_average_output
