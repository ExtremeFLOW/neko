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
!
!> A simulation component that streams data using ADIOS2.
module data_streamer_simcomp
  use num_types, only : rp
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get, &
       json_get_or_lookup_or_default
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use field, only : field_t
  !use field_array, only : field_array_t
  use field_list, only : field_list_t
  use case, only : case_t
  use point_zone, only : point_zone_t
  use time_state, only : time_state_t
  use data_streamer, only : data_streamer_t
  use time_based_controller, only : time_based_controller_t
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: data_streamer_simcomp_t

     !> Name of the fields to be streamed (must exist in the registry).
     character(len=20), allocatable :: field_names(:)

     !> Time after which to start streaming.
     real(kind=rp) :: start_time

     !> List of fields pointing to the fields to stream
     type(field_list_t) :: sampled_fields

     !> Data streamer instance
     type(data_streamer_t) :: dstream

     !> Array of fields to be streamed, in case they need to be subsampled by
     !! point zone.
     !type(field_array_t) :: subsampled_fields

     !> Point zone to use for subsampling
     !class(point_zone_t), pointer :: point_zone => null()

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => data_streamer_simcomp_init_from_json
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          data_streamer_simcomp_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => data_streamer_simcomp_free
     !> Compute the data_streamer_simcomp field
     procedure, pass(this) :: compute_ => data_streamer_simcomp_compute
  end type data_streamer_simcomp_t

contains

  !> Constructor from json.
  subroutine data_streamer_simcomp_init_from_json(this, json, case)
    class(data_streamer_simcomp_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=20), allocatable :: which_fields(:)
    character(len=:), allocatable :: name
    real(kind=rp) :: start_time

    call json_get_or_default(json, "name", name, "data_streamer_simcomp")
    call json_get(json, 'fields', which_fields)

    call json_get_or_lookup_or_default(json, 'start_time', start_time, &
         -1.0_rp)

    call this%init_base(json, case)
    call this%init_from_components(name, which_fields, start_time)

  end subroutine data_streamer_simcomp_init_from_json

  !> Common part of constructors.
  !! @param name The unique name of the simcomp
  subroutine data_streamer_simcomp_init_from_components(this, name, &
       which_fields, start_time)
    class(data_streamer_simcomp_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=20), intent(in) :: which_fields(:)
    real(kind=rp), intent(in) :: start_time

    this%name = name
    this%field_names = which_fields
    this%start_time = start_time

    ! Allocate the field list with the number of fields to sample
    call this%sampled_fields%init(size(which_fields))

    ! Assign fields from the registry
    block
      integer :: i
      do i = 1, size(which_fields)
         call this%sampled_fields%assign_to_field(i, &
              neko_registry%get_field(which_fields(i)))
      end do
    end block

    ! Initialize the data streamer
    call this%dstream%init(this%case%fluid%c_Xh)

  end subroutine data_streamer_simcomp_init_from_components

  !> Destructor.
  subroutine data_streamer_simcomp_free(this)
    class(data_streamer_simcomp_t), intent(inout) :: this
    call this%free_base()

    if (allocated(this%field_names)) deallocate(this%field_names)
    this%start_time = -1.0_rp

    call this%sampled_fields%free()

    call this%dstream%free()

  end subroutine data_streamer_simcomp_free

  !> Compute the data_streamer_simcomp field.
  !! @param time The time state.
  subroutine data_streamer_simcomp_compute(this, time)
    class(data_streamer_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    integer :: i

    if (time%t >= this%start_time) then
       do i = 1, this%sampled_fields%size()

          ! Sync from GPU to CPU
          call this%sampled_fields%items(i)%ptr%copy_from(DEVICE_TO_HOST, &
               .true.)

          call this%dstream%stream(this%sampled_fields%x(i))
       end do
    end if

  end subroutine data_streamer_simcomp_compute

end module data_streamer_simcomp
