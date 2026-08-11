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
  use field, only : field_t
  use case, only : case_t
  use comm, only : pe_rank
  use time_state, only : time_state_t
  use data_streamer, only : data_streamer_t
  use logger, only : neko_log, NEKO_LOG_DEBUG
  use time_based_controller, only : time_based_controller_t
  use registry, only : neko_registry
  use device, only : DEVICE_TO_HOST
  use utils, only : NEKO_VARNAME_LEN, neko_error
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public, extends(simulation_component_t) :: data_streamer_simcomp_t

     !> Name of the fields to be streamed (must exist in the registry).
     character(len=NEKO_VARNAME_LEN), allocatable :: field_names(:)

     !> Time after which to start streaming.
     real(kind=rp) :: start_time

     !> Data streamer instance
     type(data_streamer_t) :: dstream

     !> Whether this simcomp has been checked
     logical :: checked = .false.

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
     !> Check the simcomp
     procedure, pass(this) :: check => data_streamer_simcomp_check
  end type data_streamer_simcomp_t

contains

  !> Checks the validity of the simcomp. Currently, makes sure that:
  !!   - All fields have the same number of elements
  !!   - All fields have the same size (dof%size())
  !!   - ADIOS2 has been installed properly
  subroutine data_streamer_simcomp_check(this)
    class(data_streamer_simcomp_t), intent(inout) :: this

    type(field_t), pointer :: f
    integer :: i, num_elements
    integer :: dof_size

    if (size(this%field_names) .eq. 0) call neko_error("Empty list of fields")

    f => neko_registry%get_field_by_name(this%field_names(1))
    num_elements = f%msh%nelv
    dof_size = f%size()

    ! Loop through all fields and make sure they are "the same".
    do i = 2, size(this%field_names)
       f => neko_registry%get_field_by_name(this%field_names(i))

       if (f%msh%nelv .ne. num_elements .or. f%size() .ne. dof_size) then
          call neko_error("Invalid field :" // trim(f%name) // &
               "! All streamed fields must have the same size.")
       end if

    end do

#ifdef HAVE_ADIOS2
#else
    call neko_error('Neko has not been configured with ADIOS2 support! ' // &
         'Use --with-adios2=your/adios2/installation')
#endif

    this%checked = .true.

  end subroutine data_streamer_simcomp_check

  !> Constructor from json.
  subroutine data_streamer_simcomp_init_from_json(this, json, case)
    class(data_streamer_simcomp_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=NEKO_VARNAME_LEN), allocatable :: which_fields(:)
    character(len=:), allocatable :: name
    real(kind=rp) :: start_time
    logical :: stream_mesh

    call json_get_or_default(json, "name", name, "data_streamer_simcomp")
    call json_get(json, 'fields', which_fields)

    call json_get_or_lookup_or_default(json, 'start_time', start_time, &
         -1.0_rp)

    call json_get_or_default(json, 'stream_mesh', stream_mesh, .false.)

    call this%init_base(json, case)
    call this%init_from_components(name, which_fields, start_time, stream_mesh)

  end subroutine data_streamer_simcomp_init_from_json

  !> Common part of constructors.
  !! @param name The unique name of the simcomp.
  !! @param which_fields The names of the fields to be streamed.
  !! @param start_time Time after which to start streaming.
  subroutine data_streamer_simcomp_init_from_components(this, name, &
       which_fields, start_time, stream_mesh)
    class(data_streamer_simcomp_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: which_fields(:)
    real(kind=rp), intent(in) :: start_time
    logical, intent(in) :: stream_mesh
    type(field_t), pointer :: f

    this%name = name
    this%field_names = which_fields
    this%start_time = start_time

    call this%check()

    ! Use the mesh/element information in the first field to initialize the
    ! streamer. We assume that all fields share the same mesh.
    f => neko_registry%get_field(which_fields(1))
    call this%dstream%init(f%msh, f%xh)

    if (stream_mesh) then
       call neko_log%message("Using field " // trim(f%name) // &
            " for streaming mesh", lvl = NEKO_LOG_DEBUG)
       call neko_log%message("Streaming mesh: x-coordinates", &
            lvl = NEKO_LOG_DEBUG)
       call this%dstream%stream(f%dof%x)
       call neko_log%message("Streaming mesh: y-coordinates", &
            lvl = NEKO_LOG_DEBUG)
       call this%dstream%stream(f%dof%y)
       call neko_log%message("Streaming mesh: z-coordinates", &
            lvl = NEKO_LOG_DEBUG)
       call this%dstream%stream(f%dof%z)
    end if

    nullify(f)

  end subroutine data_streamer_simcomp_init_from_components

  !> Destructor.
  subroutine data_streamer_simcomp_free(this)
    class(data_streamer_simcomp_t), intent(inout) :: this
    call this%free_base()

    if (allocated(this%field_names)) deallocate(this%field_names)
    this%start_time = -1.0_rp
    call this%dstream%free()

  end subroutine data_streamer_simcomp_free

  !> Compute the data_streamer_simcomp field.
  !! @param time The time state.
  subroutine data_streamer_simcomp_compute(this, time)
    class(data_streamer_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    integer :: i
    type(field_t), pointer :: f

    if (.not. this%checked) call neko_error("Simcomp not checked!")

    if (time%t .ge. this%start_time) then
       do i = 1, size(this%field_names)

          ! Sync from GPU to CPU
          f => neko_registry%get_field_by_name(this%field_names(i))
          call f%copy_from(DEVICE_TO_HOST, .true.)

          call neko_log%message("Streaming field: " // this%field_names(i))

          ! Stream the field
          call this%dstream%stream(f%x)
       end do
    end if

  end subroutine data_streamer_simcomp_compute

end module data_streamer_simcomp
