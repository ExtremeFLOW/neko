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
!
!> Implements the `fluid_stats_simcomp_t` type.

module fluid_stats_simcomp
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use fluid_stats
  use fluid_stats_output, only : fluid_stats_output_t
  use case, only : case_t
  use coefs, only : coef_t
  use comm
  use logger, only : LOG_SIZE, neko_log
  use json_utils, only : json_get, json_get_or_default
  implicit none
  private

  !> A simulation component that writes a 3d field to a file.
  type, public, extends(simulation_component_t) :: fluid_stats_simcomp_t
     !> Output writer.
     type(fluid_stats_t) :: stats              !< Fluid statistics
     type(fluid_stats_output_t) :: stats_output
     real(kind=rp) :: start_time
     real(kind=rp) :: time
   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => fluid_stats_simcomp_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        fluid_stats_simcomp_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => fluid_stats_simcomp_free
     !> Does sampling for statistics
     procedure, pass(this) :: compute_ => fluid_stats_simcomp_compute
     procedure, pass(this) :: output_ => fluid_stats_simcomp_compute
     procedure, pass(this) :: restart_ => fluid_stats_simcomp_restart
  end type fluid_stats_simcomp_t

contains

  !> Constructor from json.
  !> @param json JSON object with the parameters.
  !! @param case The case object.
  subroutine fluid_stats_simcomp_init_from_json(this, json, case)
    class(fluid_stats_simcomp_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision
    character(len=20), allocatable :: fields(:)
    character(len=:), allocatable :: hom_dir
    real(kind=rp) :: start_time
    type(field_t), pointer :: u, v, w, p
    type(coef_t), pointer :: coef

    call this%init_base(json, case)
    call json_get_or_default(json, 'avg_direction', &
         hom_dir, 'none')
    call json_get_or_default(json, 'start_time', &
         start_time, 0.0_rp)
    
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")
    p => neko_field_registry%get_field("p")
    coef => case%fluid%c_Xh
    call fluid_stats_simcomp_init_from_attributes(this, u, v, w, p, coef, start_time, hom_dir)

  end subroutine fluid_stats_simcomp_init_from_json

  !> Actual constructor.
  !! @param u x-velocity
  !! @param v x-velocity
  !! @param w x-velocity
  !! @param coef sem coefs
  !! @param start_time time to start sampling stats
  !! @param hom_dir directions to average in
  subroutine fluid_stats_simcomp_init_from_attributes(this, u, v, w, p, coef, start_time, hom_dir)
    class(fluid_stats_simcomp_t), intent(inout) :: this
    character(len=*), intent(in) :: hom_dir 
    real(kind=rp), intent(in) :: start_time
    type(field_t), intent(inout) :: u, v, w, p !>Should really be intent in I think
    type(coef_t), intent(in) :: coef

    call this%stats%init(coef, u, &
         v, w, p)
    this%start_time = start_time
    this%time = start_time

    this%stats_output = fluid_stats_output_t(this%stats, &
            this%start_time, hom_dir=hom_dir, path=this%case%output_directory)
    call this%case%s%add(this%stats_output,&
                        this%output_controller%control_value,&
                        this%output_controller%control_mode)
  end subroutine fluid_stats_simcomp_init_from_attributes

  !> Destructor.
  subroutine fluid_stats_simcomp_free(this)
    class(fluid_stats_simcomp_t), intent(inout) :: this
    call this%free_base()
  end subroutine fluid_stats_simcomp_free

  subroutine fluid_stats_simcomp_restart(this, t)
    class(fluid_stats_simcomp_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    if(t .gt. this%time) this%time = t
  end subroutine fluid_stats_simcomp_restart

  !> fluid_stats, called depending on compute_control and compute_value
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine fluid_stats_simcomp_compute(this, t, tstep)
    class(fluid_stats_simcomp_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp) :: deltaT
    real(kind=rp) :: sample_start_time, sample_time
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr

    if (t .ge. this%start_time) then
       deltaT = t - this%time

       call MPI_Barrier(NEKO_COMM, ierr)

       sample_start_time = MPI_WTIME()

       call this%stats%update(deltaT)
       call MPI_Barrier(NEKO_COMM, ierr)
       this%time = t

       sample_time = MPI_WTIME() - sample_start_time

       call neko_log%section('Fluid stats')
       write(log_buf,'(A,E15.7)') 'Sampling at time:', t
       call neko_log%message(log_buf)
       write(log_buf,'(A33,E15.7)') 'Simulationtime since last sample:', deltaT
       call neko_log%message(log_buf)
       write(log_buf,'(A,E15.7)') 'Sampling time (s):', sample_time
       call neko_log%message(log_buf)
       call neko_log%end_section()
    end if

  end subroutine fluid_stats_simcomp_compute

end module fluid_stats_simcomp
