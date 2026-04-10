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
!> Implements the `scalar_sgs_stats_simcomp_t` type.
module scalar_sgs_stats_simcomp
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use time_state, only : time_state_t
  use field, only : field_t
  use scalar_sgs_stats, only : scalar_sgs_stats_t
  use scalar_sgs_stats_output, only : scalar_sgs_stats_output_t
  use case, only : case_t
  use coefs, only : coef_t
  use utils, only : NEKO_FNAME_LEN, filename_suffix, filename_tslash_pos, &
       neko_error ! just for now
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use json_utils, only : json_get, json_get_or_default
  use comm, only : NEKO_COMM
  use amr_reconstruct, only : amr_reconstruct_t
  use mpi_f08, only : MPI_WTIME, MPI_Barrier
  implicit none
  private

  !> A simulation component that computes the subgrid-scale contributions
  !! to the Reynolds stresses in LES.
  !!
  !! For further details see the Neko documentation.
  type, public, extends(simulation_component_t) :: scalar_sgs_stats_simcomp_t
     !> Backbone object computing the satistics
     type(scalar_sgs_stats_t) :: stats
     !> Output writer.
     type(scalar_sgs_stats_output_t) :: stats_output
     !> Time value at which the sampling of statistics is initiated.
     real(kind=rp) :: start_time
     real(kind=rp) :: time
     logical :: default_fname = .true.

   contains
     generic :: init_from_components => &
          init_from_components_alphat, &
          init_from_components_nut
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => scalar_sgs_stats_simcomp_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_components_alphat => &
          scalar_sgs_stats_simcomp_init_from_components_alphat
     !> Actual constructor.
     procedure, pass(this) :: init_from_components_nut => &
          scalar_sgs_stats_simcomp_init_from_components_nut
     !> Destructor.
     procedure, pass(this) :: free => scalar_sgs_stats_simcomp_free
     !> Does sampling for statistics.
     procedure, pass(this) :: compute_ => scalar_sgs_stats_simcomp_compute
     !> Write the statistics to disk.
     procedure, pass(this) :: output_ => scalar_sgs_stats_simcomp_compute
     !> Restart the simcomp.
     procedure, pass(this) :: restart_ => scalar_sgs_stats_simcomp_restart
     !> AMR restart
     procedure, pass(this) :: amr_restart => &
          scalar_sgs_stats_simcomp_amr_restart
  end type scalar_sgs_stats_simcomp_t

contains

  !> Constructor from json.
  !> @param json JSON object with the parameters.
  !! @param case The case object.
  subroutine scalar_sgs_stats_simcomp_init_from_json(this, json, case)
    class(scalar_sgs_stats_simcomp_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    type(json_file) :: json_subdict
    character(len=:), allocatable :: filename
    character(len=20), allocatable :: fields(:)
    character(len=:), allocatable :: hom_dir
    character(len=:), allocatable :: sname
    character(len=:), allocatable :: name
    real(kind=rp) :: start_time
    type(field_t), pointer :: s
    type(coef_t), pointer :: coef
    character(len=:), allocatable :: alphat_field, nut_field
    real(kind=rp) :: pr_turb
    logical :: nut_dependency

    call json_get_or_default(json, "name", name, "scalar_sgs_stats")
    call this%init_base(json, case)
    call json_get_or_default(json, 'avg_direction', &
         hom_dir, 'none')
    call json_get_or_default(json, 'start_time', &
         start_time, 0.0_rp)
    call json_get_or_default(json, 'field', &
         sname, 's')

    call json_get(json, 'alphat', json_subdict)
    call json_get(json_subdict, 'nut_dependency', nut_dependency)
    if (nut_dependency) then
       call json_get(json_subdict, 'Pr_t', pr_turb)
       call json_get(json_subdict, 'nut_field', nut_field)
    else
       call json_get(json_subdict, 'alphat_field', alphat_field)
    end if

    s => neko_registry%get_field(sname)
    coef => case%fluid%c_Xh
    this%name = name

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", filename)
       if (nut_dependency) then
          call this%init_from_components(s, coef, &
               start_time, hom_dir, nut_field, pr_turb, filename)
       else
          call this%init_from_components(s, coef, &
               start_time, hom_dir, alphat_field, filename)
       end if
    else
       if (nut_dependency) then
          call this%init_from_components(s, coef, &
               start_time, hom_dir, nut_field, pr_turb)
       else
          call this%init_from_components(s, coef, &
               start_time, hom_dir, alphat_field)
       end if
    end if

  end subroutine scalar_sgs_stats_simcomp_init_from_json

  !> Actual constructor using directly the alphat field.
  !! @param s scalar
  !! @param coef sem coefs
  !! @param start_time time to start sampling stats
  !! @param hom_dir directions to average in
  !! @param alphat_field name of the eddy diffusivity field
  !! @param fname name of the outut file
  subroutine scalar_sgs_stats_simcomp_init_from_components_alphat(this, &
       s, coef, start_time, hom_dir, alphat_field, fname)
    class(scalar_sgs_stats_simcomp_t), target, intent(inout) :: this
    character(len=*), intent(in) :: hom_dir
    real(kind=rp), intent(in) :: start_time
    type(field_t), intent(in), target :: s
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: alphat_field
    character(len=*), intent(in), optional :: fname
    character(len=NEKO_FNAME_LEN) :: stats_fname
    character(len=LOG_SIZE) :: log_buf
    character(len=5) :: prefix

    call neko_log%section('scalar SGS stats')
    write(log_buf, '(A,A)') 'Scalar field: ', trim(s%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') 'Eddy diffusivity field: ', trim(alphat_field)
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7)') 'Start time: ', start_time
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') 'Averaging in direction: ', trim(hom_dir)
    call neko_log%message(log_buf)


    call this%stats%init(coef, s, alphat_field)

    this%start_time = start_time
    this%time = start_time
    if (present(fname)) then
       this%default_fname = .false.
       stats_fname = fname
    else
       stats_fname = "scalar_sgs_stats0"
       this%default_fname = .true.
    end if

    call this%stats_output%init(this%stats, this%start_time, &
         hom_dir = hom_dir, name = stats_fname, &
         path = this%case%output_directory)

    call this%case%output_controller%add(this%stats_output, &
         this%output_controller%control_value, &
         this%output_controller%control_mode)

    call neko_log%end_section()

  end subroutine scalar_sgs_stats_simcomp_init_from_components_alphat

  !> Actual constructor using directly the nut field and the turbulent
  !> Prandtl number.
  !! @param s scalar
  !! @param coef sem coefs
  !! @param start_time time to start sampling stats
  !! @param hom_dir directions to average in
  !! @param nut_field name of the eddy diffusivity field
  !! @param pr_turb turbulent Prandtl number
  !! @param fname name of the output file
  subroutine scalar_sgs_stats_simcomp_init_from_components_nut(this, &
       s, coef, start_time, hom_dir, nut_field, pr_turb, fname)
    class(scalar_sgs_stats_simcomp_t), target, intent(inout) :: this
    character(len=*), intent(in) :: hom_dir
    real(kind=rp), intent(in) :: start_time
    type(field_t), intent(in), target :: s
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: nut_field
    real(kind=rp), intent(in) :: pr_turb
    character(len=*), intent(in), optional :: fname
    character(len=NEKO_FNAME_LEN) :: stats_fname
    character(len=LOG_SIZE) :: log_buf
    character(len=5) :: prefix

    call neko_log%section('scalar stats')
    write(log_buf, '(A,A)') 'Scalar field: ', trim(s%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') 'Eddy viscosity field: ', trim(nut_field)
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7)') 'Turbulent Prandtl number: ', pr_turb
    call neko_log%message(log_buf)
    write(log_buf, '(A,E15.7)') 'Start time: ', start_time
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') 'Averaging in direction: ', trim(hom_dir)
    call neko_log%message(log_buf)


    call this%stats%init(coef, s, nut_field, pr_turb)

    this%start_time = start_time
    this%time = start_time
    if (present(fname)) then
       this%default_fname = .false.
       stats_fname = fname
    else
       stats_fname = "scalar_sgs_stats0"
       this%default_fname = .true.
    end if

    call this%stats_output%init(this%stats, this%start_time, &
         hom_dir = hom_dir, name = stats_fname, &
         path = this%case%output_directory)

    call this%case%output_controller%add(this%stats_output, &
         this%output_controller%control_value, &
         this%output_controller%control_mode)

    call neko_log%end_section()

  end subroutine scalar_sgs_stats_simcomp_init_from_components_nut

  !> Destructor.
  subroutine scalar_sgs_stats_simcomp_free(this)
    class(scalar_sgs_stats_simcomp_t), intent(inout) :: this
    call this%free_base()
    call this%stats%free()
    call this%stats_output%free()

    call this%free_amr_base()

  end subroutine scalar_sgs_stats_simcomp_free

  subroutine scalar_sgs_stats_simcomp_restart(this, time)
    class(scalar_sgs_stats_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    character(len=NEKO_FNAME_LEN) :: fname
    character(len=5) :: prefix, suffix
    integer :: last_slash_pos
    real(kind=rp) :: t
    t = time%t
    if (t .gt. this%time) this%time = t
    if (this%default_fname) then
       fname = this%stats_output%file_%get_base_fname()
       write (prefix, '(I5)') this%stats_output%file_%get_counter()
       call filename_suffix(fname, suffix)
       last_slash_pos = &
            filename_tslash_pos(fname)
       if (last_slash_pos .ne. 0) then
          fname = &
               trim(fname(1:last_slash_pos))// &
               "scalar_sgs_stats"//trim(adjustl(prefix))//"."//suffix
       else
          fname = "scalar_sgs_stats"// &
               trim(adjustl(prefix))//"."//suffix
       end if
       call this%stats_output%init_base(fname)
    end if
  end subroutine scalar_sgs_stats_simcomp_restart

  !> scalar_sgs_stats, called depending on compute_control and compute_value
  !! @param time The current time info
  subroutine scalar_sgs_stats_simcomp_compute(this, time)
    class(scalar_sgs_stats_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: delta_t, t
    real(kind=rp) :: sample_start_time, sample_time
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr

    if (time%start_time .gt. this%start_time) then
       write(log_buf, '(A)') 'Simulation start time is later than the ' &
            // 'scalar stats start time.'
       call neko_log%warning(log_buf)
       write(log_buf, '(A,E15.7)') 'Simulation start time:', time%start_time
       call neko_log%warning(log_buf)
       write(log_buf, '(A,E15.7)') 'scalar stats start time:', this%start_time
       call neko_log%warning(log_buf)
       write(log_buf, '(A)') 'Assigning the statistics start time to ' &
            // 'the simulation start time.'
       call neko_log%warning(log_buf)
       this%start_time = time%start_time
       this%time = time%start_time
    end if

    t = time%t

    if (t .ge. this%start_time) then
       delta_t = t - this%time !This is only a real number

       call MPI_Barrier(NEKO_COMM, ierr)

       sample_start_time = MPI_WTIME()

       call this%stats%update(delta_t)
       call MPI_Barrier(NEKO_COMM, ierr)
       this%time = t

       sample_time = MPI_WTIME() - sample_start_time

       call neko_log%section('scalar stats')
       write(log_buf, '(A,E15.7)') 'Sampling at time:', t
       call neko_log%message(log_buf)
       write(log_buf, '(A33,E15.7)') 'Simulation time since last sample:', &
            delta_t
       call neko_log%message(log_buf)
       write(log_buf, '(A,E15.7)') 'Sampling time (s):', sample_time
       call neko_log%message(log_buf)
       call neko_log%end_section()
    end if

  end subroutine scalar_sgs_stats_simcomp_compute

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine scalar_sgs_stats_simcomp_amr_restart(this, reconstruct, counter, &
       tstep)
    class(scalar_sgs_stats_simcomp_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep
    character(len=LOG_SIZE) :: log_buf

    call neko_error('Nothing done for AMR reconstruction')

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    log_buf = 'Scalar sgs stats'
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)
!    call neko_log%section(log_buf, NEKO_LOG_VERBOSE)
!    call neko_log%end_section(lvl = NEKO_LOG_VERBOSE)

  end subroutine scalar_sgs_stats_simcomp_amr_restart

end module scalar_sgs_stats_simcomp
