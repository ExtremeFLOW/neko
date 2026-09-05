! Copyright (c) 2020-2026, The Neko Authors
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
!> Implements `output_controller_t`
module output_controller
  use output, only : output_t, output_ptr_t
  use fld_file, only : fld_file_t
  use comm
  use time_state, only : time_state_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use profiler, only : profiler_start_region, profiler_end_region
  use num_types, only : dp
  use time_based_controller, only : time_based_controller_t
  use mpi_f08, only : MPI_WTIME, MPI_Barrier
  implicit none
  private


  !> Centralized controller for a list of outputs.
  !! @details Holds a list of `output_t` and corresponding
  !! `time_based_controller_t`s. Uses the latter to determine, which outputs
  !! need to be sampled and written to disk at a given time step.
  type, public :: output_controller_t
     !> List of outputs.
     type(output_ptr_t), allocatable :: output_list(:)
     !> List of controllers for determining wether we should write.
     type(time_based_controller_t), allocatable :: controllers(:)
     !> Number of outputs.
     integer :: n
     !> Number of entries in the list.
     integer :: size
     !> Start time of the simulation, the default anchor of the schedules.
     real(kind=dp) :: time_start
     !> Whether outputs write the initial state of the simulation by default.
     logical :: write_at_start = .true.
     !> Final time of the simulation.
     real(kind=dp) :: time_end
   contains
     !> Constructor.
     procedure, pass(this) :: init => output_controller_init
     !> Destructor.
     procedure, pass(this) :: free => output_controller_free
     !> Add an output to the controller.
     procedure, pass(this) :: add => output_controller_add
     !> Sample the fields and output.
     procedure, pass(this) :: execute => output_controller_execute
     !> Set output counter based on time (after restart)
     procedure, pass(this) :: set_counter => output_controller_set_counter
  end type output_controller_t

contains

  !> Constructor.
  !! @param time_end The end time of the simulation.
  !! @param size The number of controllers to allocate for. Optional, defaults
  !! to 1.
  !! @param time_start The start time of the simulation. Outputs that do not
  !! have a start time of their own do not write before it. Optional, defaults
  !! to 0.
  !! @param write_at_start Whether outputs write the initial state of the
  !! simulation. The default for the outputs added afterwards, which each of
  !! them can override. Optional, defaults to `.true.`.
  subroutine output_controller_init(this, time_end, size, time_start, &
       write_at_start)
    class(output_controller_t), intent(inout) :: this
    integer, intent(in), optional :: size
    real(kind=dp), intent(in) :: time_end
    real(kind=dp), intent(in), optional :: time_start
    logical, intent(in), optional :: write_at_start
    character(len=LOG_SIZE) :: log_buf
    integer :: n, i

    call this%free()

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(this%output_list(n))
    allocate(this%controllers(n))

    do i = 1, n
       this%output_list(i)%ptr => null()
    end do

    this%size = n
    this%n = 0
    this%time_end = time_end

    if (present(time_start)) then
       this%time_start = time_start
    else
       this%time_start = 0.0_dp
    end if

    if (present(write_at_start)) then
       this%write_at_start = write_at_start
    else
       this%write_at_start = .true.
    end if

  end subroutine output_controller_init

  !> Destructor.
  subroutine output_controller_free(this)
    class(output_controller_t), intent(inout) :: this

    if (allocated(this%output_list)) then
       deallocate(this%output_list)
    end if
    if (allocated(this%controllers)) then
       deallocate(this%controllers)
    end if

    this%n = 0
    this%size = 0

  end subroutine output_controller_free

  !> Add an output @a out to the controller
  !! @param out The output to add.
  !! @param write_par The output frequency value, in accordance with
  !! `write_control`.
  !! @param write_control Determines the meaning of `write_par`. Accepts the
  !! usual list of control options.
  !! @param start_time When to start writing the output. Also the time the
  !! output schedule is anchored to. Optional, defaults to the start time of
  !! the simulation.
  !! @param write_at_start Whether to write at `start_time` itself. Should be
  !! `.false.` for outputs for which a write at the very first time step is
  !! meaningless, such as checkpoints and running statistics. Optional,
  !! defaults to what the controller was constructed with.
  subroutine output_controller_add(this, out, write_par, write_control, &
       start_time, write_at_start)
    class(output_controller_t), intent(inout) :: this
    class(output_t), intent(inout), target :: out
    real(kind=dp), intent(in) :: write_par
    character(len=*), intent(in) :: write_control
    real(kind=dp), optional, intent(in) :: start_time
    logical, optional, intent(in) :: write_at_start
    real(kind=dp) :: start_time_
    logical :: write_at_start_
    type(output_ptr_t), allocatable :: tmp(:)
    type(time_based_controller_t), allocatable :: tmp_ctrl(:)
    character(len=LOG_SIZE) :: log_buf
    integer :: n, nexecutions
    class(*), pointer :: ft

    if (present(start_time)) then
       start_time_ = start_time
    else
       start_time_ = this%time_start
    end if

    if (present(write_at_start)) then
       write_at_start_ = write_at_start
    else
       write_at_start_ = this%write_at_start
    end if


    if (this%n .ge. this%size) then
       allocate(tmp(this%size * 2))
       tmp(1:this%size) = this%output_list
       call move_alloc(tmp, this%output_list)

       allocate(tmp_ctrl(this%size * 2))
       tmp_ctrl(1:this%size) = this%controllers
       call move_alloc(tmp_ctrl, this%controllers)

       this%size = this%size * 2
    end if


    this%n = this%n + 1
    n = this%n
    this%output_list(this%n)%ptr => out

    if (trim(write_control) .eq. "org") then
       this%controllers(n) = this%controllers(1)
    else
       call this%controllers(n)%init(start_time_, this%time_end, &
            write_control, write_par, write_at_start_, &
            direction = this%time_end - this%time_start)
    end if

    ! The code below only prints to console
    call neko_log%section('Adding write output')
    call neko_log%message('File name        : '// &
         trim(this%output_list(this%n)%ptr%file_%file_type%get_fname()))
    call neko_log%message('Write control    : '//trim(write_control))
    if (.not. this%controllers(n)%never) then
       write(log_buf, '(A,ES13.6)') 'First write at   : ', &
            this%controllers(n)%next_time()
       if (trim(write_control) .ne. 'tsteps') call neko_log%message(log_buf)
    end if

    ! Show the output precision if we are outputting an fld file
    select type (ft => out%file_%file_type)
    type is (fld_file_t)
       if (ft%dp_precision) then
          call neko_log%message('Output precision : double')
       else
          call neko_log%message('Output precision : single')
       end if
    end select

    if (trim(write_control) .eq. 'simulationtime') then
       write(log_buf, '(A,ES13.6)') 'Writes per time unit (Freq.): ', &
            this%controllers(n)%frequency
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'Time between writes: ', &
            this%controllers(n)%time_interval
       call neko_log%message(log_buf)
    else if (trim(write_control) .eq. 'nsamples') then
       write(log_buf, '(A,I13)') 'Total samples: ', int(write_par)
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'Writes per time unit (Freq.): ', &
            this%controllers(n)%frequency
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'Time between writes: ', &
            this%controllers(n)%time_interval
       call neko_log%message(log_buf)
    else if (trim(write_control) .eq. 'tsteps') then
       write(log_buf, '(A,I13)') 'Time step interval: ', int(write_par)
       call neko_log%message(log_buf)
    else if (trim(write_control) .eq. 'org') then
       write(log_buf, '(A)') &
            'Write control not set, defaulting to first output settings'
       call neko_log%message(log_buf)
    end if

    call neko_log%end_section()
  end subroutine output_controller_add

  !> Query each of the `controllers` whether it is time to write, and if so,
  !! do so for the corresponding output.
  !! @param t The time value.
  !! @param tstep The current time-stepper iteration.
  !! @param ifforce Whether to force a write. Optional, defaults to 0.
  subroutine output_controller_execute(this, time, ifforce)
    class(output_controller_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    logical, intent(in), optional :: ifforce
    real(kind=dp) :: sample_start_time, sample_end_time
    real(kind=dp) :: sample_time
    character(len=LOG_SIZE) :: log_buf
    character(len=1024) :: output_fname
    integer :: i, ierr
    logical :: force, write_output, write_output_test

    if (present(ifforce)) then
       force = ifforce
    else
       force = .false.
    end if

    call profiler_start_region('Output controller', 22)
    !Do we need this Barrier?
    call MPI_Barrier(NEKO_COMM, ierr)
    sample_start_time = MPI_WTIME()

    write_output = .false.
    ! Determine if at least one output needs to be written
    ! We should not need this extra select block, and it works great
    ! without it for GNU, Intel and NEC, but breaks horribly on Cray
    ! (>11.0.x) when using high opt. levels.
    select type (samp => this)
    type is (output_controller_t)
       do i = 1, samp%n
          if (this%controllers(i)%check(time, force)) then
             write_output = .true.
             exit
          end if
       end do
    end select

    if (write_output) then
       call neko_log%section('Writer output ')
    end if

    ! Loop through the outputs and write if necessary.
    ! We should not need this extra select block, and it works great
    ! without it for GNU, Intel and NEC, but breaks horribly on Cray
    ! (>11.0.x) when using high opt. levels.
    select type (samp => this)
    type is (output_controller_t)
       do i = 1, this%n
          if (this%controllers(i)%check(time, force)) then
             output_fname = &
                  samp%output_list(i)%ptr%file_%file_type% &
                  get_next_output_fname()
             call neko_log%message('File name     : '//trim(output_fname))

             write(log_buf, '(A,I6)') 'Output number :', &
                  int(this%controllers(i)%nexecutions)
             call neko_log%message(log_buf)

             call samp%output_list(i)%ptr%sample(time%t)

             call this%controllers(i)%register_execution(time)
          end if
       end do
    class default
       call neko_error('Invalid output_controller output list')
    end select

    call MPI_Barrier(NEKO_COMM, ierr)
    sample_end_time = MPI_WTIME()

    sample_time = sample_end_time - sample_start_time
    if (write_output) then
       write(log_buf, '(A16,1x,F12.6,A,F9.6)') 'Writing at time:', time%t, &
            ' Output time (s): ', sample_time
       call neko_log%message(log_buf)
       call neko_log%end_section()
    end if
    call profiler_end_region('Output controller', 22)
  end subroutine output_controller_execute

  !> Set write counter based on time (after restart)
  !> @param time Current time info.
  subroutine output_controller_set_counter(this, time)
    class(output_controller_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    character(len=LOG_SIZE) :: log_buf
    character(len=1024) :: output_fname
    integer :: i, nexecutions
    logical :: file_exists

    do i = 1, this%n
       call this%controllers(i)%set_counter(time)
       ! A step based schedule cannot be reconstructed from the time alone,
       ! so the file counter of such an output is left alone.
       if (this%controllers(i)%nsteps .eq. 0) then
          nexecutions = this%controllers(i)%nexecutions
          call this%output_list(i)%ptr%set_counter(-1)
          call this%output_list(i)%ptr%set_start_counter(nexecutions)
       end if

       ! The file counter is where the schedule says the run has got to, so
       ! it points past the files the previous run wrote. It does not when
       ! the run is repeating an interval it has already covered, or when the
       ! output frequency was changed, and the files of the previous run are
       ! then written over.
       if (this%controllers(i)%never) cycle
       file_exists = .false.
       output_fname = &
            this%output_list(i)%ptr%file_%file_type%get_next_output_fname()
       if (pe_rank .eq. 0) then
          inquire(file = trim(output_fname), exist = file_exists)
       end if
       if (file_exists) then
          write(log_buf, '(A,A)') 'Overwriting from: ', trim(output_fname)
          call neko_log%warning(log_buf)
       end if
    end do

  end subroutine output_controller_set_counter

  !> Set write counter (after restart) explicitly
  !> @param counter The value of the write coutner to be set.
  subroutine output_controller_set_write_count(this, counter)
    class(output_controller_t), intent(inout) :: this
    integer, intent(in) :: counter
    integer :: i

    do i = 1, this%n
       this%controllers(i)%nexecutions = counter
       call this%output_list(i)%ptr%set_counter(counter)
       call this%output_list(i)%ptr%set_start_counter(counter)
    end do

  end subroutine output_controller_set_write_count


end module output_controller
