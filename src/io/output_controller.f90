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
!> Implements `output_controller_t`
module output_controller
  use output, only: output_t, output_ptr_t
  use fld_file, only: fld_file_t
  use comm
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use profiler, only : profiler_start_region, profiler_end_region
  use num_types, only : rp, dp
  use time_based_controller, only : time_based_controller_t
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
     !> Final time of the simulation.
     real(kind=rp) :: time_end
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
  !! @param time_end The end time of thesimulation.
  !! @param size The number of controllers to allocate for. Optional, defaults
  !! to 1.
  subroutine output_controller_init(this, time_end, size)
    class(output_controller_t), intent(inout) :: this
    integer, intent(in), optional :: size
    real(kind=rp), intent(in) :: time_end
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
  subroutine output_controller_add(this, out, write_par, write_control)
    class(output_controller_t), intent(inout) :: this
    class(output_t), intent(inout), target :: out
    real(kind=rp), intent(in) :: write_par
    character(len=*), intent(in) :: write_control
    type(output_ptr_t), allocatable :: tmp(:)
    type(time_based_controller_t), allocatable :: tmp_ctrl(:)
    character(len=LOG_SIZE) :: log_buf
    integer :: n
    class(*), pointer :: ft

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
       call this%controllers(n)%init(this%time_end, write_control, write_par)
    end if

    ! The code below only prints to console
    call neko_log%section('Adding write output')
    call neko_log%message('File name        : '// &
          trim(this%output_list(this%n)%ptr%file_%file_type%fname))
    call neko_log%message('Write control    : '//trim(write_control))

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
       write(log_buf, '(A,I13)') 'Total samples: ',  int(write_par)
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'Writes per time unit (Freq.): ',  &
             this%controllers(n)%frequency
       call neko_log%message(log_buf)
       write(log_buf, '(A,ES13.6)') 'Time between writes: ', &
          this%controllers(n)%time_interval
       call neko_log%message(log_buf)
    else if (trim(write_control) .eq. 'tsteps') then
       write(log_buf, '(A,I13)') 'Time step interval: ',  int(write_par)
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
  subroutine output_controller_execute(this, t, tstep, ifforce)
    class(output_controller_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    logical, intent(in), optional :: ifforce
    real(kind=dp) :: sample_start_time, sample_end_time
    real(kind=dp) :: sample_time
    character(len=LOG_SIZE) :: log_buf
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
          if (this%controllers(i)%check(t, tstep, force)) then
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
          if (this%controllers(i)%check(t, tstep, force)) then
             call neko_log%message('File name     : '// &
                  trim(samp%output_list(i)%ptr%file_%file_type%fname))

             write(log_buf, '(A,I6)') 'Output number :', &
                  int(this%controllers(i)%nexecutions)
             call neko_log%message(log_buf)

             call samp%output_list(i)%ptr%sample(t)

             call this%controllers(i)%register_execution()
          end if
       end do
    class default
       call neko_error('Invalid output_controller output list')
    end select

    call MPI_Barrier(NEKO_COMM, ierr)
    sample_end_time = MPI_WTIME()

    sample_time = sample_end_time - sample_start_time
    if (write_output) then
       write(log_buf, '(A16,1x,F10.6,A,F9.6)') 'Writing at time:', t, &
            ' Output time (s): ', sample_time
       call neko_log%message(log_buf)
       call neko_log%end_section()
    end if
    call profiler_end_region('Output controller', 22)
  end subroutine output_controller_execute

  !> Set write counter based on time (after restart)
  !> @param t Time value.
  subroutine output_controller_set_counter(this, t)
    class(output_controller_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i, nexecutions


    do i = 1, this%n
       if (this%controllers(i)%nsteps .eq. 0) then
          nexecutions = int(t / this%controllers(i)%time_interval) + 1
          this%controllers(i)%nexecutions = nexecutions

          call this%output_list(i)%ptr%set_counter(nexecutions)
          call this%output_list(i)%ptr%set_start_counter(nexecutions)
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
