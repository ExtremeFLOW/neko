! Copyright (c) 2020-2021, The Neko Authors
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
!> Defines a sampler 
module sampler
  use num_types
  use output
  use comm
  use logger
  use utils
  use profiler
  implicit none
  private

  !> Pointer to an arbitrary output
  type, private :: outp_t
     class(output_t), pointer :: outp
  end type outp_t

  !> Sampler
  type, public :: sampler_t
     type(outp_t), allocatable :: output_list(:) !< List of outputs
     real(kind=rp), allocatable :: freq_list(:) !< Write frequency (sim time)
     real(kind=rp), allocatable :: T_list(:) !< List of times before write
     integer, allocatable :: tstep_interval_list(:) !< Sample every tstep interval
     integer, allocatable :: nsample_list(:) !< Sample number
     integer :: n !< number of outputs
     integer :: size !< number of entries in list
     real(kind=rp) :: T_end
   contains
     procedure, pass(this) :: init => sampler_init
     procedure, pass(this) :: free => sampler_free
     procedure, pass(this) :: add => sampler_add
     procedure, pass(this) :: sample => sampler_sample
     procedure, pass(this) :: set_counter => sampler_set_counter
  end type sampler_t

contains

  !> Initialize a sampler
  subroutine sampler_init(this, T_end, size)
    class(sampler_t), intent(inout) :: this
    integer, intent(in), optional :: size
    real(kind=rp), intent(in) :: T_end
    character(len=LOG_SIZE) :: log_buf
    integer :: n, i


    call this%free()

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(this%output_list(n))
    allocate(this%freq_list(n))
    allocate(this%T_list(n))
    allocate(this%tstep_interval_list(n))
    allocate(this%nsample_list(n))

    do i = 1, n
       this%output_list(i)%outp => null()
       this%nsample_list(i) = 0
       this%tstep_interval_list(i) = 0
       this%freq_list(i) = 0
       this%T_list(i) = 0
    end do

    this%size = n
    this%n = 0
    this%T_end = T_end

  end subroutine sampler_init

  !> Deallocate a sampler
  subroutine sampler_free(this)
    class(sampler_t), intent(inout) :: this

    if (allocated(this%output_list)) then
       deallocate(this%output_list)       
    end if
    if (allocated(this%freq_list)) then
       deallocate(this%freq_list)       
    end if
    if (allocated(this%T_list)) then
       deallocate(this%T_list)       
    end if
    if (allocated(this%tstep_interval_list)) then
       deallocate(this%tstep_interval_list)       
    end if
    if (allocated(this%nsample_list)) then
       deallocate(this%nsample_list)       
    end if

    this%n = 0
    this%size = 0

  end subroutine sampler_free

  !> Add an output @a out to the sampler
  subroutine sampler_add(this, out, write_par, write_control)
    class(sampler_t), intent(inout) :: this
    class(output_t), intent(inout), target :: out
    real(kind=rp), intent(in) :: write_par
    character(len=*), intent(in) :: write_control
    type(outp_t), allocatable :: tmp(:)
    integer, allocatable :: tmpi(:)
    real(kind=rp), allocatable :: tmpr(:)
    character(len=LOG_SIZE) :: log_buf
    integer :: n

    if (this%n .ge. this%size) then
       allocate(tmp(this%size * 2))
       tmp(1:this%size) = this%output_list
       call move_alloc(tmp, this%output_list)
       allocate(tmpr(this%size * 2))
       tmpr(1:this%size) = this%freq_list
       call move_alloc(tmpr, this%freq_list)
       allocate(tmpr(this%size * 2))
       tmpr(1:this%size) = this%T_list
       call move_alloc(tmpr, this%T_list)
       allocate(tmpi(this%size * 2))
       tmpi(1:this%size) = this%tstep_interval_list
       call move_alloc(tmpi, this%tstep_interval_list)
       allocate(tmpi(this%size * 2))
       tmpi(1:this%size) = this%nsample_list
       call move_alloc(tmpi, this%nsample_list)
       this%size = this%size * 2
    end if

    this%n = this%n + 1
    n = this%n
    this%output_list(this%n)%outp => out
    
    call neko_log%section('Adding write output')
    call neko_log%message('File name: '// &
          trim(this%output_list(this%n)%outp%file_%file_type%fname))
    call neko_log%message( 'Write control: '//trim(write_control))
    if (trim(write_control) .eq. 'simulationtime') then
        this%T_list(n) = write_par
        this%freq_list(n) = 1.0_rp/this%T_list(n)
        this%tstep_interval_list(n) = 0
        write(log_buf, '(A,ES13.6)') 'Writes per time unit (Freq.): ',  this%freq_list(n)
        call neko_log%message(log_buf)
        write(log_buf, '(A,ES13.6)') 'Time between writes: ',  this%T_list(n)
        call neko_log%message(log_buf)
    else if (trim(write_control) .eq. 'nsamples') then
        this%freq_list(n) = ( write_par / this%T_end )
        this%T_list(n) = real(1d0,rp) / this%freq_list(n)
        this%tstep_interval_list(n) = 0
        write(log_buf, '(A,I13)') 'Total samples: ',  int(write_par)
        call neko_log%message(log_buf)
        write(log_buf, '(A,ES13.6)') 'Writes per time unit (Freq.): ',  this%freq_list(n)
        call neko_log%message(log_buf)
        write(log_buf, '(A,ES13.6)') 'Time between writes: ',  this%T_list(n)
        call neko_log%message(log_buf)
    else if (trim(write_control) .eq. 'tsteps') then 
        this%tstep_interval_list(n) = int(write_par)
        write(log_buf, '(A,I13)') 'Time step interval: ',  int(write_par)
        call neko_log%message(log_buf)
    else if (trim(write_control) .eq. 'org') then
        this%tstep_interval_list(n) = this%tstep_interval_list(1)
        this%freq_list(n) = this%freq_list(1)
        this%T_list(n) = this%T_list(1)
        write(log_buf, '(A)') 'Write control not set, defaulting to first output settings'
        call neko_log%message(log_buf)
    end if
    this%nsample_list(n) = 0
    
    call neko_log%end_section()
  end subroutine sampler_add

  !> Sample all outputs in the sampler
  subroutine sampler_sample(this, t, tstep, ifforce)
    class(sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    logical, intent(in), optional :: ifforce
    real(kind=dp) :: sample_start_time, sample_end_time
    real(kind=dp) :: sample_time
    character(len=LOG_SIZE) :: log_buf
    integer :: i, ierr
    logical :: force = .false.

    if (present(ifforce)) then
       force = ifforce
    end if

    call profiler_start_region('Sampler')
    !Do we need this Barrier?
    call MPI_Barrier(NEKO_COMM, ierr)
    sample_start_time = MPI_WTIME()

    ! We should not need this extra select block, and it works great
    ! without it for GNU, Intel and NEC, but breaks horribly on Cray         
    call neko_log%message('Writer output: ')
    ! (>11.0.x) when using high opt. levels.  
    select type (samp => this)
    type is (sampler_t)
       do i = 1, this%n
          if (force) then
             call neko_log%message(trim(this%output_list(i)%outp%file_%file_type%fname))
             write(log_buf, '(A,I6)') 'Output number:',  int(this%nsample_list(i))
             call neko_log%message(log_buf)
             call samp%output_list(i)%outp%sample(t)
             this%nsample_list(i) = this%nsample_list(i) + 1
          else if ((this%tstep_interval_list(i) .eq. 0) .and. (t .ge. (this%nsample_list(i) * this%T_list(i)))) then
             call neko_log%message(trim(this%output_list(i)%outp%file_%file_type%fname))
             write(log_buf, '(A,I6)') 'Output number:',  int(this%nsample_list(i))
             call neko_log%message(log_buf)
             call samp%output_list(i)%outp%sample(t)
             this%nsample_list(i) = this%nsample_list(i) + 1
          else if ((this%tstep_interval_list(i) .gt. 0) .and. (mod(tstep,this%tstep_interval_list(i)) .eq. 0)) then
             call neko_log%message(trim(this%output_list(i)%outp%file_%file_type%fname))
             write(log_buf, '(A,I6)') 'Output number:',  int(this%nsample_list(i))
             call neko_log%message(log_buf)
             call samp%output_list(i)%outp%sample(t)
             this%nsample_list(i) = this%nsample_list(i) + 1
          end if
       end do
    class default
       call neko_error('Invalid sampler output list')
    end select
    
    call MPI_Barrier(NEKO_COMM, ierr)
    sample_end_time = MPI_WTIME()

    sample_time = sample_end_time - sample_start_time
    write(log_buf,'(A23,1x,F10.6,A,F9.6)') 'Output written at time:', t, &
          ' Output time (s): ', sample_time
    call neko_log%message(log_buf)
    call profiler_end_region
  end subroutine sampler_sample

  !> Set sampling counter based on time (after restart)
  subroutine sampler_set_counter(this, t)
    class(sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i


    do i = 1, this%n
       if (this%tstep_interval_list(i) .eq. 0) then
          this%nsample_list(i) = int(t / this%T_list(i)) + 1
          call this%output_list(i)%outp%set_counter(this%nsample_list(i))
          call this%output_list(i)%outp%set_start_counter(this%nsample_list(i))
       end if
    end do
    
  end subroutine sampler_set_counter
 
  !> Set sampling counter (after restart) explicitly
  subroutine sampler_set_sample_count(this, sample_number)
    class(sampler_t), intent(inout) :: this
    integer, intent(in) :: sample_number
    integer :: i

    do i = 1, this%n
       this%nsample_list(i) = sample_number
       call this%output_list(i)%outp%set_counter(this%nsample_list(i))
       call this%output_list(i)%outp%set_start_counter(this%nsample_list(i))
    end do
    
  end subroutine sampler_set_sample_count
  
 
end module sampler
