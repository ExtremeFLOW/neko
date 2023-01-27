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
     type(outp_t), allocatable :: output_list(:)
     integer :: n              
     integer :: size
     integer :: nsample
     real(kind=rp) :: freq
     real(kind=rp) :: T
   contains
     procedure, pass(this) :: init => sampler_init
     procedure, pass(this) :: free => sampler_free
     procedure, pass(this) :: add => sampler_add
     procedure, pass(this) :: sample => sampler_sample
     procedure, pass(this) :: set_counter => sampler_set_counter
     procedure, pass(this) :: set_count => sampler_set_sample_count
  end type sampler_t

contains

  !> Initialize a sampler
  subroutine sampler_init(this, nsamp, T_end, freq, T, size)
    class(sampler_t), intent(inout) :: this
    integer, intent(in) :: nsamp
    real(kind=rp), intent(in) :: T_end
    real(kind=rp), optional, intent(in) :: T
    integer, intent(in), optional :: size
    real(kind=rp), optional, intent(in) :: freq
    character(len=LOG_SIZE) :: log_buf
    integer :: n, i


    call this%free()

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(this%output_list(n))

    do i = 1, n
       this%output_list(i)%outp => null()
    end do

    this%n = 0
    this%size = n

    this%nsample = 0
    if (present(freq)) then
       this%freq = freq
       this%T = real(1d0,rp) / this%freq
    else if (present(T)) then
       this%T = T
       this%freq = 1.0_rp/this%T
    else
       this%freq = ( real(nsamp,rp) / T_end )
       this%T = real(1d0,rp) / this%freq
    end if

    call neko_log%section('Sampler')
    write(log_buf, '(A,I13)') 'Samples   :',  nsamp
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'Freq.     :',  this%freq
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'Sample t  :',  this%T
    call neko_log%message(log_buf)
    call neko_log%end_section()
    
  end subroutine sampler_init

  !> Deallocate a sampler
  subroutine sampler_free(this)
    class(sampler_t), intent(inout) :: this

    if (allocated(this%output_list)) then
       deallocate(this%output_list)       
    end if

    this%n = 0
    this%size = 0

  end subroutine sampler_free

  !> Add an output @a out to the sampler
  subroutine sampler_add(this, out)
    class(sampler_t), intent(inout) :: this
    class(output_t), intent(inout), target :: out
    type(outp_t), allocatable :: tmp(:)

    if (this%n .ge. this%size) then
       allocate(tmp(this%size * 2))
       tmp(1:this%size) = this%output_list
       call move_alloc(tmp, this%output_list)
       this%size = this%size * 2
    end if

    this%n = this%n + 1
    this%output_list(this%n)%outp => out
    
  end subroutine sampler_add

  !> Sample all outputs in the sampler
  subroutine sampler_sample(this, t, ifforce)
    class(sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    logical, intent(in), optional :: ifforce
    real(kind=dp) :: sample_start_time, sample_end_time
    real(kind=dp) :: sample_time
    character(len=LOG_SIZE) :: log_buf
    integer :: i, ierr
    logical :: force = .false.

    if (present(ifforce)) then
       force = ifforce
    end if

    if (t .ge. (this%nsample * this%T) .or. force) then
       call profiler_start_region('Sampler')
       call MPI_Barrier(NEKO_COMM, ierr)
       sample_start_time = MPI_WTIME()

       ! We should not need this extra select block, and it works great
       ! without it for GNU, Intel and NEC, but breaks horribly on Cray         
       ! (>11.0.x) when using high opt. levels.  
       select type (samp => this)
       type is (sampler_t)
          do i = 1, this%n
             call samp%output_list(i)%outp%sample(t)
          end do
       class default
          call neko_error('Invalid sampler output list')
       end select
       
       call MPI_Barrier(NEKO_COMM, ierr)
       sample_end_time = MPI_WTIME()
       this%nsample = this%nsample + 1

       sample_time = sample_end_time - sample_start_time
       write(log_buf,'(A24,1x,F10.6,A,F9.6)') 'Sampling fields at time:', t, &
             ' Sample time (s): ', sample_time
       call neko_log%message(log_buf)
       call profiler_end_region
    end if
    
  end subroutine sampler_sample

  !> Set sampling counter based on time (after restart)
  subroutine sampler_set_counter(this, t)
    class(sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    this%nsample = int(t / this%T) + 1

    do i = 1, this%n
       call this%output_list(i)%outp%set_counter(this%nsample)
       call this%output_list(i)%outp%set_start_counter(this%nsample)
    end do
    
  end subroutine sampler_set_counter
 
  !> Set sampling counter (after restart) explicitly
  subroutine sampler_set_sample_count(this, sample_number)
    class(sampler_t), intent(inout) :: this
    integer, intent(in) :: sample_number
    integer :: i

    this%nsample = sample_number
    do i = 1, this%n
       call this%output_list(i)%outp%set_counter(this%nsample)
       call this%output_list(i)%outp%set_start_counter(this%nsample)
    end do
    
  end subroutine sampler_set_sample_count
  
 
end module sampler
