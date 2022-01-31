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
  end type sampler_t

contains

  !> Initialize a sampler
  subroutine sampler_init(this, nsamp, T_end, size)
    class(sampler_t), intent(inout) :: this
    integer, intent(inout) :: nsamp
    real(kind=rp) :: T_end
    integer, intent(inout), optional :: size
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
    this%freq = ( real(nsamp,rp) / T_end )
    this%T = real(1d0,rp) / this%freq
    
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
  subroutine sampler_sample(this, t)
    class(sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    real(kind=dp) :: sample_start_time, sample_end_time
    character(len=LOG_SIZE) :: log_buf
    integer :: i

    if (t .ge. (this%nsample * this%T)) then

       sample_start_time = MPI_WTIME()
       do i = 1, this%n
          call this%output_list(i)%outp%sample(t)
       end do
       sample_end_time = MPI_WTIME()
       this%nsample = this%nsample + 1

       write(log_buf,'(a23,1x,e15.7,A,F8.4)') 'Sampling fields at time:', t, &
             ' Sample time (s): ', sample_end_time - sample_start_time
       call neko_log%message(log_buf)

    end if
    
  end subroutine sampler_sample

  !> Set sampling counter (after restart)
  subroutine sampler_set_counter(this, t)
    class(sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer :: i

    this%nsample = int(t / this%T) + 1

    do i = 1, this%n
       call this%output_list(i)%outp%set_counter(this%nsample)
    end do
    
  end subroutine sampler_set_counter
  
end module sampler
