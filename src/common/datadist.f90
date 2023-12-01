! Copyright (c) 2019-2021, The Neko Authors
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
!> Defines practical data distributions
module datadist
  use mpi_f08, only : MPI_Comm
  implicit none
  private

  type dist_t
     type(MPI_Comm) :: comm !< Communicator on which the dist. is defined
     integer :: pe_rank     !< Pe's rank in the given distribution
     integer :: pe_size     !< Size of communicator in the given dist.
     integer :: L               
     integer :: R
     integer :: M !< Total, global, size
     integer :: Ip !< Number of local values on this process
  end type dist_t

  !> Load-balanced linear distribution \f$ M = PL + R \f$
  type, extends(dist_t) :: linear_dist_t
   contains 
     procedure :: num_local => linear_dist_Ip
     procedure :: num_global => linear_dist_M
     procedure :: start_idx => linear_dist_start
     procedure :: end_idx => linear_dist_end
  end type linear_dist_t

  interface linear_dist_t
     module procedure linear_dist_init
  end interface linear_dist_t

  public :: linear_dist_t

contains
  
  function linear_dist_init(n, rank, size, comm) result(this)
    integer, intent(in) :: n    !< Total size
    integer :: rank             !< PE's rank to define the dist. over
    integer :: size             !< Size of comm where the dist. is def. on
    type(MPI_Comm) :: comm      !< comm. to define the dist. over
    type(linear_dist_t), target :: this

    this%M = n
    this%comm = comm
    this%pe_rank = rank
    this%pe_size = size
    
    this%L = floor(dble(this%M) / dble(this%pe_size))
    this%R = modulo(this%M, this%pe_size)
    this%Ip = floor((dble(this%M) + dble(this%pe_size) - &
         dble(this%pe_rank) - 1d0) / dble(this%pe_size))    
  end function linear_dist_init

  pure function linear_dist_Ip(this) result(n)
    class(linear_dist_t), intent(in) :: this
    integer :: n
    n = this%Ip
  end function linear_dist_Ip

  pure function linear_dist_M(this) result(M)
    class(linear_dist_t), intent(in) :: this
    integer :: M
    M = this%M
  end function linear_dist_M

  pure function linear_dist_start(this) result(start)
    class(linear_dist_t), intent(in) :: this
    integer :: start
    start = this%pe_rank * this%L + min(this%pe_rank, this%R)
  end function linear_dist_start

  function linear_dist_end(this) result(end)
    class(linear_dist_t), intent(inout) :: this
    integer :: end
    end = linear_dist_start(this) + (this%Ip - 1)
  end function linear_dist_end
end module datadist
