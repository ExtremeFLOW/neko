! Copyright (c) 2021, The Neko Authors
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
!> Defines a gather-scatter backend
module gs_bcknd
  use num_types
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none

  integer, public, parameter :: GS_BCKND_CPU = 1, GS_BCKND_SX = 2, &
       GS_BCKND_DEV = 3
  
  !> Gather-scatter backend
  type, public, abstract :: gs_bcknd_t
     type(c_ptr) :: gather_event = C_NULL_PTR
   contains
     procedure(gs_backend_init), pass(this), deferred :: init
     procedure(gs_backend_free), pass(this), deferred :: free
     procedure(gs_gather), pass(this), deferred :: gather
     procedure(gs_scatter), pass(this), deferred :: scatter
  end type gs_bcknd_t

  !> Abstract interface for initialising a Gather-Scatter backend
  abstract interface
     subroutine gs_backend_init(this, nlocal, nshared, nlcl_blks, nshrd_blks)
       import gs_bcknd_t
       class(gs_bcknd_t), intent(inout) :: this
       integer, intent(in) :: nlocal
       integer, intent(in) :: nshared
       integer, intent(in) :: nlcl_blks
       integer, intent(in) :: nshrd_blks
     end subroutine gs_backend_init
  end interface

  !> Abstract interface for deallocating a Gather-Scatter backend
  abstract interface
     subroutine gs_backend_free(this)
       import gs_bcknd_t
       class(gs_bcknd_t), intent(inout) :: this
     end subroutine gs_backend_free
  end interface

  !> Abstract interface for the Gather kernel
  !! \f$ v(dg(i)) = op(v(dg(i)), u(gd(i)) \f$
  abstract interface
     subroutine gs_gather(this, v, m, o, dg, u, n, gd, nb, b, op, shrd)
       import gs_bcknd_t       
       import rp
       integer, intent(in) :: m
       integer, intent(in) :: n
       integer, intent(in) :: nb
       class(gs_bcknd_t), intent(inout) :: this
       real(kind=rp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=rp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       integer, intent(in) :: o
       integer, intent(in) :: op
       logical, intent(in) :: shrd    
     end subroutine gs_gather
  end interface
  
  !> Abstract interface for the Scatter kernel
  !! \f$ u(gd(i) = v(dg(i)) \f$
  abstract interface
     subroutine gs_scatter(this, v, m, dg, u, n, gd, nb, b, shrd)
       import gs_bcknd_t       
       import rp
       integer, intent(in) :: m
       integer, intent(in) :: n
       integer, intent(in) :: nb
       class(gs_bcknd_t), intent(inout) :: this              
       real(kind=rp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=rp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       logical, intent(in) :: shrd
     end subroutine gs_scatter
  end interface

  
end module gs_bcknd
