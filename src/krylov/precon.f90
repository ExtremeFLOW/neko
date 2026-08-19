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
!> Krylov preconditioner
module precon
  use num_types, only : rp
  implicit none
  private

  !> Defines a canonical Krylov preconditioner
  type, public, abstract :: pc_t
   contains
     procedure(pc_solve), pass(this), deferred :: solve
     procedure(pc_update), pass(this), deferred :: update
  end type pc_t

  !> Abstract interface for solving \f$ M z = r \f$
  !!
  !! @param z vector of length @a n
  !! @param r vector of length @a n
  abstract interface
     subroutine pc_solve(this, z, r, n)
       import rp
       import :: pc_t
       implicit none
       integer, intent(in) :: n
       class(pc_t), intent(inout) :: this
       real(kind=rp), dimension(n), intent(inout) :: z
       real(kind=rp), dimension(n), intent(inout) :: r
     end subroutine pc_solve
     subroutine pc_update(this)
       import rp
       import :: pc_t
       implicit none
       class(pc_t), intent(inout) :: this
     end subroutine pc_update
  end interface

  interface
     !> Allocate a preconditioner
     module subroutine precon_allocator(pc, type_name)
       class(pc_t), allocatable, intent(inout) :: pc
       character(len=*), intent(in) :: type_name
     end subroutine precon_allocator

     !> Destroy a preconditioner
     module subroutine precon_destroy(pc)
       class(pc_t), allocatable, intent(inout) :: pc
     end subroutine precon_destroy
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for a preconditioner allocator.
  !! Implemented in user modules, should allocate `obj` to the custom user
  !! type.
  abstract interface
     subroutine precon_allocate(obj)
       import pc_t
       class(pc_t), allocatable, intent(inout) :: obj
     end subroutine precon_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_precon(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(precon_allocate), pointer, intent(in) :: allocator
     end subroutine register_precon
  end interface

  !> A name-allocator pair for user-defined preconditioner types.
  type precon_allocator_entry
     character(len=20) :: type_name
     procedure(precon_allocate), pointer, nopass :: allocator
  end type precon_allocator_entry

  !> Registry of allocators for user-defined preconditioner types.
  type(precon_allocator_entry), allocatable :: precon_registry(:)

  !> The size of the `precon_registry`.
  integer :: precon_registry_size = 0

  public :: precon_allocator, precon_destroy, register_precon, precon_allocate

end module precon
