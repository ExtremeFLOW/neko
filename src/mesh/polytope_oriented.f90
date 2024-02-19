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
 !> Abstract type for abstract oriented polytope class for topology
module polytope_oriented
  use num_types, only : i4
  use polytope, only : polytope_t
  use polytope_aligned, only : polytope_aligned_t
  implicit none
  private

  public :: polytope_oriented_t

  !> Base type for an abstract oriented polytope
  !! @details This is an abstract type basically equal to
  !! @ref polytope_aligned_t This type corresponds to building blocs of
  !! the higher dimension abstract objects the mesh topology consists of.
  !! @note The reason we need this type is to define a deferred
  !! constructor. We cannot define it in the base type because the type
  !! @ref polytope_actualisation_t also descends from it  but requires a
  !! constructor with a different signature.
  type, extends(polytope_aligned_t), abstract :: polytope_oriented_t
   contains
     !> Free aligned polytope
     procedure, pass(this) :: free => polytope_free
     !> Initialise an aligned polytope
     procedure(polytope_oriented_init), pass(this), deferred :: init
  end type polytope_oriented_t

  !> Initialise a polytope with alignment information for topology
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  abstract interface
     subroutine polytope_oriented_init(this, pltp, algn)
       import i4
       import polytope_t
       import polytope_oriented_t
       class(polytope_oriented_t), intent(inout) :: this
       class(polytope_t), target, intent(in) :: pltp
       integer(i4), intent(in) :: algn
     end subroutine polytope_oriented_init
  end interface

contains

  !> Free oriented polytope
  subroutine polytope_free(this)
    class(polytope_oriented_t), intent(inout) :: this

    call this%free_base()
  end subroutine polytope_free

end module polytope_oriented
