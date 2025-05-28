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

module el_finder
  use num_types, only : rp
  use stack, only: stack_i4_t
  use point, only: point_t
  implicit none
  private

  !> Base type for element finder providing element candidates
  !> for a given point in the domain.
  type, public, abstract :: el_finder_t
   contains
     procedure(el_finder_free), pass(this), deferred :: free
     procedure(el_finder_find), pass(this), deferred :: find
     procedure(el_finder_find_batch), pass(this), deferred :: find_batch
  end type el_finder_t

  abstract interface
     subroutine el_finder_free(this)
       import el_finder_t
       class(el_finder_t), intent(inout) :: this
     end subroutine el_finder_free
  end interface

  abstract interface
     subroutine el_finder_find(this, my_point, el_candidates)
       import rp
       import stack_i4_t
       import point_t
       import el_finder_t
       implicit none
       class(el_finder_t), intent(inout) :: this
       type(point_t), intent(in) :: my_point
       type(stack_i4_t), intent(inout) :: el_candidates
     end subroutine el_finder_find
  end interface

  abstract interface
     subroutine el_finder_find_batch(this, points, n_points, all_el_candidates, n_el_cands)
       import rp
       import stack_i4_t
       import point_t
       import el_finder_t
       implicit none
       class(el_finder_t), intent(inout) :: this
       integer, intent(in) :: n_points
       real(kind=rp), intent(in) :: points(3,n_points)
       type(stack_i4_t), intent(inout) :: all_el_candidates
       integer, intent(inout) :: n_el_cands(n_points)
     end subroutine el_finder_find_batch
  end interface
end module el_finder
