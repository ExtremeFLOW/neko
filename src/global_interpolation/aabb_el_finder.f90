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

module aabb_el_finder
  use num_types, only : rp, dp, xp
  use neko_config, only : NEKO_BCKND_DEVICE
  use el_finder, only : el_finder_t
  use space, only : space_t
  use stack, only : stack_i4_t
  use tuple, only : tuple_i4_t
  use point, only : point_t
  use aabb, only : aabb_t
  use aabb_tree, only : aabb_tree_t, aabb_node_t, AABB_NULL_NODE
  implicit none
  private

  !> Implements global interpolation for arbitrary points in the domain.
  type, public, extends(el_finder_t) :: aabb_el_finder_t
     real(kind=dp) :: padding
     !> Structure to find rank candidates
     type(aabb_t), allocatable :: local_aabb(:)
     type(aabb_tree_t) :: local_aabb_tree

   contains
     procedure, pass(this) :: init => aabb_el_finder_init
     procedure, pass(this) :: free => aabb_el_finder_free
     procedure, pass(this) :: find => aabb_el_finder_find_candidates
     procedure, pass(this) :: find_batch => aabb_el_finder_find_candidates_batch

  end type aabb_el_finder_t

contains

  !> Initialize the AABB element finder.
  subroutine aabb_el_finder_init(this, x, y, z, nel, Xh, padding)
    class(aabb_el_finder_t), intent(inout) :: this
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(in), target :: x(:), y(:), z(:)
    integer, intent(in) :: nel
    real(kind=dp), intent(in) :: padding
    integer :: id1, id2, i, lx, ly, lz
    this%padding = padding
    lx = Xh%lx
    ly = Xh%ly
    lz = Xh%lz
    if (allocated(this%local_aabb)) deallocate(this%local_aabb)
    allocate(this%local_aabb(nel))
    !> Create a local tree for each element at this rank
    call this%local_aabb_tree%init(nel)
    do i = 1, nel
       id1 = lx*ly*lz*(i-1)+1
       id2 = lx*ly*lz*(i)
       call this%local_aabb(i)%init( real((/minval(x(id1:id2)), &
            minval(y(id1:id2)), &
            minval(z(id1:id2))/), dp), &
            real((/maxval(x(id1:id2)), &
            maxval(y(id1:id2)), &
            maxval(z(id1:id2))/), dp))
    end do
    call this%local_aabb_tree%build_from_aabb(this%local_aabb, padding)
  end subroutine aabb_el_finder_init

  !> Free the AABB element finder.
  subroutine aabb_el_finder_free(this)
    class(aabb_el_finder_t), intent(inout) :: this

    ! Free the AABB element finder
    if (allocated(this%local_aabb)) deallocate(this%local_aabb)


  end subroutine aabb_el_finder_free


  !> It uses the AABB tree to find the elements that overlap with the point.
  subroutine aabb_el_finder_find_candidates(this, my_point, el_candidates)
    class(aabb_el_finder_t), intent(inout) :: this
    type(point_t), intent(in) :: my_point
    type(stack_i4_t), intent(inout) :: el_candidates

    ! Find the element candidates for a given point
    call this%local_aabb_tree%query_overlaps(my_point, -1, el_candidates)

  end subroutine aabb_el_finder_find_candidates

  ! Might be better to organize this slightly differently
  ! In order to get more cache hits
  subroutine aabb_el_finder_find_candidates_batch(this, points, n_points, &
       all_el_candidates, n_el_cands)
    class(aabb_el_finder_t), intent(inout) :: this
    integer, intent(in) :: n_points
    real(kind=rp), intent(in) :: points(3, n_points)
    type(stack_i4_t), intent(inout) :: all_el_candidates
    integer, intent(inout) :: n_el_cands(n_points)
    type(stack_i4_t) :: el_candidates
    type(point_t) :: my_point
    integer :: i, j
    integer, pointer :: el_cands(:)
    integer :: stupid_intent
    real(kind=dp) :: pt_xyz(3)

    call all_el_candidates%clear()
    call el_candidates%init()
    n_el_cands = 0

    do i = 1, n_points
       pt_xyz = (/ points(1,i), points(2,i), points(3,i) /)
       call my_point%init(pt_xyz)
       call el_candidates%clear()
       call this%find(my_point, el_candidates)
       el_cands => el_candidates%array()
       do j = 1, el_candidates%size()
          stupid_intent = el_cands(j) - 1
          call all_el_candidates%push(stupid_intent) !< OBS c indexing
       end do
       n_el_cands(i) = el_candidates%size()
    end do

  end subroutine aabb_el_finder_find_candidates_batch


end module aabb_el_finder
