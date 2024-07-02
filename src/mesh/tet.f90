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
!> Defines a tetrahedral element
module tet
  use num_types, only : dp
  use element, only : element_t
  use tuple
  use point, only : point_t
  implicit none
  private

  integer, public, parameter :: NEKO_TET_NPTS = 4  !< Number of points
  integer, public, parameter :: NEKO_TET_NFCS = 4  !< Number of faces
  integer, public, parameter :: NEKO_TET_NEDS = 6 !< Number of edges
  integer, public, parameter :: NEKO_TET_GDIM = 3  !< Geometric dimension


  !> Tetrahedral element
  !! @details
  !! 3D element composed of 4 points
  !! @verbatim
  !! Node numbering
  !!
  !!        3 +           ^ s
  !!         /|\          |
  !!        / | \         |
  !!     1 +..|..+ 2      +----> r
  !!        \ | /        /
  !!         \|/        /
  !!        4 +        t
  !!
  !! @endverbatim
  type, public, extends(element_t) :: tet_t
   contains
     procedure, pass(this) :: init => tet_init
     procedure, pass(this) :: facet_id => tet_facet_id
     procedure, pass(this) :: facet_order => tet_facet_order
     procedure, pass(this) :: diameter => tet_diameter
     procedure, pass(this) :: centroid => tet_centroid
     procedure, pass(this) :: edge_id => tet_edge_id
     procedure, pass(this) :: equal => tet_equal
     generic :: operator(.eq.) => equal
  end type tet_t

  !> Face node ids
  !! @details
  !! @verbatim
  !! Face numbering
  !!
  !!          +   4       ^ s
  !!         /|\ /        |
  !!        / | \         |
  !!       / 1|2 \        |
  !!      +...|...+       +----> r
  !!       \  |3 /       /
  !!        \ | /       /
  !!         \|/       /
  !!          +       t
  !!
  !! @endverbatim
  !! @note Local node numbering (points)
  integer, parameter, dimension(3, 4) :: face_nodes = reshape((/1,3,4,&
                                                                2,3,4,&
                                                                1,2,4,&
                                                                1,2,3/),&
                                                                (/3,4/))

  !> Edge node ids
  !! @details
  !! @verbatim
  !! Edge numbering
  !!
  !!      2   +   3       ^ s
  !!       \ /|\ /        |
  !!        / | \         |
  !!       /  |  \        |
  !!      +.1.|...+       +----> r
  !!       \  4  /       /
  !!   5--> \ | / <--6  /
  !!         \|/       /
  !!          +       t
  !!
  !! @endverbatim
  integer, parameter, dimension(2, 6) :: edge_nodes = reshape((/1,2,&
                                                                1,3,&
                                                                2,3,&
                                                                3,4,&
                                                                1,4,&
                                                                2,4/),&
                                                                (/2,6/))

contains

  !> Create a tetrahedral element based upon four points
  subroutine tet_init(this, id, p1, p2, p3, p4)
    class(tet_t), intent(inout) :: this
    integer, intent(inout) :: id
    type(point_t), target, intent(in) :: p1, p2, p3, p4

    call this%element(id, NEKO_TET_GDIM, NEKO_TET_NPTS)

    this%pts(1)%p => p1
    this%pts(2)%p => p2
    this%pts(3)%p => p3
    this%pts(4)%p => p4

  end subroutine tet_init

  !> Return the facet id for face @a i as a 3-tuple @a t
  subroutine tet_facet_id(this, t, side)
    class(tet_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    integer :: i, j, temp
    type(point_t), pointer :: p1,p2,p3

    p1 => this%p(face_nodes(1, side))
    p2 => this%p(face_nodes(2, side))
    p3 => this%p(face_nodes(3, side))

    select type(t)
    type is(tuple3_i4_t)
       t%x = (/ p1%id(), p2%id(), p3%id() /)
       do i = 1, 2
          do j = i+1,3
             if(t%x(j) .lt. t%x(i)) then
                temp = t%x(i)
                t%x(i) = t%x(j)
                t%x(j) = temp
             endif
          enddo
       enddo
    end select

  end subroutine tet_facet_id

  !> Return the ordered points for face @a i as a 3-tuple @a t
  subroutine tet_facet_order(this, t, side)
    class(tet_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    type(point_t), pointer :: p1,p2,p3

    p1 => this%p(face_nodes(1, side))
    p2 => this%p(face_nodes(2, side))
    p3 => this%p(face_nodes(3, side))

    select type(t)
    type is(tuple3_i4_t)
       t%x = (/ p1%id(), p2%id(), p3%id() /)
    end select

  end subroutine tet_facet_order


  !> Return the edge id for an edge @a i as a 2-tuple @a t
  subroutine tet_edge_id(this, t, side)
    class(tet_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    type(point_t), pointer :: p1,p2

    p1 => this%p(edge_nodes(1, side))
    p2 => this%p(edge_nodes(2, side))

    select type(t)
    type is(tuple_i4_t)
       if (p1%id() .lt. p2%id()) then
          t%x = (/ p1%id(), p2%id() /)
       else
          t%x = (/ p2%id(), p1%id() /)
       endif

    end select

  end subroutine tet_edge_id

  !> Compute the diameter of a tetrahedral element
  function tet_diameter(this) result(res)
    class(tet_t), intent(in) :: this
    real(kind=dp) :: d1, d2, d3, d4, d5, d6, res
    type(point_t), pointer :: p1, p2, p3, p4
    integer :: i

    d1 = 0d0
    d2 = 0d0
    d3 = 0d0
    d4 = 0d0
    d5 = 0d0
    d6 = 0d0

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)
    p4 => this%p(4)


    do i = 1, NEKO_TET_GDIM
       d1 = d1 + (p2%x(i) - p3%x(i))**2
       d2 = d2 + (p1%x(i) - p3%x(i))**2
       d3 = d3 + (p1%x(i) - p2%x(i))**2
       d4 = d4 + (p1%x(i) - p4%x(i))**2
       d5 = d5 + (p2%x(i) - p4%x(i))**2
       d6 = d6 + (p3%x(i) - p4%x(i))**2
    end do

    res = d1
    res = max(res, d2)
    res = max(res, d3)
    res = max(res, d4)
    res = max(res, d5)
    res = max(res, d6)
    res = sqrt(res)

  end function tet_diameter

  !> Compute the centroid of a tetrahedral element
  function tet_centroid(this) result(res)
    class(tet_t), intent(in) :: this
    type(point_t) :: res
    type(point_t), pointer :: p1, p2, p3, p4
    integer :: i

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)
    p4 => this%p(4)
    res%x = 0d0

    do i = 1, this%gdim()
       res%x(i) = 0.25 * (p1%x(i) + p2%x(i) + p3%x(i) + p4%x(i))
    end do

  end function tet_centroid

  !> Check if two tet elements are equal
  !! @note Based on coordinates not global ids
  pure function tet_equal(this, other) result(res)
    class(tet_t), intent(in) :: this
    class(element_t), intent(in) :: other
    integer :: i
    logical :: res

    res = .false.
    select type(other)
    type is (tet_t)
       if ((this%gdim() .eq. other%gdim()) .and. &
            (this%npts() .eq. other%npts())) then
          do i = 1, this%npts()
             if (this%pts(i)%p .ne. other%pts(i)%p) then
                return
             end if
          end do
          res = .true.
       end if
    end select

  end function tet_equal

end module tet
