! Copyright (c) 2021-2023, The Neko Authors
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
!> Defines a triangular element
module tri
  use num_types, only : dp
  use element, only : element_t
  use tuple, only : tuple_t, tuple_i4_t
  use point, only : point_t
  implicit none
  private

  integer, public, parameter :: NEKO_TRI_NPTS = 3 !< Number of points
  integer, public, parameter :: NEKO_TRI_NEDS = 3 !< Number of edges
  integer, public, parameter :: NEKO_TRI_GDIM = 2 !< Geometric dimension

  !> Triangular element
  !! @details
  !! 2D element composed of 3 points
  !! @verbatim
  !! Node numbering
  !!
  !!    3+
  !!     |\
  !!     | \        ^ s
  !!     |  \       |
  !!     |   \      |
  !!    1+----+2    +---> r
  !!
  !! @endverbatim
  type, public, extends(element_t) :: tri_t
   contains
     procedure, pass(this) :: init => tri_init
     procedure, pass(this) :: facet_id => tri_facet_id
     procedure, pass(this) :: facet_order => tri_facet_order
     procedure, pass(this) :: diameter => tri_diameter
     procedure, pass(this) :: centroid => tri_centroid
     procedure, pass(this) :: equal => tri_equal
     generic :: operator(.eq.) => equal
  end type tri_t

  !> Edge node ids
  !! @details
  !! @verbatim
  !! Edge numbering
  !!
  !!     +
  !!     |\
  !!     | \       ^ s
  !!   1 |  \ 2    |
  !!     |   \     |
  !!     +----+    +---> r
  !!       3
  !!
  !! @endverbatim
  integer, parameter, dimension(2, 3) :: edge_nodes = reshape((/1,3,&
                                                                2,3,&
                                                                1,2 /),&
                                                                (/2,3/))

contains

  !> Create a trinagular element based upon three points
  subroutine tri_init(this, id, p1, p2, p3)
    class(tri_t), intent(inout) :: this
    integer, intent(inout) :: id
    type(point_t), target, intent(in) :: p1, p2, p3

    call this%element(id, NEKO_TRI_GDIM, NEKO_TRI_NPTS)

    this%pts(1)%p => p1
    this%pts(2)%p => p2
    this%pts(3)%p => p3

  end subroutine tri_init

  !> Return the edge id for face @a i as a 2-tuple @a t
  !! @todo sort this
  subroutine tri_facet_id(this, t, side)
    class(tri_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    type(point_t), pointer :: p1, p2

    p1 => this%p(edge_nodes(1, side))
    p2 => this%p(edge_nodes(2, side))

    select type(t)
    type is(tuple_i4_t)
       if (p1 .lt. p2) then
          t%x = (/ p1%id(), p2%id() /)
       else
          t%x = (/ p2%id(), p1%id() /)
       end if
    end select

  end subroutine tri_facet_id

  !> Return the ordered edge for face @a i as a 2-tuple @a t
  subroutine tri_facet_order(this, t, side)
    class(tri_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    type(point_t), pointer :: p1, p2

    p1 => this%p(edge_nodes(1, side))
    p2 => this%p(edge_nodes(2, side))

    select type(t)
    type is(tuple_i4_t)
       t%x = (/ p1%id(), p2%id() /)
    end select

  end subroutine tri_facet_order

  !> Compute the diameter of a triangular element
  function tri_diameter(this) result(res)
    class(tri_t), intent(in) :: this
    real(kind=dp) :: d1, d2, d3, res
    type(point_t), pointer :: p1, p2, p3
    integer :: i

    d1 = 0d0
    d2 = 0d0
    d3 = 0d0

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)

    do i = 1, NEKO_TRI_GDIM
       d1 = d1 + (p2%x(i) - p1%x(i))**2
       d2 = d2 + (p3%x(i) - p2%x(i))**2
       d3 = d3 + (p1%x(i) - p3%x(i))**2
    end do

    res = sqrt(max(max(d1, d2), d3))

  end function tri_diameter

  !> Compute the centroid of a triangular element
  function tri_centroid(this) result(res)
    class(tri_t), intent(in) :: this
    type(point_t) :: res
    type(point_t), pointer :: p1, p2, p3
    integer :: i

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)
    res%x = 0d0

    do i = 1, this%gdim()
       res%x(i) = 1d0/3d0 * (p1%x(i) + p2%x(i) + p3%x(i))
    end do
  end function tri_centroid

  !> Check if two triangle elements are equal
  !! @note Based on coordinates not global ids
  pure function tri_equal(this, other) result(res)
    class(tri_t), intent(in) :: this
    class(element_t), intent(in) :: other
    integer :: i
    logical :: res

    res = .false.
    select type(other)
    class is (tri_t)
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

  end function tri_equal

end module tri
