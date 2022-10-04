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
!> Defines a hexahedron element
module hex
  use num_types
  use element
  use tuple
  use point
  implicit none
  private

  integer, public, parameter :: NEKO_HEX_NPTS = 8  !< Number of points
  integer, public, parameter :: NEKO_HEX_NFCS = 6  !< Number of faces
  integer, public, parameter :: NEKO_HEX_NEDS = 12 !< Number of edges
  integer, public, parameter :: NEKO_HEX_GDIM = 3  !< Geometric dimension


  !> Hexahedron element
  !! @details
  !! 3D element composed of 8 points
  !! @verbatim
  !! Node numbering (NEKTON preprocessor notation)
  !!
  !!          3+-----+4    ^ s                 
  !!          /     /|     |                   
  !!         /     / |     |                   
  !!       7+-----+8 +2    +----> r            
  !!        |     | /     /                    
  !!        |     |/     /                     
  !!       5+-----+6    t
  !!
  !! @endverbatim
  type, public, extends(element_t) :: hex_t
   contains
     procedure, pass(this) :: init => hex_init
     procedure, pass(this) :: facet_id => hex_facet_id
     procedure, pass(this) :: facet_order => hex_facet_order
     procedure, pass(this) :: diameter => hex_diameter
     procedure, pass(this) :: centroid => hex_centroid
     procedure, pass(this) :: edge_id => hex_edge_id          
     procedure, pass(this) :: equal => hex_equal
     generic :: operator(.eq.) => equal
  end type hex_t

  !> Face node ids
  !! @details
  !! @verbatim
  !! Face numbering (NEKTON symmetric notation)
  !!                     
  !!          +--------+     ^ S
  !!         /        /|     |
  !!        /    4   / |     |
  !!  1--> /        /  |     |
  !!      +--------+ 2 +     +----> R
  !!      |        |  /     /
  !!      |    6   | /     /
  !!      |        |/     /
  !!      +--------+     T
  !!          3
  !!
  !! @endverbatim
  !! @note Local node numbering (points)
  integer, parameter, dimension(4, 6) :: face_nodes = reshape((/1,5,7,3,&
                                                                2,6,8,4,&
                                                                1,2,6,5,&
                                                                3,4,8,7,&
                                                                1,2,4,3,&
                                                                5,6,8,7/),&
                                                                (/4,6/))
  
  !> Edge node ids
  !! @details
  !! @verbatim
  !! Edge numbering (similar to NEKTON symmetric notation)
  !!
  !!              2      
  !!          +--------+        ^ S
  !!         /        /|        |
  !!  11--> /   12-->/ | <--6   |
  !!       /   4    /  |        |
  !!      +--------+   +        +----> R
  !!      |        |  /        /
  !!  7-->|    8-->| /<--10   /
  !!      |        |/        /
  !!      +--------+        T
  !!           3
  !!
  !! @endverbatim
  integer, parameter, dimension(2, 12) :: edge_nodes = reshape((/1,2,&
                                                                3,4,&
                                                                5,6,&
                                                                7,8,&
                                                                1,3,&
                                                                2,4,&
                                                                5,7,&
                                                                6,8,&
                                                                1,5,&
                                                                2,6,&
                                                                3,7,&
                                                                4,8/),&
                                                                (/2,12/))
  
contains
  
  !> Create a hexahedron element based upon eight points
  subroutine hex_init(this, id, p1, p2, p3, p4, p5, p6, p7, p8)
    class(hex_t), intent(inout) :: this
    integer, intent(inout) :: id
    type(point_t), target, intent(in) :: p1, p2, p3, p4, p5, p6, p7, p8

    call this%element(id, NEKO_HEX_GDIM, NEKO_HEX_NPTS)
    
    this%pts(1)%p => p1
    this%pts(2)%p => p2
    this%pts(3)%p => p3
    this%pts(4)%p => p4
    this%pts(5)%p => p5
    this%pts(6)%p => p6
    this%pts(7)%p => p7
    this%pts(8)%p => p8

  end subroutine hex_init

  !> Return the facet id for face @a i as a 4-tuple @a t
  subroutine hex_facet_id(this, t, side) 
    class(hex_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    integer :: i, j, temp
    type(point_t), pointer :: p1,p2,p3,p4

    p1 => this%p(face_nodes(1, side))
    p2 => this%p(face_nodes(2, side))
    p3 => this%p(face_nodes(3, side))
    p4 => this%p(face_nodes(4, side))

    select type(t)
    type is(tuple4_i4_t)
       t%x = (/ p1%id(), p2%id(), p3%id(), p4%id() /)
       do i = 1, 3 
          do j = i+1,4
             if(t%x(j) .lt. t%x(i)) then
                temp = t%x(i)
                t%x(i) = t%x(j)
                t%x(j) = temp
             endif
          enddo
       enddo
    end select

  end subroutine hex_facet_id

  !> Return the ordered points for face @a i as a 4-tuple @a t
  subroutine hex_facet_order(this, t, side) 
    class(hex_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    type(point_t), pointer :: p1,p2,p3,p4

    p1 => this%p(face_nodes(1, side))
    p2 => this%p(face_nodes(2, side))
    p3 => this%p(face_nodes(3, side))
    p4 => this%p(face_nodes(4, side))

    select type(t)
    type is(tuple4_i4_t)
       t%x = (/ p1%id(), p2%id(), p3%id(), p4%id() /)
    end select

  end subroutine hex_facet_order


  !> Return the edge id for an edge @a i as a 2-tuple @a t
  subroutine hex_edge_id(this, t, side) 
    class(hex_t), intent(in) :: this
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

  end subroutine hex_edge_id
  
  !> Compute the diameter of a hexahedron element
  function hex_diameter(this) result(res)
    class(hex_t), intent(in) :: this
    real(kind=dp) :: d1, d2, d3, d4, res
    type(point_t), pointer :: p1, p2, p3, p4, p5, p6, p7, p8
    integer :: i

    d1 = 0d0
    d2 = 0d0
    d3 = 0d0
    d4 = 0d0

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)
    p4 => this%p(4)
    p5 => this%p(5)
    p6 => this%p(6)
    p7 => this%p(7)
    p8 => this%p(8)

    do i = 1, NEKO_HEX_GDIM
       d1 = d1 + (p8%x(i) - p1%x(i))**2
       d2 = d2 + (p7%x(i) - p2%x(i))**2
       d3 = d3 + (p5%x(i) - p4%x(i))**2
       d4 = d4 + (p6%x(i) - p3%x(i))**2
    end do

    res = sqrt(max(max(d1, d2), max(d3, d4)))

  end function hex_diameter

  !> Compute the centroid of a hexahedron element
  function hex_centroid(this) result(res)
    class(hex_t), intent(in) :: this
    type(point_t) :: res
    type(point_t), pointer :: p1, p2, p3, p4, p5, p6, p7, p8
    integer :: i

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)
    p4 => this%p(4)
    p5 => this%p(5)
    p6 => this%p(6)
    p7 => this%p(7)
    p8 => this%p(8)
    res%x = 0d0

    do i = 1, this%gdim()
       res%x(i) = 0.125 * (p1%x(i) + p2%x(i) + p3%x(i) + p4%x(i) + &
            p5%x(i) + p6%x(i) + p7%x(i) + p8%x(i))
    end do
    
  end function hex_centroid

  !> Check if two hex elements are equal
  !! @note Based on coordinates not global ids
  pure function hex_equal(this, other) result(res)
    class(hex_t), intent(in) :: this
    class(element_t), intent(in) :: other
    integer :: i
    logical :: res

    res = .false.
    select type(other)
    type is (hex_t)
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

  end function hex_equal
  
end module hex
