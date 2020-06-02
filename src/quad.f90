!> Defines a quadrilateral element
module quad
  use num_types
  use element
  use tuple
  use point
  implicit none
  private

  integer, public, parameter :: NEKO_QUAD_NPTS = 4 !< Number of points
  integer, public, parameter :: NEKO_QUAD_GDIM = 2 !< Geometric dimension

  !> Quadrilateral element
  !! @details
  !! 2D element composed of 4 points
  !! @verbatim
  !! Node numbering (NEKTON preprocessor notation)
  !!
  !!      4+-----+3    ^ s                 
  !!       |     |     |                   
  !!       |     |     |                   
  !!      1+-----+2    +----> r            
  !!
  !! @endverbatim
  type, public, extends(element_t) :: quad_t
   contains
     procedure, pass(this) :: init => quad_init
     procedure, pass(this) :: facet_id => quad_facet_id
     procedure, pass(this) :: diameter => quad_diameter
     procedure, pass(this) :: centroid => quad_centroid
     procedure, pass(this) :: equal => quad_equal
     generic :: operator(.eq.) => equal
  end type quad_t

  !> Edge node ids
  !! @details
  !! @verbatim
  !! Edge numbering (similar to NEKTON symmetric notation)
  !!          4
  !!       +------+      ^ s                 
  !!       |      |      |
  !!     1 |      | 2    |
  !!       |      |      |                   
  !!       +------+      +-----> r            
  !!          3
  !! @endverbatim
  integer, parameter, dimension(2, 4) :: edge_nodes = reshape((/1,4,&
                                                                2,3,&
                                                                1,2,&
                                                                4,3 /),&
                                                                (/2,4/))
  
contains

  !> Create a quadrilateral element based upon four points
  subroutine quad_init(this, id, p1, p2, p3, p4)
    class(quad_t), intent(inout) :: this
    integer, intent(inout) :: id
    type(point_t), target, intent(in) :: p1, p2, p3, p4

    call this%element(id, NEKO_QUAD_GDIM, NEKO_QUAD_NPTS)

    this%pts(1)%p => p1
    this%pts(2)%p => p2
    this%pts(3)%p => p3
    this%pts(4)%p => p4

  end subroutine quad_init

  !> Return the edge id for face @a i as a 2-tuple @a t
  subroutine quad_facet_id(this, t, side)
    class(quad_t), intent(in) :: this
    class(tuple_t), intent(inout) :: t
    integer, intent(in) :: side
    type(point_t), pointer :: p1, p2

    p1 => this%p(edge_nodes(1, side))
    p2 => this%p(edge_nodes(2, side))

    select type(t)
    type is(tuple_i4_t)
       t = (/ p1%id(), p2%id() /)
    end select
    
  end subroutine quad_facet_id

  !> Compute the diameter of a quadrilateral element
  function quad_diameter(this) result(res)
    class(quad_t), intent(in) :: this
    real(kind=dp) :: d1, d2, res
    type(point_t), pointer :: p1, p2, p3, p4
    integer :: i

    d1 = 0d0
    d2 = 0d0

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)
    p4 => this%p(4)

    do i = 1, NEKO_QUAD_GDIM
       d1 = d1 + (p4%x(i) - p1%x(i))**2
       d2 = d2 + (p3%x(i) - p2%x(i))**2
    end do

    res = sqrt(max(d1, d2))

  end function quad_diameter

  !> Compute the centroid of a quadrilateral element
  function quad_centroid(this) result(res)
    class(quad_t), intent(in) :: this
    type(point_t) :: res
    type(point_t), pointer :: p1, p2, p3, p4
    integer :: i

    p1 => this%p(1)
    p2 => this%p(2)
    p3 => this%p(3)
    p4 => this%p(4)
    res%x = 0d0

    do i = 1, this%gdim()
       res%x(i) = 0.25d0 * p1%x(i) + p2%x(i) + p3%x(i) + p4%x(i)
    end do
  end function quad_centroid

  !> Check if two quad elements are equal
  !! @note Based on coordinates not global ids
  pure function quad_equal(this, other) result(res)
    class(quad_t), intent(in) :: this
    class(element_t), intent(in) :: other
    integer :: i    
    logical :: res

    res = .false.
    select type(other)
    class is (quad_t)
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
    
  end function quad_equal
  
end module quad
