!> Defines a quadrilateral element
!
module quad
  use num_types
  use element
  use point
  implicit none
  private

  integer, public, parameter :: NEKO_QUAD_NPTS = 4 !< Number of points
  integer, public, parameter :: NEKO_QUAD_GDIM = 2 !< Geometric dimension

  type, public, extends(element_t) :: quad_t
   contains
     procedure, pass(this) :: quad_init
     procedure, pass(this) :: diameter => quad_diameter
     procedure, pass(this) :: equal => quad_equal
     generic, public :: init => quad_init
     generic :: operator(.eq.) => equal
  end type quad_t

contains
  
  subroutine quad_init(this, id, p1, p2, p3, p4)
    class(quad_t), intent(inout) :: this
    integer, intent(inout) :: id
    type(point_t), target, intent(in) :: p1, p2, p3, p4

    call this%init(id, NEKO_QUAD_GDIM, NEKO_QUAD_NPTS)

    this%pts(1)%p => p1
    this%pts(2)%p => p2
    this%pts(3)%p => p3
    this%pts(4)%p => p4

  end subroutine quad_init

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
