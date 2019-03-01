!> Defines a hexahedron element
!
module hex
  use num_types
  use element
  use point
  implicit none
  private

  integer, public, parameter :: NEKO_HEX_NPTS = 8 !< Number of points
  integer, public, parameter :: NEKO_HEX_GDIM = 3 !< Geometric dimension

  type, public, extends(element_t) :: hex_t
   contains
     procedure, pass(this) :: init => hex_init
     procedure, pass(this) :: diameter => hex_diameter
     procedure, pass(this) :: equal => hex_equal
     generic :: operator(.eq.) => equal
  end type hex_t

contains
  
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
       d1 = d1 + (p7%x(i) - p1%x(i))**2
       d2 = d2 + (p8%x(i) - p2%x(i))**2
       d3 = d3 + (p5%x(i) - p3%x(i))**2
       d4 = d4 + (p6%x(i) - p4%x(i))**2
    end do

    res = sqrt(max(max(d1, d2), max(d3, d4)))

  end function hex_diameter

  pure function hex_equal(this, other) result(res)
    class(hex_t), intent(in) :: this
    class(element_t), intent(in) :: other
    integer :: i
    logical :: res

    res = .false.
    select type(other)
    class is (hex_t)
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
