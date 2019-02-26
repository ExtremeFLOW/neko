!> Defines a hexahedron element
!
module hex
  use num_types
  use element
  use point
  implicit none
  private

  integer, parameter :: NEKO_HEX_NPTS = 8 !< Number of points
  integer, parameter :: NEKO_HEX_GDIM = 3 !< Geometric dimension

  type, public, extends(element_t) :: hex_t
   contains
     procedure, pass(this) :: hex_init
     procedure, pass(this) :: diameter => hex_diameter
     procedure, pass(this) :: equal => hex_equal
     generic, public :: init => hex_init
     generic :: operator(.eq.) => equal
  end type hex_t

contains
  
  subroutine hex_init(this, id, p1, p2, p3, p4, p5, p6, p7, p8)
    class(hex_t), intent(inout) :: this
    integer, intent(inout) :: id
    type(point_t), target, intent(in) :: p1, p2, p3, p4, p5, p6, p7, p8
    integer :: i

    call this%init(id, NEKO_HEX_GDIM, NEKO_HEX_NPTS)
    
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
    real(kind=dp) :: x(3, NEKO_HEX_NPTS)
    integer :: i

    d1 = 0d0
    d2 = 0d0
    d3 = 0d0
    d4 = 0d0

    do i = 1, NEKO_HEX_GDIM
       x(:, i) = this%pts(i)%p%x
    end do

    do i = 1, NEKO_HEX_GDIM
       d1 = d1 + (x(i, 7) - x(i, 1))**2
       d2 = d2 + (x(i, 8) - x(i, 2))**2
       d3 = d3 + (x(i, 5) - x(i, 3))**2
       d4 = d4 + (x(i, 6) - x(i, 4))**2
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
