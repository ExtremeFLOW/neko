!> Defines a hexahedron element
!
module hex
  use num_types
  use element
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
  
  subroutine hex_init(this, id, pts)
    class(hex_t), intent(inout) :: this
    integer, intent(inout) :: id
    integer, dimension(NEKO_HEX_NPTS), intent(in) :: pts
    integer :: i

    call this%init(id, NEKO_HEX_GDIM, NEKO_HEX_NPTS)
    
    do i = 1, NEKO_HEX_NPTS
       this%pts(i) = pts(i)
    end do

  end subroutine hex_init

  function hex_diameter(this) result(res)
    class(hex_t), intent(in) :: this
    real(kind=dp) :: res
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
             if (this%pts(i) .ne. other%pts(i)) then
                return
             end if
          end do
          res = .true.
       end if
    end select

  end function hex_equal
  
end module hex
