!> Defines a quadrilateral element
!
module quad
  use num_types
  use element
  implicit none
  private

  integer, parameter :: NEKO_QUAD_NPTS = 4 !< Number of points
  integer, parameter :: NEKO_QUAD_GDIM = 2 !< Geometric dimension

  type, public, extends(element_t) :: quad_t
   contains
     procedure, pass(this) :: quad_init
     procedure, pass(this) :: diameter => quad_diameter
     procedure, pass(this) :: equal => quad_equal
     generic, public :: init => quad_init
     generic :: operator(.eq.) => equal
  end type quad_t

contains
  
  subroutine quad_init(this, id, pts)
    class(quad_t), intent(inout) :: this
    integer, intent(inout) :: id
    integer, dimension(NEKO_QUAD_NPTS), intent(in) :: pts
    integer :: i

    call this%init(id, NEKO_QUAD_GDIM, NEKO_QUAD_NPTS)

    do i = 1, NEKO_QUAD_NPTS
       this%pts(i) = pts(i)
    end do

  end subroutine quad_init

  function quad_diameter(this) result(res)
    class(quad_t), intent(in) :: this
    real(kind=dp) :: res
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
             if (this%pts(i) .ne. other%pts(i)) then
                return
             end if
          end do
          res = .true.
       end if
    end select
    
  end function quad_equal
  
end module quad
