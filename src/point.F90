!> Implements a point
!
module point
  use num_types
  implicit none
  private

  type, public ::  point_t
     real(kind=dp), private :: x(3)
   contains
     procedure :: point_eq
     procedure :: point_lt
     procedure :: point_gt
     procedure :: point_assign
     procedure :: point_write
     generic :: operator(.eq.) => point_eq
     generic :: operator(.lt.) => point_lt
     generic :: operator(.gt.) => point_gt
     generic :: assignment(=) => point_assign
#ifdef HAVE_DERIVED_TYPE_IO
     generic :: write(formatted) => point_write
#endif
  end type point_t

  interface point_t
     module procedure point_init
  end interface point_t

contains
  
  pure function point_init(x) result(this)
    real(kind=dp), dimension(3), intent(in) :: x
    type(point_t) :: this

    this%x = x
  end function point_init
  
  subroutine point_assign(this, x)
    class(point_t), intent(inout) :: this
    real(kind=dp), dimension(3), intent(in) :: x
  end subroutine point_assign

  !> Write point information
  subroutine point_write(this, unit, iotype, v_list, iostat, iomsg)
    class(point_t), intent(in) :: this
    integer, intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg
    
    write(*,*) this%x

    iostat = 0
    
  end subroutine point_write

  !> Check if \f$ p_{1} = p_{2} \f$
  pure function point_eq(p1, p2) result(res)
    class(point_t), intent(in) :: p1
    class(point_t), intent(in) :: p2
    logical :: res

    if (p1%x(1) .eq. p2%x(1) .and. &
         p1%x(2) .eq. p2%x(2) .and. &
         p1%x(3) .eq. p2%x(3)) then
       res = .true.
    else
       res = .false.
    end if
    
  end function point_eq
  
  !> Check if \f$ p_{1} < p_{2} \f$
  pure function point_lt(p1, p2) result(res)
    class(point_t), intent(in) :: p1
    class(point_t), intent(in) :: p2
    logical :: res

    if (p1%x(1) .lt. p2%x(1) .or. &
         (p1%x(1) .eq. p2%x(1) .and. &
         (p1%x(2) .lt. p2%x(2) .or. &
         (p1%x(2) .eq. p2%x(2) .and. p1%x(3) .lt. p2%x(3))))) then
       res = .true.
    else
       res = .false.
    end if
    
  end function point_lt

  !> Check if \f$ p_{1} > p_{2} \f$
  pure function point_gt(p1, p2) result(res)
    class(point_t), intent(in) :: p1
    class(point_t), intent(in) :: p2
    logical :: res
    
    if (point_lt(p1, p2)) then
       res = .false.
    else
       res = .true.
    end if

  end function point_gt  
  
end module point
