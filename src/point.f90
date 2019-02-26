!> Implements a point
!
module point
  use num_types
  use entity
  implicit none
  private

  type, extends(entity_t), public ::  point_t
     real(kind=dp), private :: x(3)
   contains
     procedure :: point_eq
     procedure :: point_lt
     procedure :: point_gt
     procedure :: point_assign
     generic :: operator(.eq.) => point_eq
     generic :: operator(.lt.) => point_lt
     generic :: operator(.gt.) => point_gt
     generic :: assignment(=) => point_assign
  end type point_t

  interface point_t
     module procedure point_init
  end interface point_t

contains
  
  function point_init(x, id) result(this)
    real(kind=dp), dimension(3), intent(in) :: x
    integer, intent(inout) :: id
    type(point_t) :: this

    call this%init(id)

    this%x = x

  end function point_init
  
  subroutine point_assign(this, x)
    class(point_t), intent(inout) :: this
    real(kind=dp), dimension(3), intent(in) :: x
  end subroutine point_assign

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
