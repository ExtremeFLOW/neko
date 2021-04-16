!> Implements a point
!
module point
  use num_types
  use math
  use entity
  implicit none
  private

  !> A point in \f$ R^d \f$ with coordinates \f$ (x,y,z)\f$
  type, extends(entity_t), public ::  point_t
     real(kind=dp), dimension(3) :: x
   contains
     procedure :: point_eq
     procedure :: point_ne
     procedure :: point_lt
     procedure :: point_gt
     procedure :: point_assign
     generic :: operator(.eq.) => point_eq
     generic :: operator(.ne.) => point_ne
     generic :: operator(.lt.) => point_lt
     generic :: operator(.gt.) => point_gt
     generic :: assignment(=) => point_assign
  end type point_t

  !> Defines a pointer to a point type
  type, public ::  point_ptr
     type(point_t), pointer :: p
  end type point_ptr

  interface point_t
     module procedure point_init, point_init_xyz
  end interface point_t

contains
  
  !> Initialize a point from an array @a x of \f$ (x,y,z) \f$ coordinates
  function point_init(x, id) result(this)
    real(kind=dp), dimension(3), intent(in) :: x
    integer, optional, intent(inout) :: id
    type(point_t) :: this

    if (present(id)) then
       call this%set_id(id)
    else
       call this%set_id(-1)
    end if
    
    this%x = x

  end function point_init

  !> Initialize a point from \f$ (x,y,z) \f$ coordinates
  function point_init_xyz(x, y, z, id) result(this)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: y
    real(kind=dp), intent(in) :: z
    integer, optional, intent(inout) :: id
    type(point_t) :: this

    if (present(id)) then
       call this%set_id(id)
    else
       call this%set_id(-1)
    end if

    this%x(1) = x
    this%x(2) = y
    this%x(3) = z

  end function point_init_xyz
  
  !> Assigns coordinates @a x to a point
  subroutine point_assign(this, x)
    class(point_t), intent(inout) :: this
    real(kind=dp), dimension(3), intent(in) :: x

    this%x = x

  end subroutine point_assign

  !> Check if \f$ p_{1} = p_{2} \f$
  !! @note this only checks coordinates
  pure function point_eq(p1, p2) result(res)
    class(point_t), intent(in) :: p1
    class(point_t), intent(in) :: p2
    logical :: res

    if (abscmp(p1%x(1), p2%x(1)) .and. &
         abscmp(p1%x(2), p2%x(2)) .and. &
         abscmp(p1%x(3), p2%x(3))) then
       res = .true.
    else
       res = .false.
    end if
    
  end function point_eq

  !> Check if \f$ p_{1} \neq p_{2} \f$
  !! @note this only checks coordinates
  pure function point_ne(p1, p2) result(res)
    class(point_t), intent(in) :: p1
    class(point_t), intent(in) :: p2
    logical :: res

    if (.not. abscmp(p1%x(1), p2%x(1)) .and. &
         .not. abscmp(p1%x(2), p2%x(2)) .and. &
         .not. abscmp(p1%x(3), p2%x(3))) then
       res = .true.
    else
       res = .false.
    end if
    
  end function point_ne
  
  !> Check if \f$ p_{1} < p_{2} \f$
  !! @note this only checks coordinates
  pure function point_lt(p1, p2) result(res)
    class(point_t), intent(in) :: p1
    class(point_t), intent(in) :: p2
    logical :: res

    if (p1%x(1) .lt. p2%x(1) .or. &
         (abscmp(p1%x(1), p2%x(1)) .and. &
         (p1%x(2) .lt. p2%x(2) .or. &
         (abscmp(p1%x(2), p2%x(2)) .and. p1%x(3) .lt. p2%x(3))))) then
       res = .true.
    else
       res = .false.
    end if
    
  end function point_lt

  !> Check if \f$ p_{1} > p_{2} \f$
  !! @note this only checks coordinates
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
