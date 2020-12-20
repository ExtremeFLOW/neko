module element
  use num_types
  use entity
  use tuple
  use point
  implicit none
  private 

  !> Base type for an element
  !! @details An element is a collection of @a npts_ points forming an
  !! element of dimension @a gdim_
  type, public, extends(entity_t), abstract :: element_t
     integer, private :: gdim_              !< Geometric dimension
     integer, private :: npts_              !< number of points
     type(point_ptr), allocatable :: pts(:) !< Points of an element 
   contains
     procedure, pass(this) :: element => element_init
     procedure, pass(this) :: free => element_free
     procedure, pass(this) :: gdim => element_gdim
     procedure, pass(this) :: npts => element_npts
     procedure, pass(this) :: p => element_point 
     procedure, pass(this) :: n_points => element_npts
     procedure, pass(this), non_overridable :: element_point
     procedure(element_equal), pass(this), deferred :: equal
     procedure(element_diameter), pass(this), deferred :: diameter
     procedure(element_centroid), pass(this), deferred :: centroid
     procedure(element_facet_id), pass(this), deferred :: facet_id
     procedure(element_facet_order), pass(this), deferred :: facet_order
  end type element_t

  abstract interface
     function element_diameter(this) result(res)
       import :: element_t
       import :: dp
       class(element_t), intent(in) :: this
       real(kind=dp) :: res
     end function element_diameter
  end interface

  abstract interface
     function element_centroid(this) result(res)
       import :: element_t
       import :: point_t
       class(element_t), intent(in) :: this
       type(point_t) :: res
     end function element_centroid
  end interface

  abstract interface
     pure function element_equal(this, other) result(res)
       import :: element_t
       class(element_t), intent(in) :: this
       class(element_t), intent(in) :: other
       logical :: res
     end function element_equal
  end interface

  abstract interface
     subroutine element_facet_id(this, t, side) 
       import :: element_t
       import :: tuple_t
       class(element_t), intent(in) :: this
       class(tuple_t), intent(inout) :: t
       integer, intent(in) :: side
     end subroutine element_facet_id
  end interface

  abstract interface
     subroutine element_facet_order(this, t, side) 
       import :: element_t
       import :: tuple_t
       class(element_t), intent(in) :: this
       class(tuple_t), intent(inout) :: t
       integer, intent(in) :: side
     end subroutine element_facet_order
  end interface
contains

  !> Create an element with @a npts
  subroutine element_init(this, id, gdim, npts)
    class(element_t), intent(inout)  :: this
    integer, intent(inout) :: id
    integer, intent(in) :: gdim
    integer, intent(in) :: npts

    call this%free()

    call this%set_id(id)
    
    this%gdim_ = gdim
    this%npts_ = npts

    allocate(this%pts(this%npts_))

  end subroutine element_init

  !> Deallocate an element
  subroutine element_free(this)
    class(element_t), intent(inout) :: this

    if (allocated(this%pts)) then
       deallocate(this%pts)
    end if
    
  end subroutine element_free
  
  !> Get the geometric dimension of an element
  pure function element_gdim(this) result(gdim)
    class(element_t), intent(in) :: this
    integer :: gdim
    gdim = this%gdim_
  end function element_gdim

  !> Get the number of points in an element
  pure function element_npts(this) result(npts)
    class(element_t), intent(in) :: this
    integer :: npts
    npts = this%npts_
  end function element_npts

  !> Return a pointer to point @a i of the element
  function element_point(this, i) result(pt)
    class(element_t), intent(in) :: this
    integer, intent(in) :: i
    type(point_t), pointer :: pt
    pt => this%pts(i)%p
  end function element_point

end module element
