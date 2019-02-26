module entity
  implicit none
  private
  
  !> Base type for an entity
  type, public, abstract :: entity_t 
     integer, private :: id_ = -1 !< Entity index
   contains
     procedure, pass(this) :: entity_init
     procedure, public :: id => entity_id
     generic, public :: init => entity_init
  end type entity_t

contains

  subroutine entity_init(this, id)
    class(entity_t), intent(inout) :: this
    integer, intent(in) :: id
    this%id_ = id
  end subroutine entity_init

  !> Return id if an entity
  pure function entity_id(this) result(id)
    class(entity_t), intent(in) :: this
    integer :: id
    id = this%id_
  end function entity_id
    
end module entity
