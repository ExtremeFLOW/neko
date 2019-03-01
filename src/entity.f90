module entity
  implicit none
  private
  
  !> Base type for an entity
  type, public, abstract :: entity_t 
     integer, private :: id_ = -1 !< Entity index
   contains
     procedure, pass(this) :: id => entity_id
     procedure, pass(this) :: set_id => entity_set_id     
  end type entity_t

contains

  !> Return @a id_ of an entity
  pure function entity_id(this) result(id)
    class(entity_t), intent(in) :: this
    integer :: id
    id = this%id_
  end function entity_id

  !> Update @a id_ of an entity
  subroutine entity_set_id(this, id)
    class(entity_t), intent(inout) :: this
    integer, intent(in) :: id
    this%id_ = id
  end subroutine entity_set_id
    
end module entity
