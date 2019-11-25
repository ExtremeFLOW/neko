module entity
  implicit none
  private
  
  !> Base type for an entity
  type, public, abstract :: entity_t 
     integer, private :: id_ = -1 !< Entity index
     integer, private :: glb_id_ = -1 !< Global entitiy index
   contains
     procedure, pass(this) :: id => entity_id
     procedure, pass(this) :: glb_id => entity_glb_id
     procedure, pass(this) :: set_id => entity_set_id     
     procedure, pass(this) :: set_glb_id => entity_set_glb_id
  end type entity_t

contains

  !> Return the index of an entity
  pure function entity_id(this) result(id)
    class(entity_t), intent(in) :: this
    integer :: id
    id = this%id_
  end function entity_id

  !> Return the global index of an entity
  pure function entity_glb_id(this) result(glb_id)
    class(entity_t), intent(in) :: this
    integer :: glb_id
    glb_id = this%glb_id_
  end function entity_glb_id

  !> Update the index of an entity
  subroutine entity_set_id(this, id)
    class(entity_t), intent(inout) :: this
    integer, intent(in) :: id
    this%id_ = id
  end subroutine entity_set_id

  !> Update the global index of an entitiy
  subroutine entity_set_glb_id(this, glb_id)
    class(entity_t), intent(inout) :: this
    integer, intent(in) :: glb_id
    this%glb_id_ = glb_id
  end subroutine entity_set_glb_id
    
end module entity
