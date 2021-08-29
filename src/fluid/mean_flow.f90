!> Defines a mean flow field
module mean_flow
  use mean_field
  implicit none

  type :: mean_flow_t
     type(mean_field_t) :: u
     type(mean_field_t) :: v
     type(mean_field_t) :: w
     type(mean_field_t) :: p
   contains
     procedure, pass(this) :: init => mean_flow_init
     procedure, pass(this) :: free => mean_flow_free
  end type mean_flow_t

contains

  !> Initialize a mean flow field
  subroutine mean_flow_init(this, u, v, w, p)
    class(mean_flow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p

    call this%free()

    call this%u%init(u)
    call this%v%init(v)
    call this%w%init(w)
    call this%p%init(p)
    
  end subroutine mean_flow_init


  !> Deallocates a mean flow field
  subroutine mean_flow_free(this)
    class(mean_flow_t), intent(inout) :: this

    call this%u%free()
    call this%v%free()
    call this%w%free()
    call this%p%free()

  end subroutine mean_flow_free
  
end module mean_flow
