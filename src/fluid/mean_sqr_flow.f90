!> Defines a mean squared flow field
module mean_sqr_flow
  use mean_sqr_field
  implicit none

  type :: mean_sqr_flow_t
     type(mean_sqr_field_t) :: uu
     type(mean_sqr_field_t) :: vv
     type(mean_sqr_field_t) :: ww
     type(mean_sqr_field_t) :: pp
   contains
     procedure, pass(this) :: init => mean_sqr_flow_init
     procedure, pass(this) :: free => mean_sqr_flow_free
  end type mean_sqr_flow_t
  
contains
  
  !> Initialize a mean squared flow field
  subroutine mean_sqr_flow_init(this, u, v, w, p)
    class(mean_sqr_flow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p

    call this%free()

    call this%uu%init(u)
    call this%vv%init(v)
    call this%ww%init(w)
    call this%pp%init(p)
    
  end subroutine mean_sqr_flow_init


  !> Deallocates a mean squared flow field
  subroutine mean_sqr_flow_free(this)
    class(mean_sqr_flow_t), intent(inout) :: this

    call this%uu%free()
    call this%vv%free()
    call this%ww%free()
    call this%pp%free()

  end subroutine mean_sqr_flow_free
  
end module mean_sqr_flow
