!> Defines a mean square field
!
module mean_sqr_field
  use mean_field
  use num_types
  use field
  use math
  implicit none

  type, extends(mean_field_t) :: mean_sqr_field_t
   contains
     procedure, pass(this) :: update => mean_sqr_field_update
  end type mean_sqr_field_t

contains    

  !> Update a mean sqr field
  subroutine mean_sqr_field_update(this, k)
    class(mean_sqr_field_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k !< Time since last sample

    this%mf%x = this%mf%x * this%time
    call addsqr2s2(this%mf%x, this%f%x, k, this%mf%dof%n_dofs)
    this%time = this%time + k
    this%mf%x = this%mf%x / this%time

  end subroutine mean_sqr_field_update

end module mean_sqr_field
