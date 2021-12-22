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

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_cmult(this%mf%x_d, this%time, size(this%mf%x))
       call device_addsqr2s2(this%mf%x_d, this%f%x_d, k, size(this%mf%x))
       this%time = this%time + k
       call device_cmult(this%mf%x_d, 1.0_rp / this%time, size(this%mf%x))
    else
       this%mf%x = this%mf%x * this%time
       call addsqr2s2(this%mf%x, this%f%x, k, this%mf%dof%n_dofs)
       this%time = this%time + k
       this%mf%x = this%mf%x / this%time
    end if

  end subroutine mean_sqr_field_update

end module mean_sqr_field
