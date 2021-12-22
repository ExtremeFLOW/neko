!> Defines a mean field
!
module mean_field
  use stats_quant
  use num_types
  use field
  use math
  implicit none
  
  type, extends(stats_quant_t) ::  mean_field_t
     type(field_t), pointer :: f => null()
     type(field_t) :: mf
     real(kind=rp) :: time
   contains
     procedure, pass(this) :: init => mean_field_init
     procedure, pass(this) :: free => mean_field_free
     procedure, pass(this) :: update => mean_field_update
  end type mean_field_t

contains

  !> Initialize a mean field for a field @a f
  subroutine mean_field_init(this, f)
    class(mean_field_t), intent(inout) :: this
    type(field_t), intent(inout), target :: f
    character(len=80) :: name
    
    call this%free()

    this%f => f
    this%time = 0.0_rp

    name = 'mean_'//trim(f%name)

    call field_init(this%mf, f%dof, name)

  end subroutine mean_field_init

  !> Deallocates a mean field
  subroutine mean_field_free(this)
    class(mean_field_t), intent(inout) :: this

    if (associated(this%f)) then
       nullify(this%f)
    end if

    call field_free(this%mf)

  end subroutine mean_field_free

  !> Update a mean field
  subroutine mean_field_update(this, k)
    class(mean_field_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k !< Time since last sample

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_cmult(this%mf%x_d, this%time, size(this%mf%x))
       call device_add2s2(this%mf%x_d, this%f%x_d, k, size(this%mf%x))
       this%time = this%time + k
       call device_cmult(this%mf%x_d, 1.0_rp / this%time, size(this%mf%x))
    else
       this%mf%x = this%mf%x * this%time
       call add2s2(this%mf%x, this%f%x, k, this%mf%dof%n_dofs)
       this%time = this%time + k
       this%mf%x = this%mf%x / this%time
    end if
       
  end subroutine mean_field_update
  
end module mean_field

