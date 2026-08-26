!> Template for a user-defined wall model.
module wall_model_template
  use coefs, only : coef_t
  use json_module, only : json_file
  use num_types, only : rp
  use user_intf, only : user_t
  use wall_model, only : wall_model_t, wall_model_allocate, &
       register_wall_model
  implicit none
  private

  type, extends(wall_model_t) :: wall_model_template_t
   contains
     procedure :: init => wall_model_template_init
     procedure :: partial_init => wall_model_template_partial_init
     procedure :: finalize => wall_model_template_finalize
     procedure :: free => wall_model_template_free
     procedure :: compute => wall_model_template_compute
  end type wall_model_template_t

  public :: wall_model_template_register_types

contains

  subroutine wall_model_template_init(this, scheme_name, coef, msk, facet, json)
    class(wall_model_template_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:), facet(:)
    type(json_file), intent(inout) :: json

    call this%partial_init(coef, scheme_name, json)
    call this%finalize(msk, facet)
  end subroutine wall_model_template_init

  subroutine wall_model_template_partial_init(this, coef, scheme_name, json)
    class(wall_model_template_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: scheme_name
    type(json_file), intent(inout) :: json

    call this%partial_init_base(coef, scheme_name, json)
    ! TODO: Read wall-model-specific parameters from json.
  end subroutine wall_model_template_partial_init

  subroutine wall_model_template_finalize(this, msk, facet, bc_name, user)
    class(wall_model_template_t), intent(inout) :: this
    integer, intent(in) :: msk(:), facet(:)
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user

    if (present(bc_name) .and. present(user)) then
       call this%finalize_base(msk, facet, bc_name, user)
    else if (present(bc_name)) then
       call this%finalize_base(msk, facet, bc_name)
    else
       call this%finalize_base(msk, facet)
    end if
    ! TODO: Initialize arrays which depend on this%n_nodes.
  end subroutine wall_model_template_finalize

  subroutine wall_model_template_free(this)
    class(wall_model_template_t), intent(inout) :: this
    ! TODO: Release fields owned by this derived type.
    call this%free_base()
  end subroutine wall_model_template_free

  subroutine wall_model_template_compute(this, t, tstep)
    class(wall_model_template_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    ! TODO: Fill this%tau_x, this%tau_y, and this%tau_z.
    call this%compute_mag_field()
  end subroutine wall_model_template_compute

  subroutine wall_model_template_register_types()
    procedure(wall_model_allocate), pointer :: allocator
    allocator => wall_model_template_allocate
    call register_wall_model('wall_model_template', allocator)
  end subroutine wall_model_template_register_types

  subroutine wall_model_template_allocate(obj)
    class(wall_model_t), allocatable, intent(inout) :: obj
    allocate(wall_model_template_t :: obj)
  end subroutine wall_model_template_allocate
end module wall_model_template
