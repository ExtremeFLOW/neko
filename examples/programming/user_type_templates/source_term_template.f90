!> Template for a user-defined source term.
module source_term_template
  use coefs, only : coef_t
  use field_list, only : field_list_t
  use json_module, only : json_file
  use num_types, only : rp
  use source_term, only : source_term_t, source_term_allocate, &
       register_source_term
  use time_state, only : time_state_t
  implicit none
  private

  type, extends(source_term_t) :: source_term_template_t
   contains
     procedure :: init => source_term_template_init
     procedure :: free => source_term_template_free
     procedure :: compute_ => source_term_template_compute
  end type source_term_template_t

  public :: source_term_template_register_types

contains

  subroutine source_term_template_init(this, json, fields, coef, variable_name)
    class(source_term_template_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), target, intent(in) :: fields
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: variable_name

    call this%free()
    ! TODO: Read source-term-specific parameters from json.
    call this%init_base(fields, coef, 0.0_rp, huge(0.0_rp))
  end subroutine source_term_template_init

  subroutine source_term_template_free(this)
    class(source_term_template_t), intent(inout) :: this
    ! TODO: Release fields owned by this derived type.
    call this%free_base()
  end subroutine source_term_template_free

  subroutine source_term_template_compute(this, time)
    class(source_term_template_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    ! TODO: Add this source term to this%fields.
  end subroutine source_term_template_compute

  subroutine source_term_template_register_types()
    procedure(source_term_allocate), pointer :: allocator
    allocator => source_term_template_allocate
    call register_source_term('source_term_template', allocator)
  end subroutine source_term_template_register_types

  subroutine source_term_template_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(source_term_template_t :: obj)
  end subroutine source_term_template_allocate
end module source_term_template
