!> Template for a user-defined LES model.
module les_model_template
  use fluid_scheme_base, only : fluid_scheme_base_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use les_model, only : les_model_t, les_model_allocate, register_les_model
  use num_types, only : rp
  implicit none
  private

  type, extends(les_model_t) :: les_model_template_t
   contains
     procedure :: init => les_model_template_init
     procedure :: free => les_model_template_free
     procedure :: compute => les_model_template_compute
  end type les_model_template_t

  public :: les_model_template_register_types

contains

  subroutine les_model_template_init(this, fluid, json)
    class(les_model_template_t), intent(inout) :: this
    class(fluid_scheme_base_t), target, intent(inout) :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name, delta_type
    logical :: extrapolation

    call json_get_or_default(json, 'nut_field', nut_name, 'nut')
    call json_get_or_default(json, 'delta_type', delta_type, 'pointwise')
    call json_get_or_default(json, 'extrapolation', extrapolation, .false.)
    call this%free()
    call this%init_base(fluid, nut_name, delta_type, extrapolation)
    ! TODO: Read model-specific parameters from json.
  end subroutine les_model_template_init

  subroutine les_model_template_free(this)
    class(les_model_template_t), intent(inout) :: this
    ! TODO: Release fields owned by this derived type.
    call this%free_base()
  end subroutine les_model_template_free

  subroutine les_model_template_compute(this, t, tstep)
    class(les_model_template_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    ! TODO: Compute this%nut from the resolved flow field.
  end subroutine les_model_template_compute

  subroutine les_model_template_register_types()
    procedure(les_model_allocate), pointer :: allocator
    allocator => les_model_template_allocate
    call register_les_model('les_model_template', allocator)
  end subroutine les_model_template_register_types

  subroutine les_model_template_allocate(obj)
    class(les_model_t), allocatable, intent(inout) :: obj
    allocate(les_model_template_t :: obj)
  end subroutine les_model_template_allocate
end module les_model_template
