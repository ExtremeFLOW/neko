!> Template for a user-defined point zone.
module point_zone_template
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use num_types, only : rp
  use point_zone, only : point_zone_t, point_zone_allocate, &
       register_point_zone
  implicit none
  private

  type, extends(point_zone_t) :: point_zone_template_t
   contains
     procedure :: init => point_zone_template_init
     procedure :: free => point_zone_template_free
     procedure :: criterion => point_zone_template_criterion
  end type point_zone_template_t

  public :: point_zone_template_register_types

contains

  subroutine point_zone_template_init(this, json, size)
    class(point_zone_template_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    integer, intent(in) :: size
    character(len=:), allocatable :: name
    logical :: invert, full_elements

    call json_get(json, 'name', name)
    call json_get_or_default(json, 'invert', invert, .false.)
    call json_get_or_default(json, 'full_elements', full_elements, .false.)
    call this%init_base(size, name, invert, full_elements)
    ! TODO: Read criterion-specific parameters from json.
  end subroutine point_zone_template_init

  subroutine point_zone_template_free(this)
    class(point_zone_template_t), intent(inout) :: this
    ! TODO: Release fields owned by this derived type.
    call this%free_base()
  end subroutine point_zone_template_free

  pure function point_zone_template_criterion(this, x, y, z, j, k, l, e) &
       result(is_inside)
    class(point_zone_template_t), intent(in) :: this
    real(kind=rp), intent(in) :: x, y, z
    integer, intent(in) :: j, k, l, e
    logical :: is_inside
    ! TODO: Replace with the point-selection criterion.
    is_inside = .false.
  end function point_zone_template_criterion

  subroutine point_zone_template_register_types()
    procedure(point_zone_allocate), pointer :: allocator
    allocator => point_zone_template_allocate
    call register_point_zone('point_zone_template', allocator)
  end subroutine point_zone_template_register_types

  subroutine point_zone_template_allocate(obj)
    class(point_zone_t), allocatable, intent(inout) :: obj
    allocate(point_zone_template_t :: obj)
  end subroutine point_zone_template_allocate
end module point_zone_template
