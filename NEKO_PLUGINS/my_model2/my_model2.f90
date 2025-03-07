module my_model2
  use num_types, only : rp
  use les_model, only : les_model_t, register_les_model, les_model_allocate
  use field, only : field_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use json_utils, only : json_get_or_default
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use logger, only : LOG_SIZE, neko_log
  use field_math, only : field_cfill

  implicit none

  !> Implements a dummy user model
  type, public, extends(les_model_t) :: my_model2_t
     !> Model constant, defaults to 0.07.
     real(kind=rp) :: c
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => my_model2_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          my_model2_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => my_model2_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => my_model2_compute
  end type my_model2_t

contains
  !> Constructor.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param json A dictionary with parameters.
  subroutine my_model2_init(this, fluid, json)
    class(my_model2_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name
    real(kind=rp) :: c
    character(len=:), allocatable :: delta_type
    logical :: if_ext
    character(len=LOG_SIZE) :: log_buf

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")
    call json_get_or_default(json, "extrapolation", if_ext, .true.)

    call neko_log%section('LES model')
    write(log_buf, '(A)') 'Model : CUSTOM MODEL 2 TEST'
    call neko_log%message(log_buf)
    call neko_log%end_section()

    call my_model2_init_from_components(this, fluid, c, nut_name, &
         delta_type, if_ext)
  end subroutine my_model2_init

  !> Constructor from components.
  !! @param fluid The fluid_scheme_base_t object.
  !! @param c The model constant.
  !! @param nut_name The name of the SGS viscosity field.
  !! @param delta_type The type of filter size.
  !! @param if_ext Whether trapolate the velocity.
  subroutine my_model2_init_from_components(this, fluid, c, nut_name, &
       delta_type, if_ext)
    class(my_model2_t), intent(inout) :: this
    class(fluid_scheme_base_t), intent(inout), target :: fluid
    real(kind=rp) :: c
    character(len=*), intent(in) :: nut_name
    character(len=*), intent(in) :: delta_type
    logical, intent(in) :: if_ext

    call this%free()

    call this%init_base(fluid, nut_name, delta_type, if_ext)
    this%c = c

  end subroutine my_model2_init_from_components

  !> Destructor for the les_model_t (base) class.
  subroutine my_model2_free(this)
    class(my_model2_t), intent(inout) :: this

    call this%free_base()
  end subroutine my_model2_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine my_model2_compute(this, t, tstep)
    class(my_model2_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    write(*,*) "Executing custom SGS model 2!!"
    ! Just a value to test
    call field_cfill(this%nut, 5e-3_rp)

  end subroutine my_model2_compute

  !> The allocator for my_model2_t
  subroutine custom_allocator(obj)
    class(les_model_t), allocatable, intent(inout) :: obj

    allocate(my_model2_t::obj)
  end subroutine custom_allocator

  !> module name + register_types routine
  !! Can register all the custom types from the module here!
  subroutine my_model2_register_types()
    procedure(les_model_allocate), pointer :: allocator_ptr

    allocator_ptr => custom_allocator  ! Assign procedure pointer

    write(*,*) "Registering the my_model2_t SGS model"
    call register_les_model("my_model2", allocator_ptr)
    write(*,*) "Done"

  end subroutine my_model2_register_types


end module my_model2
