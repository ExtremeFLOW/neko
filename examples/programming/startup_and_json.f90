! This module demonstrates how to define user-specific routines in Neko. It
! provides an example of how to interact with the JSON parameter dictionary
! used for simulation configuration. The key points covered in this tutorial
! are:
!
! - Registering user-defined functions using the `user_setup` routine.
! - Using the `user_startup` routine to inspect and manipulate the JSON
!   parameter dictionary before the simulation starts.
! - Extracting parameters from the JSON file using `json_get` and
!   `json_get_or_default`.
! - Printing JSON objects using the `print` method.
! - Adding or modifying parameters in the JSON file using the `add` method.
! - Saving the JSON parameter dictionary for later use in the module.
!
! This tutorial highlights the flexibility of the JSON-based configuration
! system in Neko and demonstrates how to customize simulations by interacting
! with the JSON parameter dictionary.


! The user module always needs to be named "user"!
module user
  ! This use statement populates the scope of the module with all the types and
  ! public procedures defined in Neko. So, this is convenient, but your scope
  ! does get very polluted. It is generallys a good idea to `use` only the
  ! modules you need, and moreover specify what exactly you want from those
  ! modules using `use ..., only :`.
  use neko
  implicit none

  ! A module-scope variable to hold a copy of the case file JSON parameter
  ! dictionary.
  type(json_file) :: case_params

  ! Some variable that we want to use in the user code.
  real(kind=rp) :: some_variable

contains

  ! This is a special routine that registers the user-defined functions with
  ! Neko. Based on the inerface defined in user_intf.f90, we can register our
  ! user-defined implementations of the various routines. You do this by
  ! assigning procedure pointers to subroutines in the user module. Here, we
  ! register only the startup routine.
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup

  end subroutine user_setup

  ! The startup routine, provides the possibility to inspect and manipulate
  ! the JSON parameter dictionary before the simulation starts. The routine is
  ! called very early, before any of the solvers or simulation components are
  ! initialized. This is also a good place to set up some constants that are
  ! needed in the user code.
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    ! Some auxillary variables to extract various type of parameters from the
    ! JSON
    real(kind=rp) :: some_real
    ! Note that we need an allocatable when we extract strings or arrays.
    ! The json_fortran library will allocate the string or array for us.
    character(len=:), allocatable :: some_string
    real(kind=rp), allocatable :: some_real_array(:)
    ! We can extract an entire JSON object into a new json_file.
    type(json_file) :: some_json_object

    ! Assign our constant to something. Can be based on JSON parameters
    ! extracted below. A common use case is to set material properties based
    ! on custom JSON entries in the case file, e.g. some non-dimensional number.
    some_variable = 1.0_rp

    ! Neko relies on json_fortran to parse and manipulate the JSON file.
    ! In addition, there are some convenience subroutines in the json_utils
    ! module.

    ! Extract a string parameter. This will through an error if
    ! case.fluid.scheme does not exist in the JSON file.
    call json_get(params, "case.fluid.scheme", some_string)


    ! Extract a real parameter, providing a default if it does not exist.
    call json_get_or_default(params, "case.fluid.Re", some_real, 100.0_rp)

    ! Extract the object with the velocity initial conditions.
    call json_get(params, "case.fluid.initial_condition", &
         some_json_object)

    ! We can print the contents to the console.
    call some_json_object%print()

    ! We can use the json_fotran methods to manipulate the JSON.
    ! Here we add an additional parameter.
    call params%add("case.fluid.initial_condition.my_param", 3.0_rp)

    ! Add can also be used to modify existing parameters.
    some_real_array = [1.0_rp, 2.0_rp, 3.0_rp]
    call params%add("case.fluid.initial_condition.value", some_real_array)
    call params%add("case.end_time", 0.0_rp)

    ! Show the updated part of the JSON file.
    call json_get(params, "case.fluid.initial_condition", &
         some_json_object)
    call some_json_object%print()

    ! There are a few other json utilities in neko, as well as more methods in
    ! the json_fortran library. So, you can pretty much do anything you want
    ! with the case file, if you so need.

    ! We can save the params into our module-scope variable for later use.
    case_params = params

  end subroutine startup


end module user
