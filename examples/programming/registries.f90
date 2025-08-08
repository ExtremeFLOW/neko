! This tutorial demonstrates how to use registries in Neko to manage fields and
! temporary data. It covers the following key points:
!
! - Registering user-defined routines using the `user_setup` subroutine.
! - Accessing and managing fields using the `neko_field_registry`:
!   - Retrieving existing fields registered by solvers (e.g., velocity fields).
!   - Registering new fields for use across modules or for output purposes.
! - Using the `neko_scratch_registry` to manage temporary fields:
!   - Requesting temporary fields for intermediate calculations.
!   - Reliquishing temporary fields to free resources.
!
module user
  use neko
  implicit none

contains


  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup
    user%initialize => initialize
    user%compute => compute

  end subroutine user_setup

  ! We will use the startup routine to manipulate the end time.
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call params%add("case.end_time", 0.0_rp)
  end subroutine startup

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: my_field_ptr

    ! Sometimes we need fields to be accessible across different modules. This
    ! can include fields that we create in the user module. To that end, Neko
    ! provides a special singleton object of type `field_registry_t`. The name
    ! of the singleton is `neko_field_registry`. Solvers register their solution
    ! fields there, for example.

    ! Here we can grab the x-velocity from the fluid_scheme_t, which is
    ! registered as "u" in the registry.
    my_field_ptr => neko_field_registry%get_field("u")

    ! We can also register new fields. For that we need a dofmap_t object, just
    ! like when we initialize a field. Once in the registry, the field is
    ! available for access and manipulation in other modules.
    ! Practical example: you may want to run some calculations on this field and
    ! then use the field_writer simcomp to output it to disk during the
    ! simulation. The field_writer looks for fields in the registry, given their
    ! names.
    call neko_field_registry%add_field(my_field_ptr%dof, "my_field")

  end subroutine initialize

  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: temp_index !<- For the scratch registry
    type(field_t), pointer :: temp_field_ptr !<- Will be our temporary field
    type(field_t), pointer :: u, v

    ! Sometimes we need some temporary fields to perform certain calculations.
    ! Instead of declaring them inside a subroutine, Neko provides yet another
    ! singleton: `neko_scratch_registry`, of type `scratch_registry_t`. The
    ! advantage of using it is saving memory and time to allocate things.

    ! We get a temporary by using a field_t pointer and an index.
    call neko_scratch_registry%request_field(temp_field_ptr, temp_index)

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")

    ! Perform some operations and use the temporary for something
    call field_cfill(temp_field_ptr, 1.0_rp)
    call field_add2(temp_field_ptr, u)
    v = temp_field_ptr

    ! At the end of the subroutine, the temporary field needs to be reliquished.
    ! This is where the integer index comes in.
    call neko_scratch_registry%relinquish_field(temp_index)

  end subroutine compute

end module user
