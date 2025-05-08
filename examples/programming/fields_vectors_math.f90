! A template user file containing the user-defined functions
!
module user
  use neko
  implicit none


  ! One of the most important types in Neko is the field_t type. It is used to
  ! represent the solution fields in the simulation.
  type(field_t) :: vel_mag

  ! A vector is essentially a 1D array that handles the storage and association
  ! between data on the host and the device.
  type(vector_t) :: vec

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%user_startup => user_startup
    user%user_init_modules => user_init_modules
    user%user_check => user_check
    user%user_finalize_modules => user_finalize_modules
  end subroutine user_setup

  ! We will use the user_startup routine to inspect and manipulate the end time.
  subroutine user_startup(params)
    type(json_file), intent(inout) :: params

    call params%add("case.end_time", 0.0_rp)
  end subroutine user_startup

  ! The user_init_modules routine is called after all the objects used to run
  ! the simulation are created. Thus, we can initialize user variables that,
  ! for example, depend on the mesh or the fluid solver.
  subroutine user_init_modules(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u, v, w, p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    ! Like almost all other Neko types, the field_t has an init method that is
    ! effectively a constructor. Typically, we use a dofmap_t object and a
    ! a name to call the init method. The dofmap_t is usually fetched from the
    ! solver, like a `fluid_scheme_t`.
    call vel_mag%init(coef%dof, "vel_mag")

    ! The actual values of the field are stored in the x array or the x_d
    ! pointer. The x_d pointer is used to access the data on the GPU. The x
    ! array has a rank of 4, corresponding to GLL points in 3 dimensions, and
    ! the number of elements in the mesh (i,j,k,e). 

    ! Some operators are overloaded for the field_t type. For example, we can
    ! assign it to a scalar. At construction, the values are set to zero.
    vel_mag = 1.0_rp


    ! The situation for vector is similar, but we only need to provide the size
    ! of the array we want in the constructor.
    call vec%init(50)

    vec = 1.0_rp

  end subroutine user_init_modules

  ! This routine
  subroutine user_check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), intent(inout) :: u, v, w, p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

  end subroutine user_check

  ! If you declare objects at module scope that need to be destroyed, you can
  ! do it here. In our case it is the vel_mag field.
  subroutine user_finalize_modules(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    ! All types the allocate memory in Neko have a free method, which is a 
    ! destructor.
    call vel_mag%free()

  end subroutine user_finalize_modules

end module user
