! This tutorial demonstrates how to work with fields and vectors in Neko. It
! covers the following key points:
!
! - Initializing and working with `field_t` and `vector_t` types.
! - Performing mathematical operations on fields and vectors, including:
!   - Using overloaded operators for basic operations.
!   - Accessing and modifying field data directly on the CPU.
!   - Using low-level math routines for device and CPU compatibility.
!   - Leveraging the `field_math` module for backend-agnostic computations.
! - Demonstrating the lifecycle of user-defined objects:
!   - Initializing fields and vectors in `user_init_modules`.
!   - Performing custom calculations in `user_check`.
!   - Cleaning up allocated resources in `user_finalize_modules`.
!
! This tutorial highlights the flexibility of Neko's field and vector types,
! as well as the tools provided for efficient and portable numerical
! computations across different backends.
module user
  use neko
  implicit none


  ! One of the most important types in Neko is the field_t type. It is used to
  ! represent the solution fields in the simulation.
  type(field_t) :: my_field

  ! A vector is essentially a 1D array that handles the storage and association
  ! between data on the host and the device.
  type(vector_t) :: vec

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup
    user%initialize => initialize
    user%compute => compute
    user%finalize => finalize
  end subroutine user_setup

  ! We will use the user_startup routine to manipulate the end time.
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call params%add("case.end_time", 0.0_rp)
  end subroutine startup

  ! The user_init_modules routine is called after all the objects used to run
  ! the simulation are created. Thus, we can initialize user variables that,
  ! for example, depend on the mesh or the fluid solver.
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    ! Like almost all other Neko types, the field_t has an init method that is
    ! effectively a constructor. Typically, we use a dofmap_t object and a
    ! a name to call the init method. The dofmap_t is usually fetched from the
    ! solver, like a `fluid_scheme_t`. To do this, we make use of the special
    ! neko_user_access object, which is a singleton that provides access to
    ! various objects used in the simulation.
    call my_field%init(neko_user_access%case%fluid%dm_Xh, "my_field")

    ! The actual values of the field are stored in the x array or the x_d
    ! pointer. The x_d pointer is used to access the data on the GPU. The x
    ! array has a rank of 4, corresponding to GLL points in 3 dimensions, and
    ! the number of elements in the mesh (i,j,k,e).

    ! Some operators are overloaded for the field_t type. For example, we can
    ! assign it to a scalar. At construction, the values are set to zero.
    my_field = 1.0_rp


    ! The situation for vector is similar, but we only need to provide the size
    ! of the array we want in the constructor.
    call vec%init(50)

    vec = 1.0_rp

  end subroutine initialize

  ! This routine is run at the end of each time step. It is mean to hold
  ! arbitrary calculations that the user wants to do.
  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    integer :: i

    ! Let us consider doing some simple math on our field. Here we just add 1
    ! to the value at each GLL point. Note that we use a linear index to loop
    ! over the entire field, and use %x(i,1,1,1). This is very common in Neko,
    ! and relies on contiguous memory layout. Sometimes you will see 4D arrays
    ! passed to dummy arguments that are 1D, 2D, etc. This is possible for the
    ! same reason.
    do i = 1, my_field%size() ! <- Number of GLL points.
       my_field%x(i,1,1,1) = my_field%x(i,1,1,1) * 2.0_rp
    end do

    ! There is an issue with the above code: it only works on the CPU. You can
    ! not write GPU kernels as a user. Therefore, Neko has a whole range of
    ! low-level math functions that can be used. You can find them in math.f90.
    ! Corresponding routines for the GPU are in device_math.f90, and
    ! field_math.f90 contains wrappers that dispatch correctly to either a CPU
    ! or device routine based on the current backend. The names of the routines
    ! are taken from Nek5000, so if you are familiar with that code, you will
    ! find it easy to use.

    ! A generic way to perform the computation above would be
    if (NEKO_BCKND_DEVICE .eq. 1) then ! <- Check if we are on the device.
       call device_cmult(my_field%x_d, 2.0_rp, my_field%size())
    else
       call cmult(my_field%x, 2.0_rp, my_field%size())
    end if

    ! Alternatively, (and more concisely) we can use the field_math module,
    call field_cmult(my_field, 2.0_rp)

    ! For vector_t there are no math wrappers, so the latter approach is the
    ! one to go for. The vector_t does overload some operators though, so use
    ! those when possible.

  end subroutine compute

  ! If you declare objects at module scope that need to be destroyed, you can
  ! do it here. In our case it is the my_field field.
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    ! All types the allocate memory in Neko have a free method, which is a
    ! destructor.
    call my_field%free()

  end subroutine finalize

end module user
