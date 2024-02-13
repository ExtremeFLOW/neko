# User File {#user-file}

## High-level structure

The current high-level structure of the user file is shown below.

```fortran
module user
  use neko
  implicit none

contains

  ! Register user defined functions here (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u

    
  end subroutine user_setup
 
end module user

```

The user file implements the `user` module. The `user` modules contains a subroutine named `user_setup`, which we use to interface the internal procedures defined in `src/common/user_intf.f90` with the subroutines that you will implement in your user file. Each user subroutine should be implemented under the `contains` statement, below `user_setup`.

A basic working example is demonstrated [below](#basic-working-example).

@note The above code snippet is the most basic code structure for the user file. Compiling it and running it would be equivalent to running the "vanilla" neko executable `bin/neko` in your local neko installation folder.

## Default user functions

The following user functions, if defined in the user file, will always be executed, regardless of what is set in the case file: 

- [user_init_modules](#user-file_init-and-final): For initializing user variables and objects
- [user_finalize_modules](#user-file_init-and-final): For finalizing, e.g freeing variables and terminating processes
- [user_check](#user-file_user-check): Executed at the end of every time step, for e.g. computing and/or outputting user defined quantities.
- [material_properties](#user-file_mat-prop): For computing and setting material properties such as `rho`, `mu`, `cp` and `lambda`.
- [user_mesh_setup](#user-file_user-mesh-setup): For applying any deformation to the element nodes of the mesh specified in the case file.
- [scalar_user_bc](#user-file_scalar-bc): For applying boundary conditions to the scalar, on all zones that are not already specified with uniform dirichlet values e.g. `d=1`. For more information on the scalar, see the [relevant section of the case file](#case-file_scalar).

### Initializing and finalizing {#user-file_init-and-final}

These two subroutines may be used to initialize/finalize any user defined variables, external objects, or processes. They are respectively executed right before/after the simulation time loop.

```.f90

  ! Initialize user variables or external objects
  subroutine initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    ! insert your initialization code here

  end subroutine initialize

  ! Finalize user variables or external objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    ! insert your code here

  end subroutine initialize
```

In the example above, the subroutines `initialize` and `finalize` contains the actual implementations, and both need to be interfaced to the internal procedures `user_init_modules` and `user_finalize_modules` in `user_setup`:

```.f90

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u

    u%user_init_modules => initialize 
    u%user_finalize_modules => finalize 

  end subroutine user_setup

```

@note `user_init_modules` and `user_finalize_modules` are independent of each other. Using one does not require the use of the other.

### Computing at every time step {#user-file_user-check}

This subroutine gets its name from its Nek5000 counterpart, `usrcheck`. It is executed at the end of every time step. It can be used for computing and/or outputting your own variables/quantities at every time step.
```.f90
  ! This is called at the end of every time step
  subroutine usercheck(t, tstep, u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param

    ! insert code below

  end subroutine usercheck

```

In the example above, the subroutine `usercheck` contains the actual implementation, and need to be registered by adding: 

```.f90
u%user_check => usercheck
```

to our `user_setup`.

### Setting material properties {#user-file_mat-prop}

This subroutine allow for more complex computations and setting of various material properties, such as `rho`, `mu` for the fluid and `cp`, `lambda` for the scalar. The example below is taken from the [rayleigh-benard-cylinder example](https://github.com/ExtremeFLOW/neko/blob/564686b127ff75a362a06126c6b23e9b4e21879e/examples/rayleigh-benard-cylinder/rayleigh.f90#L22C1-L38C41). 

```.f90
  
  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params
    real(kind=rp) :: Re

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.scalar.Pr", Pr)

    Re = sqrt(Ra / Pr)
    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

```

And of course not forgetting to register our function in `user_setup` by adding the following line:

```.f90
u%material_properties => set_material_properties
```

### Runtime mesh deformation {#user-file_user-mesh-setup}

This user function allows for the modification of the mesh at runtime, by acting on the element nodes of the mesh specified in the case file. The example below is taken from the [compression example](https://github.com/ExtremeFLOW/neko/blob/a0613606360240e5059e65d6d98f4a57cf73e237/examples/tgv/tgv.f90#L27).

```.f90
  ! Rescale mesh
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert
    real(kind=rp) :: d
    d = 4._rp

    ! original mesh has size 0..8 to be mapped onto -pi..pi
    ! will be updated later to a method giving back the vertices of the mesh
    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = (msh%points(i)%x(1) - d) / d * pi
       msh%points(i)%x(2) = (msh%points(i)%x(2) - d) / d * pi
       msh%points(i)%x(3) = (msh%points(i)%x(3) - d) / d * pi
    end do

  end subroutine user_mesh_scale

```

The registering of the above function in `user_setup` should then be done as follows:

```.f90
u%user_mesh_setup => user_mesh_scale
```

### Scalar boundary conditions {#user-file_scalar-bc}

This user functions allows for the speci

## Case-specific user functions 

As explained in the [case file](case-file.md) page, certain components of the simulation can be set to be user defined. These components are:

- The initial condition defined in `case.fluid.initial_condition`,
- The inflow boundary condition defined in `case.fluid.inflow_condition`,
- The fluid forcing/source term defined in `case.fluid.source_terms`,
- The scalar source term defined in `case.scalar.source_terms`,

### Initial condition



### Inflow condition

### Fluid source term

### Scalar source term

### Scalar boundary condition

## Basic working example {#basic-working-example}

## Compiling and running

The user file is a regular Fortran `.f90` file that needs to be compiled with `makeneko`, located in the `bin` folder of your neko installation. To compile a user file `user.f90`, run:

```bash
makeneko user.f90
```

If everything goes well, you should observe the following output:

```bash
N E K O build tool, Version 0.7.99
(build: 2024-02-13 on x86_64-pc-linux-gnu using gnu)

Building user NEKO ... done!
```

Compiling your user file with `makeneko` will create a `neko` executable, which you will need to execute with your case file as an argument. For example, if your case file is called `user.case`:

```bash
./neko user.case
```

Or in parallel using MPI:

```bash
mpirun -n 8 ./neko user.case
```
