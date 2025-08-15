# User File {#user-file}

\tableofcontents

The user file is a fortran file where the user can implement their own functions
to extend the capabilities of the default Neko executable. The user file can be
used for setting advanced initial/boundary conditions, source terms, I/O
operations, and interactions with the Neko framework.

This section will provide a written explanation of the user file, its
structure, and how to compile and run it. For a more hands-on approach, see the
[examples](@ref programming-examples) section, which contains several user file
examples that demonstrate various aspects of programming with Neko.

## Compiling and running

The user file is a regular Fortran `.f90` file that needs to be compiled with
`makeneko`, located in the `bin` folder of your neko installation. To compile a
user file `user.f90`, run:

```bash
makeneko user.f90
```

If everything goes well, you should observe the following output:

```bash
N E K O build tool, Version 1.99.1
(build: 2025-08-01 on x86_64-pc-linux-gnu using gnu)
Building user NEKO ...
Detected the module named 'user' in user.f90
No custom modules detected.
No custom modules register types.
Done!
```

Compiling your user file with `makeneko` will create a `neko` executable, which
you will need to execute with your case file as an argument. For example, if
your case file is called `user.case`:

```bash
./neko user.case
```

Or in parallel using MPI:

```bash
mpirun -n 8 ./neko user.case
```


## High-level structure

The current high-level structure of the user file is shown below.

```fortran
module user
  use neko
  implicit none

contains

  ! Register user defined functions here (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user


  end subroutine user_setup

end module user

```

The user file implements the `user` module. This modules contains a subroutine
named `user_setup`, which we use to interface the internal procedures defined in
`src/common/user_intf.f90` with the subroutines that you will implement in your
user file. Each user subroutine should be implemented under the `contains`
statement, below `user_setup`.

@note The above code snippet is the most basic code structure for the user file.
Compiling it and running it would be equivalent to running the "vanilla" neko
executable `bin/neko` in your local neko installation folder.

## Default user functions

The following user functions, if defined in the user file, will always be
executed, regardless of what is set in the case file:

- [startup](@ref user-file_init-and-final): For early access to the case
  file, its manipulation or initializing simple user parameter variables.
- [initialize](@ref user-file_init-and-final): For initializing most user
  variables and objects.
- [finilaze](@ref user-file_init-and-final): For finalizing, e.g
  freeing variables and terminating processes
- [compute](@ref user-file_user-check): Executed at the end of every time
  step, for e.g. computing and/or outputting user defined quantities.
- [preprocess](@ref user-file_user-check): Similar to `compute` but executed at
  the beginning of the time step.
- [material_properties](@ref user-file_mat-prop): For computing and setting
  material properties such as `rho`, `mu`, `cp` and `lambda`.
- [mesh_setup](@ref user-file_user-mesh-setup): For applying a deformation
  to the mesh element nodes, before the simulation time loop.

### Initializing and finalizing {#user-file_init-and-final}

Three subroutines `startup`, `initialize` and `finalize` may be used to
initialize/finalize any user defined variables, external objects, or processes.
The `startup` routine is called immediately after the case file is read in,
meaning that no important objects are initialized yet (e.g. the mesh, fluid,
etc.). The `initialize` and `finalize` are respectively executed right
before/after the simulation time loop.

In most cases, one can use `initialize` routine and not the `startup`. The
latter is only necessary when you want to either manipulate the case file
programmatically before it is used in the simulation, or define some variables
that will be used in the constructors of some of the types. An example of the
latter is defining some material property constants that can be used in the user
`material_properties` routine, which is run by the constructor of the
`case%fluid` object.


```fortran
  ! Manipulate the case file and define simple user variables
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    ! insert your initialization code here

  end subroutine startup

  ! Initialize user variables or external objects
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    ! insert your initialization code here

  end subroutine initialize

  ! Finalize user variables or external objects
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    ! insert your code here

  end subroutine finalize
```

In the example above, the subroutines `startup`, `initialize` and `finalize`
contain the actual implementations. They must also be associated to the internal
procedures inside, inside `user_setup`:

```fortran

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup
    user%initialize => initialize
    user%finalize => finalize

  end subroutine user_setup

```

@note All three routines are independent of each other. Using one does not
require the use of the other. Moreover, the names of the routines implemented in
the user file can be chosen freely, as long as they implement the correct
interface.

### Computing at every time step {#user-file_user-check}

The subroutine `compute` is executed at the end of every time step. It can be
used for computing and/or outputting your own variables/quantities at every time
step.
```fortran
  ! This is called at the end of every time step
  subroutine user_check(time)
    type(time_state_t), intent(in) :: time

    ! insert code below

  end subroutine user_check

```

In the example above, the subroutine `user_check` contains the actual
implementation, as a homage to the Nek5000 routine, which served a similar
purpose. As usual, the routine needs to be registered by adding:

```fortran
user%compute => user_check
```

to our `user_setup`. The `preprocess` routine can be implemented similarly.

### Setting material properties {#user-file_mat-prop}

The `material_properties` routine allows for more complex computations and
setting of various material properties, such as `rho`, `mu` for the fluid and
`cp`, `lambda` for the scalar. The example below is taken from the
[rayleigh_benard_cylinder
example](https://github.com/ExtremeFLOW/neko/blob/564686b127ff75a362a06126c6b23e9b4e21879e/examples/rayleigh_benard_cylinder/rayleigh.f90#L22C1-L38C41).

```fortran
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    if (scheme_name .eq. "fluid") then
       call field_cfill(properties%get("fluid_rho"), 1.0_rp)
       call field_cfill(properties%get("fluid_mu"), mu)
    else if (scheme_name .eq. "temperature") then
       call field_cfill(properties%get("temperature_cp"), 1.0_rp)
       call field_cfill(properties%get("temperature_lambda"), mu / Pr)
    end if
  end subroutine material_properties
```

And of course not forgetting to register our function in `user_setup` by adding
the following line:

```fortran
u%material_properties => material_properties
```

Note the usage of the `scheme_name` argument. This is used in several user
routines, to provide the information about what solver (scheme) called the
routine. In this example, we make use of it to distinguish the fluid and scalar
solvers. The default scheme names are `fluid` and `scalar`, but can be
controlled by the `name` entry in the JSON configuration of the solver. Here,
the name for the scalar was set to `temperature`.

### Runtime mesh deformation {#user-file_user-mesh-setup}

This user function allows for the modification of the mesh at runtime, by acting
on the element nodes of the mesh specified in the case file. This function is
only called once before the simulation time loop. The example below is taken
from the [tgv
example](https://github.com/ExtremeFLOW/neko/blob/a0613606360240e5059e65d6d98f4a57cf73e237/examples/tgv/tgv.f90#L27-L42).

```fortran
  ! Rescale mesh
  subroutine user_mesh_scale(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time

    integer :: i, p, nvert
    real(kind=rp) :: d
    d = 4._rp

    ! The original mesh has size 0..8 to be mapped onto -pi..pi
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

```fortran
    u%mesh_setup => user_mesh_scale
```

## Case-specific user functions

As explained in the [case file](case-file.md) page, certain components of the
simulation can be set to be user defined. These components and their associated
user functions are:

| Description                   | User function                                                   | JSON Object in the case file                                      |
| ----------------------------- | --------------------------------------------------------------- | ------------------------------------------------------------------|
| Initial conditions            | [initial_conditions](@ref user-file_user-ic)                    | `case.fluid.initial_condition` or `case.scalar.initial_condition` |
| Source terms                  | [source_term](@ref user-file_user-f)                            | `case.fluid.source_terms` or  `case.scalar.source_terms`          |
| Dirichlet boundary conditions | [dirichlet_conditions](@ref user-file_field-dirichlet-update)   | `case.fluid.boundary_types` or `case.scalar.boundary_types`       |

@note For the sake of simplicity, we refer to the setup with one scalar, i.e.
`case.scalar` in the JSON. For multiple scalars, the same things apply, but the
configuration is inside each individual element of `case.scalars`. In the case of
multiple scalars, users should use `"checkpoint_format": "hdf5"` when enabling
checkpointing.

### Fluid and Scalar initial conditions {#user-file_user-ic}

Enabling user defined initial conditions for the fluid and/or scalar is done by
setting the `initial_condition.type` to `"user"` in the relevant sections of the
case file, `case.fluid` and/or `case.scalar`.

```.json

"case": {
    "fluid": {
        "initial_condition": {
            "type": "user"
        }
    }
}
```

See the relevant sections on the [fluid](@ref case-file_fluid-ic) and
[scalar](@ref case-file_scalar) initial conditions in the
[case file page](@ref case-file) for more details.

The associated user function for initial conditions must then be added to the
user file. An example for the fluid taken from the [advecting cone
example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/scalar_mms/scalar_mms.f90#L55-L79),
is shown below.

```fortran
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    integer :: i, e, k, j
    real(kind=rp) :: x, y
    type (field_t), pointer :: u, v, w, s
    type(dofmap_t), pointer :: dof

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       call field_rone(u)
       call field_rzero(v)
       call field_rzero(w)
    else !scalar
       s => fields%get(scheme_name)
       if (scheme_name .eq. 's1') then
          do i = 1, s%dof%size()
             x = s%dof%x(i,1,1,1)
             y = s%dof%y(i,1,1,1)
             s%x(i,1,1,1) = sin(x)
          end do
       end if
    end if
  end subroutine initial_conditions

```

Not again the usage of `scheme_name` to distinguish between the fluid and the
scalar. Depending on the scheme, the contents of the field list `fields`
changes, and we extract individual fields via field pointers accordingly.
The incompressible fluid solver always generates solution fields, `u`, `v` and
`w`. For the scalar, the name of the field coincides with `scheme_name`. For single
scalar cases, this defaults to `s`. For multiple scalar cases, the field name is
set to the scalar name specified in the JSON configuration (e.g., "s1", "s2", etc.).

@note Notice that the code for the scalar runs on the CPU. There is no need to
add the transfer to GPU memory in this user routine, it will be done under the
hood afterwards.

We should also add of the following lines in `user_setup`, as usual.

```fortran
user%initial_condtions => initial_conditions.
```

### Fluid and scalar source terms {#user-file_user-f}

Enabling user defined source terms for the fluid and/or scalar is done by adding
JSON Objects to the `case.fluid.source_terms` and/or `case.scalar.source_terms`
lists.

```.json

"case": {
    "fluid": {
        "source_terms":
        [
            {
                "type": "user"
            }
        ]
    }
}
```

See the relevant sections on the [fluid](@ref case-file_fluid-source-term) and
[scalar](@ref case-file_scalar) source terms in the [case file page](@ref
case-file) for more details.

The associated user functions for the fluid and/or scalar source terms can then
be added to the user file. An example for the fluid, taken from the `scalar_mms`
example, is shown below.

```fortran
  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    real(kind=rp) :: x
    type(field_t), pointer :: f
    integer :: i

    if (scheme_name .eq. 'fluid') return

    f => rhs%items(1)%ptr
    do i = 1, f%size()
       x = f%dof%x(i,1,1,1)

       ! 0.01 is the viscosity
       f%x(i,1,1,1) = cos(x) - 0.01 * sin(x) - 1.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f%x, f%x_d, f%size(), &
            HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine source_term
```

@note Notice the use of the `NEKO_BCKND_DEVICE` flag, which will be set to 1 if
running on GPUs, and the use of `device_` functions. See [Running on GPUs](@ref
user-file_tips_running-on-gpus) for more information on how this works.

As usual, we point the user routine to our implementation in `user_setup`.

```fortran
user%source_term => source_term
```

### Dirichlet boundary conditions {#user-file_field-dirichlet-update}

This user function can be used to specify Dirichlet boundary values for velocity
components `u,v,w`, the pressure `p`, and/or the scalar(s). This type of
boundary condition allows for time-dependent velocity profiles or non-uniform
 pressure profiles to e.g. impose an outlet pressure computed from another
simulation.

The user routine is called by the `user_velocity` and `user_pressure` boundary
conditions for the fluid, and the `user` boundary condtition for the scalar.
Once the appropriate boundaries have been identified, the user function
`dirichlet_conditions` should be used to compute and apply the desired values to
our velocity/pressure/scalar field(s).

The structure of the interface is very similar to e.g. the initial conditions.
One gets a list of solution fields, the contents of which depends on what scheme
owns the boundary condition. For `user_velocity`, it will be a list of 3 fields
with names `u`, `v`, `w`; For `user_pressure`, a list with a single field `p`;
For the scalar, also a single field with the same name as the solution field for
the scalar (`s` by default).

It is crucial to understand that all three boundary condition will call the same
routine! So, if one has, for example, both `user_velocity` for the fluid and
`user` for the scalar, it is necessary to have an `if` statement in the user
routine to distinguish between the two cases. The convenient way to do that is
to check the size of the passed field list and the names of the fields inside.
For example, if there is one field and it is called `s`, one executes the code
setting the boundary values for the scalar `s`.

Note that the fields that one manipulates in the user routine are not the actual
velocity fields, etc. Instead, the code does a masked copy from the dummy fields
manipulated in the user routine to the actual fields of unknowns, with the mask
corresponding to the boundary faces. So, even if you somehow manipulate the
fields elsewhere in the domain inside the user routine, that will not affect the
solution.

In the following example, we indicate in `case.fluid.boundary_conditions` that
we would like to apply a velocity profile on the boundary number 1 (in this
case, the inlet boundary). On boundary number 2 (the outlet boundary), we change
the pressure. In `case.scalar.boundary_conditions`, we indicate the same for the
scalar on boundaries 1 and 2 (inlet and outlet).

The header of the user function is given in the code snippet below.

```fortran
  subroutine dirichlet_update(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(inout) :: bc
    type(time_state_t) :: time
```

The arguments and their purpose are as follows:

* `fields` is the list of the fields that can be edited. It is a list of
`field_t` objects.
  * The field `i` contained in `fields` is accessed using
  `fields%%items(i)%%ptr` and will refer to a `field_t` object. Alternatively,
  one can use the `get` method to retrieve a field by name or index, as done in
  the examples above for other routines.
* `bc` contains a `field_dirichlet_t` object to help access the boundary indices
  through the boundary mask, `msk`.
  * The boundary mask of the `bc `object is accessed via `bc%%msk`. It contains
  the linear indices of each GLL point on the boundary facets. @note
  `msk(0)` contains the size of the array. The first boundary index is `msk(1)`.
* `time`, is a simple structure that contains various info on time stepping,
  notably, the current time iteration and time value.

Links to the documentation to learn more about what the types mentioned above
contain and how to use them: `src/field/field.f90`, `src/bc/bc.f90`.

The user function should be registered in `user_setup` with the following line:

```fortran
u%dirichlet_conditions => dirichlet_update
```

A very simple example illustrating the above is shown below, which is taken from the
[cyl_boundary_layer example](https://github.com/ExtremeFLOW/neko/blob/feature/field_bcs/examples/cyl_boundary_layer/cyl_bl.f90)

```fortran
  subroutine dirichlet_update(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

    integer :: i
    real(kind=rp) :: y,z

    ! Only do this at the first time step since our BCs are constant in time.
    if (tstep .ne. 1) return

    ! Check that we are being called by the fluid via the name of the field
    if (fields%items(1)%ptr%name .eq. "u") then

       associate(u => field_bc_list%items(1)%ptr, &
                 v => field_bc_list%items(2)%ptr, &
                 w => field_bc_list%items(3)%ptr)
         !
         ! Perform operations on u%x, v%x, w%x here
         ! Here we are applying very simple uniform boundaries (u,v,w) = (1,0,0)
         ! Technically the values are put in the interior as well, but this
         ! does not matter, only the boundary values will be copied to the
         ! actual fields
         u = 1.0_rp
         v = 0.0_rp
         w = 0.0_rp

       end associate

    ! Check that we are being called by the user_pressure bc via the name
    ! of the field
    else if (fields%items(1)%ptr%name .eq. "p") then
       associate( p => fields%items(1)%ptr)
         !
         ! Perform operations on the pressure field here
         !

         do i = 1, bc%msk(0)
            p%x(bc%msk(i), 1, 1, 1) = -1
         end do

       end associate

    ! Check that we are being called by the scalar via the name of the field
    else if (fields%items(1)%ptr%name .eq. "s") then

       associate( s => field_bc_list%items(1)%ptr)
         !
         ! Perform operations on the scalar field here
         !

         do i = 1, bc%msk(0)
            y = bc%dof%y(bc%msk(i), 1, 1, 1)
            z = bc%dof%z(bc%msk(i), 1, 1, 1)
            s%x(bc%msk(i), 1, 1, 1) = sin(y)*sin(z)
         end do

       end associate

    end if
  end subroutine dirichlet_update
```

This example is applying constant dirichlet values at the selected boundaries
for the velocity components and presure. The scalar is applied a function
`s(y,z) = sin(y)*sin(z)` to demonstrate the usage of boundary masks.

@attention The notation `u = 1.0_rp` is only possible because of the overloading
of the assignement operator `=` in `field_t`. In general, a field's array should
be accessed and modified with `u%%x`.

Note that we are only applying our boundary values at the first timestep, which
is done simply with the line `if (time%tstep .ne. 1) return`. This is a trick
that can be used for time-independent boundary profiles that require some kind
of time consuming operation like interpolation or reading from a file, which
would add overhead if executed at every time step.

@attention All the rules for [Running on GPUs](@ref
user-file_tips_running-on-gpus) apply when working on field arrays. Use
`device_memcpy` to make sure the device arrays are also updated.

## Additional remarks and tips

### Running on GPUs {#user-file_tips_running-on-gpus}

When running on GPUs, special care must be taken when using certain user
functions. The short explanation is that the device (GPU) has its own memory and
cannot directly access the memory on the host (CPU). This means that data and
more specifically arrays must be copied manually from the host to the device
(see device::device_memcpy).

@attention In some cases, data transfer via `device_memcpy` is avoidable. Neko
has some device math functions implemented that operate directly on device
arrays. If you can decompose whatever operations you are performing in a user
function into a set of instructions from the `math` module (e.g. `cadd`,
`cfill`, `sub2`, ...), you may use the corresponding `device_math` functions to
[offload work to the GPU](@ref accelerators_offload-work).
Additionally, when you operate on `field_t` objects (which is very common), you
can use routines in the `field_math` module, that will automatically run on the
right backend. This is highly recommended and avoids `if` statements.

An example where manual transfer to the GPU was necessary has already been
given for the source term. It is reiterated below.

```fortran
  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    real(kind=rp) :: x
    type(field_t), pointer :: f
    integer :: i

    if (scheme_name .eq. 'fluid') return

    f => rhs%items(1)%ptr
    do i = 1, f%size()
       x = f%dof%x(i,1,1,1)

       ! 0.01 is the viscosity
       f%x(i,1,1,1) = cos(x) - 0.01 * sin(x) - 1.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f%x, f%x_d, f%size(), &
            HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine source_term
```

Note that `f%%x` is the field's data on the CPU. Therefore, to populate the
correspoding array on GPU, we need to call to `device_memcpy`.
Looking at the detailof the `device_memcpy` call, we note the following:
- Device arrays are refered to by appending the suffix `_d` to the host array
  variable name (e.g. `f%%x` and `f%%x_d`). This is the standard in Neko.
- We specify the direction of the data movement with the flag `HOST_TO_DEVICE`.
  Other flags can also be used to move data from device to host
  (`DEVICE_TO_HOST`) or device to device (`DEVICE_TO_DEVICE`). See the
  [accelerators page](@ref accelerators_data-transfer) for more details on this.
- The `sync` argument is a non-optional argument which dictates wether or not to
  perform the data transfer synchronously.

@attention Use asynchronous data transfers at your own risk! If you are unsure,
use `sync = .true.` as a starting point.

Finally, observe that we use the flag `NEKO_BCKND_DEVICE` to check if we are
indeed running on GPUs. In that case, `NEKO_BCKND_DEVICE` would be equal to 1.

#### Custom GPU kernels {#user-file_tips_running-on-gpus-custom-kernels}

When running on GPUs it is possible to call an own custom kernel. This
could be more performant for more complex source terms or boundary
conditions.

Assume we have a custom kernel called `mydevice_kernel` in a CUDA file
called `mykernel.cu`, to set the initial velicity field. To call the
custom kernel, the user file must define a C interface to the routine
(inside the `.cu` file) that launches the kernel

```fortran
interface
  subroutine mydevice_kernel(u_d, v_d, w_d, n) &
        bind(c, name = 'mydevice_kernel')
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr
    type(c_ptr), value :: u_d, v_d, w_d
    integer(c_int) :: n
  end subroutine mydevice_kernel
end interface
```

Furthermore, the CUDA/HIP file must allow for C linkage, hence the
routine `mydevice_kernel` must be inside an `extern "C"` block.

```cpp
extern "C" {
  void mydevice_kernel(void *u, void *v, void *w, int *n) {
    /* Launch the device kernel here */
  }
}
```
The user defined routine for initial conditions calls the kernel in
the following way.

```fortran

  !> Set the advecting velocity field.
  subroutine set_velocity(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call mydevice_kernel(u%x_d, v%x_d, w%x_d, u%dof%size())
    endif

  end subroutine set_velocity

```

Finally, compile using `makeneko` and provide both the `user.f90` and
`mykernel.cu` file.

```bash
makeneko user.f90 mykernel.cu
```

@note `makeneko` does currently only support custom kernels written in CUDA or HIP.

### Registries {#user-file_tips_registries}

Neko uses the concept of `registry` as a practical way to retrieve fields and
point zones anywhere in the user file.

The field registry `neko_field_registry` is often used in user functions where
certain fields are not directly accessible as arguments. One can retrieve any
field in the registry by its `name` with `neko_field_registry%%get_field(name)`.
Default fields that are added to the registry are `u,v,w,p` and `s` if running
with the scalar enabled. For a practical example of usage, see the
[rayleigh benard example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/rayleigh_benard/rayleigh.f90#L102-L105)

Other fields may be added to the registry by various simulation components. For
example:
- If running with `simulation_components.vorticity` enabled, the fields
  `omega_x, omega_y, omega_z` will be accessible in the registry.
- If running with `simulation_components.lambda2` enabled, the field `lambda2`
  will be accessible in the registry.

@note You can add your own fields to the registry with `neko_field_registry%%add_field` (see field_registry::field_add).

The point zone registry, `neko_point_zone_registry`, can be used to retrieve
pointers to `point_zone_t` objects defined in the case file. See
[using point zones](#point-zones_using-point-zones) for detailed instructions.

### User access to solver internals {#user-file_access}

Neko proives a special variable called `neko_user_access` that can be used in
the user file to access the internal state of the entire case and the underlying
solvers.

In particular the case object (`case_t`) is reached via
`neko_user_access%%case`. There-in the `fluid` and `scalars` components refer
to the asscociated schemes. You are encouraged to look a bit at `case.f90` to
see the overall structure.

A very common use case is to get access to various SEM-related objects, as well
as the mesh object.

* `neko_user_access%%case%%fluid%%msh` -- the mesh.
* `neko_user_access%%case%%fluid%%dm_Xh` -- the function space used in the SEM.
* `neko_user_access%%case%%fluid%%dm_Xh` -- the map of degrees of freedom,
  contains the locations of the GLL nodes.
* `neko_user_access%%case%%fluid%%c_Xh` -- the coefficients of the SEM,
* `neko_user_access%%case%%fluid%%gs_Xh` -- the gather-scatter kernels used for
  direct stiffness summation.
