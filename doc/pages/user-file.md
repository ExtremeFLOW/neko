# User File {#user-file}

The user file is a fortran file where the user can implement their own functions
to extend the capabilities of the default Neko executable. The user file can be
used for setting advanced initial/boundary conditions, source terms, I/O
operations, and interactions with the Neko framework.

## Compiling and running

The user file is a regular Fortran `.f90` file that needs to be compiled with
`makeneko`, located in the `bin` folder of your neko installation. To compile a
user file `user.f90`, run:

```bash
makeneko user.f90
```

If everything goes well, you should observe the following output:

```bash
N E K O build tool, Version 0.7.99
(build: 2024-02-13 on x86_64-pc-linux-gnu using gnu)

Building user NEKO ... done!
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
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u


  end subroutine user_setup

end module user

```

The user file implements the `user` module. The `user` modules contains a
subroutine named `user_setup`, which we use to interface the internal procedures
defined in `src/common/user_intf.f90` with the subroutines that you will
implement in your user file. Each user subroutine should be implemented under
the `contains` statement, below `user_setup`.

@note The above code snippet is the most basic code structure for the user file.
Compiling it and running it would be equivalent to running the "vanilla" neko
executable `bin/neko` in your local neko installation folder.

## Default user functions

The following user functions, if defined in the user file, will always be
executed, regardless of what is set in the case file:

- [user_init_modules](#user-file_init-and-final): For initializing user
  variables and objects
- [user_finalize_modules](#user-file_init-and-final): For finalizing, e.g
  freeing variables and terminating processes
- [user_check](#user-file_user-check): Executed at the end of every time step,
  for e.g. computing and/or outputting user defined quantities.
- [material_properties](#user-file_mat-prop): For computing and setting material
  properties such as `rho`, `mu`, `cp` and `lambda`.
- [user_mesh_setup](#user-file_user-mesh-setup): For applying a deformation to
  the mesh element nodes, before the simulation time loop.
- [scalar_user_bc](#user-file_scalar-bc): For applying boundary conditions to
  the scalar, on all zones that are not already specified with uniform dirichlet
  values e.g. `d=1`. For more information on the scalar, see the [relevant
  section of the case file](#case-file_scalar).

### Initializing and finalizing {#user-file_init-and-final}

The two subroutines `user_init_modules` and `user_finalize_modules` may be used
to initialize/finalize any user defined variables, external objects, or
processes. They are respectively executed right before/after the simulation time
loop.

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

In the example above, the subroutines `initialize` and `finalize` contain the
actual implementations. They must also be interfaced to the internal procedures
`user_init_modules` and `user_finalize_modules` in `user_setup`:

```.f90

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u

    u%user_init_modules => initialize
    u%user_finalize_modules => finalize

  end subroutine user_setup

```

@note `user_init_modules` and `user_finalize_modules` are independent of each
other. Using one does not require the use of the other.

### Computing at every time step {#user-file_user-check}

The subroutine `user_check` is executed at the end of every time step. It can be
used for computing and/or outputting your own variables/quantities at every time
step.
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

In the example above, the subroutine `usercheck` contains the actual
implementation, and needs to be registered by adding:

```.f90
u%user_check => usercheck
```

to our `user_setup`.

### Setting material properties {#user-file_mat-prop}

`material_properties` allows for more complex computations and setting of
various material properties, such as `rho`, `mu` for the fluid and `cp`,
`lambda` for the scalar. The example below is taken from the
[rayleigh-benard-cylinder
example](https://github.com/ExtremeFLOW/neko/blob/564686b127ff75a362a06126c6b23e9b4e21879e/examples/rayleigh-benard-cylinder/rayleigh.f90#L22C1-L38C41).

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

And of course not forgetting to register our function in `user_setup` by adding
the following line:

```.f90
u%material_properties => set_material_properties
```

### Runtime mesh deformation {#user-file_user-mesh-setup}

This user function allows for the modification of the mesh at runtime, by acting
on the element nodes of the mesh specified in the case file. This function is
only called once before the simulation time loop. The example below is taken
from the [tgv
example](https://github.com/ExtremeFLOW/neko/blob/a0613606360240e5059e65d6d98f4a57cf73e237/examples/tgv/tgv.f90#L27-L42).

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

This user function can be used to specify the scalar boundary values, on all
zones that are not already set to uniform Dirichlet or Neumann values e.g. `d=1`
or `n=0`. For more information on the scalar, see the [relevant section of the
case file](#case-file_scalar). The example below sets the scalar boundary
condition values to be a linear function of the `z` coordinate (taken from the
[rayleigh-benard
example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/rayleigh-benard-cylinder/rayleigh.f90#L41-L63)).

```.f90

  subroutine set_scalar_boundary_conditions(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! This will be used on all zones without labels
    s = 1.0_rp - z

  end subroutine set_scalar_boundary_conditions

```

This function will be called on all the points on the relevant boundaries. The
registering of the above function in `user_setup` should be done as follows:

```.f90
u%scalar_user_bc => set_scalar_boundary_conditions
```


## Case-specific user functions

As explained in the [case file](case-file.md) page, certain components of the
simulation can be set to be user defined. These components and their associated
user functions are:

| Description                     | User function                                              | JSON Object in the case file     |
|---------------------------------|------------------------------------------------------------|----------------------------------|
| Fluid initial condition         | [fluid_user_ic](#user-file_user-ic)                        | `case.fluid.initial_condition`   |
| Scalar initial condition        | [scalar_user_ic](#user-file_user-ic)                       | `case.scalar.initial_condition`  |
| Fluid inflow boundary condition | [fluid_user_if](#user-file_fluid-user-if)                  | `case.fluid.inflow_condition`    |
| Scalar boundary conditions      | [scalar_user_bc](#user-file_scalar-bc)                     | (user function is always called) |
| Fluid source term               | [fluid_user_f_vector or fluid_user_f](#user-file_user-f)   | `case.fluid.source_terms`        |
| Scalar source term              | [scalar_user_f_vector or scalar_user_f](#user-file_user-f) | `case.scalar.source_terms`       |

Note that `scalar_user_bc` is included for completeness but is technically not case-specific.

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

See the relevant sections on the [fluid](#case-file_fluid-ic) and
[scalar](#case-file_scalar) initial conditions in the [case file
page](#case-file) for more details.

The associated user functions for the fluid and/or scalar initial conditions can
then be added to the user file. An example for the fluid taken from the
[advecting cone
example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/advecting_cone/advecting_cone.f90#L48-L75),
is shown below.

```.f90

  !> Set the advecting velocity field.
  subroutine set_velocity(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y

    do i = 1, u%dof%size()
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)

       ! Angular velocity is pi, giving a full rotation in 2 sec
       u%x(i,1,1,1) = -y*pi
       v%x(i,1,1,1) = x*pi
       w%x(i,1,1,1) = 0
    end do

     if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_velocity

```

@note Notice the use of the `NEKO_BCKND_DEVICE` flag, which will be set to 1 if
running on GPUs, and the calls to `device_memcpy` to transfer data between the
host and the device. See [Running on GPUs](#user-file_tips_running-on-gpus) for
more information on how this works.

The same can be done for the scalar, with the example below also inspired from
the [advecting cone
example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/advecting_cone/advecting_cone.f90#L14-L45):

```.f90

  !> User initial condition for the scalar
  subroutine set_s_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: cone_radius, mux, muy, x, y, r, theta

    ! Center of the cone
    mux = 1
    muy = 0

    cone_radius = 0.5

    do i = 1, s%dof%size()
       x = s%dof%x(i,1,1,1) - mux
       y = s%dof%y(i,1,1,1) - muy

       r = sqrt(x**2 + y**2)
       theta = atan2(y, x)

       ! Check if the point is inside the cone's base
       if (r > cone_radius) then
         s%x(i,1,1,1) = 0.0
       else
         s%x(i,1,1,1) = 1.0 - r / cone_radius
       endif
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_s_ic

```

@note Notice the use of the `NEKO_BCKND_DEVICE` flag, which will be set to 1 if
running on GPUs, and the calls to `device_memcpy` to transfer data between the
host and the device. See [Running on GPUs](#user-file_tips_running-on-gpus) for
more information on how this works.

We should also add of the following lines in `user_setup`, registering our user
functions `set_velocity` and `set_s_ic` to be used as the fluid and scalar
initial conditions:

```.f90
u%fluid_user_ic => set_velocity
u%scalar_user_ic => set_s_ic
```

### Fluid inflow condition {#user-file_fluid-user-if}

Enabling user defined inflow condition for the fluid is done by setting
the `case.fluid.inflow_condition.type` to `"user"`:

```.json

"case": {
    "fluid": {
        "inflow_condition": {
            "type": "user"
        }
    }
}
```

See the [the relevant section](#case-file_fluid-if) in the [case file
page](#case-file) for more details. The associated user function for the fluid
inflow condition can then be added to the user file. An example inspired from
the [lid-driven cavity
example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/lid/lid.f90#L29-L53)
is shown below.

```.f90
 ! user-defined boundary condition
  subroutine user_bc(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    real(kind=rp) lsmoothing
    lsmoothing = 0.05_rp    ! length scale of smoothing at the edges

    u = step( x/lsmoothing ) * step( (1._rp-x)/lsmoothing )
    v = 0._rp
    w = 0._rp

  end subroutine user_bc
```

We should also add of the following line in `user_setup`, registering our user
function `user_bc` to be used as the fluid inflow conditions:

```.f90
u%fluid_user_if => user_bc
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
                "type": "user_vector"
            }
        ]
    }
}
```

See the relevant sections on the [fluid](#case-file_fluid-source-term) and
[scalar](#case-file_scalar) source terms in the [case file page](#case-file) for
more details.

@attention There are two variants of the source term user functions: `_user_f`
and `_user_f_vector`. The former is called when setting `"user_pointwise"` as
the source term type, while the latter requires the use of the `"user_vector"`
keyword in the case file. The pointwise variant, `fluid_user_f` is not supported
on GPUs. In general, `fluid_user_f_vector` is the prefered variant, and is the
one which will be use in our examples below. The same applies for the scalar
source term user functions.

The associated user functions for the fluid and/or scalar source terms can then
be added to the user file. An example for the fluid, taken from the
[rayleigh-benard-cylinder
example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/rayleigh-benard-cylinder/rayleigh.f90#L101C1-L121C44),
is shown below.

```.f90
  ! Sets the z-component of the fluid forcing term = scalar
  subroutine set_bousinesq_forcing_term(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

    ! Retrieve u,v,w,s fields from the field registry
    type(field_t), pointer :: u, v, w, s
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term
```

@note Notice the use of the `neko_field_registry` to retrieve the velocity and
scalar fields. See [Registries](#user-file_tips_registries) for more information
about registries in neko. @note Notice the use of the `NEKO_BCKND_DEVICE` flag,
which will be set to 1 if running on GPUs, and the use of `device_` functions.
See [Running on GPUs](#user-file_tips_running-on-gpus) for more information on
how this works.

The same can be done for the scalar, with the example below also taken from the
[scalar_mms
example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/scalar_mms/scalar_mms.f90#L28-L47):

```.f90

  !> Set source term
  subroutine set_source(f, t)
    class(scalar_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: x, y
    integer :: i

    do i = 1, f%dm%size()
       x = f%dm%x(i,1,1,1)
       y = f%dm%y(i,1,1,1)

       ! 0.01 is the viscosity
       f%s(i,1,1,1) = cos(x) - 0.01 * sin(x) - 1.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f%s, f%s_d, f%dm%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_source

```

@note Notice the use of the `NEKO_BCKND_DEVICE` flag, which will be set to 1 if
running on GPUs, and the call to `device_memcpy` to transfer data between the
host and the device. See [Running on GPUs](#user-file_tips_running-on-gpus) for
more information on how this works.

We should also add of the following lines in `user_setup`, registering our user
functions `set_boussinesq_forcing_term` and `set_source` to be used as the fluid
and scalar source terms:

```.f90
u%fluid_user_f_vector => set_boussinesq_forcing_term
u%scalar_user_f_vector => set_source
```

## Additional remarks and tips

### Running on GPUs {#user-file_tips_running-on-gpus}

When running on GPUs, special care must be taken when using certain user
functions. The short explanation is that the device (GPU) has its own memory and
cannot directly access the memory on the host (CPU). This means that data and
more specifically arrays must be copied manually from the host to the device
(see
[device_memcpy](https://neko.cfd/docs/d6/dac/interfacedevice_1_1device__memcpy.html)).

@attention In some cases, data transfer via `device_memcpy` is avoidable. Neko
has some device math functions implemented that operate directly on device
arrays. If you can decompose whatever operations you are performing in a user
function into a set of instructions from the `math` module (e.g. `cadd`,
`cfill`, `sub2`, ...), you may use the corresponding `device_math` functions to
[offload work to the GPU](#accelerators_offload-work). See the [fluid forcing
code snippet](#user-file_user-f) for a simple example. For more advanced
examples, see the [rayleigh-benard
example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/rayleigh-benard/rayleigh.f90#L96-119)
or the [tgv
example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/tgv/tgv.f90#L146-172).

To illustrate this, let us have a look at the [fluid initial condition code
snippet](#user-file_user-ic):

```.f90

  !> Set the advecting velocity field.
  subroutine set_velocity(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y

    !
    ! 1. Set the initial condition in fields u%x, v%x, w%x
    !
    do i = 1, u%dof%size()
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)

       ! Angular velocity is pi, giving a full rotation in 2 sec
       u%x(i,1,1,1) = -y*pi
       v%x(i,1,1,1) = x*pi
       w%x(i,1,1,1) = 0
    end do

    !
    ! 2. Copy the data set in u%x, v%x, w%x to the device arrays
    ! u%x_d, v%x_d, w%x_d.
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_velocity

```

The code above is used to set the fluid initial condition, by specifying the
values of fields `u,v,w` (and `p`) at all points in the domain. Notice that we
have divided the above code into two parts.

In the first part, we set the velocity components `u=-y*pi*`, `v=x*pi*`, and
`w=0`, which updates the velocity field arrays `u%%x, v%%x, w%%x` allocated on
the **host (CPU)**. If we were to run on GPUs, these lines of code would only
act on the velocity arrays on the host (CPU), leaving the device (GPU) arrays
untouched.

We take care of this in the second part, for all three velocity arrays. To
update the device (GPU) arrays, we use `device_memcpy` to copy the data
contained in a host (CPU) array to a device (GPU) array. Looking at the details
of the `device_memcpy` calls, we note the following:
- Device arrays are refered to by appending the suffix `_d` to the host array
  variable name (e.g. `u%%x` and `u%%x_d`). This is the standard in Neko.
- We specify the direction of the data movement with the flag `HOST_TO_DEVICE`.
  Other flags can also be used to move data from device to host
  (`DEVICE_TO_HOST`) or device to device (`DEVICE_TO_DEVICE`). See the
  [accelerators page](#accelerators_data-transfer) for more details on this.
- The `sync` argument is a non-optional argument which dictates wether or not to
  perform the data transfer synchronously.

@attention Use asynchronous data transfers at your own risk! If you are unsure,
use `sync = .true.` as a starting point.

Finally, observe that we use the flag `NEKO_BCKND_DEVICE` to check if we are
indeed running on GPUs. In that case, `NEKO_BCKND_DEVICE` would be equal to 1.

### Registries {#user-file_tips_registries}

Neko uses the concept of `registry` as a practical way to retrieve fields and
point zones anywhere in the user file.

The field registry `neko_field_registry` is often used in user functions where
certain fields are not directly accessible as arguments. One can retrieve any
field in the registry by its `name` with `neko_field_registry%%get_field(name)`.
Default fields that are added to the registry are `u,v,w,p` and `s` if running
with the scalar enabled. For a practical example of usage, see the [rayleigh
benard
example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/rayleigh-benard/rayleigh.f90#L102-L105)

Other fields may be added to the registry by various simulation components. For
example:
- If running with `simulation_components.vorticity` enabled, the fields
  `omega_x, omega_y, omega_z` will be accessible in the registry.
- If running with `simulation_components.lambda2` enabled, the field `lambda2`
  will be accessible in the registry.

@note You can add your own fields to the registry with [neko_field_registry%%add_field](https://neko.cfd/docs/d1/dc4/namespacefield__registry.html#a10e4570d0cae94f110c9b0eacbb56769).

The point zone registry, `neko_point_zone_registry`, can be used to retrieve
pointers to `point_zone_t` objects defined in the case file. See [using point
zones](#point-zones_using-point-zones) for detailed instructions.
