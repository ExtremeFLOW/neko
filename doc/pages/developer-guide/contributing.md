# Contributing to Neko {#contributing} 

\tableofcontents

Please read the following guide before contributing new code or a bug fix to Neko.

All contributions to Neko must be made under the 3-Clause BSD license. Please refer to the `COPYING` file.

## Git branches
Neko follows the Git branching model described in https://nvie.com/posts/a-successful-git-branching-model, where `develop` contains the latest contributions, and all pull requests should start from `develop` and be merged back into `develop`.  New branches should be named `feature/<description>` for new features or `fix/<description>` for bug fixes.

When a pull request is submitted, a series of continuous integration tests will be run. A pull request will not be accepted nor merged into `develop` until it passes the test suite.

## Code style
Fortran code should conform to the Fortran 2008 standard and should use an indentation level of 2, except for the extra indentation within `do` `if`, `select` or `where` statements and for each level inside a structure e.g. `type`, `interface`, where the indentation level is 3, and a 0 indentation is used for `module` or `contains` (except for `contains` inside a derived type, where a single indentation level is used). These are the default rules in Emacs' Fortran mode, an example is given below,
```fortran
module example
  use mod
  implicit none

  type :: derived_t
     integer :: x
   contains
     procedure, pass(this) :: bar     
  end type derived_t

contains

  subroutine foo(x, y)
    integer, intent(in) :: x
    integer, intent(inout) :: y
    real(kind=rp) :: value

    do i = 1, 10
       ...
    end do

    if (x .lt. y) then
       ...
    end if
    
  end subroutine foo
end module example
```
Please note that the maximum line length in Neko should not exceed 80 columns.

Additional details on the code style, as well as tools used to verify it, can be
found under the [code style](@ref code-style).

### Data types
For portability reasons, it is essential to only use data type kinds defined in `src/config/num_types.f90` and avoid legacy constructs like `real*8` or `integer(kind=8)`

Floating-point numbers should be declared using `real` with the kind `rp`, which is the configured working precision for Neko (defaults to double). If single, double or quad precision is explicitly needed, use the kinds `sp`, `dp` or `qp`, respectively. For 16, 32, or 64-bit integers, Neko has defined the kinds ` i2`, `i4` or `i8`, respectively; however, for standard integers, it is perfectly fine to use a plain `integer`.

## Build system and code organization
This section contains information on how to add new source files to the build
system. It also give an overview of the common folder structures you will
encounter in the code.

Neko uses Autotools for building all sources. You will need to have at least
`autoconf`, `libtool`, and `automake` installed for development work. It is also
highly recommended to have `makedepf90` installed to avoid error-prone manual
dependency tracking. Since Neko uses `submodules`, a recent version of
`makedepf90` from  https://salsa.debian.org/science-team/makedepf90 is needed.

The vast majority of the code is located inside the `src` folder, where it is
organized into subfolders, roughly based on the functionality, which is
implemented inside the them. A common folder structure is used to organize the
implementation of the same functionality on various compute backends. This is
discussed further below. Outside, `src`, the `contrib` folder is home for the
implementations of various executable utilities.

### Building CPU Fortran code
The following steps describe how to add a new Fortran file to Neko's build
system.
1. Place the file in an appropriate subdirectory under `src/`. Either create a
   new subdirectory or place the file in `common` if none of the existing
   directories is a good match. Avoid placing the file directly under `src`
2. Add the file to the `neko_fortran_SOURCES` list in `src/Makefile.am`,
   following the pattern of `<subdir under src>/newfile.f90`.
3. Ensure correct dependency tracking
    * Using `makedepf90`
      - Regenerate build system by running `./regen.sh` at the top level, and
        reconfigure using `configure`
      - Regenerate dependency file (`src/.depends`) by issuing `make depend`
        under `src`
      - Do a final regeneration of the build system (`./regen.sh`)
    * Manually 
      - Manually add the new file and all its dependencies to `src/.depends`,
        using the following pattern, 
      ```Make
      directory/newfile.o : directory/newfile.f90 directory/dependency.o anotherdirectory/anotherdependency.o
      ```
      - Regenerate build system by running `./regen.sh` at the top level.
4. Finally reconfigure using `configure` and rebuild Neko with your new
   contribution!

### Device code
Here we focus on three particular backends: CUDA, HIP, and OpenCL. In Neko these
three are collectively referred to as the "device", distinguishing them from
both ordinary CPU code, but also other special backends like the SX. Whenever a
type needs to implement dedicated compute kernels, Neko often follows a
standardized directory and file structure. 

As an example, we will consider the folder `src/bc`, please open it in a
terminal alongside reading this text. Inside, you can see a bunch of `.f90`
files, corresponding to different boundary conditions, for example,
`symmtery.f90` and `dirichlet.f90`. These files define the corresponding types,
like `dirichlet_t` and are included in the build system as per the instructions
above. So, these files are just pure Fortran.

Observe that there is a directory called `bcknd`. This is where the low-level
implementation of compute kernels. Sometimes this folder contains a `cpu`
folder, with modules implementing CPU kernels. In the case of boundary
conditions, the CPU kernels are directly implemented in the high-level `.f90`
files mentioned above. Here, we focus on the folder `device`, which contains the
implementation of kernels for the 3 "device" backends. Note the contents: a
bunch of `.F90` files (note the big F in the extension!) and a folder for each
backend (`cuda`, `hip`, `opencl`).

Let us start with the `.F90`s. The reason for the big F is that these files
contain preprocessor directives (like `#ifdef`). This is just a general
convention, not invented by Neko authors. The `.F90` files act as a device 
abstraction layer. What that means in practice is that they contain:
- Wrapper routines, calling the correct backend kernel, depending on the device
  Neko has been compiled for.
- Several `interface` blocks, for the backend kernels, which are essentially
  glue between Fortran and the actual backend code, which is called from a C
  routine. This is standard "C interoperability" stuff, which you can read about
  in books or online. 
  
Let's take a look at `device_symmetry.F90`. It contains one subroutine,
`device_symmetry_apply_vector`, that just dispatches the call to one of the
backends, via the corresponding interfaces defined in the beginning of the file.
If we now go back to `symmetry.f90`, we see that it imports this routine, and
calls it inside `symmetry_apply_vector_dev`. 

Since the `.F90`s are ordinary Fortran files, they are also added to the build
system in the already specified way.

Now we leave the familiar Fortran territory and look at individual backends.
Generally, CUDA and HIP are very similar (in fact, even the code of the kernels
is often pretty much a copy-paste), so we will look at CUDA in detail and then
quickly mentioned adjustments for HIP. OpenCL is a different beast, but the 
pattern is still quite similar.

#### CUDA

To follow along, navigate to `src/bc/bcknd/device/cuda`. We see three types of
files: `.cu`, `*_kernel.h` and potentially other, auxillary, `.h` files.

- The `.cu` files contain a C routine, which corresponds to the `interface`
  defined in the `.F90` file. We need this because we cannot just launch a CUDA
  kernel from Fortran, so instead we created this binding between a Fortran and
  a C routine, and the latter launches the CUDA kernel for us.
- The kernels themselves are not in the `.cu` files, but rather in the
  `*_kernel.h` files. The `.cu` files simply `#include` them. For example,
  `symmtery.cu` includes `symmetry_kernel.h`,  which has the
  `symmetry_apply_vector_kernel`, which is then launched in the `.cu` file.
- The other `.h` files just contain implementations of functions and kernels
used in several of the `*_kernel.h` files.

Let's zoom out again, because now we have the full picture for CUDA, which is
similar for other device backends.
- The `.f90` file.
  - Defines the type, heavy on OOP. 
  - With a `use` statement, imports kernel subroutines from modules in the
    `bcknd` subdirectory.
  - Sometimes contains the CPU implementation instead of using a dedicated
    kernel for it.
- The `.F90` files inside `bcknd/device`.
  - Contains wrapper routines that are called from the `.f90` and dispatch to
    the implementation on the correct device backend.
  - Contains `interface`s for the implementation on each device, binding Fortran
    that we can call from the wrappers in the `.F90` with C code that launches
    actual device kernels.
- Thn `.cu` files in `bcknd/device/cuda` that contain C routines adhering to the
  interface defined in the `.F90`. They `#include` the kernel code, and launch 
  the kernel.
- The `_kernel.h` files in `bcknd/device/cuda` that contain implementations of
  CUDA kernels.


Now, let's see how the `.cu` and `.h` files are properly added to the build
system. First, consider the `src/Makefile.am` file. The `.cu`s are added under
the
```bash
if ENABLE_CUDA
libneko_la_SOURCES += \
```
list. On the other hand, all the `.h` files are added to the 
```bash
EXTRA_DIST = \
```
list.

Unfortunately, we also have to manually provide the dependencies between the
files. That is done in the `src/.depends_device` file. A typical line will bind
a `.lo` file with the corresponding `.cu` and `.h` files. For example:
```bash
bc/bcknd/device/cuda/symmetry.lo : bc/bcknd/device/cuda/symmetry.cu bc/bcknd/device/cuda/symmetry_kernel.h
```
So, one should make similar entries for each new `.cu` file one adds.

#### HIP
Thankfully, HIP works in precisely the same way as CUDA, with the following
differenece: instead of `.cu` we have `.hip` files. These are to be added under
```bash
if ENABLE_HIP
libneko_la_SOURCES += \
```
for the build system to see them. The `.h` files go to `EXTRA_DIST` as with
CUDA. The dependencies are also specified in the exact same way.

#### OpenCL
To follow along, navigate to `src/bc/bcknd/device/opencl`. Quite similar to
CUDA, there are up to 3 types of files that can be present here.
- The `.c` files, which are just like `.hip` or `.cu`. They contain a C routine,
  which launches the OpenCL kernel.
- The `_kernel.cl` files, which are like the `_kernel.h` files for CUDA and HIP.
  You will notice that the `.c` files `#include` the `.cl` files, but with an
  additional `.h` in the end like

  ```c
  #include "symmetry_kernel.cl.h"
  ```
  Don't let that lead you astray with your file naming: the `.h` is added
automatically, and the files are just `.cl`!
- There may possibly be auxillary `.h` files for common functionlity across
  kernels.

Note that the `.c` files contain also the following `include` statement.
```c
#include <device/opencl/prgm_lib.h>
```
This is an important file, which we have to modify to register our kernel. Each
OpenCL kernel has an associated `_program`, for example, if you look at
`symmetry.c`, you will see the variable `symmetry_program` used in a few places.
This is not really important for the implementation per se, but the point is
that twe have to define these `*_program` variables somewhere, and in Neko this
is done in two files. One is the above mentioned `src/device/opencl/prgm_lib.h`,
where you can find
```c
/** Device Symmetry kernels */
extern void *symmetry_program;
```
One needs to add a similar line for the `_program` corresponding to the new
kernel. The other file is very similar, but written in Fortran:
`src/device/opencl/prh_lib.F90`. Here you can e.g. find

```fortran
  !> Device Symmetry kernels
  type(c_ptr), public, bind(c) :: symmetry_program = C_NULL_PTR
```
and just like with the `.h` file a corresponding entry for the new kernel should
be added.

This finally brings us to the build system. The `.c` files are to be added under
```bash
if ENABLE_OPENCL
libneko_la_SOURCES += \
```
in `src/Makefile.am`, and the `.cl` and `.h` files go under the `EXTRA_DIST`
list, just like the CUDA and HIP kernel files did. 

The entry to `src/.depends_device` also follows a familiar pattern, e.g.
```bash
bc/bcknd/device/opencl/symmetry.lo : bc/bcknd/device/opencl/symmetry.c bc/bcknd/device/opencl/symmetry_kernel.cl.h
```

#### Device-based type polymorphism

Neko is developing, and not all the conventions active now have been applied
from the beginning, and even in new code there are sometimes exceptions. For
example, it has already been mentioned tha the CPU kernels sometimes reside in
a seprate module, like the device ones, and sometimes just implemented inside 
the main module defining the type.

Some important types in Neko use polymorphism in addition to the `#ifdef`
statements we saw in the `.F90` files in order define the behavior on the
device. If we look at the `pnpn_res.f90` file, we see that an abstract type is
defined, `pnpn_pres_res_t`. Then in `pnpn_res_device.F90` the type is extended:
```fortran
  type, public, extends(pnpn_prs_res_t) :: pnpn_prs_res_device_t
   contains
     procedure, nopass :: compute => pnpn_prs_res_device_compute
  end type pnpn_prs_res_device_t
``` 
An inspection of the `compute` routine reveals a complex code that both
launches kernels and also makes use of the operators and math routines that also
dispatch to the correct backend. Similar types are defined for the CPU and also
other backends (SX) in this case. The correct type is then instantiated at run
time in a factory routine located in `pnpn_res_fctry.f90`.

What this achieves is that you don't need the extra `if` statement to
distinguish between the CPU, the device, and the SX backends each time you want
`pnpn_res` to do its job. Instead, you have that `if` statement in the factory
routine, and you get the correct type to use for your backend for the rest of
the run time.

In the majority of situations, one can take the hit of having that extra `if`
statement and not complicate the code with the extra types and the factory.


#### Summary of build system entires
Here is concise summary of how to add device kernels to the build system:

In `src/Makefile.am`
- Add `.f90` and `.F90` files under `neko_fortran_SOURCES = \`
- Add `.c` files under `neko_c_SOURCES = \`
- Add all `.h` and all `.cl` files under `EXTRA_DIST = \`
- Add `.cu` files under 
```make
if ENABLE_CUDA
libneko_la_SOURCES += \
```
- Add `.hip` files under 
```make
if ENABLE_HIP
libneko_la_SOURCES += \
```
- Add all OpenCL `.c` files under
```make
if ENABLE_OPENCL
libneko_la_SOURCES += \
```
- Add associated code entries to `prgm_lib.h` and `prgm_lib.F90`, located in
  `/src/device/opencl/`.
- Add appropriate entries to `src/.depends_device`.

**Happy hacking!** 🍻
