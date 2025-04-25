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

## Build system
This section contains information on how to add new source files to the build
system. 

Neko uses Autotools for building all sources. You will need to have at least
`autoconf`, `libtool`, and `automake` installed for development work. It is also
highly recommended to have `makedepf90` installed to avoid error-prone manual
dependency tracking. Since Neko uses `submodules`, a recent version of
`makedepf90` from  https://salsa.debian.org/science-team/makedepf90 is needed.


### CPU Fortran code
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
Here we focus on three particular backends: CUDA, HIP, and OpenCL. Whenever a
type needs to implement dedicated kernels, Neko often follows a standardized
directory and file structure. 

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
files mentioned above. Here, we focus on the folder `device`, which contain the
implementation of kernels for various backends. Note the contents: a bunch of
`.F90` files (note the big F in the extension!) and a folder for each backend
(`cuda`, `hip`, `opencl`).

Let us start with the `.F90`s. The reason for the big F is that these files
contain preprocessor directives (like `#ifdef`). This is just a general
convention, not invented by Neko authors. The `.F90` files act as a device 
abstraction layer. What that means in practice is that they contain:
- Wrapper routines, calling the correct backend kernel, depending on the backend Neko has been compiled for.
- Several `interface` blocks, for the backend kernels, which are essentially
  glue between Fortran and the actual backend code, which is written in C/C++.
  This is standard "C interoperability" stuff, which you can read about in books
  or online. Let's take a look at `device_symmetry.F90`. It contains one
subroutine, `device_symmetry_apply_vector`, that just dispatches the call to one
of the backends, via the corresponding interfaces defined in the beginning of
the file. If we now go back to `symmetry.f90`, we see that it imports this
routine, and calls it inide `symmetry_apply_vector_dev`. Since the `.F90`s are
ordinary Fortran files, they are also added to the build system in the already
familiar way.

Now we leave the familiar Fortran territory and look at individual backends.
Generally, CUDA and HIP are very similar (in fact, even the code of the kernels
is often pretty much a copy-paste), so we will look at CUDA in detail and then
quickly mentioned adjustmentst for HIP. OpenCL is a different beast and we
discuss it separately.

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

Let's zoom out again, because now we have the full picture for CUDA. So, we have
an `.f90` file, which defines the type, and via `use` statements imports
subroutines that are low level kernels. The `.F90` file, that further dispatches
the call to a routine for the correct backend. The routines are actually C
routines, and the `.F90` defines the corresponding `interface`s for Fortran-C
interoperability. The C routine is in the `.cu` file, and it launches the CUDA
kernel located in the corresponding `_kernel.h` file, which it `#include`s.

Now, let's see how the `.cu` and `.h` files are properly added to the build
system. First, consider the `src/Makefile.am` file. The `.cu`s are added under 
```make
if ENABLE_CUDA
libneko_la_SOURCES += \
```
list. On the other hand, all the `.h` files are added to the 
```
EXTRA_DIST = \
```
list.

Unfortunately, we also have to manually provide the dependencies between the
files. That is done in the `src/.depends-device` file. A typical line will bind
a `.lo` file with the corresponding `.cu` and `.h` files. For example:
```
bc/bcknd/device/cuda/symmetry.lo : bc/bcknd/device/cuda/symmetry.cu bc/bcknd/device/cuda/symmetry_kernel.h
```
So, one should make similar entries for each new `.cu` file one adds.

#### HIP
Thankfully, HIP works in precisely the same way as CUDA, with the following
differenece: instead of `.cu` we have `.hip` files. These are to be added under
```make
if ENABLE_HIP
libneko_la_SOURCES += \
```
for the build system to see them. The `.h` files go to `EXTRA_DIST` as with
CUDA. The dependencies are also specified in the exact same way.

#### OpenCL
To follow along, navigate to `src/bc/bcknd/device/opencl`. Quite similar to CUDA, there are up to 3 types of files that can be present here.
- The `.c` files, which are just like `.hip` or `.cu`. They contain a C routine, which launches the OpenCL kernel.
- The `_kernel.cl` files, which are like the `_kernel.h` files for CUDA and HIP. You will notice that the `.c` files `#include` the `.cl` files, but with an additional `.h` in the end like

```c
#include "symmetry_kernel.cl.h"
```
Don't let that lead you astray with your file naming: the `.h` is added automatically, and the files are just `.cl`!
- There may possibly be auxillary `.h` files for common functionlity across kernels.

Note that the `.c` files contain also the following include.
```c
#include <device/opencl/prgm_lib.h>
```
This is an important file, which we have to modify to register our kernel. Each
OpenCL kernel as an associated `program`, for example if you look at
`symmetry.c`, you will see `symmetry_program` mentioned. This is not really important for the implementation per se, but the point is that twe have to define these `*_program` variables somewhere, and in Neko this is done in two files.
One is the above mentioned `src/device/opencl/prgm_lib.h`, where you can find
```c
/** Device Symmetry kernels */
extern void *symmetry_program;
```
One needs to add a similar line for the `_program` corresponding to the new
kernel. The other file is very similar, but written in Fortran: `src/device/opencl/prh_lib.F90`. Here you can e.g. find
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

#### Summary
Here is concise summary of how to add device kernels to the build system:

In `src/Makefile.am`
- Add `.f90` and `.F90` files under `neko_fortran_SOURCES = \`
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
