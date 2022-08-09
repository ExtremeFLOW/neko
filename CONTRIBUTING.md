# Contributing to Neko 
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

### Data types
For portability reasons, it is essential to only use data type kinds defined in `src/config/num_types.f90` and avoid legacy constructs like `real*8` or `integer(kind=8)`

Floating-point numbers should be declared using `real` with the kind `rp`, which is the configured working precision for Neko (defaults to double). If single, double or quad precision is explicitly needed, use the kinds `sp`, `dp` or `qp`, respectively. For 16, 32, or 64-bit integers, Neko has defined the kinds ` i2`, `i4` or `i8`, respectively; however, for standard integers, it is perfectly fine to use a plain `integer`.

## Build system
This section contains information on how to add new source files to the build system. _Note that this section currently only covers Fortran code. It will be updated with information on how to add accelerator code and unit tests in a near future._

Neko uses Autotools for building all sources. You will need to have at least `autoconf` and `automake` installed for development work. It is also highly recommended to have `makedepf90` installed to avoid error-prone manual dependency tracking.

The following steps describe how to add a new Fortran file to Neko`s build system
1. Place the file in an appropriate subdirectory under `src/`. Either create a new subdirectory or place the file in `common` if none of the existing directories is a good match. Avoid placing the file directly under `src`
2. Add the file to the `neko_fortran_SOURCES` list in `src/Makefile.am`, following the pattern of `<subdir under src>/newfile.f90`.
3. Ensure correct dependency tracking
    * Using `makedepf90`
      - Regenerate build system by running `./regen.sh` at the top level, and reconfigure using `configure`
      - Regenerate dependency file (`src/.depends`) by issuing `make depend` under `src`
      - Do a final regeneration of the build system (`./regen.sh`)
    * Manually 
      - Manually add the new file and all its dependencies to `src/.depends`, using the following pattern, 
      ```Make
      directory/newfile.o : directory/newfile.f90 directory/dependency.o anotherdirectory/anotherdependency.o
      ```
      - Regenerate build system by running `./regen.sh` at the top level.
4. Finally reconfigure using `configure` and rebuild Neko with your new contribution!

For more information, please refer to the documentation https://extremeflow.github.io/neko

**Happy hacking!** 🍻
