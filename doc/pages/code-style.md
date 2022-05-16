# Code style {#code-style}
Fortran code should conform to the Fortran 2008 standard and should use an indentation level of 2, except for the extra indentation within `do` `if`, `select` or `where` statements and for each level inside a structure e.g. `type`, `interface`, where the indentation level is 3, and a 0 indentation is used for `module` or `contains` (except for `contains` inside a derived type, where a single indentation level is used).

These are the default rules in Emacs' Fortran mode, an example is given below,

~~~~~~~~~~~~~~~{.f90}
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
~~~~~~~~~~~~~~~

Please note that the maximum line length in Neko should not exceed 80 columns.

## Data types
For portability reasons, it is essential to only use data type kinds defined in .num_types.f90 and avoid legacy constructs like `real*8` or `integer(kind=8)`

Floating-point numbers should be declared using `real` with the kind `rp`, which is the configured working precision for Neko (defaults to double). If single, double or quad precision is explicitly needed, use the kinds `sp`, `dp` or `qp`, respectively. For 16, 32, or 64-bit integers, Neko has defined the kinds ` i2`, `i4` or `i8`, respectively; however, for standard integers, it is perfectly fine to use a plain `integer`.
