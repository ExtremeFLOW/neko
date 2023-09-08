# Programming patterns and conventions {#dev_patterns}

This section aims to summarize the programming conventions in Neko to guide all
developers into writing code with good style and facilitate the reuse of the
same programming patterns throughout the code base. It is a good idea to check
your code against the guidelines here before opening a PR.

It should be noted that these are note followed universally throughout the code
base as of now, but we should at least make all new code as clean as possible.

## A. Naming 

1. Generally all names are lowercase, connected with `_` when needed. Some
   exceptions take place, e.g. `Re` for the Reynolds number and `Xh` for the
   function space, but these should be kept to a minimum.

2. Follow established naming conventions for certain types. For example, `Xh` for the function space.  

3. Use easy-to-remember unambiguous variable names. Avoid abbreviations unless
   they are universally understandable (e.g. RANS).
   
   The basic rationale is that code is read much more often than it is written,
   so clarity is very important. Moreover, it reduces the cognitive load of the
   developers since one does not have to remember the abbreviations, which often
   follow no logic.  Longer variable names are also no longer a practical issue
   due to assistance from IDEs. Moreover, there is always the option to use
   `associate`, if contracting a variable name is really necessary.

   *Antipattern*: `usr`, `msh` , `intrp`, `tnsr`.

   *Pattern*: `user`, `mesh` , `interpolate`, `tensor`.

   This rule is
   * Critical to observe for all keywords added to the JSON case file without
     exceptions.
   * Very important to observe in any code, which is intended as `public`,
     particularly interfaces.
   * Desirable to observe even for code not derictly exposed to the outside world.
   
4. All derived types should end with `_t` in the name.

5. For a module implementing type `mytype_t`, the module should reside in
   `mytype.f90` and the name of the module should be `mytype`.

6. Prepend the implementations of type-bound procedures (TBPs) with the name of
   the type, sans the `_t`.

   This may lead to long names when both the name of the type and of the
   procedure by  themselves have long names. But since this name is fully
   internal, it is not a major concern. By following the pattern, we never have
   to think how to choose the name.

   *Antipattern*:

   ```fortran
   type, public :: interpolator_t
    contains
      procedure, pass(this) :: init => interp_init
      procedure, pass(this) :: free => interp_free
      procedure, pass(this) :: map => interpolate
   end type
   ```

   *Pattern*:

   ```fortran
   type, public :: interpolator_t
    contains
      procedure, pass(this) :: init => interpolator_init
      procedure, pass(this) :: free => interpolator_free
      procedure, pass(this) :: map => interpolator_map
   end type
   ```

7. For TBPs or type components that are `private` or just generally intended as
   internal to the type, add an underscore to the end of name. Whe you wish to
   add a getter for the component, it can be called the same name as the
   component, sans the underscore.

8. Use `init` for constructor TBPs.

9. Use `free` for destructor TBPs.

## B. Scope

1. Always use `only` when `using` something from another module.  The `neko`
   module is an exception and imports everything. This is done so that the
   `.usr` files need only `use` the `neko` module to get access to everything.

2. All modules must have `implicit none` at module level.

3. All modules should be `private` at module level.

4. Apply own taste and judgement regarding which procedures or components
   should be `private`.

   We are generally quite permissive in terms of access and rely on the
   developers' good judgement to not mess with the types in the wrong way. That
   being said, using `private` is not discouraged in any way.

## C. Documentation

1. All procedures and interfaces must be documented, with the *minimal*
   requirement being a one-liner telling what the procedure does followed by
   `@param`s describing all the dummy  arguments. 

2. Feel free to add more documentation under `@details` or other Doxygen
   decorators. It is highly encouraged.

3. The one-liner describing the TBP should also be in the type definition. It is
   best to match the text in both of these places, since it is unclear which one
   will be picked up by Doxygen.

4. When a module implements one type, the documentation should be for the type.
   The module documentation can simply be "implements type `mytype_t`".

5. All components of the type should be documented. Put the docstring above the component.

   The alternative is to put it on the same line as the component, but it only
   works for short docstrings and cannot be multiline. Since a future component
   may require a long docstring, and mixing the styles does not look nice, it is
   better to default to putting the string above.














