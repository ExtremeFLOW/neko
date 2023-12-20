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

2. Follow established naming conventions for certain types. For example, `Xh`
   for the function space.  

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
   * Desirable to observe even for code not derictly exposed to the outside
     world.
   
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

10. The dummy argument, which is `pass`ed to TBPs should be called  `this`.

## B. Scope

1. Always use `only` when `using` something from another module.  The `neko`
   module is an exception and imports everything. 
   
   This latter is done so that the user `.f90` files need only `use` the `neko`
   module to get access to everything.

2. All modules must have `implicit none` at module level.

3. All modules should be `private` at module level.

4. Apply own taste and judgement regarding which procedures or components
   should be `private`.

   We are generally quite permissive in terms of access and rely on the
   developers' good judgement to not mess with the types in the wrong way. That
   being said, using `private` is not discouraged in any way.

## C. Constructors and destructors. 

1.  If the type has a destructor, the constructor must begin with calling the
    destructor.

2.  Implement constructors and destructors as TBPs. This may change in the
    future, but is currently the safest and most OOP-friendly way to do this.

3.  Do not use the `final` attribute for destructors, the compilers don't
    support it well.

4. If a base abstract class needs a constructor, call it `init_base`. Call the
   destructor `free_base`.

5. For types initializing from JSON, let `init` be the constructor from the
   `json_file` type and `init_from_components` a constructor directly for
   component types. The `init` constructor should parse the JSON and the call
   `init_from_components`.

   This way we always have constructors for concrete types that are independent
   of JSON.

6. In the type hierarchy, `init` should be the name of the constructor, which is
   introduced as a deferred procedure in the base abstract class. Same for
   `free` for the destructor.

7. All type components should be fully initialized in the constructor, e.g.
   pointers should be associated, allocatable types allocated, etc. This implies
   that the constructor interface should be sufficient for the type to fully
   take care of its own initialization.


## D. Documentation

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

5. All types must have a docstring describing its purpose. The more detailed the
   better.

6. All components of the type should be documented. Put the docstring above the
   component.

   The alternative is to put it on the same line as the component, but it only
   works for short docstrings and cannot be multiline. Since a future component
   may require a long docstring, and mixing the styles does not look nice, it is
   better to default to putting the string above.

## E. Design

The points here are mostly things to consider rather than strict rules.
Hopefully, they can guide towards better implementation.

1. For types doing have computations, consider creating separate subroutines for
   each compute backend, and putting them in separate files in the `bcknd`.

   This way nitty-gritty optimisation can happen in `bcknd`, whereas the main
   code for the type is kept clean.

   Try to add support for as many backends as possible. CPU-only will be merged,
   but should throw a clear error message when run on an accelerator.
   

2. Wrap legacy code with wrappers following correct naming conventions.

   We will surely leverage the fact that Neko is written in Fortran, and lift
   some old Nek5000 code to quickly add some functionalities. This old code
   should be in an internal subroutine and called from a wrapper.


3. If you find yourself using `optional` dummy arguments, consider whether it
   makes more sense to split the code into separate procedures instead. However,
   `optional` is not as such disapproved of.

4. Follow the single responsibility principle for types.

   For example, don't mix I\O and other functionality into one type. It is very
   common and completely fine to create helper types encapsulating some
   necessary functionality. Don't try to map types to "real-world concepts",
   that doesn't work.

   This point represents the S in SOLID, which is a set of design principles for
   types. Highly recommended to get acquainted with all of them. 

5. When encountering code repetition, consider whether these two pieces of code
will likely change at different rates and different reasons in the future. If
so, code repetition is fine, it will eventually cease to exist naturally. 

   True code duplication is when a change in one place will also necessary a
   change in the other.









