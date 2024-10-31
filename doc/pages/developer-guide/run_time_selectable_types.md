# Run-time selectable types {#rts_types}

\tableofcontents

It is often necessary to determine the type of objects at run time, for example,
based on user input. Neko implements the factory pattern to that end. The
associated code structure can be found in many components of the code, so it
worth becoming familiar with it. Here we will use the LES models as an example.

The creation of objects is organized as follows. A base abstract type is
defined, that will serve as the parent of all the types that can be constructed
by the factory. In our example, this is `les_model_t`. A type is abstract if it
defines at least one `deferred` procedure. This means that the procedure is not
implemented in the type itself, but all descendants of the type *must* implement
it. Crucial to our setup here is that the base abstract type defines a deferred
constructor, always called `init` (with a few exceptions). This means that all
descendants can be constructed using the same subroutine dummy arguments. These
dummy arguments are defined in an `abstract interface`. The example below shows
the interface for the constructor of all LES models. It is quite common that a
deferred destructor is defined, called `free`. In the case of LES models, the
routine `compute` is also deferred; this is where the models compute the subgrid
viscosity.

~~~~~~~~~~~~~~~{.f90}
  abstract interface
     subroutine les_model_init(this, dofmap, coef, json)
       import les_model_t, json_file, dofmap_t, coef_t
       class(les_model_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       type(dofmap_t), intent(in) :: dofmap
       type(json_file), intent(inout) :: json
     end subroutine les_model_init
  end interface
~~~~~~~~~~~~~~~

As discussed above, the descendants of the abstract type implement different
variations of the functionality in question, in our case LES models. One such
type is `smagorinsky.f90`, which implements the Smagorinsky model. Of course, an
implementation of all deferred routines is present in the type.

The factory iteself is just a subroutine. It lives in a separate file, which
ends with `_fctry.f90`, so `les_model_fctry.f90` in our example. The factory
accepts a pointer to the base abstract type as one of its dummy arguments,
usually called `object`. This pointer is to be allocated to one of the
descendants of the abstract type. The factory chooses what type to allocate
simply based on a provided name. In the top of the factory module, type names
known to the factory are listed in a variable ending with `KNOWN_TYPES`. The
name is typically passed to the dummy argument called `type_name`. This name
could, for example, be read from the case file. For LES models, one provides the
`model` keyword in the case file, and this gets passed to the factory. So if
`"model": "smagorinsky"` is in the case file, the factory will allocate the
pointer to a `smagorinsky_t` and initialize it using the common `init`
constructor. Calling `compute` on the allocated pointer will then run the
routine defined for the `smagorinsky_t` type.

For the particular case of LES models, the pointer to `les_model_t` is stored in
the simcomp `les_simcomp`. Note that it is defined as `class`, not `type`. This
is necessary to make use of polymorphism, i.e. that the pointer can "act" as any
of its descendants after being allocated to one. The factory is then called in
the `init` routine of the simcomp.

~~~~~~~~~~~~~~~{.f90}
  type, public, extends(simulation_component_t) :: les_simcomp_t
    class(les_model_t), allocatable :: les_model
~~~~~~~~~~~~~~~

So, what needs to be done to add a new LES model? Simple, one creates a new type
descending from `les_model_t`, implements all the deferred routines, adds the
type name to the known types of the factory, and finally expands the `if`
statement in the factory to allocate to the new type given the corresponding
type name. The same steps are followed for a lot of other functionality! There
are some minor variations, for example, not all factories initialize the object,
some just allocate it and leave the construction to the callee. Eventually these
differences may be ironed out.
