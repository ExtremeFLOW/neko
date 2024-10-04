# Important types {#important_types}

\tableofcontents

This section is meant to help new developers get a quick introduction to some of
the most important types in Neko, to get an idea of where to look for what. In
this spirit, the descriptions here are not meant to be exhaustive. In many
cases, we will just point to top type in a hierarchy, and leave the reader to
explore the descendants.


## SEM foundation types

- [mesh_t](#mesh::mesh_t): The computational mesh.
- [space_t](#space::space_t): The function space in which the solution is
  sought. Contains things related to the polynomial basis within each element.
- [dofmap_t](#dofmap::dofmap_t): Map of degrees of freedom. Most importantly, it
  holds all the GLL nodes locations.
- [coef_t](#coefs::coef_t): Stores coefficients for transformation to and from
  the reference element, along with some other auxillary data. 
- [gs_t](#gather_scatter::gs_t): Gather-scatter kernels used to make the
  solution continuous, i.e. perform direct stiffness summation. 
- [field_t](#field::field_t): The main type for storing the unknowns, and
  essentially everything else that lives on the mesh. 

## Basic math routines
Here, we also list file names rather than types, since the basic math is implemented
as subroutines.

- `math.f90`: Basic math operations on raw arrays.
- `device_math.F90`: Basic math operations on device arrays.
- `field_math.f90`: Basic math operations on [field_t](#field::field_t).
- `operators.f90`: Various explicit operators, including derivatives, etc.
- [vector_t](#vector::vector_t) and [matrix_t](#matrix::matrix_t): 1D and 2D arrays
  with support for computation on [accelerators](#accelerators).

## Governing equation solvers and related types

- [case_t](#case::case_t): An abstraction for the simulation case. Stores the
  fluid and scalar solver as components, handles IO via
  [sampler_t](#sampler::sampler_t).
- [fluid_scheme_t](#fluid_scheme::fluid_scheme_t): Navier-Stokes solvers.
- [scalar_scheme_t](#scalar_scheme::scalar_scheme_t): Scalar
  advection-diffusion-reaction solvers.
- [bc_t](#bc::bc_t): Boundary conditions.
- [time_scheme_t](#time_scheme::time_scheme_t): Time integration schemes.
- [source_term_t](#source_term::source_term_t): Source terms.

## Singletons

Singleton types are meant to only have a single object of their kind to be
created. These objects are declared in the same module where the type resides,
and all have their name starting with `neko_`.

- [field_registry_t](#field_registry::field_registry_t): A registry of
  [field_t](#field::field_t), retrievable by name or index. This is the main
  object used to access the fields of unknowns for any place in the code.
- [scratch_registry_t](#scratch_registry::scratch_registry_t): Provides a
  mechanism to get a temporary [field_t](#field::field_t) for doing some work.
  Use this instead of creating temporary fields inside a subroutine. 
- [simcomp_executor_t](#simcomp_executor::simcomp_executor_t): Driver for
  simulation components. The object is called `neko_simcomps`.
- [log_t](#logger::log_t): Used to write to the simulation log.

## Linear algebra

- [ksp_t](#krylov::ksp_t): Krylov solvers.
- [pc_t](#precon::pc_t): Preconditioners.


