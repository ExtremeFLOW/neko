# Important types {#important_types}

\tableofcontents

This section is meant to help new developers get a quick introduction to some of
the most important types in Neko, to get an idea of where to look for what. In
this spirit, the descriptions here are not meant to be exhaustive. In many
cases, we will just point to top type in a hierarchy, and leave the reader to
explore the descendants.


## SEM foundation types {#important_types_sem-foundation}

- [mesh_t](#mesh::mesh_t): The computational mesh.
- [space_t](#space::space_t): The function space in which the solution is
  sought. Contains things related to the polynomial basis within each element.
- [dofmap_t](#dofmap::dofmap_t): Map of degrees of freedom. Most importantly, it
  holds all the GLL nodes locations.
- [coef_t](#coefs::coef_t): Stores coefficients for transformation to and from
  the reference element, along with some other auxiliary data.
- [gs_t](#gather_scatter::gs_t): Gather-scatter kernels used to make the
  solution continuous, i.e. perform direct stiffness summation.
- [field_t](#field::field_t): The main type for storing the unknowns, and
  essentially everything else that lives on the mesh.

## Basic math routines {#important_types_basic-math}
Here, we also list file names rather than types, since the basic math is
implemented as subroutines.

- `math.f90`: Basic math operations on raw arrays.
- `device_math.F90`: Basic math operations on device arrays.
- `field_math.f90`: Basic math operations on [field_t](#field::field_t).
- `vector_math.f90` and `matrix_math.f90`: Basic math operations on
  [vector_t](#vector::vector_t) and [matrix_t](#matrix::matrix_t).
- `operators.f90`: Various explicit operators, including derivatives, etc.

## Array container types {#important_types_array-containers}

Neko wraps raw Fortran arrays in a handful of container types. The main
motivation is running on [accelerators](#accelerators): each container holds
both a host array `x` and a device pointer `x_d`, and manages the allocation of
both, so that code higher up does not have to deal with device memory
explicitly. The containers also carry a `name`, which is what the registries
below use to look them up. All containers share a common lifecycle: `init` to
allocate, `free` to deallocate, `size` to query the number of entries and
`copy_from` to move data between the host and the device. The assignment
operator `=` is overloaded to copy from another container of the same kind or
to set all entries to a scalar.

- [vector_t](#vector::vector_t): A rank-1 array.
- [matrix_t](#matrix::matrix_t): A rank-2 array, with `get_nrows` and
  `get_ncols` for the dimensions. Also provides `inverse`.
- [tensor3_t](#tensor3::tensor3_t) and [tensor4_t](#tensor4::tensor4_t):
  Rank-3 and rank-4 arrays, with `get_n1`, `get_n2`, etc. for the dimensions.
  Not to be confused with the `tensor.f90` file, which implements tensor
  product operations on raw arrays.

Fortran does not allow arrays of pointers, so for every container there is a
corresponding `*_ptr_t` type holding a single pointer component `ptr`.

## Generic containers {#important_types_generic-containers}

The `adt` directory holds abstract data types that are used mostly in the
mesh and gather-scatter setup, rather than for the solution itself. As
Fortran has no templates, each comes as a family of types specialised to a
data type, with the suffix indicating the payload: `i4` for 32-bit integers,
`i8` for 64-bit integers, `r8` for double precision reals, and so on.

- [stack_t](#stack::stack_t): A dynamically growing stack, with `push`, `pop`
  and `array` to get the contents as a plain array. This is the go-to
  container when the final number of entries is not known in advance.
- [htable_t](#htable::htable_t): A hash table keyed by an integer, real,
  tuple or point, storing arbitrary data. Each variant comes with an iterator
  type.
- [uset_t](#uset::uset_t): An unordered set, built on top of the hash table.
- [tuple_t](#tuple::tuple_t): Small fixed-size tuples of integers and reals,
  used as keys and stack entries in the types above.

## Governing equation solvers and related types {#important_types_solvers}

- [case_t](#case::case_t): An abstraction for the simulation case. Stores the
  fluid and scalar solver as components, handles IO via
  [sampler_t](#sampler::sampler_t).
- [fluid_scheme_incompressible_t](#fluid_scheme_incompressible::fluid_scheme_incompressible_t):
  incompressible Navier-Stokes solvers.
- [fluid_scheme_compressible_t](#fluid_scheme_compressible::fluid_scheme_compressible_t):
  compressible Navier-Stokes solvers.
- [scalar_scheme_t](#scalar_scheme::scalar_scheme_t): Scalar
  advection-diffusion-reaction solvers.
- [bc_t](#bc::bc_t): Boundary conditions.
- [time_scheme_t](#time_scheme::time_scheme_t): Time integration schemes.
- [source_term_t](#source_term::source_term_t): Source terms.

## Singletons {#important_types_singletons}

Singleton types are meant to only have a single object of their kind to be
created. These objects are declared in the same module where the type resides,
and all have their name starting with `neko_`.

- [registry_t](#registry::registry_t): A registry of the
  [array container types](#important_types_array-containers), i.e.
  [field_t](#field::field_t), [vector_t](#vector::vector_t),
  [matrix_t](#matrix::matrix_t), [tensor3_t](#tensor3::tensor3_t) and
  [tensor4_t](#tensor4::tensor4_t), retrievable by name or index. This is the
  main object used to access the array types of unknowns for any place in the
  code.
- [scratch_registry_t](#scratch_registry::scratch_registry_t): Provides a
  mechanism to get a temporary [field_t](#field::field_t),
  [vector_t](#vector::vector_t), [matrix_t](#matrix::matrix_t),
  [tensor3_t](#tensor3::tensor3_t), [tensor4_t](#tensor4::tensor4_t), or raw
  host or device array for doing some work. Use this instead of creating
  temporary fields inside a subroutine.
- [simcomp_executor_t](#simcomp_executor::simcomp_executor_t): Driver for
  simulation components. The object is called `neko_simcomps`.
- [log_t](#logger::log_t): Used to write to the simulation log.
- [user_access_t](#user_access_singleton::usr_access_t): Used to access the
  internals of the [case_t](#case::case_t) object from within the user file.

## Linear algebra {#important_types_linear-algebra}

- [ksp_t](#krylov::ksp_t): Krylov solvers.
- [pc_t](#precon::pc_t): Preconditioners.


