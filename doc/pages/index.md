# Neko Manual {#index}

Neko is a portable framework for high-order spectral element flow simulations.
Written in modern Fortran, Neko adopts an object-oriented approach, allowing
multi-tier abstractions of the solver stack and facilitating various hardware
backends ranging from general-purpose processors, CUDA and HIP enabled
accelerators to SX-Aurora vector processors. Neko has its roots in the spectral
element code Nek5000 from UChicago/ANL, from where many of the namings, code
structure and numerical methods are adopted.

## Structure of the Manual

In order to facilitate reading of the documentation. The manual is divided into
several sections. The [user guides](@ref user-guide) are intended for users of the code, while the
[developer guides](@ref developer-guide) are intended for developers of the code. The [appendix](@ref appendices) contains
additional information that is not directly related to the usage of the code.

- \subpage user-guide
  - [Installation](@ref installation) explains how to download and compile Neko
  on your platform.
  - [Case File](@ref case-file) discusses the various parameters and options 
  for your case setup such as boundary conditions, output control, etc.
  - [User File](@ref user-file) explains all the user functions and how they can
  be used to run more advanced simulations.
  - [Simulation Components](@ref simcomps) presents some extra functionalities 
  and tools such as computation and output of additional fields, in-situ 
post-processing operations, data sampling, etc.
  - [Point Zones](@ref point-zones) allow you to select zones in the mesh for
application of source terms, initial conditions, etc.
  - [Statistics Guide](@ref statistics-guide) outlines the steps to generate
3D and 2D field statistics.
  - [Input/Output](@ref io) explains how to read and write various types of data
using the Neko framework.
- \subpage developer-guide
  - [Contributing](@ref contributing) presents basic instructions to add
your contributions to Neko.
  - [Development Patterns](@ref dev_patterns) outlines the standards to be used
when developing code in the Neko framework, such as naming conventions and
documentation.
  - [Code Style](@ref code-style) introduces some extra programming conventions
related to coding style and IDEs.
  - [Testing](@ref testing) outlines the steps to run and add unit tests to the
code base with pFUnit.
  - [Accelerators](@ref accelerators) discusses important concepts and 
conventions related to GPU programming in the Neko framework.
  - [Run-time selectable types](@ref rts_types) presents the standard programming
pattern used to select object types at run time.
- \subpage appendices
  - [Governing Equations](@ref governing-equations) used in our solvers.
  - [Publications](@ref publications)
