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
several sections. The user guides are intended for users of the code, while the
developer guides are intended for developers of the code. The appendix contains
additional information that is not directly related to the usage of the code.

- \subpage user-guide
  - [Installation](@ref installation)
  - [Case File](@ref case-file)
  - [User File](@ref user-file)
  - [Simulation Components](@ref simcomps)
  - [Point Zones](@ref point-zones)
  - [Statistics Guide](@ref statistics-guide)
  - [Input/Output](@ref io)
- \subpage developer-guide
  - [Contributing](@ref contributing)
  - [Development Patterns](@ref dev_patterns)
  - [Code Style](@ref code-style)
  - [Testing](@ref testing)
  - [Accelerators](@ref accelerators)
- \subpage appendices
  - [Governing Equations](@ref governing-equations)
  - [Publications](@ref publications)
