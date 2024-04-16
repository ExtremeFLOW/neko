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
  - [Installation](#installation)
  - [Case File](#case-file)
  - [User File](#user-file)
  - [Simulation Components](#simcomps)
  - [Point Zones](#point-zones)
  - [Statistics Guide](#statistics-guide)
  - [Input/Output](#io)
- \subpage developer-guide
  - [Contributing](#contributing)
  - [Development Patterns](#dev_patterns)
  - [Code Style](#code-style)
  - [Testing](#testing)
  - [Accelerators](#accelerators)
- \subpage appendices
  - [Governing Equations](#governing-equations)
  - [Publications](#publications)
