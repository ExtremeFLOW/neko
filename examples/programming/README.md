# Examples: Programming the user file {#programming-examples}
\tableofcontents

This directory contains tutorial files that demonstrate various aspects of
programming with Neko. Each file provides examples and explanations for specific
features, concepts, or workflows in Neko, helping users to customize and extend
their simulations effectively.

The files don't necessarily need to be run, but executing `prepare.sh` will
copy over a mesh and case file from the cylinder example, so you can use
`makeneko` on a `.f90` and run it.

### `user_file_template.f90`
A user file template with stubs for every possible user subroutine.

\include user_file_template.f90

### `startup_and_json.f90`
Demonstrates how to define user-specific routines and interact with the JSON
parameter dictionary used for simulation configuration.

\include startup_and_json.f90

### `fields_vectors_math.f90`
Explains how to work with the `field_t` and `vector_t` types in Neko, including
mathematical operations and lifecycle management.

\include fields_vectors_math.f90

### `registries.f90`
Shows how to use registries in Neko to manage fields and temporary data
efficiently.

\include registries.f90

### `output.f90`
Shows how to output simulation data to .fld and .csv files.

\include output.f90

### `custom_types.f90`
Shows how to add a new run-time-selectable type to Neko, using a source term
as an example.

\include custom_types.f90
