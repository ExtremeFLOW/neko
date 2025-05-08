# Examples: Programming the user file 

This directory contains tutorial files that demonstrate various aspects of
programming with Neko. Each file provides examples and explanations for specific
features, concepts, or workflows in Neko, helping users to customize and extend
their simulations effectively.

The files don't necessarily need to be run, but executing `prepare.sh` will
copy over a mesh and case file from the cylinder example, so you can use
`makeneko` on a `.f90` and run it.

## Tutorials

### `startup_and_json.f90`
Demonstrates how to define user-specific routines and interact with the JSON
  parameter dictionary used for simulation configuration.

### `fields_vectors_math.f90`
Explains how to work with `field_t` and `vector_t` types in Neko, including
  mathematical operations and lifecycle management.

### `registries.f90`
Shows how to use registries in Neko to manage fields and temporary data
  efficiently.