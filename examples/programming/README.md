# Examples: Programming the user file {#programming-examples}
\tableofcontents

This directory contains tutorial files that demonstrate various aspects of
programming with Neko. Each file provides examples and explanations for specific
features, concepts, or workflows in Neko, helping users to customize and extend
their simulations effectively.

The files don't necessarily need to be run, but executing `prepare.sh` will
copy over a mesh and case file from the cylinder example, so you can use
`makeneko` on a `.f90` and run it.

### A user file template

A user file template with stubs for every possible user subroutine.

<details>
   <summary><b><u>Show code example</u></b></summary>
    \include user_file_template.f90
</details>

**Source file:** `user_file_template.f90`

### Runtime interaction with the JSON case file

Demonstrates how to define user-specific routines and interact with the JSON
parameter dictionary used for simulation configuration.

<details>
    <summary><b><u>Show code example</u></b></summary>
    \include startup_and_json.f90
</details>

**Source file:** `startup_and_json.f90`

### Fields and vectors

Explains how to work with the `field_t` and `vector_t` types in Neko, including
mathematical operations and lifecycle management.
<details>
    <summary><b><u>Show code example</u></b></summary>
    \include fields_vectors_math.f90
</details>

**Source file:** `fields_vectors_math.f90`

### Using registries

Shows how to use registries in Neko to manage fields and temporary data
efficiently.
<details>
    <summary><b><u>Show code example</u></b></summary>
    \include registries.f90
</details>

**Source file:** `registries.f90`

### Outputting data

Shows how to output simulation data to .fld and .csv files.
<details>
    <summary><b><u>Show code example</u></b></summary>
    \include output.f90
</details>

**Source file:** `output.f90`

### Adding custom run-time-selectable types 

Shows how to add a new run-time-selectable type to Neko, using a source term
as an example.
<details>
    <summary><b><u>Show code example</u></b></summary>
    \include custom_types.f90
</details>

**Source file:** `custom_types.f90`
