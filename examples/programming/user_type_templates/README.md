# User-type templates

Each file is a standalone module that demonstrates one of Neko's supported
user-type injection points.  Copy the file for the type you need, rename the
file, module, derived type, JSON type name, and registration routine together,
then add it to your `makeneko` command alongside the user file.

For example:

```sh
makeneko user.f90 my_source_term.f90
```

`makeneko` discovers each `module_name_register_types` routine and calls it at
startup.  The registration routine in every template is deliberately kept at
the end of its module so it is easy to find.

The templates cover every type currently supported by factory-based injection:

- `les_model_template.f90` — `les_model_t`
- `wall_model_template.f90` — `wall_model_t`
- `simulation_component_template.f90` — `simulation_component_t`
- `source_term_template.f90` — `source_term_t`
- `point_zone_template.f90` — `point_zone_t`
- `preconditioner_template.f90` — `pc_t`
- `krylov_solver_template.f90` — `ksp_t`
- `ax_helm_template.f90` — `ax_t` Helmholtz matrix-vector product
- `scalar_scheme_template.f90` — `scalar_scheme_t`

The executable bodies are safe, minimal starting points rather than useful
models: the LES and wall-model templates produce zero additional effects; the
preconditioner and Ax templates are identities; and the Krylov template stops
with an explicit error until its solve routine is implemented. Replace the
marked `TODO` sections with the actual algorithm.
