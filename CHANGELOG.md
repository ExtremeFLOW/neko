# Changelog

## Develop

- *BREAKING*, normal_outflow conditions now require specifying `value`, which
  is used to set the value of the tangential components of velocity.
- Added ALE framework.
- Added masked I/O capabilities for the field_writer via the optional 
  `point_zone` JSON keyword.
- Added the user-defined Neumann boundary conditions for the scalar solver.
- *BREAKING* Changed the user-defined scalar Dirichlet boundary conditions
  keyword from `user` to `user_dirichlet`.
- Add possibility to create mesh and dofmap objects from masked entries.
- Enabled 1D stats files in csv format as a possible input to `average_fields_in_time`.
- Added the possibility to configure interpolation parameters for `probes`.
- *BREAKING* Changed the user interface of fluid/scalar initial condition
  to read interpolation parameters from the `interpolation` JSON subdict
  instead of individual parameters.
- Added public variables `GLOBAL_INTERP_PAD` and `GLOBAL_INTERP_TOL`
  in `global_interpolation` as default values for `tolerance` and
  `padding` parameters.
- Added the possibility to initialize `global_interpolation` from a JSON subdict.
- Added the `json_get_subdict_or_empty` which seeks a JSON subdict and returns
  an empty object if not found.
- If the "field" is provided for `scalar_stats` in the case file, append this
  name to the `scalar_stats` registry prefix and the default output filename.
- Added a script to add new unit tests under `contrib/add_unit_test`. Added
  templates for serial and parallel unit tests.
- Added optional log output from the flow_rate_force, controlled by the `log`
  parameter.
- Increased precision of the time value in the log.
- Added a script to add new unit tests under `contrib/add_unit_test`. The same
  script can add a .pf file to an existing suite.
- Bugfix: Fixed a bug in the `unmap` subroutine, where the device pointer was
  used to check if the field was mapped, which lead to a crash when trying to
  unmap an array that was not associated with a device. Correctly does nothing
  now.
- Added an AI policy to the contribution guidelines.
- Added simple support for VTKHDF. For now it can be used for fluid outputs.
  Simple restarts are supported with fixed mesh and MPI configuration.
  The VTKHDF output format is still experimental and will change in the future.
- Added templates for serial and parallel unit tests.
- Added code review instructions for LLMs in a copilot-friendly location.
- Improved pixi installation. Added support to create a Python environment
  inside the pixi shell. Added support to choose real precision.
- Added the Deardorff SGS model.
- Added the optional `expected_size` argument to `json_get_*_array`
  to throw an error if the parsed array size is incorrect.
- Fixed checkpoint JSON parameter parsing and their documentation. The
  `output_checkpoints` parameter no longer has a default value.
- Added runtime statistics for subgrid-scale contribution to the anisotropic part
  of the residual stresses.
- Introduced `import_fields`: a subroutine to read and import fld data,
  with interpolation capabilities.
- Added `vector_list_t` and `name` to `vector_t`.
- Rework hash table iterators, significantly faster (O(tsize) => O(entries)
- Remove redundant directory in `site-packages` when installing pyneko
- Added options to used masked parts of the domain when performing interpolation
- Update the simcomp wrappers to better handle allocation and deallocation.
- Added a factory subroutine for scalar schemes, allowing for more flexible
  creation of scalar scheme objects based on JSON input.
- Fixed a bug in the scalar scheme handler where polymorphic objects were not
  being handled correctly.
- Support for user-defined scalar schemes are now added.
- Added source term for direct forcing from a field defined in the registry.
- Add description of the `fld` file format to the documentation.
- Added possibility to assign names to boundary conditions in the case file. The
  `bc_list_t` now supports item retrieval by name or zone_index.
- Runtime statistics fields are now retrievable from the registry, for both
  fluid_stats and user_stats. The naming convention of the fields in the
  registry is `name_of_simcomp + "/mean_" + name_of_field`.
- Updated field types with a wrapper and ensure lifetime management of field
  data in field lists and arrays.
- Updated Developer Patterns documentation with new information on how to manage
  pointers and lists of complex objects.
- Added cache cleanup job for CI workflows upon PR closure.
- Updated compiler check workflows to run on release branches and master.
- Removed commented-out workflow sections.
- Added compiler support section to README.
- Restrict the `setuptools` version to be less than 81, due to a breaking change
  in that version for flinter.
- Added the `full_elements` option to point_zones. Allows including all points
  in an element in the mask.
- *BREAKING* The sign of the Boussinesq source term is fixed such that the input
  gravity vector could be prescribed correctly.
- Added an option for writing the mesh in every output field file.
- *BREAKING* All simcomps now have a `name` keyword in the case file. A default
  name is assigned, but all `name`s must be unique. If you have two or more
  simcomps of the same `type`, you must manually provide each a unique `name`.
- *BREAKING* JSON case file parsing now uses strict type checking. This means,
  for example, that providing an integer like 2 for a real entry will throw an
  error, one should set 2.0. Descriptive error and warning messages are issued.
- Added the possibilty to provide global constants in the case file under the
  `constants` object.
  - Added real scalar entries to `registry_t`.
  - Added `neko_const_registry` to store global constants defined in the case
    file.
  - Added submodule `case_file_utils` to `json_utils` for extracting JSON
    entry values from either the JSON itself or the `neko_const_regitry`.
- Introduce deprecation warning functionality. Allowing marking functions
  and classes as deprecated, with optional custom messages.
- Add missing free operators for `output_t` class.
- Add min/max operations when applying strong boundary conditions for the
  scalar, mimicing the procedure for the fluid. Needed with meshes where an
  element touches the boundary with only an edge.
- Fix `user` scalar boundary conditions only being applied once in the beginning
  of the simulation.
- Fix `mean_field_output_t` initialization, causing `start_time` to not be
  respected by the `user_stats` simulation component.
- Fix field assignment operator to correctly handle name assignment only when
  the current field's name is empty. Caused HDF5 bugs when writing fields with
  pre-existing names.
- Fix cyclic boundary rotation device bug, which tried to launch kernels
  with zero threads for ranks not containing cyclic boundaries.
- Change default parameters for tamg and phmg to be less expensive.

### Deprecated features
- `operator::ortho` calls with implicit device arrays are deprecated. Please use
  `device_ortho` instead.
