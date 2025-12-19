# Changelog

## Develop
- Added the possibilty to provide global constants in the case file under the
  `constants` object.
  - Added real scalar entries to `registry_t`.
  - Added `neko_const_registy` to store global constants defined in the case
    file.
  - Added submodule `case_file_utils` to `json_utils` for extracting JSON
    entry values from either the JSON itself or the `neko_const_regitry`.
- Fix `mean_field_output_t` initialization, causing `start_time` to not be
  respected by the `user_stats` simulation component.
- Fix field assignment operator to correctly handle name assignment only when
  the current field's name is empty. Caused HDF5 bugs when writing fields with
  pre-existing names.
- Fix cyclic boundary rotation device bug, which tried to launch kernels
  with zero threads for ranks not containing cyclic boundaries.
