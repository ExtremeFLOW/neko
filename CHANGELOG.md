# Changelog

## Develop

- Rework hash table iterators, significantly faster (O(tsize) => O(entries)
- Remove redundant directory in `site-packages` when installing pyneko
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
