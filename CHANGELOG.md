# Changelog

## Develop

- The CUDA and HIP auto-tuners for `ax_helm` and the SEM operators (`opgrad`,
  `dudxyz`, `cdtp`, `conv1`, `convect_scalar`, `lambda2`) now also sweep the
  thread block geometry, not just the kernel formulation: chunk size for the
  1d variants and elements per block for the kstep variants. Candidates are
  timed in interleaved rounds and min-reduced, and all of them are logged.
  New environment variables `NEKO_EB_TUNE`, `NEKO_EB`, `NEKO_CHUNKS`,
  `NEKO_TUNE_ROUNDS` and `NEKO_TUNE_ITERS` control the search; the
  elements-per-block sweep is on by default for CUDA and off for HIP, where
  the blocked kernels spill.
- Fixed threaded OpenSHMEM startup, which called `shmem_init_thread` twice
  when the first call provided less than `SHMEM_THREAD_MULTIPLE` (Cray LIBSMA
  aborts on a second initialisation), and then rejected any level below
  `SHMEM_THREAD_FUNNELED` even though the host backend funnels its calls.
- Fixed material property handling in `scalar_pnpn`: user-defined properties
  are now updated before right-hand-side assembly, Neumann values are treated
  as physical fluxes, and spatially varying `rho * cp` is applied consistently.
- Added OpenMP threading to the Coarray Fortran and host OpenSHMEM
  gather-scatter backends, which duplicated the one-sided puts and the
  signalling on every thread in hybrid MPI+OpenMP runs.
- Added `boundary_data_t`, `wall_shear_stress_simcomp`, and `boundary_data_writer_simcomp`.
- Fixed a segfault at initialisation in the OpenSHMEM gather-scatter backend
  (`NEKO_GS_COMM=SHMEM`), which freed its own dof lists at the start of
  `init` and then read from them to size the symmetric buffers.
- Added runtime autotuning of the host gather-scatter communication backend,
  mirroring the device MPI strategy tuning: with `NEKO_GS_COMM` unset, each
  `gs_t` benchmarks the host backends the build supports and keeps the
  fastest one. Setting `NEKO_GS_COMM` pins a backend and skips the benchmark.
- Added `NEKO_GS_CAF_SIGNALING=auto`, selecting the fastest coarray
  signalling mode by benchmarking. Bound once per program, and only in
  effect when the comm. backend is autotuned as well.
- Fixed a data race in openMP block in `adv_dealias` for scalar and ALE.
- Added `phmg_update` to propagate mesh change to coarse level grids.
- Fixed extrusion of curved edges when reading a 2D .nmsh file.
- Added mathematical expressions as case file values, available as the
  `expression` initial condition for the fluid and the scalar, and as the
  `expression_velocity`, `expression_pressure` and `expression_dirichlet`
  boundary conditions. Expressions can use `x`, `y`, `z`, `t`, `dt`, `pi`,
  elementary functions, and any scalar declared under `case.constants`, so
  simple spatially varying conditions no longer require a user file.
- Added integration test for ALE (test_ale).
- Updated simulation_components documentation to mirror the latest codebase.
- Fixed stale accumulator in the SX gather-scatter backend (min/max/mul).
- *BREAKING* Renamed the allocation-only `precon_factory` API to
  `precon_allocator`. Added runtime registration of user-defined
  preconditioner and Krylov solver types.
- Added opt-in zero-copy unified memory mapping for the HIP backend on AMD
  MI300A APUs: with `NEKO_HIP_ZEROCOPY=1` (and `HSA_XNACK=1`), mapped arrays
  alias their host allocation instead of being replicated on the device,
  roughly halving the memory footprint of mapped data. Off by default
  (currently slower per step than replicated buffers on MI300A); discrete GPUs
  always use replicated buffers.
- Added `only_facets = .true.` to `bc%finalize` in force_torque.
- Removed the hardcoded `CUDA_ARCH="-arch sm_60"` fallback used for
  `--enable-device-mpi` CUDA builds, which CUDA >= 13 toolkits can no longer
  compile (Pascal support was dropped). `configure` now requires `CUDA_ARCH`
  to be set explicitly for such builds and errors out otherwise, rather than
  silently guessing an architecture that may not match the target GPU.
- Add physical viscous and conductive flux support to the compressible solver.
- Fixed stale facet-normals in fluid_pnpn pressure surface terms.
  This was only affecting ALE simulations containing rotations.
- *BREAKING* Changed the size of mesh velocity lag arrays in ALE to 2. 
  Old ALE restart files (before this commit) should not be used anymore. 
  Non-ALE restart files are totally unaffected.
- Added CSV/HDF5 trajectory output options for the `lagrangian_particles`
  simcomp, including `output_format`, `snapshots_per_file`, documentation, and
  the `rebounding_particles` example.
- Fixed NVSHMEM support for NVSHMEM 3.x, which no longer ships the combined
  `libnvshmem`; `configure` now detects the split `libnvshmem_host` and
  `libnvshmem_device` libraries (falling back to `-lnvshmem` for older
  releases). The device link is now performed at build time via
  `nvcc -dlink`, so final binaries link with the Fortran compiler as usual:
  NVSHMEM builds no longer require a manual `nvcc` link or `-rdc=true` in
  `CUDA_CFLAGS`, and no longer disable the `neko` binary, contrib tools,
  or the unit tests.
- Modified the additive reduction routines in `device_math` (`glsum`, `glsc2`,
  `glsc3`, `glsubnorm`, and `glsc3_many`) to do reductions in extended
  precision instead of `rp`. This improves numerical robustness when running
  in single precision.
- Fixed a linking failure with CUDA 13, which no longer implicitly links the
  host C++ runtime (`undefined reference to __cxa_guard_acquire`); `configure`
  now links `libstdc++` explicitly for CUDA >= 13, except with the Cray
  compiler (CCE/PrgEnv-cray) which provides its own C++ runtime.
- Added a native Tofu (uTofu) gather-scatter backend for Fugaku and other
  Tofu-D systems (`NEKO_GS_COMM=UTOFU`, `--with-utofu`), using one-sided
  RDMA puts with fused arrival signalling and thread-parallel injection;
  tunable via the `NEKO_GS_UTOFU_*` environment variables.
- Added a gather-scatter backend using MPI-3 neighbourhood collectives
  (`NEKO_GS_COMM=NEIGHBOUR`), one `MPI_Ineighbor_alltoallv` per operation
  instead of per-peer point-to-point messages.
- Added a Coarray Fortran gather-scatter backend (`NEKO_GS_COMM=CAF`) using
  one-sided coarray puts, with the signalling mechanism selectable at
  runtime via `NEKO_GS_CAF_SIGNALING` (`sync`, `atomic`, or F2018 `event`).
- Added an OpenSHMEM gather-scatter backend for CPU builds
  (`NEKO_GS_COMM=SHMEM`, `--with-openshmem`) using put-with-signal, for
  systems with a native OpenSHMEM such as Cray OpenSHMEMX.
- Added a fused multi-component gather-scatter, `gs%op(u1, u2, u3, n, op)`,
  exchanging three-component halos in one communication round instead of
  three; used by curl, the fluid residuals, and the coupled Krylov solvers.
- Changed the gather-scatter setup and the mesh's external point
  connectivity to an owner-rendezvous scheme via a new crystal router,
  enabling initialisation beyond ~100k MPI ranks.
- Added OpenMP threading of the CPU backend's critical path (math and
  solver kernels, Krylov solvers and preconditioners, boundary conditions,
  dealiasing, and the host gather-scatter), making hybrid MPI+OpenMP a
  supported execution mode on CPUs rather than accelerator-only.
- Added HIP and CUDA support for ALE.
- Added `spatial_average` simcomp for spatially averaging a list of registered
  fields.
- Changed the normal vectors argument type in `setup_normals` to `vector_t` and added copy to device in the routine.
- Added new math operator for device. device_masked_copy_aligned, which performs
  a masked copy of data from one field to another, for a point zone mask.
- Job control time limits can now be specified by a flexible string format, e.g.
  "30:00" for 30 minutes, "1-00:00:00" for 24 hours, "3600" for 3600 seconds,
  etc.  
- Added a Valgrind-based regression test suite under `tests/regression/valgrind`
  for detecting memory leaks in a minimal TGV run.
- Added `contrib/icem2re2`, a user-facing utility that converts ICEM/ANSYS
  Fluent `.msh` meshes to Neko `.re2` meshes. Requires Python 3 with `numpy`,
  `scipy`, and `pymech`. Supports translational periodic boundaries via
  `--periodic periodic.json`.
- Removed restart limitations on load balancing. We now cache the the balanced
  mesh and read the cache if available. Load balance name pattern updated to
  include partition number.
- Added `host_array_t` and `device_array_t` temporary array types and support
  for requesting these through `scratch_registry_t`.
- Added the `cai_sagaut_model_ii` wall model with CPU, CUDA, HIP, and OpenCL.
- Added the `create_periodic_zones` contrib utility for converting pairs of
  labeled zones in an existing `.nmsh` mesh into periodic zones.
  backends. This model is based on the work of [Cai and Sagaut (PoF,
  2021)](https://doi.org/10.1063/5.0048563).
- Added hdf5 support for probes and added hdf5 I/O helper routines
- Added `device_coef_generate_mass`and `device_coef_generate_area_and_normal`
  for hip and cuda.
- Added the `hpfrt` source term for high-pass filter-based stabilization.
- Added the `data_streamer` simulation component, allowing data streaming
- Added `device_coef_generate_mass`and `device_coef_generate_area_and_normal`
  for hip and cuda.
- Added the Richardson wall model.
- Added the variable NEKO_VARNAME_LEN in `common/utils.f90` to set a fixed
  size for `name` attributes in e.g. `field_t` and `vector_t`.
- Added the `field_subsampler` simulation component, allowing sampling of
  fields at a different polynomial order and with masking by `point_zones`.
- Basic implementation of an overset interface boundary condition is added
- Modify field_writer and probes to by default output in
  `case.output_directory`
- Added MOST wall model and added diagnostics for the wall models.
- Added the `data_streamer` simulation component, allowing data streaming
  with ADIOS2.
- Fixed a bug (mu_msk) in `device_calc_force_array` in `force_torque.f90`.
- Added MOST wall model and added diagnostics for the wall models.
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
- `operator::dudxyz` calls with implicit device arrays are deprecated. Please
  use `opr_device_dudxyz` instead.
- `operator::div` calls with implicit device arrays are deprecated. Please use
  `div_d` instead.
- `operator::ortho` calls with implicit device arrays are deprecated. Please use
  `device_ortho` instead.
- `operator::cfl_r4` calls with implicit device arrays are deprecated. Please
  use `cfl_d` instead.
- `operator::strain_rate_r4` calls with implicit device arrays are deprecated.
  Please use `strain_rate_d` instead.
- `operator::rotate_cyc_r1` and `operator::rotate_cyc_r4` calls with implicit
  device arrays are deprecated. Please use `rotate_cyc_d` instead.
