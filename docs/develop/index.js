var index =
[
    [ "Structure of the Manual", "index.html#autotoc_md68", null ],
    [ "User guide", "dd/d04/user-guide.html", [
      [ "Installing Neko", "d5/dfc/installation.html", [
        [ "Building from source", "d5/dfc/installation.html#autotoc_md142", [
          [ "Dependencies", "d5/dfc/installation.html#autotoc_md143", [
            [ "Building JSON Fortran", "d5/dfc/installation.html#autotoc_md144", null ],
            [ "Building HDF5 (optional, but highly recommended)", "d5/dfc/installation.html#autotoc_md145", null ],
            [ "Building ParMETIS (optional)", "d5/dfc/installation.html#autotoc_md146", null ],
            [ "Bulding PFunit (optional)", "d5/dfc/installation.html#autotoc_md147", null ],
            [ "Python (optional)", "d5/dfc/installation.html#deps-python", null ]
          ] ],
          [ "Building Neko", "d5/dfc/installation.html#building-neko", [
            [ "Compiling Neko for CPU or SX-Aurora", "d5/dfc/installation.html#autotoc_md148", null ],
            [ "Compiling Neko for NVIDIA GPUs", "d5/dfc/installation.html#autotoc_md149", null ],
            [ "Compiling Neko for AMD GPUs", "d5/dfc/installation.html#autotoc_md150", null ],
            [ "Compiling Neko with a collective communications library", "d5/dfc/installation.html#autotoc_md151", null ]
          ] ]
        ] ],
        [ "Installing via Spack", "d5/dfc/installation.html#autotoc_md152", [
          [ "Quick start guide with Spack", "d5/dfc/installation.html#autotoc_md153", null ]
        ] ],
        [ "Installing using pixi", "d5/dfc/installation.html#autotoc_md154", null ],
        [ "Using a Docker container", "d5/dfc/installation.html#autotoc_md155", null ],
        [ "Testing", "d5/d75/testing.html", [
          [ "pFUnit", "d5/d75/testing.html#autotoc_md64", null ],
          [ "Configuring Neko", "d5/d75/testing.html#autotoc_md65", null ],
          [ "Running the tests", "d5/d75/testing.html#autotoc_md66", null ],
          [ "Adding a new test", "d5/d75/testing.html#autotoc_md67", null ]
        ] ]
      ] ],
      [ "Meshing", "d9/df2/meshing.html", [
        [ "General considerations", "d9/df2/meshing.html#general-considerations", null ],
        [ "Constructing meshes", "d9/df2/meshing.html#constructing-meshes", null ],
        [ "Adding periodicity to an existing <tt>.nmsh</tt>", "d9/df2/meshing.html#adding-periodicity-to-an-existing-nmsh", null ]
      ] ],
      [ "Case File", "dd/d33/case-file.html", [
        [ "High-level structure", "dd/d33/case-file.html#autotoc_md72", null ],
        [ "Output frequency control", "dd/d33/case-file.html#autotoc_md73", null ],
        [ "The case object", "dd/d33/case-file.html#autotoc_md74", [
          [ "Constants", "dd/d33/case-file.html#autotoc_md75", null ],
          [ "Time control", "dd/d33/case-file.html#autotoc_md76", null ],
          [ "Restarts and joblimit", "dd/d33/case-file.html#autotoc_md77", null ],
          [ "Boundary type numbering in the \"output_boundary\" field", "dd/d33/case-file.html#autotoc_md78", null ]
        ] ],
        [ "Numerics", "dd/d33/case-file.html#autotoc_md79", null ],
        [ "Fluid", "dd/d33/case-file.html#case-file_fluid", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md80", null ],
          [ "Compressible flows", "dd/d33/case-file.html#autotoc_md81", [
            [ "Compressible boundary conditions", "dd/d33/case-file.html#autotoc_md82", null ]
          ] ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md83", null ],
          [ "Schwarz iterations", "dd/d33/case-file.html#autotoc_md84", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#case-file_fluid-boundary-conditions", [
            [ "Specifying the boundaries", "dd/d33/case-file.html#autotoc_md85", null ],
            [ "Available conditions", "dd/d33/case-file.html#autotoc_md86", null ],
            [ "MOST wall model", "dd/d33/case-file.html#most-wall-model", null ]
          ] ],
          [ "Richardson wall model", "dd/d33/case-file.html#richardson-wall-model", null ],
          [ "Initial conditions", "dd/d33/case-file.html#case-file_fluid-ic", null ],
          [ "Source terms", "dd/d33/case-file.html#case-file_fluid-source-term", [
            [ "Brinkman", "dd/d33/case-file.html#autotoc_md87", null ],
            [ "Gradient Jump Penalty", "dd/d33/case-file.html#autotoc_md88", null ],
            [ "Sponge", "dd/d33/case-file.html#autotoc_md89", null ]
          ] ],
          [ "Arbitrary Lagrangian-Eulerian Framework", "dd/d33/case-file.html#case-file_fluid-ale", [
            [ "Solver", "dd/d33/case-file.html#case-file_fluid-ale-solver", [
              [ "Output Files and Diagnostics", "dd/d33/case-file.html#autotoc_md90", null ]
            ] ],
            [ "Mesh preview", "dd/d33/case-file.html#autotoc_md91", null ],
            [ "Bodies", "dd/d33/case-file.html#autotoc_md92", [
              [ "Oscillation", "dd/d33/case-file.html#autotoc_md93", null ],
              [ "Rotation", "dd/d33/case-file.html#autotoc_md94", null ],
              [ "Pivot", "dd/d33/case-file.html#case-file_fluid-ale-pivot", null ],
              [ "Mesh Stiffness", "dd/d33/case-file.html#case-file_fluid_ale_stiff_geom", null ]
            ] ],
            [ "Restarting ALE simulations", "dd/d33/case-file.html#autotoc_md95", null ]
          ] ]
        ] ],
        [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md96", [
          [ "Multilevel preconditioners", "dd/d33/case-file.html#autotoc_md97", null ],
          [ "Flow rate forcing", "dd/d33/case-file.html#autotoc_md98", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md99", null ]
        ] ],
        [ "Scalar", "dd/d33/case-file.html#case-file_scalar", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md100", null ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md101", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#autotoc_md102", null ],
          [ "Initial conditions", "dd/d33/case-file.html#autotoc_md103", null ],
          [ "Source terms", "dd/d33/case-file.html#autotoc_md104", null ],
          [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md105", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md106", null ]
        ] ],
        [ "Simulation components", "dd/d33/case-file.html#autotoc_md107", null ],
        [ "Point zones", "dd/d33/case-file.html#autotoc_md108", null ],
        [ "Runtime statistics", "dd/d33/case-file.html#autotoc_md109", null ]
      ] ],
      [ "User File", "d6/def/user-file.html", [
        [ "Compiling and running", "d6/def/user-file.html#autotoc_md195", null ],
        [ "High-level structure", "d6/def/user-file.html#autotoc_md196", null ],
        [ "Default user functions", "d6/def/user-file.html#autotoc_md197", [
          [ "Initializing and finalizing", "d6/def/user-file.html#user-file_init-and-final", null ],
          [ "Computing at every time step", "d6/def/user-file.html#user-file_user-check", null ],
          [ "Setting material properties", "d6/def/user-file.html#user-file_mat-prop", null ],
          [ "Runtime mesh deformation", "d6/def/user-file.html#user-file_user-mesh-setup", null ]
        ] ],
        [ "Case-specific user functions", "d6/def/user-file.html#autotoc_md198", [
          [ "Fluid and Scalar initial conditions", "d6/def/user-file.html#user-file_user-ic", null ],
          [ "Fluid and scalar source terms", "d6/def/user-file.html#user-file_user-f", null ],
          [ "Dirichlet boundary conditions", "d6/def/user-file.html#user-file_field-dirichlet-update", null ],
          [ "Neumann boundary conditions", "d6/def/user-file.html#user-file_field-neumann-update", null ]
        ] ],
        [ "Arbitrary Lagrangian-Eulerian (ALE) user functions", "d6/def/user-file.html#user-file_ale", [
          [ "Custom rigid body motion", "d6/def/user-file.html#user-file_ale-rigid_motion", null ],
          [ "Mesh Velocity", "d6/def/user-file.html#user-file_ale-mesh-velocity", null ],
          [ "Custom Base Shapes", "d6/def/user-file.html#user-file_ale-base-shapes", null ]
        ] ],
        [ "Additional remarks and tips", "d6/def/user-file.html#autotoc_md199", [
          [ "Running on GPUs", "d6/def/user-file.html#user-file_tips_running-on-gpus", null ],
          [ "Running in multiple program multiple data (MPMD) mode", "d6/def/user-file.html#user-file_tips_mpmd", [
            [ "Custom GPU kernels", "d6/def/user-file.html#user-file_tips_running-on-gpus-custom-kernels", null ]
          ] ],
          [ "Registries", "d6/def/user-file.html#user-file_tips_registries", null ],
          [ "User access to solver internals", "d6/def/user-file.html#user-file_access", null ]
        ] ]
      ] ],
      [ "Simulation components", "d3/d84/simcomps.html", [
        [ "What are simulation components?", "d3/d84/simcomps.html#autotoc_md178", null ],
        [ "Adding simulation components to the case", "d3/d84/simcomps.html#autotoc_md179", null ],
        [ "List of simulation components", "d3/d84/simcomps.html#autotoc_md180", null ],
        [ "Controlling execution and file output", "d3/d84/simcomps.html#autotoc_md181", [
          [ "Differential operators", "d3/d84/simcomps.html#autotoc_md182", [
            [ "derivative", "d3/d84/simcomps.html#simcomp_derivative", null ],
            [ "curl", "d3/d84/simcomps.html#simcomp_curl", null ],
            [ "divergence", "d3/d84/simcomps.html#simcomp_divergence", null ]
          ] ],
          [ "gradient", "d3/d84/simcomps.html#simcomp_gradient", null ],
          [ "weak_gradient", "d3/d84/simcomps.html#simcomp_weak_gradient", null ],
          [ "lambda2", "d3/d84/simcomps.html#simcomp_lambda2", null ],
          [ "boundary_operation", "d3/d84/simcomps.html#simcomp_boundary_operation", null ],
          [ "boundary_flux", "d3/d84/simcomps.html#simcomp_boundary_flux", null ],
          [ "probes", "d3/d84/simcomps.html#simcomp_probes", [
            [ "Supported types", "d3/d84/simcomps.html#autotoc_md183", null ],
            [ "Example usage", "d3/d84/simcomps.html#autotoc_md184", null ]
          ] ],
          [ "field_writer", "d3/d84/simcomps.html#simcomp_field_writer", null ],
          [ "force_torque", "d3/d84/simcomps.html#simcomp_force_torque", [
            [ "Torque calculation for moving bodies", "d3/d84/simcomps.html#autotoc_md185", null ]
          ] ],
          [ "les_model", "d3/d84/simcomps.html#simcomp_les_model", null ],
          [ "User statistics", "d3/d84/simcomps.html#user_stats", null ],
          [ "Spectral error indicator", "d3/d84/simcomps.html#simcomp_speri", null ],
          [ "Data streamer", "d3/d84/simcomps.html#simcomp_data_streamer", null ],
          [ "Field subsampler", "d3/d84/simcomps.html#simcomp_field_subsampler", null ]
        ] ]
      ] ],
      [ "Point zones", "da/dd0/point-zones.html", [
        [ "What are point zones?", "da/dd0/point-zones.html#autotoc_md168", null ],
        [ "Predefined geometrical shapes", "da/dd0/point-zones.html#autotoc_md169", [
          [ "Box", "da/dd0/point-zones.html#autotoc_md170", null ],
          [ "Sphere", "da/dd0/point-zones.html#autotoc_md171", null ],
          [ "Cylinder", "da/dd0/point-zones.html#autotoc_md172", null ]
        ] ],
        [ "Operations on point zones", "da/dd0/point-zones.html#autotoc_md173", [
          [ "Inversion", "da/dd0/point-zones.html#autotoc_md174", null ],
          [ "Combination", "da/dd0/point-zones.html#autotoc_md175", null ],
          [ "Including full element data", "da/dd0/point-zones.html#autotoc_md176", null ]
        ] ],
        [ "User-defined geometrical shapes", "da/dd0/point-zones.html#autotoc_md177", null ],
        [ "Using point zones", "da/dd0/point-zones.html#point-zones_using-point-zones", null ]
      ] ],
      [ "Statistics guide", "df/d8f/statistics-guide.html", [
        [ "Fluid Statistics", "df/d8f/statistics-guide.html#fluid-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md186", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md187", null ]
        ] ],
        [ "Postprocessing", "df/d8f/statistics-guide.html#autotoc_md188", null ],
        [ "Scalar Statistics", "df/d8f/statistics-guide.html#scalar-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md189", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md190", null ]
        ] ],
        [ "Fluid Subgrid-Scale (SGS) Statistics", "df/d8f/statistics-guide.html#fluid-sgs-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md191", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md192", null ]
        ] ],
        [ "Scalar Subgrid-Scale (SGS) Statistics", "df/d8f/statistics-guide.html#scalar-sgs-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md193", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md194", null ]
        ] ],
        [ "Note to users", "df/d8f/statistics-guide.html#note-to-users", null ]
      ] ],
      [ "Input-output", "d7/d7f/io.html", [
        [ "Mesh", "d7/d7f/io.html#autotoc_md156", null ],
        [ "Three-dimensional field output", "d7/d7f/io.html#autotoc_md157", [
          [ "Compression of field output", "d7/d7f/io.html#autotoc_md158", null ]
        ] ],
        [ "Checkpoint files", "d7/d7f/io.html#autotoc_md159", null ],
        [ "VTKHDF output", "d7/d7f/io.html#vtkhdf-output", [
          [ "Prerequisites", "d7/d7f/io.html#vtkhdf-prerequisites", null ],
          [ "Enabling VTKHDF output", "d7/d7f/io.html#vtkhdf-enabling", null ],
          [ "File structure", "d7/d7f/io.html#vtkhdf-file-structure", null ],
          [ "Cell representation", "d7/d7f/io.html#vtkhdf-cell-representation", null ],
          [ "Temporal vs non-temporal output", "d7/d7f/io.html#vtkhdf-temporal-vs-non-temporal", null ],
          [ "Limitations", "d7/d7f/io.html#vtkhdf-limitations", null ]
        ] ]
      ] ],
      [ "Extending neko", "d4/d1b/extending.html", null ],
      [ "Neko API", "d8/daa/api.html", [
        [ "Installation", "d8/daa/api.html#autotoc_md69", [
          [ "Python (pyneko)", "d8/daa/api.html#autotoc_md70", null ],
          [ "Julia (Neko.jl)", "d8/daa/api.html#autotoc_md71", null ]
        ] ]
      ] ],
      [ "Performance guidelines", "dc/d3c/performance.html", [
        [ "Installation", "dc/d3c/performance.html#autotoc_md160", [
          [ "Accelerator specific options", "dc/d3c/performance.html#autotoc_md161", null ]
        ] ],
        [ "Simulation setup", "dc/d3c/performance.html#autotoc_md162", [
          [ "Load balancing", "dc/d3c/performance.html#autotoc_md163", null ],
          [ "Parameters", "dc/d3c/performance.html#autotoc_md164", null ]
        ] ],
        [ "Running a simulation", "dc/d3c/performance.html#autotoc_md165", [
          [ "Gather-scatter communication backends", "dc/d3c/performance.html#autotoc_md166", [
            [ "NVSHMEM backend", "dc/d3c/performance.html#performance-nvshmem-backend", null ],
            [ "OpenSHMEM backend", "dc/d3c/performance.html#performance-openshmem-backend", null ],
            [ "Coarray Fortran backend", "dc/d3c/performance.html#autotoc_md167", null ]
          ] ]
        ] ]
      ] ],
      [ "Global Interpolation", "dd/d61/global-interpolation.html", [
        [ "Overview", "dd/d61/global-interpolation.html#autotoc_md110", null ],
        [ "The global interpolation module", "dd/d61/global-interpolation.html#autotoc_md112", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md113", null ],
          [ "Features", "dd/d61/global-interpolation.html#autotoc_md114", null ]
        ] ],
        [ "The global interpolation type", "dd/d61/global-interpolation.html#autotoc_md116", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md117", null ],
          [ "Attributes", "dd/d61/global-interpolation.html#autotoc_md118", [
            [ "Domain and Space Information", "dd/d61/global-interpolation.html#autotoc_md119", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md120", null ],
            [ "Local Points", "dd/d61/global-interpolation.html#autotoc_md121", null ],
            [ "Interpolation Tools", "dd/d61/global-interpolation.html#autotoc_md122", null ],
            [ "Parallelism", "dd/d61/global-interpolation.html#autotoc_md123", null ],
            [ "Configuration", "dd/d61/global-interpolation.html#autotoc_md124", null ]
          ] ],
          [ "Methods", "dd/d61/global-interpolation.html#autotoc_md126", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md127", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md128", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md129", null ],
            [ "Validation", "dd/d61/global-interpolation.html#autotoc_md130", null ],
            [ "Memory Management", "dd/d61/global-interpolation.html#autotoc_md131", null ]
          ] ],
          [ "Example Usage", "dd/d61/global-interpolation.html#autotoc_md133", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md134", null ],
            [ "Finding Points", "dd/d61/global-interpolation.html#autotoc_md135", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md136", null ]
          ] ],
          [ "Notes", "dd/d61/global-interpolation.html#autotoc_md138", null ],
          [ "Environment variables", "dd/d61/global-interpolation.html#autotoc_md139", null ],
          [ "Related Modules", "dd/d61/global-interpolation.html#autotoc_md141", null ]
        ] ]
      ] ],
      [ "Filtering", "df/d4a/filter.html", [
        [ "Filter base type", "df/d4a/filter.html#filter_base_type", null ],
        [ "PDE filter", "df/d4a/filter.html#filter_pde", null ],
        [ "Elementwise filter", "df/d4a/filter.html#filter_elementwise", null ],
        [ "High-pass filter relaxation source term", "df/d4a/filter.html#filter_hpfrt", null ]
      ] ],
      [ "Examples: Programming the user file", "d5/db5/programming-examples.html", null ]
    ] ],
    [ "Developer guide", "dc/d70/developer-guide.html", [
      [ "Contributing to Neko", "d1/d5a/contributing.html", [
        [ "Git branches", "d1/d5a/contributing.html#autotoc_md40", null ],
        [ "Code style", "d1/d5a/contributing.html#autotoc_md41", [
          [ "Data types", "d1/d5a/contributing.html#autotoc_md42", null ]
        ] ],
        [ "Build system and code organization", "d1/d5a/contributing.html#autotoc_md43", [
          [ "Building CPU Fortran code", "d1/d5a/contributing.html#autotoc_md44", null ],
          [ "Device code", "d1/d5a/contributing.html#autotoc_md45", [
            [ "CUDA", "d1/d5a/contributing.html#autotoc_md46", null ],
            [ "HIP", "d1/d5a/contributing.html#autotoc_md47", null ],
            [ "OpenCL", "d1/d5a/contributing.html#autotoc_md48", null ],
            [ "Device-based type polymorphism", "d1/d5a/contributing.html#autotoc_md49", null ],
            [ "Summary of build system entires", "d1/d5a/contributing.html#autotoc_md50", null ]
          ] ]
        ] ]
      ] ],
      [ "Programming patterns and conventions", "d0/d47/dev_patterns.html", [
        [ "A. Naming", "d0/d47/dev_patterns.html#autotoc_md51", null ],
        [ "B. Scope", "d0/d47/dev_patterns.html#autotoc_md52", null ],
        [ "C. Constructors and destructors.", "d0/d47/dev_patterns.html#autotoc_md53", null ],
        [ "D. Documentation", "d0/d47/dev_patterns.html#autotoc_md54", null ],
        [ "E. Design", "d0/d47/dev_patterns.html#autotoc_md55", null ],
        [ "F. Complex and Polymorphic arrays.", "d0/d47/dev_patterns.html#autotoc_md56", [
          [ "F.1. Object wrapper pattern", "d0/d47/dev_patterns.html#autotoc_md57", null ],
          [ "F.2. Object pointer pattern", "d0/d47/dev_patterns.html#autotoc_md58", null ]
        ] ]
      ] ],
      [ "Code style", "da/db6/code-style.html", [
        [ "Data types", "da/db6/code-style.html#autotoc_md35", null ],
        [ "Enforcing rules", "da/db6/code-style.html#autotoc_md36", [
          [ "Tools", "da/db6/code-style.html#autotoc_md37", [
            [ "flint", "da/db6/code-style.html#autotoc_md38", null ],
            [ "findent", "da/db6/code-style.html#autotoc_md39", null ]
          ] ]
        ] ]
      ] ],
      [ "Testing", "d5/d75/testing.html", [
        [ "pFUnit", "d5/d75/testing.html#autotoc_md64", null ],
        [ "Configuring Neko", "d5/d75/testing.html#autotoc_md65", null ],
        [ "Running the tests", "d5/d75/testing.html#autotoc_md66", null ],
        [ "Adding a new test", "d5/d75/testing.html#autotoc_md67", null ]
      ] ],
      [ "Accelerators", "de/d06/accelerators.html", [
        [ "Device abstraction layer", "de/d06/accelerators.html#autotoc_md30", [
          [ "Memory management", "de/d06/accelerators.html#autotoc_md31", [
            [ "Allocation/deallocation", "de/d06/accelerators.html#autotoc_md32", null ],
            [ "Associate data on host and device", "de/d06/accelerators.html#autotoc_md33", null ],
            [ "Map a host array to a device", "de/d06/accelerators.html#autotoc_md34", null ],
            [ "Data transfer", "de/d06/accelerators.html#accelerators_data-transfer", null ]
          ] ],
          [ "Offload work", "de/d06/accelerators.html#accelerators_offload-work", null ]
        ] ]
      ] ],
      [ "Run-time selectable types", "d5/d5f/rts_types.html", null ],
      [ "Important types", "d3/d40/important_types.html", [
        [ "SEM foundation types", "d3/d40/important_types.html#autotoc_md59", null ],
        [ "Basic math routines", "d3/d40/important_types.html#autotoc_md60", null ],
        [ "Governing equation solvers and related types", "d3/d40/important_types.html#autotoc_md61", null ],
        [ "Singletons", "d3/d40/important_types.html#autotoc_md62", null ],
        [ "Linear algebra", "d3/d40/important_types.html#autotoc_md63", null ]
      ] ]
    ] ],
    [ "Appendices", "da/dd6/appendices.html", [
      [ "Environmental variable reference", "da/dd6/appendices.html#appendices_env-var", [
        [ "Logging level details", "da/dd6/appendices.html#autotoc_md0", null ],
        [ "Gather-scatter communication backend details", "da/dd6/appendices.html#autotoc_md1", null ],
        [ "Coarray Fortran signaling mode details", "da/dd6/appendices.html#autotoc_md2", null ]
      ] ],
      [ "Governing equations", "db/d27/governing-equations.html", [
        [ "Momentum (Fluid)", "db/d27/governing-equations.html#autotoc_md16", null ],
        [ "Scalar", "db/d27/governing-equations.html#autotoc_md17", null ],
        [ "Compressible Flows", "db/d27/governing-equations.html#autotoc_md18", null ],
        [ "Non-dimensionalisation", "db/d27/governing-equations.html#autotoc_md19", null ]
      ] ],
      [ "Neko .nmsh mesh format", "d0/d5b/nmsh-format.html", [
        [ "File overview", "d0/d5b/nmsh-format.html#autotoc_md20", null ],
        [ "Data types and portability", "d0/d5b/nmsh-format.html#autotoc_md21", null ],
        [ "Element records", "d0/d5b/nmsh-format.html#autotoc_md22", [
          [ "Vertex record (<tt>nmsh_vertex_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md23", null ],
          [ "Quad record (<tt>nmsh_quad_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md24", null ],
          [ "Hex record (<tt>nmsh_hex_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md25", null ],
          [ "Vertex ordering", "d0/d5b/nmsh-format.html#autotoc_md26", null ]
        ] ],
        [ "Zone records (boundary and periodic data)", "d0/d5b/nmsh-format.html#autotoc_md27", null ],
        [ "Curve records (curved edges/faces)", "d0/d5b/nmsh-format.html#autotoc_md28", null ],
        [ "2D meshes", "d0/d5b/nmsh-format.html#autotoc_md29", null ]
      ] ],
      [ "Neko .fld field format", "d8/d00/fld-format.html", [
        [ "File naming and series", "d8/d00/fld-format.html#autotoc_md3", null ],
        [ "File overview", "d8/d00/fld-format.html#autotoc_md4", null ],
        [ "Header", "d8/d00/fld-format.html#autotoc_md5", [
          [ "rdcode", "d8/d00/fld-format.html#autotoc_md6", null ]
        ] ],
        [ "Test pattern", "d8/d00/fld-format.html#autotoc_md7", null ],
        [ "Element index list", "d8/d00/fld-format.html#autotoc_md8", null ],
        [ "Field blocks", "d8/d00/fld-format.html#autotoc_md9", [
          [ "Vector fields (X and U)", "d8/d00/fld-format.html#autotoc_md10", null ],
          [ "Scalar fields (P, T, S)", "d8/d00/fld-format.html#autotoc_md11", null ],
          [ "Precision", "d8/d00/fld-format.html#autotoc_md12", null ]
        ] ],
        [ "Metadata blocks (3D only)", "d8/d00/fld-format.html#autotoc_md13", null ],
        [ "2D fields and Z-velocity", "d8/d00/fld-format.html#autotoc_md14", null ],
        [ "Portability notes", "d8/d00/fld-format.html#autotoc_md15", null ]
      ] ],
      [ "Publications", "de/d26/publications.html", null ]
    ] ]
];