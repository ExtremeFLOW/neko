var index =
[
    [ "Structure of the Manual", "index.html#autotoc_md72", null ],
    [ "User guide", "dd/d04/user-guide.html", [
      [ "Installing Neko", "d5/dfc/installation.html", [
        [ "Building from source", "d5/dfc/installation.html#autotoc_md146", [
          [ "Dependencies", "d5/dfc/installation.html#autotoc_md147", [
            [ "Building JSON Fortran", "d5/dfc/installation.html#autotoc_md148", null ],
            [ "Building HDF5 (optional, but highly recommended)", "d5/dfc/installation.html#autotoc_md149", null ],
            [ "Building ParMETIS (optional)", "d5/dfc/installation.html#autotoc_md150", null ],
            [ "Bulding PFunit (optional)", "d5/dfc/installation.html#autotoc_md151", null ],
            [ "Python (optional)", "d5/dfc/installation.html#deps-python", null ]
          ] ],
          [ "Building Neko", "d5/dfc/installation.html#building-neko", [
            [ "Compute backends and support tiers", "d5/dfc/installation.html#compute-backends", null ],
            [ "Communication backends", "d5/dfc/installation.html#communication-backends", null ],
            [ "Compiling Neko for CPU or SX-Aurora", "d5/dfc/installation.html#autotoc_md152", null ],
            [ "Compiling Neko for NVIDIA GPUs", "d5/dfc/installation.html#autotoc_md153", null ],
            [ "Compiling Neko for AMD GPUs", "d5/dfc/installation.html#autotoc_md154", null ],
            [ "Compiling Neko for Apple Silicon GPUs", "d5/dfc/installation.html#autotoc_md155", null ],
            [ "Compiling Neko with a collective communications library", "d5/dfc/installation.html#autotoc_md156", null ]
          ] ]
        ] ],
        [ "Installing via Spack", "d5/dfc/installation.html#autotoc_md157", [
          [ "Quick start guide with Spack", "d5/dfc/installation.html#autotoc_md158", null ]
        ] ],
        [ "Installing using pixi", "d5/dfc/installation.html#autotoc_md159", null ],
        [ "Installing via FreeBSD ports", "d5/dfc/installation.html#autotoc_md160", null ],
        [ "Using a Docker container", "d5/dfc/installation.html#autotoc_md161", null ],
        [ "Testing", "d5/d75/testing.html", [
          [ "pFUnit", "d5/d75/testing.html#autotoc_md68", null ],
          [ "Configuring Neko", "d5/d75/testing.html#autotoc_md69", null ],
          [ "Running the tests", "d5/d75/testing.html#autotoc_md70", null ],
          [ "Adding a new test", "d5/d75/testing.html#autotoc_md71", null ]
        ] ]
      ] ],
      [ "Meshing", "d9/df2/meshing.html", [
        [ "General considerations", "d9/df2/meshing.html#general-considerations", null ],
        [ "Constructing meshes", "d9/df2/meshing.html#constructing-meshes", [
          [ "Box meshes", "d9/df2/meshing.html#constructing-meshes-box", null ],
          [ "Converting a Nek5000 mesh", "d9/df2/meshing.html#constructing-meshes-nek5000", null ],
          [ "Experimental workflows", "d9/df2/meshing.html#constructing-meshes-experimental", null ]
        ] ],
        [ "High-order boundary representation", "d9/df2/meshing.html#high-order-boundary-representation", [
          [ "Curves stored in the mesh file", "d9/df2/meshing.html#autotoc_md166", null ],
          [ "User-defined geometry deformation", "d9/df2/meshing.html#autotoc_md167", null ]
        ] ],
        [ "Adding periodicity to an existing <tt>.nmsh</tt>", "d9/df2/meshing.html#adding-periodicity-to-an-existing-nmsh", null ]
      ] ],
      [ "Case File", "dd/d33/case-file.html", [
        [ "High-level structure", "dd/d33/case-file.html#autotoc_md76", null ],
        [ "Output frequency control", "dd/d33/case-file.html#autotoc_md77", null ],
        [ "The case object", "dd/d33/case-file.html#autotoc_md78", [
          [ "Constants", "dd/d33/case-file.html#autotoc_md79", null ],
          [ "Time control", "dd/d33/case-file.html#autotoc_md80", null ],
          [ "Restarts and joblimit", "dd/d33/case-file.html#autotoc_md81", null ],
          [ "Boundary type numbering in the \"output_boundary\" field", "dd/d33/case-file.html#autotoc_md82", null ]
        ] ],
        [ "Numerics", "dd/d33/case-file.html#autotoc_md83", null ],
        [ "Fluid", "dd/d33/case-file.html#case-file_fluid", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md84", null ],
          [ "Compressible flows", "dd/d33/case-file.html#autotoc_md85", [
            [ "Compressible boundary conditions", "dd/d33/case-file.html#autotoc_md86", null ]
          ] ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md87", null ],
          [ "Schwarz iterations", "dd/d33/case-file.html#autotoc_md88", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#case-file_fluid-boundary-conditions", [
            [ "Specifying the boundaries", "dd/d33/case-file.html#autotoc_md89", null ],
            [ "Available conditions", "dd/d33/case-file.html#autotoc_md90", null ],
            [ "MOST wall model", "dd/d33/case-file.html#most-wall-model", null ]
          ] ],
          [ "Richardson wall model", "dd/d33/case-file.html#richardson-wall-model", null ],
          [ "Initial conditions", "dd/d33/case-file.html#case-file_fluid-ic", null ],
          [ "Source terms", "dd/d33/case-file.html#case-file_fluid-source-term", [
            [ "Brinkman", "dd/d33/case-file.html#autotoc_md91", null ],
            [ "Gradient Jump Penalty", "dd/d33/case-file.html#autotoc_md92", null ],
            [ "Sponge", "dd/d33/case-file.html#autotoc_md93", null ]
          ] ],
          [ "Arbitrary Lagrangian-Eulerian Framework", "dd/d33/case-file.html#case-file_fluid-ale", [
            [ "Solver", "dd/d33/case-file.html#case-file_fluid-ale-solver", [
              [ "Output Files and Diagnostics", "dd/d33/case-file.html#autotoc_md94", null ]
            ] ],
            [ "Mesh preview", "dd/d33/case-file.html#autotoc_md95", null ],
            [ "Bodies", "dd/d33/case-file.html#autotoc_md96", [
              [ "Oscillation", "dd/d33/case-file.html#autotoc_md97", null ],
              [ "Rotation", "dd/d33/case-file.html#autotoc_md98", null ],
              [ "Pivot", "dd/d33/case-file.html#case-file_fluid-ale-pivot", null ],
              [ "Mesh Stiffness", "dd/d33/case-file.html#case-file_fluid_ale_stiff_geom", null ]
            ] ],
            [ "Restarting ALE simulations", "dd/d33/case-file.html#autotoc_md99", null ]
          ] ]
        ] ],
        [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md100", [
          [ "Multilevel preconditioners", "dd/d33/case-file.html#autotoc_md101", null ],
          [ "Flow rate forcing", "dd/d33/case-file.html#autotoc_md102", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md103", null ]
        ] ],
        [ "Scalar", "dd/d33/case-file.html#case-file_scalar", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md104", null ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md105", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#autotoc_md106", null ],
          [ "Initial conditions", "dd/d33/case-file.html#autotoc_md107", null ],
          [ "Source terms", "dd/d33/case-file.html#autotoc_md108", null ],
          [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md109", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md110", null ]
        ] ],
        [ "Simulation components", "dd/d33/case-file.html#autotoc_md111", null ],
        [ "Point zones", "dd/d33/case-file.html#autotoc_md112", null ],
        [ "Runtime statistics", "dd/d33/case-file.html#autotoc_md113", null ]
      ] ],
      [ "User File", "d6/def/user-file.html", [
        [ "Compiling and running", "d6/def/user-file.html#autotoc_md204", null ],
        [ "High-level structure", "d6/def/user-file.html#autotoc_md205", null ],
        [ "Default user functions", "d6/def/user-file.html#autotoc_md206", [
          [ "Initializing and finalizing", "d6/def/user-file.html#user-file_init-and-final", null ],
          [ "Computing at every time step", "d6/def/user-file.html#user-file_user-check", null ],
          [ "Setting material properties", "d6/def/user-file.html#user-file_mat-prop", null ],
          [ "Runtime mesh deformation", "d6/def/user-file.html#user-file_user-mesh-setup", null ]
        ] ],
        [ "Case-specific user functions", "d6/def/user-file.html#autotoc_md207", [
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
        [ "Additional remarks and tips", "d6/def/user-file.html#autotoc_md208", [
          [ "Running on GPUs", "d6/def/user-file.html#user-file_tips_running-on-gpus", null ],
          [ "Running in multiple program multiple data (MPMD) mode", "d6/def/user-file.html#user-file_tips_mpmd", [
            [ "Custom GPU kernels", "d6/def/user-file.html#user-file_tips_running-on-gpus-custom-kernels", null ]
          ] ],
          [ "Registries", "d6/def/user-file.html#user-file_tips_registries", null ],
          [ "User access to solver internals", "d6/def/user-file.html#user-file_access", null ]
        ] ]
      ] ],
      [ "Simulation components", "d3/d84/simcomps.html", [
        [ "What are simulation components?", "d3/d84/simcomps.html#autotoc_md187", null ],
        [ "Adding simulation components to the case", "d3/d84/simcomps.html#autotoc_md188", null ],
        [ "List of simulation components", "d3/d84/simcomps.html#autotoc_md189", null ],
        [ "Controlling execution and file output", "d3/d84/simcomps.html#autotoc_md190", [
          [ "Differential operators", "d3/d84/simcomps.html#autotoc_md191", [
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
            [ "Supported types", "d3/d84/simcomps.html#autotoc_md192", null ],
            [ "Example usage", "d3/d84/simcomps.html#autotoc_md193", null ]
          ] ],
          [ "field_writer", "d3/d84/simcomps.html#simcomp_field_writer", null ],
          [ "lagrangian_particles", "d3/d84/simcomps.html#simcomp_lagrangian_particles", null ],
          [ "force_torque", "d3/d84/simcomps.html#simcomp_force_torque", [
            [ "Torque calculation for moving bodies", "d3/d84/simcomps.html#autotoc_md194", null ]
          ] ],
          [ "les_model", "d3/d84/simcomps.html#simcomp_les_model", null ],
          [ "User statistics", "d3/d84/simcomps.html#user_stats", null ],
          [ "Spatial average", "d3/d84/simcomps.html#simcomp_spatial_average", null ],
          [ "Spectral error indicator", "d3/d84/simcomps.html#simcomp_speri", null ],
          [ "Data streamer", "d3/d84/simcomps.html#simcomp_data_streamer", null ],
          [ "Field subsampler", "d3/d84/simcomps.html#simcomp_field_subsampler", null ]
        ] ]
      ] ],
      [ "Point zones", "da/dd0/point-zones.html", [
        [ "What are point zones?", "da/dd0/point-zones.html#autotoc_md177", null ],
        [ "Predefined geometrical shapes", "da/dd0/point-zones.html#autotoc_md178", [
          [ "Box", "da/dd0/point-zones.html#autotoc_md179", null ],
          [ "Sphere", "da/dd0/point-zones.html#autotoc_md180", null ],
          [ "Cylinder", "da/dd0/point-zones.html#autotoc_md181", null ]
        ] ],
        [ "Operations on point zones", "da/dd0/point-zones.html#autotoc_md182", [
          [ "Inversion", "da/dd0/point-zones.html#autotoc_md183", null ],
          [ "Combination", "da/dd0/point-zones.html#autotoc_md184", null ],
          [ "Including full element data", "da/dd0/point-zones.html#autotoc_md185", null ]
        ] ],
        [ "User-defined geometrical shapes", "da/dd0/point-zones.html#autotoc_md186", null ],
        [ "Using point zones", "da/dd0/point-zones.html#point-zones_using-point-zones", null ]
      ] ],
      [ "Statistics guide", "df/d8f/statistics-guide.html", [
        [ "Fluid Statistics", "df/d8f/statistics-guide.html#fluid-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md195", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md196", null ]
        ] ],
        [ "Postprocessing", "df/d8f/statistics-guide.html#autotoc_md197", null ],
        [ "Scalar Statistics", "df/d8f/statistics-guide.html#scalar-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md198", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md199", null ]
        ] ],
        [ "Fluid Subgrid-Scale (SGS) Statistics", "df/d8f/statistics-guide.html#fluid-sgs-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md200", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md201", null ]
        ] ],
        [ "Scalar Subgrid-Scale (SGS) Statistics", "df/d8f/statistics-guide.html#scalar-sgs-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md202", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md203", null ]
        ] ],
        [ "Note to users", "df/d8f/statistics-guide.html#note-to-users", null ]
      ] ],
      [ "Input-output", "d7/d7f/io.html", [
        [ "Mesh", "d7/d7f/io.html#autotoc_md162", null ],
        [ "Three-dimensional field output", "d7/d7f/io.html#autotoc_md163", [
          [ "Compression of field output", "d7/d7f/io.html#autotoc_md164", null ]
        ] ],
        [ "Checkpoint files", "d7/d7f/io.html#autotoc_md165", null ],
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
        [ "Installation", "d8/daa/api.html#autotoc_md73", [
          [ "Python (pyneko)", "d8/daa/api.html#autotoc_md74", null ],
          [ "Julia (Neko.jl)", "d8/daa/api.html#autotoc_md75", null ]
        ] ]
      ] ],
      [ "Performance guidelines", "dc/d3c/performance.html", [
        [ "Installation", "dc/d3c/performance.html#autotoc_md168", [
          [ "Accelerator specific options", "dc/d3c/performance.html#autotoc_md169", null ],
          [ "CPU specific options", "dc/d3c/performance.html#autotoc_md170", null ]
        ] ],
        [ "Simulation setup", "dc/d3c/performance.html#autotoc_md171", [
          [ "Load balancing", "dc/d3c/performance.html#autotoc_md172", null ],
          [ "Parameters", "dc/d3c/performance.html#autotoc_md173", null ]
        ] ],
        [ "Running a simulation", "dc/d3c/performance.html#autotoc_md174", [
          [ "Gather-scatter communication backends", "dc/d3c/performance.html#autotoc_md175", [
            [ "MPI neighbourhood collective backend", "dc/d3c/performance.html#performance-neighbour-backend", null ],
            [ "NVSHMEM backend", "dc/d3c/performance.html#performance-nvshmem-backend", null ],
            [ "OpenSHMEM backend", "dc/d3c/performance.html#performance-openshmem-backend", null ],
            [ "Coarray Fortran backend", "dc/d3c/performance.html#autotoc_md176", null ],
            [ "uTofu backend", "dc/d3c/performance.html#performance-utofu-backend", null ]
          ] ]
        ] ]
      ] ],
      [ "Global Interpolation", "dd/d61/global-interpolation.html", [
        [ "Overview", "dd/d61/global-interpolation.html#autotoc_md114", null ],
        [ "The global interpolation module", "dd/d61/global-interpolation.html#autotoc_md116", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md117", null ],
          [ "Features", "dd/d61/global-interpolation.html#autotoc_md118", null ]
        ] ],
        [ "The global interpolation type", "dd/d61/global-interpolation.html#autotoc_md120", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md121", null ],
          [ "Attributes", "dd/d61/global-interpolation.html#autotoc_md122", [
            [ "Domain and Space Information", "dd/d61/global-interpolation.html#autotoc_md123", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md124", null ],
            [ "Local Points", "dd/d61/global-interpolation.html#autotoc_md125", null ],
            [ "Interpolation Tools", "dd/d61/global-interpolation.html#autotoc_md126", null ],
            [ "Parallelism", "dd/d61/global-interpolation.html#autotoc_md127", null ],
            [ "Configuration", "dd/d61/global-interpolation.html#autotoc_md128", null ]
          ] ],
          [ "Methods", "dd/d61/global-interpolation.html#autotoc_md130", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md131", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md132", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md133", null ],
            [ "Validation", "dd/d61/global-interpolation.html#autotoc_md134", null ],
            [ "Memory Management", "dd/d61/global-interpolation.html#autotoc_md135", null ]
          ] ],
          [ "Example Usage", "dd/d61/global-interpolation.html#autotoc_md137", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md138", null ],
            [ "Finding Points", "dd/d61/global-interpolation.html#autotoc_md139", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md140", null ]
          ] ],
          [ "Notes", "dd/d61/global-interpolation.html#autotoc_md142", null ],
          [ "Environment variables", "dd/d61/global-interpolation.html#autotoc_md143", null ],
          [ "Related Modules", "dd/d61/global-interpolation.html#autotoc_md145", null ]
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
        [ "Git branches", "d1/d5a/contributing.html#autotoc_md44", null ],
        [ "Code style", "d1/d5a/contributing.html#autotoc_md45", [
          [ "Data types", "d1/d5a/contributing.html#autotoc_md46", null ]
        ] ],
        [ "Build system and code organization", "d1/d5a/contributing.html#autotoc_md47", [
          [ "Building CPU Fortran code", "d1/d5a/contributing.html#autotoc_md48", null ],
          [ "Device code", "d1/d5a/contributing.html#autotoc_md49", [
            [ "CUDA", "d1/d5a/contributing.html#autotoc_md50", null ],
            [ "HIP", "d1/d5a/contributing.html#autotoc_md51", null ],
            [ "OpenCL", "d1/d5a/contributing.html#autotoc_md52", null ],
            [ "Device-based type polymorphism", "d1/d5a/contributing.html#autotoc_md53", null ],
            [ "Summary of build system entires", "d1/d5a/contributing.html#autotoc_md54", null ]
          ] ]
        ] ]
      ] ],
      [ "Programming patterns and conventions", "d0/d47/dev_patterns.html", [
        [ "A. Naming", "d0/d47/dev_patterns.html#autotoc_md55", null ],
        [ "B. Scope", "d0/d47/dev_patterns.html#autotoc_md56", null ],
        [ "C. Constructors and destructors.", "d0/d47/dev_patterns.html#autotoc_md57", null ],
        [ "D. Documentation", "d0/d47/dev_patterns.html#autotoc_md58", null ],
        [ "E. Design", "d0/d47/dev_patterns.html#autotoc_md59", null ],
        [ "F. Complex and Polymorphic arrays.", "d0/d47/dev_patterns.html#autotoc_md60", [
          [ "F.1. Object wrapper pattern", "d0/d47/dev_patterns.html#autotoc_md61", null ],
          [ "F.2. Object pointer pattern", "d0/d47/dev_patterns.html#autotoc_md62", null ]
        ] ]
      ] ],
      [ "Code style", "da/db6/code-style.html", [
        [ "Data types", "da/db6/code-style.html#autotoc_md39", null ],
        [ "Enforcing rules", "da/db6/code-style.html#autotoc_md40", [
          [ "Tools", "da/db6/code-style.html#autotoc_md41", [
            [ "flint", "da/db6/code-style.html#autotoc_md42", null ],
            [ "findent", "da/db6/code-style.html#autotoc_md43", null ]
          ] ]
        ] ]
      ] ],
      [ "Testing", "d5/d75/testing.html", [
        [ "pFUnit", "d5/d75/testing.html#autotoc_md68", null ],
        [ "Configuring Neko", "d5/d75/testing.html#autotoc_md69", null ],
        [ "Running the tests", "d5/d75/testing.html#autotoc_md70", null ],
        [ "Adding a new test", "d5/d75/testing.html#autotoc_md71", null ]
      ] ],
      [ "Accelerators", "de/d06/accelerators.html", [
        [ "Device abstraction layer", "de/d06/accelerators.html#autotoc_md34", [
          [ "Memory management", "de/d06/accelerators.html#autotoc_md35", [
            [ "Allocation/deallocation", "de/d06/accelerators.html#autotoc_md36", null ],
            [ "Associate data on host and device", "de/d06/accelerators.html#autotoc_md37", null ],
            [ "Map a host array to a device", "de/d06/accelerators.html#autotoc_md38", null ],
            [ "Data transfer", "de/d06/accelerators.html#accelerators_data-transfer", null ]
          ] ],
          [ "Offload work", "de/d06/accelerators.html#accelerators_offload-work", null ]
        ] ]
      ] ],
      [ "Run-time selectable types", "d5/d5f/rts_types.html", null ],
      [ "Important types", "d3/d40/important_types.html", [
        [ "SEM foundation types", "d3/d40/important_types.html#autotoc_md63", null ],
        [ "Basic math routines", "d3/d40/important_types.html#autotoc_md64", null ],
        [ "Governing equation solvers and related types", "d3/d40/important_types.html#autotoc_md65", null ],
        [ "Singletons", "d3/d40/important_types.html#autotoc_md66", null ],
        [ "Linear algebra", "d3/d40/important_types.html#autotoc_md67", null ]
      ] ]
    ] ],
    [ "Appendices", "da/dd6/appendices.html", [
      [ "Environmental variable reference", "da/dd6/appendices.html#appendices_env-var", [
        [ "Logging level details", "da/dd6/appendices.html#autotoc_md0", null ],
        [ "Gather-scatter communication backend details", "da/dd6/appendices.html#autotoc_md1", null ],
        [ "MPI thread level details", "da/dd6/appendices.html#autotoc_md2", null ],
        [ "Coarray Fortran signaling mode details", "da/dd6/appendices.html#autotoc_md3", null ],
        [ "uTofu injection details", "da/dd6/appendices.html#autotoc_md4", null ]
      ] ],
      [ "Governing equations", "db/d27/governing-equations.html", [
        [ "Momentum (Fluid)", "db/d27/governing-equations.html#autotoc_md18", null ],
        [ "Scalar", "db/d27/governing-equations.html#autotoc_md19", null ],
        [ "Compressible Flows", "db/d27/governing-equations.html#autotoc_md20", null ],
        [ "Non-dimensionalisation", "db/d27/governing-equations.html#autotoc_md21", null ]
      ] ],
      [ "Neko .nmsh mesh format", "d0/d5b/nmsh-format.html", [
        [ "File overview", "d0/d5b/nmsh-format.html#autotoc_md22", null ],
        [ "Data types and portability", "d0/d5b/nmsh-format.html#autotoc_md23", null ],
        [ "Element records", "d0/d5b/nmsh-format.html#autotoc_md24", [
          [ "Vertex record (<tt>nmsh_vertex_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md25", null ],
          [ "Quad record (<tt>nmsh_quad_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md26", null ],
          [ "Hex record (<tt>nmsh_hex_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md27", null ],
          [ "Vertex ordering", "d0/d5b/nmsh-format.html#autotoc_md28", null ]
        ] ],
        [ "Zone records (boundary and periodic data)", "d0/d5b/nmsh-format.html#autotoc_md29", null ],
        [ "Curve records (curved edges/faces)", "d0/d5b/nmsh-format.html#autotoc_md30", null ],
        [ "2D meshes", "d0/d5b/nmsh-format.html#autotoc_md31", null ]
      ] ],
      [ "Neko .fld field format", "d8/d00/fld-format.html", [
        [ "File naming and series", "d8/d00/fld-format.html#autotoc_md5", null ],
        [ "File overview", "d8/d00/fld-format.html#autotoc_md6", null ],
        [ "Header", "d8/d00/fld-format.html#autotoc_md7", [
          [ "rdcode", "d8/d00/fld-format.html#autotoc_md8", null ]
        ] ],
        [ "Test pattern", "d8/d00/fld-format.html#autotoc_md9", null ],
        [ "Element index list", "d8/d00/fld-format.html#autotoc_md10", null ],
        [ "Field blocks", "d8/d00/fld-format.html#autotoc_md11", [
          [ "Vector fields (X and U)", "d8/d00/fld-format.html#autotoc_md12", null ],
          [ "Scalar fields (P, T, S)", "d8/d00/fld-format.html#autotoc_md13", null ],
          [ "Precision", "d8/d00/fld-format.html#autotoc_md14", null ]
        ] ],
        [ "Metadata blocks (3D only)", "d8/d00/fld-format.html#autotoc_md15", null ],
        [ "2D fields and Z-velocity", "d8/d00/fld-format.html#autotoc_md16", null ],
        [ "Portability notes", "d8/d00/fld-format.html#autotoc_md17", null ]
      ] ],
      [ "Publications", "de/d26/publications.html", [
        [ "Computer science", "de/d26/publications.html#autotoc_md32", null ],
        [ "Flow physics", "de/d26/publications.html#autotoc_md33", null ]
      ] ]
    ] ]
];