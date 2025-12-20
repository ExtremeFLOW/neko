var index =
[
    [ "Structure of the Manual", "index.html#autotoc_md41", null ],
    [ "User guide", "dd/d04/user-guide.html", [
      [ "Installing Neko", "d5/dfc/installation.html", [
        [ "Building from source", "d5/dfc/installation.html#autotoc_md111", [
          [ "Dependencies", "d5/dfc/installation.html#autotoc_md112", [
            [ "Building JSON Fortran", "d5/dfc/installation.html#autotoc_md113", null ],
            [ "Building HDF5 (optional, but highly recommended)", "d5/dfc/installation.html#autotoc_md114", null ],
            [ "Building ParMETIS (optional)", "d5/dfc/installation.html#autotoc_md115", null ],
            [ "Bulding PFunit (optional)", "d5/dfc/installation.html#autotoc_md116", null ]
          ] ],
          [ "Building Neko", "d5/dfc/installation.html#autotoc_md117", [
            [ "Compiling Neko for CPU or SX-Aurora", "d5/dfc/installation.html#autotoc_md118", null ],
            [ "Compiling Neko for NVIDIA GPUs", "d5/dfc/installation.html#autotoc_md119", null ],
            [ "Compiling Neko for AMD GPUs", "d5/dfc/installation.html#autotoc_md120", null ],
            [ "Compiling Neko with a collective communications library", "d5/dfc/installation.html#autotoc_md121", null ]
          ] ]
        ] ],
        [ "Installing via Spack", "d5/dfc/installation.html#autotoc_md122", [
          [ "Quick start guide with Spack", "d5/dfc/installation.html#autotoc_md123", null ]
        ] ],
        [ "Installing using pixi", "d5/dfc/installation.html#autotoc_md124", null ],
        [ "Using a Docker container", "d5/dfc/installation.html#autotoc_md125", null ],
        [ "Testing", "d5/d75/testing.html", [
          [ "pFUnit", "d5/d75/testing.html#autotoc_md37", null ],
          [ "Configuring Neko", "d5/d75/testing.html#autotoc_md38", null ],
          [ "Running the tests", "d5/d75/testing.html#autotoc_md39", null ],
          [ "Adding a new test", "d5/d75/testing.html#autotoc_md40", null ]
        ] ]
      ] ],
      [ "Meshing", "d9/df2/meshing.html", [
        [ "General considerations", "d9/df2/meshing.html#general-considerations", null ],
        [ "Constructing meshes", "d9/df2/meshing.html#autotoc_md130", null ]
      ] ],
      [ "Case File", "dd/d33/case-file.html", [
        [ "High-level structure", "dd/d33/case-file.html#autotoc_md45", null ],
        [ "Output frequency control", "dd/d33/case-file.html#autotoc_md46", null ],
        [ "The case object", "dd/d33/case-file.html#autotoc_md47", [
          [ "Time control", "dd/d33/case-file.html#autotoc_md48", null ],
          [ "Restarts and joblimit", "dd/d33/case-file.html#autotoc_md49", null ],
          [ "Boundary type numbering in the \"output_boundary\" field", "dd/d33/case-file.html#autotoc_md50", null ]
        ] ],
        [ "Numerics", "dd/d33/case-file.html#autotoc_md51", null ],
        [ "Fluid", "dd/d33/case-file.html#autotoc_md52", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md53", null ],
          [ "Compressible flows", "dd/d33/case-file.html#autotoc_md54", [
            [ "Compressible boundary conditions", "dd/d33/case-file.html#autotoc_md55", null ]
          ] ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md56", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#case-file_fluid-boundary-conditions", [
            [ "Specifying the boundaries", "dd/d33/case-file.html#autotoc_md57", null ],
            [ "Available conditions", "dd/d33/case-file.html#autotoc_md58", null ]
          ] ],
          [ "Initial conditions", "dd/d33/case-file.html#case-file_fluid-ic", null ],
          [ "Source terms", "dd/d33/case-file.html#case-file_fluid-source-term", [
            [ "Brinkman", "dd/d33/case-file.html#autotoc_md59", null ],
            [ "Gradient Jump Penalty", "dd/d33/case-file.html#autotoc_md60", null ],
            [ "Sponge", "dd/d33/case-file.html#autotoc_md61", null ]
          ] ]
        ] ],
        [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md62", [
          [ "Multilevel preconditioners", "dd/d33/case-file.html#autotoc_md63", null ],
          [ "Flow rate forcing", "dd/d33/case-file.html#autotoc_md64", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md65", null ]
        ] ],
        [ "Scalar", "dd/d33/case-file.html#case-file_scalar", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md66", null ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md67", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#autotoc_md68", null ],
          [ "Initial conditions", "dd/d33/case-file.html#autotoc_md69", null ],
          [ "Source terms", "dd/d33/case-file.html#autotoc_md70", null ],
          [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md71", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md72", null ]
        ] ],
        [ "Simulation components", "dd/d33/case-file.html#autotoc_md73", null ],
        [ "Point zones", "dd/d33/case-file.html#autotoc_md74", null ],
        [ "Runtime statistics", "dd/d33/case-file.html#autotoc_md75", null ]
      ] ],
      [ "User File", "d6/def/user-file.html", [
        [ "Compiling and running", "d6/def/user-file.html#autotoc_md160", null ],
        [ "High-level structure", "d6/def/user-file.html#autotoc_md161", null ],
        [ "Default user functions", "d6/def/user-file.html#autotoc_md162", [
          [ "Initializing and finalizing", "d6/def/user-file.html#user-file_init-and-final", null ],
          [ "Computing at every time step", "d6/def/user-file.html#user-file_user-check", null ],
          [ "Setting material properties", "d6/def/user-file.html#user-file_mat-prop", null ],
          [ "Runtime mesh deformation", "d6/def/user-file.html#user-file_user-mesh-setup", null ]
        ] ],
        [ "Case-specific user functions", "d6/def/user-file.html#autotoc_md163", [
          [ "Fluid and Scalar initial conditions", "d6/def/user-file.html#user-file_user-ic", null ],
          [ "Fluid and scalar source terms", "d6/def/user-file.html#user-file_user-f", null ],
          [ "Dirichlet boundary conditions", "d6/def/user-file.html#user-file_field-dirichlet-update", null ]
        ] ],
        [ "Additional remarks and tips", "d6/def/user-file.html#autotoc_md164", [
          [ "Running on GPUs", "d6/def/user-file.html#user-file_tips_running-on-gpus", [
            [ "Custom GPU kernels", "d6/def/user-file.html#user-file_tips_running-on-gpus-custom-kernels", null ]
          ] ],
          [ "Registries", "d6/def/user-file.html#user-file_tips_registries", null ],
          [ "User access to solver internals", "d6/def/user-file.html#user-file_access", null ]
        ] ]
      ] ],
      [ "Simulation components", "d3/d84/simcomps.html", [
        [ "What are simulation components?", "d3/d84/simcomps.html#autotoc_md146", null ],
        [ "Adding simulation components to the case", "d3/d84/simcomps.html#autotoc_md147", null ],
        [ "List of simulation components", "d3/d84/simcomps.html#autotoc_md148", null ],
        [ "Controling execution and file output", "d3/d84/simcomps.html#autotoc_md149", [
          [ "Differential operators", "d3/d84/simcomps.html#autotoc_md150", [
            [ "derivative", "d3/d84/simcomps.html#simcomp_derivative", null ],
            [ "curl", "d3/d84/simcomps.html#simcomp_curl", null ],
            [ "divergence", "d3/d84/simcomps.html#simcomp_divergence", null ]
          ] ],
          [ "gradient", "d3/d84/simcomps.html#simcomp_gradient", null ],
          [ "weak_gradient", "d3/d84/simcomps.html#simcomp_weak_gradient", null ],
          [ "lambda2", "d3/d84/simcomps.html#simcomp_lambda2", null ],
          [ "probes", "d3/d84/simcomps.html#simcomp_probes", [
            [ "Supported types", "d3/d84/simcomps.html#autotoc_md151", null ],
            [ "Example usage", "d3/d84/simcomps.html#autotoc_md152", null ]
          ] ],
          [ "field_writer", "d3/d84/simcomps.html#simcomp_field_writer", null ],
          [ "force_torque", "d3/d84/simcomps.html#simcomp_force_torque", null ],
          [ "les_model", "d3/d84/simcomps.html#simcomp_les_model", null ],
          [ "User statistics", "d3/d84/simcomps.html#user_stats", null ],
          [ "Spectral error indicator", "d3/d84/simcomps.html#simcomp_speri", null ]
        ] ]
      ] ],
      [ "Point zones", "da/dd0/point-zones.html", [
        [ "What are point zones?", "da/dd0/point-zones.html#autotoc_md137", null ],
        [ "Predefined geometrical shapes", "da/dd0/point-zones.html#autotoc_md138", [
          [ "Box", "da/dd0/point-zones.html#autotoc_md139", null ],
          [ "Sphere", "da/dd0/point-zones.html#autotoc_md140", null ],
          [ "Cylinder", "da/dd0/point-zones.html#autotoc_md141", null ]
        ] ],
        [ "Operations on point zones", "da/dd0/point-zones.html#autotoc_md142", [
          [ "Inversion", "da/dd0/point-zones.html#autotoc_md143", null ],
          [ "Combination", "da/dd0/point-zones.html#autotoc_md144", null ]
        ] ],
        [ "User-defined geometrical shapes", "da/dd0/point-zones.html#autotoc_md145", null ],
        [ "Using point zones", "da/dd0/point-zones.html#point-zones_using-point-zones", null ]
      ] ],
      [ "Statistics guide", "df/d8f/statistics-guide.html", [
        [ "Fluid Statistics", "df/d8f/statistics-guide.html#autotoc_md153", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md154", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md155", null ]
        ] ],
        [ "Postprocessing", "df/d8f/statistics-guide.html#autotoc_md156", null ],
        [ "Scalar Statistics", "df/d8f/statistics-guide.html#autotoc_md157", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md158", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md159", null ]
        ] ]
      ] ],
      [ "Global Interpolation", "dd/d61/global-interpolation.html", [
        [ "Overview", "dd/d61/global-interpolation.html#autotoc_md79", null ],
        [ "Module: <tt>global_interpolation</tt>", "dd/d61/global-interpolation.html#autotoc_md81", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md82", null ],
          [ "Features", "dd/d61/global-interpolation.html#autotoc_md83", null ]
        ] ],
        [ "Type: <tt>global_interpolation_t</tt>", "dd/d61/global-interpolation.html#autotoc_md85", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md86", null ],
          [ "Attributes", "dd/d61/global-interpolation.html#autotoc_md87", [
            [ "Domain and Space Information", "dd/d61/global-interpolation.html#autotoc_md88", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md89", null ],
            [ "Local Points", "dd/d61/global-interpolation.html#autotoc_md90", null ],
            [ "Interpolation Tools", "dd/d61/global-interpolation.html#autotoc_md91", null ],
            [ "Parallelism", "dd/d61/global-interpolation.html#autotoc_md92", null ],
            [ "Configuration", "dd/d61/global-interpolation.html#autotoc_md93", null ]
          ] ],
          [ "Methods", "dd/d61/global-interpolation.html#autotoc_md95", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md96", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md97", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md98", null ],
            [ "Validation", "dd/d61/global-interpolation.html#autotoc_md99", null ],
            [ "Memory Management", "dd/d61/global-interpolation.html#autotoc_md100", null ]
          ] ],
          [ "Example Usage", "dd/d61/global-interpolation.html#autotoc_md102", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md103", null ],
            [ "Finding Points", "dd/d61/global-interpolation.html#autotoc_md104", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md105", null ]
          ] ],
          [ "Notes", "dd/d61/global-interpolation.html#autotoc_md107", null ],
          [ "Environment variables", "dd/d61/global-interpolation.html#autotoc_md108", null ],
          [ "Related Modules", "dd/d61/global-interpolation.html#autotoc_md110", null ]
        ] ]
      ] ],
      [ "Filtering", "df/d4a/filter.html", [
        [ "Filter base type", "df/d4a/filter.html#autotoc_md76", null ],
        [ "PDE filter", "df/d4a/filter.html#autotoc_md77", null ],
        [ "Elementwise filter", "df/d4a/filter.html#autotoc_md78", null ]
      ] ],
      [ "Input-output", "d7/d7f/io.html", [
        [ "Mesh", "d7/d7f/io.html#autotoc_md126", null ],
        [ "Three-dimensional field output", "d7/d7f/io.html#autotoc_md127", [
          [ "Compression of field output", "d7/d7f/io.html#autotoc_md128", null ]
        ] ],
        [ "Checkpoint files", "d7/d7f/io.html#autotoc_md129", null ]
      ] ],
      [ "Examples: Programming the user file", "d5/db5/programming-examples.html", null ],
      [ "Extending neko", "d4/d1b/extending.html", null ],
      [ "Neko API", "d8/daa/api.html", [
        [ "Installation", "d8/daa/api.html#autotoc_md42", [
          [ "Python (pyneko)", "d8/daa/api.html#autotoc_md43", null ],
          [ "Julia (Neko.jl)", "d8/daa/api.html#autotoc_md44", null ]
        ] ]
      ] ],
      [ "Performance guidelines", "dc/d3c/performance.html", [
        [ "Installation", "dc/d3c/performance.html#autotoc_md131", [
          [ "Accelerator specific options", "dc/d3c/performance.html#autotoc_md132", null ]
        ] ],
        [ "Simulation setup", "dc/d3c/performance.html#autotoc_md133", [
          [ "Load balancing", "dc/d3c/performance.html#autotoc_md134", null ],
          [ "Parameters", "dc/d3c/performance.html#autotoc_md135", null ]
        ] ],
        [ "Running a simulation", "dc/d3c/performance.html#autotoc_md136", null ]
      ] ]
    ] ],
    [ "Developer guide", "dc/d70/developer-guide.html", [
      [ "Contributing to Neko", "d1/d5a/contributing.html", [
        [ "Git branches", "d1/d5a/contributing.html#autotoc_md16", null ],
        [ "Code style", "d1/d5a/contributing.html#autotoc_md17", [
          [ "Data types", "d1/d5a/contributing.html#autotoc_md18", null ]
        ] ],
        [ "Build system and code organization", "d1/d5a/contributing.html#autotoc_md19", [
          [ "Building CPU Fortran code", "d1/d5a/contributing.html#autotoc_md20", null ],
          [ "Device code", "d1/d5a/contributing.html#autotoc_md21", [
            [ "CUDA", "d1/d5a/contributing.html#autotoc_md22", null ],
            [ "HIP", "d1/d5a/contributing.html#autotoc_md23", null ],
            [ "OpenCL", "d1/d5a/contributing.html#autotoc_md24", null ],
            [ "Device-based type polymorphism", "d1/d5a/contributing.html#autotoc_md25", null ],
            [ "Summary of build system entires", "d1/d5a/contributing.html#autotoc_md26", null ]
          ] ]
        ] ]
      ] ],
      [ "Programming patterns and conventions", "d0/d47/dev_patterns.html", [
        [ "A. Naming", "d0/d47/dev_patterns.html#autotoc_md27", null ],
        [ "B. Scope", "d0/d47/dev_patterns.html#autotoc_md28", null ],
        [ "C. Constructors and destructors.", "d0/d47/dev_patterns.html#autotoc_md29", null ],
        [ "D. Documentation", "d0/d47/dev_patterns.html#autotoc_md30", null ],
        [ "E. Design", "d0/d47/dev_patterns.html#autotoc_md31", null ]
      ] ],
      [ "Code style", "da/db6/code-style.html", [
        [ "Data types", "da/db6/code-style.html#autotoc_md11", null ],
        [ "Enforcing rules", "da/db6/code-style.html#autotoc_md12", [
          [ "Tools", "da/db6/code-style.html#autotoc_md13", [
            [ "flint", "da/db6/code-style.html#autotoc_md14", null ],
            [ "findent", "da/db6/code-style.html#autotoc_md15", null ]
          ] ]
        ] ]
      ] ],
      [ "Testing", "d5/d75/testing.html", [
        [ "pFUnit", "d5/d75/testing.html#autotoc_md37", null ],
        [ "Configuring Neko", "d5/d75/testing.html#autotoc_md38", null ],
        [ "Running the tests", "d5/d75/testing.html#autotoc_md39", null ],
        [ "Adding a new test", "d5/d75/testing.html#autotoc_md40", null ]
      ] ],
      [ "Accelerators", "de/d06/accelerators.html", [
        [ "Device abstraction layer", "de/d06/accelerators.html#autotoc_md6", [
          [ "Memory management", "de/d06/accelerators.html#autotoc_md7", [
            [ "Allocation/deallocation", "de/d06/accelerators.html#autotoc_md8", null ],
            [ "Associate data on host and device", "de/d06/accelerators.html#autotoc_md9", null ],
            [ "Map a host array to a device", "de/d06/accelerators.html#autotoc_md10", null ],
            [ "Data transfer", "de/d06/accelerators.html#accelerators_data-transfer", null ]
          ] ],
          [ "Offload work", "de/d06/accelerators.html#accelerators_offload-work", null ]
        ] ]
      ] ],
      [ "Run-time selectable types", "d5/d5f/rts_types.html", null ],
      [ "Important types", "d3/d40/important_types.html", [
        [ "SEM foundation types", "d3/d40/important_types.html#autotoc_md32", null ],
        [ "Basic math routines", "d3/d40/important_types.html#autotoc_md33", null ],
        [ "Governing equation solvers and related types", "d3/d40/important_types.html#autotoc_md34", null ],
        [ "Singletons", "d3/d40/important_types.html#autotoc_md35", null ],
        [ "Linear algebra", "d3/d40/important_types.html#autotoc_md36", null ]
      ] ]
    ] ],
    [ "Appendices", "da/dd6/appendices.html", [
      [ "Environmental variable reference", "da/dd6/appendices.html#appendices_env-var", [
        [ "Logging level details", "da/dd6/appendices.html#autotoc_md0", null ],
        [ "Gather-scatter communication backend details", "da/dd6/appendices.html#autotoc_md1", null ]
      ] ],
      [ "Governing equations", "db/d27/governing-equations.html", [
        [ "Momentum (Fluid)", "db/d27/governing-equations.html#autotoc_md2", null ],
        [ "Scalar", "db/d27/governing-equations.html#autotoc_md3", null ],
        [ "Compressible Flows", "db/d27/governing-equations.html#autotoc_md4", null ],
        [ "Non-dimensionalisation", "db/d27/governing-equations.html#autotoc_md5", null ]
      ] ],
      [ "Publications", "de/d26/publications.html", null ]
    ] ]
];