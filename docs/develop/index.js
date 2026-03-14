var index =
[
    [ "Structure of the Manual", "index.html#autotoc_md67", null ],
    [ "User guide", "dd/d04/user-guide.html", [
      [ "Installing Neko", "d5/dfc/installation.html", [
        [ "Building from source", "d5/dfc/installation.html#autotoc_md137", [
          [ "Dependencies", "d5/dfc/installation.html#autotoc_md138", [
            [ "Building JSON Fortran", "d5/dfc/installation.html#autotoc_md139", null ],
            [ "Building HDF5 (optional, but highly recommended)", "d5/dfc/installation.html#autotoc_md140", null ],
            [ "Building ParMETIS (optional)", "d5/dfc/installation.html#autotoc_md141", null ],
            [ "Bulding PFunit (optional)", "d5/dfc/installation.html#autotoc_md142", null ],
            [ "Python (optional)", "d5/dfc/installation.html#deps-python", null ]
          ] ],
          [ "Building Neko", "d5/dfc/installation.html#building-neko", [
            [ "Compiling Neko for CPU or SX-Aurora", "d5/dfc/installation.html#autotoc_md143", null ],
            [ "Compiling Neko for NVIDIA GPUs", "d5/dfc/installation.html#autotoc_md144", null ],
            [ "Compiling Neko for AMD GPUs", "d5/dfc/installation.html#autotoc_md145", null ],
            [ "Compiling Neko with a collective communications library", "d5/dfc/installation.html#autotoc_md146", null ]
          ] ]
        ] ],
        [ "Installing via Spack", "d5/dfc/installation.html#autotoc_md147", [
          [ "Quick start guide with Spack", "d5/dfc/installation.html#autotoc_md148", null ]
        ] ],
        [ "Installing using pixi", "d5/dfc/installation.html#autotoc_md149", null ],
        [ "Using a Docker container", "d5/dfc/installation.html#autotoc_md150", null ],
        [ "Testing", "d5/d75/testing.html", [
          [ "pFUnit", "d5/d75/testing.html#autotoc_md63", null ],
          [ "Configuring Neko", "d5/d75/testing.html#autotoc_md64", null ],
          [ "Running the tests", "d5/d75/testing.html#autotoc_md65", null ],
          [ "Adding a new test", "d5/d75/testing.html#autotoc_md66", null ]
        ] ]
      ] ],
      [ "Meshing", "d9/df2/meshing.html", [
        [ "General considerations", "d9/df2/meshing.html#general-considerations", null ],
        [ "Constructing meshes", "d9/df2/meshing.html#autotoc_md155", null ]
      ] ],
      [ "Case File", "dd/d33/case-file.html", [
        [ "High-level structure", "dd/d33/case-file.html#autotoc_md71", null ],
        [ "Output frequency control", "dd/d33/case-file.html#autotoc_md72", null ],
        [ "The case object", "dd/d33/case-file.html#autotoc_md73", [
          [ "Constants", "dd/d33/case-file.html#autotoc_md74", null ],
          [ "Time control", "dd/d33/case-file.html#autotoc_md75", null ],
          [ "Restarts and joblimit", "dd/d33/case-file.html#autotoc_md76", null ],
          [ "Boundary type numbering in the \"output_boundary\" field", "dd/d33/case-file.html#autotoc_md77", null ]
        ] ],
        [ "Numerics", "dd/d33/case-file.html#autotoc_md78", null ],
        [ "Fluid", "dd/d33/case-file.html#case-file_fluid", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md79", null ],
          [ "Compressible flows", "dd/d33/case-file.html#autotoc_md80", [
            [ "Compressible boundary conditions", "dd/d33/case-file.html#autotoc_md81", null ]
          ] ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md82", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#case-file_fluid-boundary-conditions", [
            [ "Specifying the boundaries", "dd/d33/case-file.html#autotoc_md83", null ],
            [ "Available conditions", "dd/d33/case-file.html#autotoc_md84", null ]
          ] ],
          [ "Initial conditions", "dd/d33/case-file.html#case-file_fluid-ic", null ],
          [ "Source terms", "dd/d33/case-file.html#case-file_fluid-source-term", [
            [ "Brinkman", "dd/d33/case-file.html#autotoc_md85", null ],
            [ "Gradient Jump Penalty", "dd/d33/case-file.html#autotoc_md86", null ],
            [ "Sponge", "dd/d33/case-file.html#autotoc_md87", null ]
          ] ]
        ] ],
        [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md88", [
          [ "Multilevel preconditioners", "dd/d33/case-file.html#autotoc_md89", null ],
          [ "Flow rate forcing", "dd/d33/case-file.html#autotoc_md90", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md91", null ]
        ] ],
        [ "Scalar", "dd/d33/case-file.html#case-file_scalar", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md92", null ],
          [ "Turbulence modelling", "dd/d33/case-file.html#autotoc_md93", null ],
          [ "Boundary conditions", "dd/d33/case-file.html#autotoc_md94", null ],
          [ "Initial conditions", "dd/d33/case-file.html#autotoc_md95", null ],
          [ "Source terms", "dd/d33/case-file.html#autotoc_md96", null ],
          [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md97", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md98", null ]
        ] ],
        [ "Simulation components", "dd/d33/case-file.html#autotoc_md99", null ],
        [ "Point zones", "dd/d33/case-file.html#autotoc_md100", null ],
        [ "Runtime statistics", "dd/d33/case-file.html#autotoc_md101", null ]
      ] ],
      [ "User File", "d6/def/user-file.html", [
        [ "Compiling and running", "d6/def/user-file.html#autotoc_md188", null ],
        [ "High-level structure", "d6/def/user-file.html#autotoc_md189", null ],
        [ "Default user functions", "d6/def/user-file.html#autotoc_md190", [
          [ "Initializing and finalizing", "d6/def/user-file.html#user-file_init-and-final", null ],
          [ "Computing at every time step", "d6/def/user-file.html#user-file_user-check", null ],
          [ "Setting material properties", "d6/def/user-file.html#user-file_mat-prop", null ],
          [ "Runtime mesh deformation", "d6/def/user-file.html#user-file_user-mesh-setup", null ]
        ] ],
        [ "Case-specific user functions", "d6/def/user-file.html#autotoc_md191", [
          [ "Fluid and Scalar initial conditions", "d6/def/user-file.html#user-file_user-ic", null ],
          [ "Fluid and scalar source terms", "d6/def/user-file.html#user-file_user-f", null ],
          [ "Dirichlet boundary conditions", "d6/def/user-file.html#user-file_field-dirichlet-update", null ]
        ] ],
        [ "Additional remarks and tips", "d6/def/user-file.html#autotoc_md192", [
          [ "Running on GPUs", "d6/def/user-file.html#user-file_tips_running-on-gpus", [
            [ "Custom GPU kernels", "d6/def/user-file.html#user-file_tips_running-on-gpus-custom-kernels", null ]
          ] ],
          [ "Registries", "d6/def/user-file.html#user-file_tips_registries", null ],
          [ "User access to solver internals", "d6/def/user-file.html#user-file_access", null ]
        ] ]
      ] ],
      [ "Simulation components", "d3/d84/simcomps.html", [
        [ "What are simulation components?", "d3/d84/simcomps.html#autotoc_md172", null ],
        [ "Adding simulation components to the case", "d3/d84/simcomps.html#autotoc_md173", null ],
        [ "List of simulation components", "d3/d84/simcomps.html#autotoc_md174", null ],
        [ "Controling execution and file output", "d3/d84/simcomps.html#autotoc_md175", [
          [ "Differential operators", "d3/d84/simcomps.html#autotoc_md176", [
            [ "derivative", "d3/d84/simcomps.html#simcomp_derivative", null ],
            [ "curl", "d3/d84/simcomps.html#simcomp_curl", null ],
            [ "divergence", "d3/d84/simcomps.html#simcomp_divergence", null ]
          ] ],
          [ "gradient", "d3/d84/simcomps.html#simcomp_gradient", null ],
          [ "weak_gradient", "d3/d84/simcomps.html#simcomp_weak_gradient", null ],
          [ "lambda2", "d3/d84/simcomps.html#simcomp_lambda2", null ],
          [ "probes", "d3/d84/simcomps.html#simcomp_probes", [
            [ "Supported types", "d3/d84/simcomps.html#autotoc_md177", null ],
            [ "Example usage", "d3/d84/simcomps.html#autotoc_md178", null ]
          ] ],
          [ "field_writer", "d3/d84/simcomps.html#simcomp_field_writer", null ],
          [ "force_torque", "d3/d84/simcomps.html#simcomp_force_torque", null ],
          [ "les_model", "d3/d84/simcomps.html#simcomp_les_model", null ],
          [ "User statistics", "d3/d84/simcomps.html#user_stats", null ],
          [ "Spectral error indicator", "d3/d84/simcomps.html#simcomp_speri", null ]
        ] ]
      ] ],
      [ "Point zones", "da/dd0/point-zones.html", [
        [ "What are point zones?", "da/dd0/point-zones.html#autotoc_md162", null ],
        [ "Predefined geometrical shapes", "da/dd0/point-zones.html#autotoc_md163", [
          [ "Box", "da/dd0/point-zones.html#autotoc_md164", null ],
          [ "Sphere", "da/dd0/point-zones.html#autotoc_md165", null ],
          [ "Cylinder", "da/dd0/point-zones.html#autotoc_md166", null ]
        ] ],
        [ "Operations on point zones", "da/dd0/point-zones.html#autotoc_md167", [
          [ "Inversion", "da/dd0/point-zones.html#autotoc_md168", null ],
          [ "Combination", "da/dd0/point-zones.html#autotoc_md169", null ],
          [ "Including full element data", "da/dd0/point-zones.html#autotoc_md170", null ]
        ] ],
        [ "User-defined geometrical shapes", "da/dd0/point-zones.html#autotoc_md171", null ],
        [ "Using point zones", "da/dd0/point-zones.html#point-zones_using-point-zones", null ]
      ] ],
      [ "Statistics guide", "df/d8f/statistics-guide.html", [
        [ "Fluid Statistics", "df/d8f/statistics-guide.html#fluid-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md179", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md180", null ]
        ] ],
        [ "Postprocessing", "df/d8f/statistics-guide.html#autotoc_md181", null ],
        [ "Scalar Statistics", "df/d8f/statistics-guide.html#scalar-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md182", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md183", null ]
        ] ],
        [ "Fluid Subgrid-Scale (SGS) Statistics", "df/d8f/statistics-guide.html#fluid-sgs-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md184", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md185", null ]
        ] ],
        [ "Scalar Subgrid-Scale (SGS) Statistics", "df/d8f/statistics-guide.html#scalar-sgs-statistics", [
          [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md186", null ],
          [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md187", null ]
        ] ]
      ] ],
      [ "Global Interpolation", "dd/d61/global-interpolation.html", [
        [ "Overview", "dd/d61/global-interpolation.html#autotoc_md105", null ],
        [ "Module: <tt>global_interpolation</tt>", "dd/d61/global-interpolation.html#autotoc_md107", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md108", null ],
          [ "Features", "dd/d61/global-interpolation.html#autotoc_md109", null ]
        ] ],
        [ "Type: <tt>global_interpolation_t</tt>", "dd/d61/global-interpolation.html#autotoc_md111", [
          [ "Description", "dd/d61/global-interpolation.html#autotoc_md112", null ],
          [ "Attributes", "dd/d61/global-interpolation.html#autotoc_md113", [
            [ "Domain and Space Information", "dd/d61/global-interpolation.html#autotoc_md114", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md115", null ],
            [ "Local Points", "dd/d61/global-interpolation.html#autotoc_md116", null ],
            [ "Interpolation Tools", "dd/d61/global-interpolation.html#autotoc_md117", null ],
            [ "Parallelism", "dd/d61/global-interpolation.html#autotoc_md118", null ],
            [ "Configuration", "dd/d61/global-interpolation.html#autotoc_md119", null ]
          ] ],
          [ "Methods", "dd/d61/global-interpolation.html#autotoc_md121", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md122", null ],
            [ "Point Management", "dd/d61/global-interpolation.html#autotoc_md123", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md124", null ],
            [ "Validation", "dd/d61/global-interpolation.html#autotoc_md125", null ],
            [ "Memory Management", "dd/d61/global-interpolation.html#autotoc_md126", null ]
          ] ],
          [ "Example Usage", "dd/d61/global-interpolation.html#autotoc_md128", [
            [ "Initialization", "dd/d61/global-interpolation.html#autotoc_md129", null ],
            [ "Finding Points", "dd/d61/global-interpolation.html#autotoc_md130", null ],
            [ "Interpolation", "dd/d61/global-interpolation.html#autotoc_md131", null ]
          ] ],
          [ "Notes", "dd/d61/global-interpolation.html#autotoc_md133", null ],
          [ "Environment variables", "dd/d61/global-interpolation.html#autotoc_md134", null ],
          [ "Related Modules", "dd/d61/global-interpolation.html#autotoc_md136", null ]
        ] ]
      ] ],
      [ "Filtering", "df/d4a/filter.html", [
        [ "Filter base type", "df/d4a/filter.html#autotoc_md102", null ],
        [ "PDE filter", "df/d4a/filter.html#autotoc_md103", null ],
        [ "Elementwise filter", "df/d4a/filter.html#autotoc_md104", null ]
      ] ],
      [ "Input-output", "d7/d7f/io.html", [
        [ "Mesh", "d7/d7f/io.html#autotoc_md151", null ],
        [ "Three-dimensional field output", "d7/d7f/io.html#autotoc_md152", [
          [ "Compression of field output", "d7/d7f/io.html#autotoc_md153", null ]
        ] ],
        [ "Checkpoint files", "d7/d7f/io.html#autotoc_md154", null ]
      ] ],
      [ "Examples: Programming the user file", "d5/db5/programming-examples.html", null ],
      [ "Extending neko", "d4/d1b/extending.html", null ],
      [ "Neko API", "d8/daa/api.html", [
        [ "Installation", "d8/daa/api.html#autotoc_md68", [
          [ "Python (pyneko)", "d8/daa/api.html#autotoc_md69", null ],
          [ "Julia (Neko.jl)", "d8/daa/api.html#autotoc_md70", null ]
        ] ]
      ] ],
      [ "Performance guidelines", "dc/d3c/performance.html", [
        [ "Installation", "dc/d3c/performance.html#autotoc_md156", [
          [ "Accelerator specific options", "dc/d3c/performance.html#autotoc_md157", null ]
        ] ],
        [ "Simulation setup", "dc/d3c/performance.html#autotoc_md158", [
          [ "Load balancing", "dc/d3c/performance.html#autotoc_md159", null ],
          [ "Parameters", "dc/d3c/performance.html#autotoc_md160", null ]
        ] ],
        [ "Running a simulation", "dc/d3c/performance.html#autotoc_md161", null ]
      ] ]
    ] ],
    [ "Developer guide", "dc/d70/developer-guide.html", [
      [ "Contributing to Neko", "d1/d5a/contributing.html", [
        [ "Git branches", "d1/d5a/contributing.html#autotoc_md39", null ],
        [ "Code style", "d1/d5a/contributing.html#autotoc_md40", [
          [ "Data types", "d1/d5a/contributing.html#autotoc_md41", null ]
        ] ],
        [ "Build system and code organization", "d1/d5a/contributing.html#autotoc_md42", [
          [ "Building CPU Fortran code", "d1/d5a/contributing.html#autotoc_md43", null ],
          [ "Device code", "d1/d5a/contributing.html#autotoc_md44", [
            [ "CUDA", "d1/d5a/contributing.html#autotoc_md45", null ],
            [ "HIP", "d1/d5a/contributing.html#autotoc_md46", null ],
            [ "OpenCL", "d1/d5a/contributing.html#autotoc_md47", null ],
            [ "Device-based type polymorphism", "d1/d5a/contributing.html#autotoc_md48", null ],
            [ "Summary of build system entires", "d1/d5a/contributing.html#autotoc_md49", null ]
          ] ]
        ] ]
      ] ],
      [ "Programming patterns and conventions", "d0/d47/dev_patterns.html", [
        [ "A. Naming", "d0/d47/dev_patterns.html#autotoc_md50", null ],
        [ "B. Scope", "d0/d47/dev_patterns.html#autotoc_md51", null ],
        [ "C. Constructors and destructors.", "d0/d47/dev_patterns.html#autotoc_md52", null ],
        [ "D. Documentation", "d0/d47/dev_patterns.html#autotoc_md53", null ],
        [ "E. Design", "d0/d47/dev_patterns.html#autotoc_md54", null ],
        [ "F. Complex and Polymorphic arrays.", "d0/d47/dev_patterns.html#autotoc_md55", [
          [ "F.1. Object wrapper pattern", "d0/d47/dev_patterns.html#autotoc_md56", null ],
          [ "F.2. Object pointer pattern", "d0/d47/dev_patterns.html#autotoc_md57", null ]
        ] ]
      ] ],
      [ "Code style", "da/db6/code-style.html", [
        [ "Data types", "da/db6/code-style.html#autotoc_md34", null ],
        [ "Enforcing rules", "da/db6/code-style.html#autotoc_md35", [
          [ "Tools", "da/db6/code-style.html#autotoc_md36", [
            [ "flint", "da/db6/code-style.html#autotoc_md37", null ],
            [ "findent", "da/db6/code-style.html#autotoc_md38", null ]
          ] ]
        ] ]
      ] ],
      [ "Testing", "d5/d75/testing.html", [
        [ "pFUnit", "d5/d75/testing.html#autotoc_md63", null ],
        [ "Configuring Neko", "d5/d75/testing.html#autotoc_md64", null ],
        [ "Running the tests", "d5/d75/testing.html#autotoc_md65", null ],
        [ "Adding a new test", "d5/d75/testing.html#autotoc_md66", null ]
      ] ],
      [ "Accelerators", "de/d06/accelerators.html", [
        [ "Device abstraction layer", "de/d06/accelerators.html#autotoc_md29", [
          [ "Memory management", "de/d06/accelerators.html#autotoc_md30", [
            [ "Allocation/deallocation", "de/d06/accelerators.html#autotoc_md31", null ],
            [ "Associate data on host and device", "de/d06/accelerators.html#autotoc_md32", null ],
            [ "Map a host array to a device", "de/d06/accelerators.html#autotoc_md33", null ],
            [ "Data transfer", "de/d06/accelerators.html#accelerators_data-transfer", null ]
          ] ],
          [ "Offload work", "de/d06/accelerators.html#accelerators_offload-work", null ]
        ] ]
      ] ],
      [ "Run-time selectable types", "d5/d5f/rts_types.html", null ],
      [ "Important types", "d3/d40/important_types.html", [
        [ "SEM foundation types", "d3/d40/important_types.html#autotoc_md58", null ],
        [ "Basic math routines", "d3/d40/important_types.html#autotoc_md59", null ],
        [ "Governing equation solvers and related types", "d3/d40/important_types.html#autotoc_md60", null ],
        [ "Singletons", "d3/d40/important_types.html#autotoc_md61", null ],
        [ "Linear algebra", "d3/d40/important_types.html#autotoc_md62", null ]
      ] ]
    ] ],
    [ "Appendices", "da/dd6/appendices.html", [
      [ "Environmental variable reference", "da/dd6/appendices.html#appendices_env-var", [
        [ "Logging level details", "da/dd6/appendices.html#autotoc_md0", null ],
        [ "Gather-scatter communication backend details", "da/dd6/appendices.html#autotoc_md1", null ]
      ] ],
      [ "Governing equations", "db/d27/governing-equations.html", [
        [ "Momentum (Fluid)", "db/d27/governing-equations.html#autotoc_md15", null ],
        [ "Scalar", "db/d27/governing-equations.html#autotoc_md16", null ],
        [ "Compressible Flows", "db/d27/governing-equations.html#autotoc_md17", null ],
        [ "Non-dimensionalisation", "db/d27/governing-equations.html#autotoc_md18", null ]
      ] ],
      [ "Neko .nmsh mesh format", "d0/d5b/nmsh-format.html", [
        [ "File overview", "d0/d5b/nmsh-format.html#autotoc_md19", null ],
        [ "Data types and portability", "d0/d5b/nmsh-format.html#autotoc_md20", null ],
        [ "Element records", "d0/d5b/nmsh-format.html#autotoc_md21", [
          [ "Vertex record (<tt>nmsh_vertex_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md22", null ],
          [ "Quad record (<tt>nmsh_quad_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md23", null ],
          [ "Hex record (<tt>nmsh_hex_t</tt>)", "d0/d5b/nmsh-format.html#autotoc_md24", null ],
          [ "Vertex ordering", "d0/d5b/nmsh-format.html#autotoc_md25", null ]
        ] ],
        [ "Zone records (boundary and periodic data)", "d0/d5b/nmsh-format.html#autotoc_md26", null ],
        [ "Curve records (curved edges/faces)", "d0/d5b/nmsh-format.html#autotoc_md27", null ],
        [ "2D meshes", "d0/d5b/nmsh-format.html#autotoc_md28", null ]
      ] ],
      [ "Neko .fld field format", "d8/d00/fld-format.html", [
        [ "File naming and series", "d8/d00/fld-format.html#autotoc_md2", null ],
        [ "File overview", "d8/d00/fld-format.html#autotoc_md3", null ],
        [ "Header", "d8/d00/fld-format.html#autotoc_md4", [
          [ "rdcode", "d8/d00/fld-format.html#autotoc_md5", null ]
        ] ],
        [ "Test pattern", "d8/d00/fld-format.html#autotoc_md6", null ],
        [ "Element index list", "d8/d00/fld-format.html#autotoc_md7", null ],
        [ "Field blocks", "d8/d00/fld-format.html#autotoc_md8", [
          [ "Vector fields (X and U)", "d8/d00/fld-format.html#autotoc_md9", null ],
          [ "Scalar fields (P, T, S)", "d8/d00/fld-format.html#autotoc_md10", null ],
          [ "Precision", "d8/d00/fld-format.html#autotoc_md11", null ]
        ] ],
        [ "Metadata blocks (3D only)", "d8/d00/fld-format.html#autotoc_md12", null ],
        [ "2D fields and Z-velocity", "d8/d00/fld-format.html#autotoc_md13", null ],
        [ "Portability notes", "d8/d00/fld-format.html#autotoc_md14", null ]
      ] ],
      [ "Publications", "de/d26/publications.html", null ]
    ] ]
];