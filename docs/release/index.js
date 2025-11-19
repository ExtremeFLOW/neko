var index =
[
    [ "Structure of the Manual", "index.html#autotoc_md29", null ],
    [ "User guide", "dd/d04/user-guide.html", [
      [ "Installing Neko", "d5/dfc/installation.html", [
        [ "Building from source", "d5/dfc/installation.html#autotoc_md51", [
          [ "Dependencies", "d5/dfc/installation.html#autotoc_md52", [
            [ "Building JSON Fortran", "d5/dfc/installation.html#autotoc_md53", null ],
            [ "Building gslib (optional)", "d5/dfc/installation.html#autotoc_md54", null ],
            [ "Building ParMETIS (optional)", "d5/dfc/installation.html#autotoc_md55", null ],
            [ "Bulding PFunit (optional)", "d5/dfc/installation.html#autotoc_md56", null ]
          ] ],
          [ "Building Neko", "d5/dfc/installation.html#autotoc_md57", [
            [ "Compiling Neko for CPU or SX-Aurora", "d5/dfc/installation.html#autotoc_md58", null ],
            [ "Compiling Neko for NVIDIA GPUs", "d5/dfc/installation.html#autotoc_md59", null ],
            [ "Compiling Neko for AMD GPUs", "d5/dfc/installation.html#autotoc_md60", null ]
          ] ]
        ] ],
        [ "Installing via Spack", "d5/dfc/installation.html#autotoc_md61", [
          [ "Quick start guide with Spack", "d5/dfc/installation.html#autotoc_md62", null ]
        ] ],
        [ "Using a Docker container", "d5/dfc/installation.html#autotoc_md63", null ],
        [ "Testing", "d5/d75/testing.html", [
          [ "pFUnit", "d5/d75/testing.html#autotoc_md25", null ],
          [ "Configuring Neko", "d5/d75/testing.html#autotoc_md26", null ],
          [ "Running the tests", "d5/d75/testing.html#autotoc_md27", null ],
          [ "Adding a new test", "d5/d75/testing.html#autotoc_md28", null ]
        ] ]
      ] ],
      [ "Case File", "dd/d33/case-file.html", [
        [ "High-level structure", "dd/d33/case-file.html#autotoc_md30", null ],
        [ "Output frequency control", "dd/d33/case-file.html#autotoc_md31", null ],
        [ "The case object", "dd/d33/case-file.html#autotoc_md32", null ],
        [ "Numerics", "dd/d33/case-file.html#autotoc_md33", null ],
        [ "Fluid", "dd/d33/case-file.html#autotoc_md34", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md35", null ],
          [ "Boundary types", "dd/d33/case-file.html#case-file_boundary-types", null ],
          [ "Inflow boundary conditions", "dd/d33/case-file.html#case-file_fluid-if", null ],
          [ "Shear stress boundary conditions", "dd/d33/case-file.html#case-file_fluid-sh", null ],
          [ "Wall model  boundary conditions", "dd/d33/case-file.html#case-file_fluid-wm", null ],
          [ "Initial conditions", "dd/d33/case-file.html#case-file_fluid-ic", null ],
          [ "Blasius profile", "dd/d33/case-file.html#autotoc_md36", null ],
          [ "Source terms", "dd/d33/case-file.html#case-file_fluid-source-term", [
            [ "Brinkman", "dd/d33/case-file.html#autotoc_md37", null ]
          ] ],
          [ "Gradient Jump Penalty", "dd/d33/case-file.html#autotoc_md38", null ]
        ] ],
        [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md39", [
          [ "Flow rate forcing", "dd/d33/case-file.html#autotoc_md40", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md41", null ]
        ] ],
        [ "Scalar", "dd/d33/case-file.html#case-file_scalar", [
          [ "Material properties", "dd/d33/case-file.html#autotoc_md42", null ],
          [ "Boundary types", "dd/d33/case-file.html#autotoc_md43", null ],
          [ "Initial conditions", "dd/d33/case-file.html#autotoc_md44", null ],
          [ "Source terms", "dd/d33/case-file.html#autotoc_md45", null ],
          [ "Full parameter table", "dd/d33/case-file.html#autotoc_md46", null ]
        ] ],
        [ "Statistics", "dd/d33/case-file.html#autotoc_md47", null ],
        [ "Simulation components", "dd/d33/case-file.html#autotoc_md48", null ],
        [ "Point zones", "dd/d33/case-file.html#autotoc_md49", null ],
        [ "Runtime statistics", "dd/d33/case-file.html#autotoc_md50", null ]
      ] ],
      [ "User File", "d6/def/user-file.html", [
        [ "Compiling and running", "d6/def/user-file.html#autotoc_md85", null ],
        [ "High-level structure", "d6/def/user-file.html#autotoc_md86", null ],
        [ "Default user functions", "d6/def/user-file.html#autotoc_md87", [
          [ "Initializing and finalizing", "d6/def/user-file.html#user-file_init-and-final", null ],
          [ "Computing at every time step", "d6/def/user-file.html#user-file_user-check", null ],
          [ "Setting material properties", "d6/def/user-file.html#user-file_mat-prop", null ],
          [ "Runtime mesh deformation", "d6/def/user-file.html#user-file_user-mesh-setup", null ],
          [ "Scalar boundary conditions", "d6/def/user-file.html#user-file_scalar-bc", null ],
          [ "User defined simulation components", "d6/def/user-file.html#user-file_simcomps", null ]
        ] ],
        [ "Case-specific user functions", "d6/def/user-file.html#autotoc_md88", [
          [ "Fluid and Scalar initial conditions", "d6/def/user-file.html#user-file_user-ic", null ],
          [ "Fluid inflow condition", "d6/def/user-file.html#user-file_fluid-user-if", null ],
          [ "Fluid and scalar source terms", "d6/def/user-file.html#user-file_user-f", null ],
          [ "Complex fluid and/or scalar boundary conditions", "d6/def/user-file.html#user-file_field-dirichlet-update", null ]
        ] ],
        [ "Additional remarks and tips", "d6/def/user-file.html#autotoc_md89", [
          [ "Running on GPUs", "d6/def/user-file.html#user-file_tips_running-on-gpus", null ],
          [ "Registries", "d6/def/user-file.html#user-file_tips_registries", null ]
        ] ]
      ] ],
      [ "Simulation components", "d3/d84/simcomps.html", [
        [ "What are simulation components?", "d3/d84/simcomps.html#autotoc_md76", null ],
        [ "Adding simulation components to the case", "d3/d84/simcomps.html#autotoc_md77", null ],
        [ "List of simulation components", "d3/d84/simcomps.html#autotoc_md78", null ],
        [ "Controling execution and file output", "d3/d84/simcomps.html#autotoc_md79", [
          [ "vorticity", "d3/d84/simcomps.html#simcomp_vorticity", null ],
          [ "lambda2", "d3/d84/simcomps.html#simcomp_lambda2", null ],
          [ "probes", "d3/d84/simcomps.html#simcomp_probes", [
            [ "Supported types:", "d3/d84/simcomps.html#autotoc_md80", null ],
            [ "Example usage:", "d3/d84/simcomps.html#autotoc_md81", null ]
          ] ],
          [ "field_writer", "d3/d84/simcomps.html#simcomp_field_writer", null ],
          [ "derivative", "d3/d84/simcomps.html#simcomp_derivative", null ],
          [ "force_torque", "d3/d84/simcomps.html#simcomp_force_torque", null ],
          [ "weak_grad", "d3/d84/simcomps.html#simcomp_weak_grad", null ],
          [ "Spectral error indicator", "d3/d84/simcomps.html#simcomp_speri", null ]
        ] ]
      ] ],
      [ "Point zones", "da/dd0/point-zones.html", [
        [ "What are point zones?", "da/dd0/point-zones.html#autotoc_md67", null ],
        [ "Predefined geometrical shapes", "da/dd0/point-zones.html#autotoc_md68", [
          [ "Box", "da/dd0/point-zones.html#autotoc_md69", null ],
          [ "Sphere", "da/dd0/point-zones.html#autotoc_md70", null ],
          [ "Cylinder", "da/dd0/point-zones.html#autotoc_md71", null ]
        ] ],
        [ "Operations on point zones", "da/dd0/point-zones.html#autotoc_md72", [
          [ "Inversion", "da/dd0/point-zones.html#autotoc_md73", null ],
          [ "Combination", "da/dd0/point-zones.html#autotoc_md74", null ]
        ] ],
        [ "User-defined geometrical shapes", "da/dd0/point-zones.html#autotoc_md75", null ],
        [ "Using point zones", "da/dd0/point-zones.html#point-zones_using-point-zones", null ]
      ] ],
      [ "Statistics guide", "df/d8f/statistics-guide.html", [
        [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md82", null ],
        [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md83", null ],
        [ "Postprocessing", "df/d8f/statistics-guide.html#autotoc_md84", null ]
      ] ],
      [ "Input-output", "d7/d7f/io.html", [
        [ "Mesh", "d7/d7f/io.html#autotoc_md64", null ],
        [ "Three-dimensional field output", "d7/d7f/io.html#autotoc_md65", null ],
        [ "Checkpoint files", "d7/d7f/io.html#autotoc_md66", null ]
      ] ]
    ] ],
    [ "Developer guide", "dc/d70/developer-guide.html", [
      [ "Contributing to Neko", "d1/d5a/contributing.html", [
        [ "Git branches", "d1/d5a/contributing.html#autotoc_md11", null ],
        [ "Code style", "d1/d5a/contributing.html#autotoc_md12", [
          [ "Data types", "d1/d5a/contributing.html#autotoc_md13", null ]
        ] ],
        [ "Build system", "d1/d5a/contributing.html#autotoc_md14", null ]
      ] ],
      [ "Programming patterns and conventions", "d0/d47/dev_patterns.html", [
        [ "A. Naming", "d0/d47/dev_patterns.html#autotoc_md15", null ],
        [ "B. Scope", "d0/d47/dev_patterns.html#autotoc_md16", null ],
        [ "C. Constructors and destructors.", "d0/d47/dev_patterns.html#autotoc_md17", null ],
        [ "D. Documentation", "d0/d47/dev_patterns.html#autotoc_md18", null ],
        [ "E. Design", "d0/d47/dev_patterns.html#autotoc_md19", null ]
      ] ],
      [ "Code style", "da/db6/code-style.html", [
        [ "Data types", "da/db6/code-style.html#autotoc_md8", null ],
        [ "Linting rules", "da/db6/code-style.html#autotoc_md9", null ],
        [ "Tools", "da/db6/code-style.html#autotoc_md10", null ]
      ] ],
      [ "Testing", "d5/d75/testing.html", [
        [ "pFUnit", "d5/d75/testing.html#autotoc_md25", null ],
        [ "Configuring Neko", "d5/d75/testing.html#autotoc_md26", null ],
        [ "Running the tests", "d5/d75/testing.html#autotoc_md27", null ],
        [ "Adding a new test", "d5/d75/testing.html#autotoc_md28", null ]
      ] ],
      [ "Accelerators", "de/d06/accelerators.html", [
        [ "Device abstraction layer", "de/d06/accelerators.html#autotoc_md3", [
          [ "Memory management", "de/d06/accelerators.html#autotoc_md4", [
            [ "Allocation/deallocation", "de/d06/accelerators.html#autotoc_md5", null ],
            [ "Associate data on host and device", "de/d06/accelerators.html#autotoc_md6", null ],
            [ "Map a host array to a device", "de/d06/accelerators.html#autotoc_md7", null ],
            [ "Data transfer", "de/d06/accelerators.html#accelerators_data-transfer", null ]
          ] ],
          [ "Offload work", "de/d06/accelerators.html#accelerators_offload-work", null ]
        ] ]
      ] ],
      [ "Run-time selectable types", "d5/d5f/rts_types.html", null ],
      [ "Important types", "d3/d40/important_types.html", [
        [ "SEM foundation types", "d3/d40/important_types.html#autotoc_md20", null ],
        [ "Basic math routines", "d3/d40/important_types.html#autotoc_md21", null ],
        [ "Governing equation solvers and related types", "d3/d40/important_types.html#autotoc_md22", null ],
        [ "Singletons", "d3/d40/important_types.html#autotoc_md23", null ],
        [ "Linear algebra", "d3/d40/important_types.html#autotoc_md24", null ]
      ] ]
    ] ],
    [ "Appendices", "da/dd6/appendices.html", [
      [ "Environmental variable reference", "da/dd6/appendices.html#appendices_env-var", [
        [ "Logging level details", "da/dd6/appendices.html#autotoc_md0", null ]
      ] ],
      [ "Governing equations", "db/d27/governing-equations.html", [
        [ "Fluid", "db/d27/governing-equations.html#autotoc_md1", null ],
        [ "Scalar", "db/d27/governing-equations.html#autotoc_md2", null ]
      ] ],
      [ "Publications", "de/d26/publications.html", null ]
    ] ]
];