var user_guide =
[
    [ "Installing Neko", "d5/dfc/installation.html", "d5/dfc/installation" ],
    [ "Meshing", "d9/df2/meshing.html", [
      [ "General considerations", "d9/df2/meshing.html#general-considerations", null ],
      [ "Constructing meshes", "d9/df2/meshing.html#constructing-meshes", [
        [ "Box meshes", "d9/df2/meshing.html#constructing-meshes-box", null ],
        [ "Converting a Gmsh mesh", "d9/df2/meshing.html#constructing-meshes-gmsh", null ],
        [ "Converting a Nek5000 mesh", "d9/df2/meshing.html#constructing-meshes-nek5000", null ],
        [ "Experimental workflows", "d9/df2/meshing.html#constructing-meshes-experimental", null ]
      ] ],
      [ "High-order boundary representation", "d9/df2/meshing.html#high-order-boundary-representation", [
        [ "Curves stored in the mesh file", "d9/df2/meshing.html#autotoc_md156", null ],
        [ "User-defined geometry deformation", "d9/df2/meshing.html#autotoc_md157", null ]
      ] ],
      [ "Adding periodicity to an existing <span class=\"tt\">.nmsh</span>", "d9/df2/meshing.html#adding-periodicity-to-an-existing-nmsh", null ]
    ] ],
    [ "Case File", "dd/d33/case-file.html", "dd/d33/case-file" ],
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
        [ "User wall sampling", "d6/def/user-file.html#user-file_wall-sampling", null ],
        [ "Fluid and Scalar initial conditions", "d6/def/user-file.html#user-file_user-ic", null ],
        [ "Fluid and scalar source terms", "d6/def/user-file.html#user-file_user-f", null ],
        [ "Dirichlet boundary conditions", "d6/def/user-file.html#user-file_field-dirichlet-update", null ],
        [ "Neumann boundary conditions", "d6/def/user-file.html#user-file_field-neumann-update", null ]
      ] ],
      [ "Arbitrary Lagrangian-Eulerian (ALE) user functions", "d6/def/user-file.html#user-file_ale", [
        [ "Custom rigid body motion", "d6/def/user-file.html#user-file_ale-rigid_motion", null ],
        [ "Registries", "d6/def/user-file.html#user-file_tips_registries", null ],
        [ "User access to solver internals", "d6/def/user-file.html#user-file_access", null ]
      ] ]
    ] ],
    [ "Simulation components", "d3/d84/simcomps.html", [
      [ "What are simulation components?", "d3/d84/simcomps.html#autotoc_md177", null ],
      [ "Adding simulation components to the case", "d3/d84/simcomps.html#autotoc_md178", null ],
      [ "List of simulation components", "d3/d84/simcomps.html#autotoc_md179", null ],
      [ "Controlling execution and file output", "d3/d84/simcomps.html#autotoc_md180", [
        [ "Differential operators", "d3/d84/simcomps.html#autotoc_md181", [
          [ "derivative", "d3/d84/simcomps.html#simcomp_derivative", null ],
          [ "curl", "d3/d84/simcomps.html#simcomp_curl", null ],
          [ "divergence", "d3/d84/simcomps.html#simcomp_divergence", null ],
          [ "gradient", "d3/d84/simcomps.html#simcomp_gradient", null ],
          [ "weak_gradient", "d3/d84/simcomps.html#simcomp_weak_gradient", null ]
        ] ],
        [ "lambda2", "d3/d84/simcomps.html#simcomp_lambda2", null ],
        [ "boundary_operation", "d3/d84/simcomps.html#simcomp_boundary_operation", null ],
        [ "boundary_flux", "d3/d84/simcomps.html#simcomp_boundary_flux", null ],
        [ "boundary_data_writer", "d3/d84/simcomps.html#simcomp_boundary_data_writer", [
          [ "Output files", "d3/d84/simcomps.html#autotoc_md182", null ]
        ] ],
        [ "probes", "d3/d84/simcomps.html#simcomp_probes", [
          [ "Supported types", "d3/d84/simcomps.html#autotoc_md183", null ],
          [ "Example usage", "d3/d84/simcomps.html#autotoc_md184", null ]
        ] ],
        [ "field_writer", "d3/d84/simcomps.html#simcomp_field_writer", null ],
        [ "lagrangian_particles", "d3/d84/simcomps.html#simcomp_lagrangian_particles", null ],
        [ "force_torque", "d3/d84/simcomps.html#simcomp_force_torque", [
          [ "Torque calculation for moving bodies", "d3/d84/simcomps.html#autotoc_md185", null ]
        ] ],
        [ "wall_shear_stress", "d3/d84/simcomps.html#simcomp_wall_shear_stress", null ],
        [ "les_model", "d3/d84/simcomps.html#simcomp_les_model", null ],
        [ "User statistics", "d3/d84/simcomps.html#user_stats", null ],
        [ "Spatial average", "d3/d84/simcomps.html#simcomp_spatial_average", null ],
        [ "Spectral error indicator", "d3/d84/simcomps.html#simcomp_speri", null ],
        [ "Data streamer", "d3/d84/simcomps.html#simcomp_data_streamer", null ],
        [ "Field subsampler", "d3/d84/simcomps.html#simcomp_field_subsampler", null ],
        [ "artificial_viscosity_model", "d3/d84/simcomps.html#simcomp_avm", null ]
      ] ]
    ] ],
    [ "Point zones", "da/dd0/point-zones.html", [
      [ "What are point zones?", "da/dd0/point-zones.html#autotoc_md167", null ],
      [ "Predefined geometrical shapes", "da/dd0/point-zones.html#autotoc_md168", [
        [ "Box", "da/dd0/point-zones.html#autotoc_md169", null ],
        [ "Sphere", "da/dd0/point-zones.html#autotoc_md170", null ],
        [ "Cylinder", "da/dd0/point-zones.html#autotoc_md171", null ]
      ] ],
      [ "Operations on point zones", "da/dd0/point-zones.html#autotoc_md172", [
        [ "Inversion", "da/dd0/point-zones.html#autotoc_md173", null ],
        [ "Combination", "da/dd0/point-zones.html#autotoc_md174", null ],
        [ "Including full element data", "da/dd0/point-zones.html#autotoc_md175", null ]
      ] ],
      [ "User-defined geometrical shapes", "da/dd0/point-zones.html#autotoc_md176", null ],
      [ "Using point zones", "da/dd0/point-zones.html#point-zones_using-point-zones", null ]
    ] ],
    [ "Statistics guide", "df/d8f/statistics-guide.html", [
      [ "Output filenames and restarts", "df/d8f/statistics-guide.html#statistics-output-filenames", null ],
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
      [ "Mesh", "d7/d7f/io.html#autotoc_md153", null ],
      [ "Three-dimensional field output", "d7/d7f/io.html#autotoc_md154", [
        [ "Compression of field output", "d7/d7f/io.html#autotoc_md155", null ]
      ] ],
      [ "Checkpoint files", "d7/d7f/io.html#checkpoint-files", null ],
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
      [ "Installation", "d8/daa/api.html#autotoc_md66", [
        [ "Python (pyneko)", "d8/daa/api.html#autotoc_md67", null ],
        [ "Julia (Neko.jl)", "d8/daa/api.html#autotoc_md68", null ]
      ] ]
    ] ],
    [ "Performance guidelines", "dc/d3c/performance.html", [
      [ "Installation", "dc/d3c/performance.html#autotoc_md158", [
        [ "Accelerator specific options", "dc/d3c/performance.html#autotoc_md159", null ],
        [ "CPU specific options", "dc/d3c/performance.html#autotoc_md160", null ]
      ] ],
      [ "Simulation setup", "dc/d3c/performance.html#autotoc_md161", [
        [ "Load balancing", "dc/d3c/performance.html#autotoc_md162", null ],
        [ "Parameters", "dc/d3c/performance.html#autotoc_md163", null ],
        [ "SEM operator auto-tuning", "dc/d3c/performance.html#performance-operator-autotuning", null ]
      ] ],
      [ "Running a simulation", "dc/d3c/performance.html#autotoc_md164", [
        [ "Gather-scatter communication backends", "dc/d3c/performance.html#autotoc_md165", [
          [ "The crystal router backends", "dc/d3c/performance.html#performance-gs-crystal", null ],
          [ "Runtime autotuning of the comm. backend", "dc/d3c/performance.html#performance-gs-autotuning", [
            [ "Choosing the candidates", "dc/d3c/performance.html#autotoc_md166", null ]
          ] ],
          [ "MPI neighbourhood collective backend", "dc/d3c/performance.html#performance-neighbour-backend", null ],
          [ "NVSHMEM backend", "dc/d3c/performance.html#performance-nvshmem-backend", null ],
          [ "OpenSHMEM backend", "dc/d3c/performance.html#performance-openshmem-backend", null ],
          [ "Coarray Fortran backend", "dc/d3c/performance.html#performance-caf-backend", null ],
          [ "uTofu backend", "dc/d3c/performance.html#performance-utofu-backend", null ]
        ] ]
      ] ]
    ] ],
    [ "Global Interpolation", "dd/d61/global-interpolation.html", [
      [ "Overview", "dd/d61/global-interpolation.html#autotoc_md105", null ],
      [ "The global interpolation module", "dd/d61/global-interpolation.html#autotoc_md107", [
        [ "Description", "dd/d61/global-interpolation.html#autotoc_md108", null ],
        [ "Features", "dd/d61/global-interpolation.html#autotoc_md109", null ]
      ] ],
      [ "The global interpolation type", "dd/d61/global-interpolation.html#autotoc_md111", [
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
      [ "Filter base type", "df/d4a/filter.html#filter_base_type", null ],
      [ "PDE filter", "df/d4a/filter.html#filter_pde", null ],
      [ "Elementwise filter", "df/d4a/filter.html#filter_elementwise", null ],
      [ "High-pass filter relaxation source term", "df/d4a/filter.html#filter_hpfrt", null ]
    ] ]
];