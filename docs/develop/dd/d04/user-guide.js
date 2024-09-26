var user_guide =
[
    [ "Installing Neko", "d5/dfc/installation.html", "d5/dfc/installation" ],
    [ "Case File", "dd/d33/case-file.html", [
      [ "High-level structure", "dd/d33/case-file.html#autotoc_md30", null ],
      [ "Output frequency control", "dd/d33/case-file.html#autotoc_md31", null ],
      [ "The case object", "dd/d33/case-file.html#autotoc_md32", null ],
      [ "Numerics", "dd/d33/case-file.html#autotoc_md33", null ],
      [ "Fluid", "dd/d33/case-file.html#autotoc_md34", [
        [ "Material properties", "dd/d33/case-file.html#autotoc_md35", null ],
        [ "Boundary types", "dd/d33/case-file.html#case-file_boundary-types", null ],
        [ "Inflow boundary conditions", "dd/d33/case-file.html#case-file_fluid-if", null ],
        [ "Initial conditions", "dd/d33/case-file.html#case-file_fluid-ic", null ],
        [ "Blasius profile", "dd/d33/case-file.html#autotoc_md36", null ],
        [ "Source terms", "dd/d33/case-file.html#case-file_fluid-source-term", [
          [ "Brinkman", "dd/d33/case-file.html#autotoc_md37", null ]
        ] ]
      ] ],
      [ "Linear solver configuration", "dd/d33/case-file.html#autotoc_md38", [
        [ "Flow rate forcing", "dd/d33/case-file.html#autotoc_md39", null ],
        [ "Full parameter table", "dd/d33/case-file.html#autotoc_md40", null ]
      ] ],
      [ "Scalar", "dd/d33/case-file.html#case-file_scalar", [
        [ "Material properties", "dd/d33/case-file.html#autotoc_md41", null ],
        [ "Boundary types", "dd/d33/case-file.html#autotoc_md42", null ],
        [ "Initial conditions", "dd/d33/case-file.html#autotoc_md43", null ],
        [ "Source terms", "dd/d33/case-file.html#autotoc_md44", null ],
        [ "Full parameter table", "dd/d33/case-file.html#autotoc_md45", null ]
      ] ],
      [ "Statistics", "dd/d33/case-file.html#autotoc_md46", null ],
      [ "Simulation components", "dd/d33/case-file.html#autotoc_md47", null ],
      [ "Point zones", "dd/d33/case-file.html#autotoc_md48", null ],
      [ "Runtime statistics", "dd/d33/case-file.html#autotoc_md49", null ]
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
    [ "Simulaticn components", "d3/d84/simcomps.html", [
      [ "What are simulation components?", "d3/d84/simcomps.html#autotoc_md75", null ],
      [ "Adding simulation components to the case", "d3/d84/simcomps.html#autotoc_md76", null ],
      [ "List of simulation components", "d3/d84/simcomps.html#autotoc_md77", null ],
      [ "Controling execution and file output", "d3/d84/simcomps.html#autotoc_md78", [
        [ "vorticity", "d3/d84/simcomps.html#simcomp_vorticity", null ],
        [ "lambda2", "d3/d84/simcomps.html#simcomp_lambda2", null ],
        [ "probes", "d3/d84/simcomps.html#simcomp_probes", [
          [ "Supported types:", "d3/d84/simcomps.html#autotoc_md79", null ],
          [ "Example usage:", "d3/d84/simcomps.html#autotoc_md80", null ]
        ] ],
        [ "field_writer", "d3/d84/simcomps.html#simcomp_field_writer", null ],
        [ "derivative", "d3/d84/simcomps.html#simcomp_derivative", null ],
        [ "force_torque", "d3/d84/simcomps.html#simcomp_force_torque", null ],
        [ "weak_grad", "d3/d84/simcomps.html#simcomp_weak_grad", null ]
      ] ]
    ] ],
    [ "Point zones", "da/dd0/point-zones.html", [
      [ "What are point zones?", "da/dd0/point-zones.html#autotoc_md66", null ],
      [ "Predefined geometrical shapes", "da/dd0/point-zones.html#autotoc_md67", [
        [ "Box", "da/dd0/point-zones.html#autotoc_md68", null ],
        [ "Sphere", "da/dd0/point-zones.html#autotoc_md69", null ],
        [ "Cylinder", "da/dd0/point-zones.html#autotoc_md70", null ]
      ] ],
      [ "Operations on point zones", "da/dd0/point-zones.html#autotoc_md71", [
        [ "Inversion", "da/dd0/point-zones.html#autotoc_md72", null ],
        [ "Combination", "da/dd0/point-zones.html#autotoc_md73", null ]
      ] ],
      [ "User-defined geometrical shapes", "da/dd0/point-zones.html#autotoc_md74", null ],
      [ "Using point zones", "da/dd0/point-zones.html#point-zones_using-point-zones", null ]
    ] ],
    [ "Statistics guide", "df/d8f/statistics-guide.html", [
      [ "Using statistics", "df/d8f/statistics-guide.html#autotoc_md81", null ],
      [ "List of fields in output files", "df/d8f/statistics-guide.html#autotoc_md82", null ],
      [ "Postprocessing", "df/d8f/statistics-guide.html#autotoc_md83", [
        [ "Notes on the statistics calculation in Neko", "df/d8f/statistics-guide.html#autotoc_md84", null ]
      ] ]
    ] ],
    [ "Input-output", "d7/d7f/io.html", [
      [ "Mesh", "d7/d7f/io.html#autotoc_md63", null ],
      [ "Three-dimensional field output", "d7/d7f/io.html#autotoc_md64", null ],
      [ "Checkpoint files", "d7/d7f/io.html#autotoc_md65", null ]
    ] ]
];