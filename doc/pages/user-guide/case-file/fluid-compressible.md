# Compressible fluid {#case-file_fluid-compressible}

Neko supports compressible flow simulations via the compressible solver.
To enable compressible flow, set `"scheme": "compressible"` in the fluid
configuration. This solver integrates the compressible Navier-Stokes equations
using a Runge-Kutta time integration scheme with artificial viscosity for
stability.

The compressible solver accepts the following parameters:

| Name    | Description                                      | Admissible values | Default value |
| ------- | ------------------------------------------------ | ----------------- | ------------- |
| `gamma` | Ratio of specific heats for ideal gas            | Positive reals    | `1.4`         |
| `mu`    | Constant physical dynamic viscosity              | Non-negative real | `0.0`         |
| `kappa` | Constant physical thermal conductivity           | Non-negative real | `0.0`         |

The compressible solver uses variable time-stepping controlled by the CFL
number. Set `variable_timestep` to `true` and specify `target_cfl` in the time
control object.

Constant physical viscosity and thermal conductivity can be specified as `mu`
and `kappa` in the fluid section of the case file. Alternatively, they can be
set via the `material_properties` user interface in the user file, which also
allows spatially or temporally varying values. When a user material-properties
routine is provided, it takes precedence over the constant JSON values.
Stabilization uses the existing Laplacian artificial viscosity, computed as the
minimum of entropy-based and low-order viscosities. When `fluid_mu` or
`fluid_kappa` are nonzero, the solver also applies the compressible
Navier-Stokes viscous stress flux and conductive energy flux. If neither
property is provided, both default to zero and no physical viscous flux is added.
See the user file documentation for details on implementing
`material_properties`.

Example configuration:
~~~~~~~~~~~~~~~{.json}
{
  "fluid": {
    "scheme": "compressible",
    "gamma": 1.4,
    "viscous_regularization":{
        "type": "artificial_viscosity",
        "reg_coeff_name": "entropy_viscosity"
    },
    "initial_condition": {
      "type": "user"
    },
    "boundary_conditions": [
      {
        "type": "velocity_value",
        "zone_indices": [1],
        "value": [3, 0, 0]
      },
      {
        "type": "density_value",
        "zone_indices": [1],
        "value": 1.4
      },
      {
        "type": "pressure_value",
        "zone_indices": [1],
        "value": 1
      }
    ],
    "output_control": "nsamples",
    "output_value": 20
  },
  "simulation_components": [
    {
      "type": "artificial_viscosity_model",
      "model": "entropy_viscosity",
      "field_name": "entropy_viscosity"
    }
  ]
}
~~~~~~~~~~~~~~~

## Compressible boundary conditions

The compressible solver supports the following boundary conditions:

| Boundary Condition | Description                               |
| ------------------ | ----------------------------------------- |
| velocity_value     | Dirichlet condition for velocity (inflow) |
| density_value      | Dirichlet condition for density           |
| pressure_value     | Dirichlet condition for pressure          |
| no_slip            | Zero velocity wall                        |
| slip               | A slip wall                               |
| symmetry           | Symmetry plane                            |
| outflow            | Pressure outlet (zero gradient)           |
| normal_outflow     | Normal outflow condition                  |

For examples of compressible flow setups, see the `euler_1d_sod`,
`euler_2d_forward_facing_step`, and `euler_tgv` examples.


## Viscous regularization {#case-file_viscous-regularization}

Users can use the viscous regularization object to enhance the smoothness or 
the numerical stability of the solution via a diffusion term. For example, 
artificial viscosity can be set up by this object to perform shock capturing.
Note that some regularization techniques do not have a diffusive mathematical 
form, for example gradient jump penalty and high-pass filter relaxation terms,
and they are included in the source terms instead. The viscous regularization
can be set up with the following options:

* `type`, the viscous regularization type.
  - `artificial_viscosity`, the standard second-order diffusion term,
    \f$\frac{\partial}{\partial x}\left(\mu_\mathrm{artificial}\frac{\partial u}{\partial x}\right)\f$.
* `reg_coeff_name`, name of the \f$\mu_\mathrm{artificial}\f$ field, usually 
computed by a simulation component.

Viscous regularization is currently supported by the compressible fluid
solver. Artificial viscosity requires both the regularization object that
consumes the coefficient and a simulation component that computes it. The
`reg_coeff_name` and `field_name` values must match:

~~~~~~~~~~~~~~~{.json}
{
  "fluid": {
    "scheme": "compressible",
    "viscous_regularization": {
      "type": "artificial_viscosity",
      "reg_coeff_name": "entropy_viscosity"
    }
  },
  "simulation_components": [
    {
      "type": "artificial_viscosity_model",
      "model": "entropy_viscosity",
      "field_name": "entropy_viscosity"
    }
  ]
}
~~~~~~~~~~~~~~~
