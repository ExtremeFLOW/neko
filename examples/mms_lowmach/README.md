# Low-Mach viscous-stress MMS verification

A method-of-manufactured-solutions (MMS) check for the viscous stress in **both**
low-Mach residuals (`lowmach_vel_res_cpu_compute` and `lowmach_prs_res_cpu_compute`),
with a spatially varying viscosity.

## What it verifies

The deviatoric viscous stress is

    tau = mu (grad u + grad u^T) - (2/3) mu (div u) I .

* **Velocity (momentum) residual** computes `div[mu(grad u+grad u^T)]` with the
  existing coupled full-stress operator `ax_helm_full` (mu in `coef%h1`, so the
  variable-mu `S^T grad mu` coupling is handled inside the operator) and adds the
  dilatation `-(2/3) grad(mu Q_T)` explicitly (`opgrad`). `ax_helm`'s interface is
  never touched.
* **Pressure residual** computes the viscous force in rotational form,
  `mu(curl curl u) - S^T(2 grad mu)`, plus the low-Mach corrections
  `-(4/3) mu grad(Q_T) + (2/3) Q_T grad(mu)`, and adds the continuity source
  `(bd/dt) Q_T B` for `div u = Q_T`.

## Manufactured solution

On the fully periodic box `[0,2pi]^3`:

    u   = (sin x, sin y, sin z)            (so div u = cos x + cos y + cos z)
    Q_T = cos x + cos y + cos z            (= div u, exactly, by construction)
    mu  = mu0 ( 1 + a sin x sin y sin z )  (variable viscosity, a = 0.5)

Because `grad u` is diagonal, the stress divergence collapses to the clean form

    (div tau)_i = -(4/3) mu sin(x_i) + d_i(mu) ( 2 cos(x_i) - (2/3) Q_T )

(reducing to `-(4/3) mu0 sin(x_i)` for constant mu). This is the analytic target.

* **Test A (momentum):** the shipped velocity residual with `p=0, f=0, bd=0`
  returns the weak `div(tau)`; gather-scatter + inverse mass gives the nodal field,
  compared to the analytic target.
* **Test B (pressure):** the shipped pressure residual with `p=0, f=0, bd=0, rho=1`
  returns the weak `div(div tau)`. It is compared to the discrete divergence
  (`cdtp`, same post-processing) of the analytic stress vector. Because
  `div u = Q_T` holds exactly, the rotational form and the direct form must agree.

## Running

    ./run.sh

## Expected output (observed)

Both errors fall ~2 decades per +2 polynomial order — spectral convergence,
confirming the signs and coefficients of the variable-mu stress and the low-Mach
corrections in both residuals.

     lx    momentum: Linf      relL2        pressure: Linf      relL2
    --------------------------------------------------------------------
      4    7.289112028E-02    2.969035996E-02    2.367338488E-01    3.783617161E-02
      6    3.819123419E-03    1.225684536E-03    6.059511539E-03    6.151142739E-04
      8    7.263377222E-05    2.005839364E-05    4.931979181E-05    3.666921495E-06
     10    7.302242137E-07    1.789389057E-07    1.922212067E-07    1.126772214E-08

## Scope / limitations

* `grad u` is diagonal for this manufactured field, so the curl-curl term is not
  exercised here (it is the pre-validated upstream incompressible-stress
  machinery); the new low-Mach terms and the variable-mu coupling are.
* The `(bd/dt) Q_T B` continuity source is a direct algebraic term and is not
  separately convergence-tested (Test B uses `bd=0`).
* CPU backend only.
