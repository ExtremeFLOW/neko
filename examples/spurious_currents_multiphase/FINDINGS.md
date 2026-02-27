# Spurious Currents: Curvature Sign Bug and Fix

## Problem

A stationary circular drop (R=0.2) simulated with the phase-field / CSF surface
tension method in Neko developed catastrophic spurious currents and blew up
within ~250 time steps.

**Parameters:** rho=300, mu=0.1, sigma=1.0, La=12000, epsilon=0.01, 10x10 SEM mesh (p=7).

## Root cause: sign error in curvature

The phase field has phi=1 inside the drop and phi=0 outside, so grad(phi) points
**inward**.  The code computed kappa = div(n_hat), giving kappa = -1/R = -5.0.
The Brackbill CSF convention requires **kappa = -div(n_hat) = +5.0** so that the
surface tension force F = sigma * kappa * grad(phi) points inward, balancing the
Laplace pressure.

With the wrong sign, the force pushed the interface **apart** instead of holding
it together, creating a positive feedback loop leading to blowup.

![Sign convention diagram](sign_convention_diagram.png)

## The fix

One line added after each `div()` call (commit `637fc27`):

```fortran
call div(temp4%x, temp1%x, temp2%x, temp3%x, coef)
call coef%gs_h%op(temp4, GS_OP_ADD)
call col2(temp4%x, coef%mult, temp4%size())
call cmult(temp4%x, -1.0_rp, temp4%size())    ! negate for Brackbill convention
```

## Verification

Rerun on Dardel with fixed code (same parameters, t=0.3):

| Metric | Buggy (t=0.25) | Fixed (t=0.3) |
|--------|---------------|--------------|
| u_max | 2.4e+01 | **4.8e-07** |
| kappa_rms | 7.9e+02 | **5.02** |
| F_ST max | 1.8e+07 | **125** |
| phi bounds | [-11.2, 9.6] | **[~0, ~1.0]** |
| Laplace Dp | diverged | **5.04** (expected: 5.0) |

The fixed simulation is stable with spurious velocities at machine-precision
level.  The Laplace pressure jump matches sigma/R = 5.0 to within 1%.
