# Progress update slides: CDI/CSF in turbulent channel flow

Each slide is marked with its layout (SINGLE or TWO-COLUMN).
All equations are in LaTeX: paste into PowerPoint via Insert > Equation.
Speaker notes are marked **[NOTES]**.

=========================================================
## SLIDE 1: Problem setup
**Layout: SINGLE COLUMN**
=========================================================

**Title:** CDI + CSF interface capturing in turbulent channel flow: method validation

**Goal**
Validate CDI (Olsson–Kreiss) + CSF (Brackbill) in a turbulent channel as a
stress-test environment. Turbulence is the environment; interface capturing is
what is being validated.

**Setup**
- Domain: $4\pi \times 2 \times 4\pi/3$, walls at $y = \pm 1$, periodic in $x,z$
- $Re_b = 2800$, $Re_\tau \approx 180$, constant bulk velocity $U_b = 1$ (body force)
- Single drop, $R = 0.3$, centred at $y = 0$, matched fluids ($\rho = \mu = \text{const}$)
- Spectral element method (Neko): $81 \times 18 \times 27$ mesh, $N = 7$

**Governing equations**

Momentum:
$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -\nabla p + \frac{1}{Re_b}\nabla^2\mathbf{u} + \underbrace{\sigma\kappa\nabla\varphi}_{\mathbf{F}_\text{ST}}$$

Phase field with CDI resharpening:
$$\frac{\partial \varphi}{\partial t} + \mathbf{u}\cdot\nabla\varphi = \nabla\cdot\!\left[\varepsilon\gamma_c\nabla\varphi - \gamma_c\,\varphi(1-\varphi)\,\hat{\mathbf{n}}\right]$$

Curvature and adaptive CDI velocity:
$$\kappa = -\nabla\cdot\hat{\mathbf{n}}, \qquad \hat{\mathbf{n}} = \frac{\nabla\varphi}{|\nabla\varphi|}, \qquad \gamma_c = \gamma_\text{param}\,u_\text{max}$$

**Workflow**
- Step 1: single-phase spin-up to $t = 20$ TU, checkpoint ($Re_\tau = 180$ confirmed)
- Step 2: inject drop via $\varphi_0 = \tfrac{1}{2}\!\left(1 + \tanh\tfrac{R - |\mathbf{x}-\mathbf{x}_c|}{2\varepsilon}\right)$, run two-phase $t = 20 \to 25$ TU

**[NOTES]**
$\gamma_c$ scales with $u_\text{max}$ so CDI keeps pace with velocity changes.
CSF sign: $\kappa < 0$ for convex drop ($\varphi = 1$ inside), so $\mathbf{F}_\text{ST}$ points inward (restoring).
Sign fixed in commit 637fc27 of this fork.

=========================================================
## SLIDE 2: Parameter selection
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: Interface and CDI parameters**

Interface width $\varepsilon = 0.07$:
$$\frac{\varepsilon}{\Delta_\text{GLL}} \geq 3, \qquad \Delta_\text{GLL,xz} = \frac{4\pi/81}{7} \approx 0.022 \quad\Rightarrow\quad \varepsilon \geq 0.066$$

Drop integrity:
$$\frac{R}{2\varepsilon} = \frac{0.3}{0.14} = 2.1 \quad\Rightarrow\quad \varphi(0) = \tfrac{1}{2}\!\left(1 + \tanh 2.1\right) = 0.986 \;\checkmark$$

CDI compression velocity $\gamma_\text{param} = 0.05$:
$$0.01 \leq \frac{\gamma_\text{param}}{u_\text{max}} \leq 0.1$$

Timescales:
$$\tau_\text{CDI} = \frac{\varepsilon^2}{\gamma_c} = 0.065\;\text{TU}, \qquad \frac{\tau_\text{conv}}{\tau_\text{CDI}} = \frac{R/U_b}{\tau_\text{CDI}} = 4.6 \;\checkmark$$

CDI resharpens $4.6\times$ faster than drop-scale convective strain.

**RIGHT COLUMN: Weber number and case matrix**

$$We = \frac{\rho U_b^2 R}{\sigma} = \frac{R}{\sigma} \quad (U_b=1,\;\rho=1) \qquad\Rightarrow\qquad \sigma = \frac{0.3}{We}$$

| $We$ | $\sigma$ | Regime |
|:---:|:---:|---|
| $1$ | $0.300$ | ST dominates; drop stays spherical |
| $10$ | $0.030$ | Inertia $\approx$ ST; moderate deformation |
| $100$ | $0.003$ | Strong deformation; CDI stress test |
| $730$ | $4\times10^{-4}$ | ST negligible; free deformation (reference) |

Primary target: $We = 1$, exact curvature $\kappa = 2/R = 6.67$ known analytically.
Any curvature error amplifies into a visible spurious velocity.

**[NOTES]**
Early runs used $\varepsilon = 0.05$ ($\varepsilon/\Delta_\text{GLL} = 2.3$, under-resolved in $x,z$) and
$\gamma_\text{param} = 0.015$ ($\gamma/u_\text{max} = 0.010$, at lower bound). Both corrected before restart cases.

=========================================================
## SLIDE 3: The capillary timestep constraint
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: Why explicit CSF is stiff**

Capillary wave stability limit:
$$\Delta t_\text{cap} = \sqrt{\frac{\rho\,\Delta x^3}{2\pi\,\sigma}}$$

SEM complication: which $\Delta x$?

GLL nodes cluster near element edges. Neko's CFL controller uses the minimum
nodal spacing, not the element average. Back-derived from observed $\Delta t$:

$$\Delta x_\text{eff} = \frac{\Delta t_\text{obs}\cdot u_\text{max}}{\text{target\_cfl}} = \frac{0.00130 \times 1.341}{0.2} = 0.0087$$

$$\Delta x_\text{eff} \approx 0.55\;\Delta_\text{GLL} \qquad \left(\Delta_\text{GLL} = 0.016\right)$$

Factor $\sim 2$ in $\Delta x$ gives factor $\sim 5$ in $\Delta t_\text{cap}$ (cubic).

The naive estimate predicted We=10 stable: $\Delta t/\Delta t_\text{cap} = 0.28$.
Corrected estimate: $\Delta t/\Delta t_\text{cap} = 0.69$; run blew up.

**RIGHT COLUMN: All cases compared** (using $\Delta x_\text{eff} = 0.0087$)

| Case | $We$ | $\sigma$ | $\Delta t_\text{cap}$ (TU) | $\Delta t$ (TU) | $\Delta t/\Delta t_\text{cap}$ | Outcome |
|---|:---:|:---:|:---:|:---:|:---:|---|
| v1, cfl=0.4 | $1$ | $0.30$ | $0.00059$ | $0.0043$ | $\mathbf{7.2}$ | Blown up |
| v2, cfl=0.2 | $1$ | $0.30$ | $0.00059$ | $0.0013$ | $\mathbf{2.2}$ | Blown up |
| We=10, cfl=0.2 | $10$ | $0.030$ | $0.00188$ | $0.0013$ | $\mathbf{0.69}$ | Blown up |
| We=100 (planned) | $100$ | $0.003$ | $0.00595$ | $0.0013$ | $\mathbf{0.22}$ | Predicted stable |
| v4, We=730 | $730$ | $4\times10^{-4}$ | $0.071$ | $0.003$ | $\mathbf{0.04}$ | Completed |

Empirical safety criterion:
$$\frac{\Delta t}{\Delta t_\text{cap}} \lesssim 0.5$$

Required target\_cfl for 50\% margin:

| $We$ | Required target\_cfl |
|:---:|:---:|
| $1$ | $0.046$ |
| $10$ | $0.144$ |
| $100$ | $\checkmark$ current $0.2$ fine |

**[NOTES]**
v4 is stable because $We = 730$ makes $\sigma$ negligible; CSF force is tiny regardless of $\Delta t$.

=========================================================
## SLIDE 4: Three blow-ups: CDI intact, CSF unstable
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: Diagnostic trace (We=10)**

| $t$ | $\kappa_\text{rms}$ | $\varphi_\text{max}$ | $u_\text{max}$ |
|:---:|:---:|:---:|:---:|
| $20.002$ | $6.1$ | $\mathbf{0.981}$ | $1.341$ |
| $20.132$ | $6.7$ | $\mathbf{0.986}$ | $1.341$ |
| $20.263$ | $16.5$ | $\mathbf{0.987}$ | $1.426$ |
| $20.317$ | $27.6$ | $\mathbf{0.988}$ | $1.711$ |
| $20.358$ | $52.0$ | $\mathbf{0.989}$ | $2.203$ |
| $20.389$ | $93.4$ | $\mathbf{0.991}$ | $3.262$ |
| $20.435$ | $158$ | $\mathbf{0.998}$ | $4.907$ |

$\varphi_\text{max} < 1$ throughout: CDI maintains the tanh profile
even as $\kappa_\text{rms}$ grows $25\times$ and $u_\text{max}$ triples.

**RIGHT COLUMN: Mechanism and all-case summary**

Explicit CSF feedback loop:
$$\kappa \uparrow \;\Rightarrow\; F_\text{ST} = \sigma\kappa\underbrace{|\nabla\varphi|}_{\approx 1/2\varepsilon} \uparrow \;\Rightarrow\; \Delta u = F_\text{ST}\,\Delta t \uparrow \;\Rightarrow\; \kappa \uparrow$$

Saturation at $\kappa_\text{rms} \approx 150$–$180$:
$$r_\text{local} = \frac{2}{\kappa_\text{rms}} \approx \frac{2}{170} = 0.012 < \varepsilon = 0.07$$

Features become sub-interface in scale; the numerical scheme saturates.

All three blow-ups show the same pattern:

| Case | $\Delta t/\Delta t_\text{cap}$ | $\kappa$ plateau | $\varphi_\text{max}$ at plateau | Stable duration |
|:---:|:---:|:---:|:---:|:---:|
| v1 | $7.2$ | $\sim 170$ | $> 1$ (late) | $1.55$ TU |
| v2 | $2.2$ | $\sim 180$ | $> 1$ (late) | $0.40$ TU |
| We=10 | $0.69$ | $\sim 155$ | $< 1$ always | $0.44$ TU |

Conclusion: instability is in the **explicit CSF integrator**, not CDI and not the turbulent flow.

**[NOTES]**
$\varphi_\text{max} > 1$ appears only after the flow is already fully diverged ($u_\text{max} \approx 26$).
CDI is working correctly throughout. This means improving CDI parameters alone
will never fix the blow-up; the issue is fundamentally in the CSF timestep.

=========================================================
## SLIDE 5: Current status and path forward
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: Laminar $We=1$: seeding test (completed)**

Same parameters ($We=1$, $\sigma=0.3$, $\varepsilon=0.07$, cfl=0.2),
Reichardt IC with no velocity perturbations.

$$\frac{\Delta t}{\Delta t_\text{cap}} \approx 2.9 \quad \text{(worse than turbulent v2!)}$$

Result: blown up at $t \approx 0.90$ TU. Same CDI-intact signature.

| $t$ | $\kappa_\text{rms}$ | $\varphi_\text{max}$ | $u_\text{max}$ |
|:---:|:---:|:---:|:---:|
| $0.000$ | $6.10$ | $0.981$ | $1.149$ |
| $0.173$ | $6.10$ | $0.985$ | $1.207$ |
| $0.520$ | $7.30$ | $0.987$ | $1.233$ |
| $0.840$ | $22.0$ | $0.990$ | $2.180$ |
| $0.893$ | $158.8$ | $\mathbf{0.992}$ | $14.70$ |

Seeding is **not necessary** — numerical round-off alone triggers the instability.
Turbulent fluctuations roughly halve the stable duration:

| Case | $\Delta t/\Delta t_\text{cap}$ | Stable duration |
|:---:|:---:|:---:|
| Turbulent We=10 | $0.69$ | $0.44$ TU |
| Laminar We=1 | $2.9$ | $0.90$ TU |

Laminar survived longer despite a worse ratio — absence of turbulent seeding
is the only explanation. But still blew up.

**RIGHT COLUMN: Cases and path forward**

All runs so far, turbulent or laminar, with $\Delta t/\Delta t_\text{cap} > 0.5$ blew up:

| Case | $We$ | $\Delta t/\Delta t_\text{cap}$ | Outcome |
|:---:|:---:|:---:|:---:|
| v1 | $1$ | $7.2$ | Blown up |
| v2 | $1$ | $2.2$ | Blown up |
| We=10 | $10$ | $0.69$ | Blown up |
| Laminar | $1$ | $2.9$ | Blown up at $0.90$ TU |

Open question: is CDI maintaining accurate normals and curvature?

$\varphi_\text{max} < 1$ only confirms interface amplitude. Inaccurate $\hat{\mathbf{n}}$
gives inaccurate $\kappa = -\nabla\cdot\hat{\mathbf{n}}$, which could contribute to CSF instability.

Immediate next step: $\sigma = 0$ turbulent restart (running)
- CDI only, no CSF force on momentum
- If $\kappa_\text{rms}$ stable near $6.67$ for 5 TU: CDI normals are accurate
- Then blow-ups are definitively a CSF timestep issue

If CDI quality confirmed, next case: We=100 ($\sigma=0.003$, $\Delta t/\Delta t_\text{cap} = 0.22$)

Open question: semi-implicit CSF
$$\left(\mathbf{I} - \Delta t\,\mathcal{L}_\text{CSF}\right)\mathbf{u}^{n+1} = \text{RHS}$$
Removes $\Delta t_\text{cap}$ constraint entirely; requires modifying Neko's momentum solver.

**[NOTES]**
Laminar survived twice as long as turbulent We=10 (0.90 vs 0.44 TU) despite worse ratio (2.9 vs 0.69).
Turbulent fluctuations seed the instability faster but are not necessary.
The sigma=0 test is the key missing piece before drawing a final conclusion.
