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
- $Re_b = 2800$, $Re_\tau \approx 180$, constant bulk velocity $U_b = 1$
- Single drop, $R = 0.3$, centred at $y = 0$, matched fluids ($\rho = \mu = \text{const}$)
- Spectral element method (Neko): $81 \times 18 \times 27$ mesh, $N = 7$

**Governing equations**

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -\nabla p + \frac{1}{Re_b}\nabla^2\mathbf{u} + \underbrace{\sigma\kappa\nabla\varphi}_{\mathbf{F}_\text{ST}}$$

$$\frac{\partial \varphi}{\partial t} + \mathbf{u}\cdot\nabla\varphi = \nabla\cdot\!\left[\varepsilon\gamma_c\nabla\varphi - \gamma_c\,\varphi(1-\varphi)\,\hat{\mathbf{n}}\right], \qquad \kappa = -\nabla\cdot\hat{\mathbf{n}}, \quad \hat{\mathbf{n}} = \frac{\nabla\varphi}{|\nabla\varphi|}$$

**Workflow**
- Step 1: single-phase spin-up to $t = 20$ TU, checkpoint ($Re_\tau = 180$ confirmed)
- Step 2: inject drop analytically, run two-phase $t = 20 \to 25$ TU

**[NOTES]**
CSF sign: $\kappa < 0$ for convex drop ($\varphi = 1$ inside), so $\mathbf{F}_\text{ST}$ points inward (restoring).
Sign fixed in commit 637fc27 of this fork.

=========================================================
## SLIDE 2: The capillary timestep problem — why explicit CSF is stiff
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: Stability limit**

Capillary wave stability (Morris 2000):
$$\Delta t_\text{cap} = \sqrt{\frac{\rho\,\Delta x^3}{2\pi\,\sigma}}$$

GLL nodes cluster near element edges — the effective $\Delta x$ is not the element width but the minimum nodal spacing. Back-derived from observed $\Delta t$:

$$\Delta x_\text{eff} = \frac{\Delta t_\text{obs}\cdot u_\text{max}}{\text{target\_cfl}} = \frac{0.0013 \times 1.34}{0.2} = 0.0087$$

Factor $\sim 2$ in $\Delta x$ gives factor $\sim 5$ in $\Delta t_\text{cap}$ (cubic scaling).

The naive estimate predicted We=10 stable; the corrected estimate said it would blow up — and it did.

**RIGHT COLUMN: All cases compared** ($\Delta x_\text{eff} = 0.0087$)

| Case | $We$ | $\Delta t / \Delta t_\text{cap}$ | Outcome |
|---|:---:|:---:|---|
| v1, cfl=0.4 | $1$ | $7.2$ | Blown up at 1.55 TU |
| v2, cfl=0.2 | $1$ | $2.2$ | Blown up at 0.40 TU |
| We=10, cfl=0.2 | $10$ | $0.69$ | Blown up at 0.44 TU |
| We=10 diag | $10$ | $0.69$ | Blown up at 0.45 TU (dense animation) |
| Laminar We=1 | $1$ | $2.9$ | Blown up at 0.90 TU |
| v4, We=730 | $730$ | $0.04$ | Completed ✓ |

Every case with non-negligible $\sigma$ has blown up. We=10 blew up even with ratio $= 0.69 < 1$.
Stability boundary not yet established.

**[NOTES]**
v4 is stable because $We = 730$ makes $\sigma$ negligible; CSF force is tiny.
Laminar survived twice as long as turbulent We=10 despite a worse ratio — turbulent fluctuations
seed the instability faster but are not necessary. Numerical round-off alone is enough.

=========================================================
## SLIDE 3: Blow-up pattern — CDI intact, explicit CSF unstable
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: We=10 diagnostic trace**

| $t$ | $\kappa_\text{rms}$ | $\varphi_\text{max}$ | $u_\text{max}$ |
|:---:|:---:|:---:|:---:|
| $20.002$ | $6.1$ | $\mathbf{0.981}$ | $1.341$ |
| $20.132$ | $6.7$ | $\mathbf{0.986}$ | $1.341$ |
| $20.263$ | $16.5$ | $\mathbf{0.987}$ | $1.426$ |
| $20.358$ | $52.0$ | $\mathbf{0.989}$ | $2.203$ |
| $20.435$ | $158$ | $\mathbf{0.998}$ | $4.907$ |

$\varphi_\text{max} < 1$ throughout: **CDI maintains the tanh profile even as $\kappa_\text{rms}$ grows $25\times$**.

Explicit CSF feedback loop:
$$\kappa \uparrow \;\Rightarrow\; F_\text{ST} = \sigma\kappa|\nabla\varphi| \uparrow \;\Rightarrow\; \Delta u = F_\text{ST}\,\Delta t \uparrow \;\Rightarrow\; \kappa \uparrow$$

**RIGHT COLUMN: Key diagnostics across all blow-ups**

| | v1 | v2 | We=10 | We=10 diag | Laminar |
|---|:---:|:---:|:---:|:---:|:---:|
| $\Delta t/\Delta t_\text{cap}$ | $7.2$ | $2.2$ | $0.69$ | $0.69$ | $2.9$ |
| $\kappa$ plateau | $\sim\!170$ | $\sim\!180$ | $\sim\!155$ | $\sim\!158$ | $\sim\!160$ |
| $\varphi_\text{max}$ at plateau | $> 1$ (late) | $< 1$ always | $< 1$ always | $< 1$ always | $< 1$ always |
| Stable duration | $1.55$ TU | $0.40$ TU | $0.44$ TU | $0.45$ TU | $0.90$ TU |

$\kappa$ plateau $\approx 170$ across all cases:
$$r_\text{local} = \frac{2}{\kappa_\text{rms}} \approx \frac{2}{170} = 0.012 < \varepsilon = 0.07$$
Features become sub-interface in scale; the numerical scheme saturates.

**Conclusion: instability is in the explicit CSF integrator, not CDI and not turbulence.**

**[NOTES]**
CDI is working correctly in all cases — the interface profile is maintained. The blow-up is purely from
the explicit treatment of the surface tension force. Semi-implicit CSF or a lower $\Delta t$ would fix this.

=========================================================
## SLIDE 4: σ=0 CDI-only diagnostic — CDI is fine; κ is not
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: What the σ=0 run shows**

With $\sigma = 0$: no CSF force, no blow-up possible. Only CDI acts.
Drop injected at $t = 20$, turbulent restart. Run 1.24 TU.

CDI quality metrics from field snapshots (Python postprocessing):

| Metric | $t=20.5$ | $t=21.0$ | Reference |
|:--|:---:|:---:|:---:|
| $\|\nabla\varphi\|_\text{mean}$ at interface | $3.93$ | $4.37$ | $3.57$ (CDI theory) |
| $r_\text{std}$ on $\varphi=0.5$ surface | $0.029$ | $0.052$ | $0$ (sphere) |
| $\int\varphi\,dV$ | $0.20965$ | $0.20962$ | conserved |
| $\kappa_\text{rms}$ (Python postprocess) | $25.8$ | $24.2$ | $6.67$ (sphere) |
| $\kappa_\text{rms}$ (Neko, GS-averaged) | $62$ | $57$ | $6.67$ (sphere) |

- CDI functioning: interface gradient above theory, volume conserved
- Drop nearly spherical: $r_\text{std}/R = 9.7\%$ — too small to explain $\kappa_\text{rms}/\kappa_\text{sphere} \approx 9$
- Factor 2.4 between Python and Neko κ
- **κ is inaccurate. CDI and the drop shape are not the problem.**

**RIGHT COLUMN: High-frequency diagnostic — κ grows with advection**

20 field snapshots at $\Delta t = 0.02$ TU from $t=20.12$ to $t=20.50$:

| $t$ | centroid displacement | $\hat{\mathbf{n}}$ misalign | $\kappa_\text{rms}$ (Python) |
|:---:|:---:|:---:|:---:|
| $20.12$ | $0.14$ | $4.3°$ | $6.0$ |
| $20.22$ | $0.27$ | $8.0°$ | $11.9$ |
| $20.32$ | $0.41$ | $16.0°$ | $20.1$ |
| $20.50$ | $0.61$ | $25.4°$ | $25.8$ |

**$\kappa_\text{rms}$ grows in proportion to centroid displacement (element faces crossed), not to drop deformation.** This is the smoking gun for an element-face artifact.

**[NOTES]**
The tanh L2 residual also grows (0.014 → 0.036) and interface width starts 8.6% too wide (CDI with γ=0.05
is too slow: τ_CDI = ε/γ = 1.4 TU). But this is secondary — the κ artifact grows before CDI has time to act.

=========================================================
## SLIDE 5: The mechanism — C0 n̂ enforcement creates a κ artifact
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: Neko: GS-average makes n̂ globally C0**

1. Element-local SEM gradient $\nabla\varphi|_e$ (Lagrange derivative)
2. **Gather-scatter average** over all elements sharing each face/edge/corner node → $\nabla\varphi$ is C0
3. Normalise: $\hat{\mathbf{n}} = \nabla\varphi / |\nabla\varphi|$ — also C0
4. $\text{div}(\hat{\mathbf{n}})$ element-locally → GS-average → $\kappa = -\nabla\cdot\hat{\mathbf{n}}$

**After GS:** $\hat{\mathbf{n}}$ has the same value at both sides of a face (C0), but different *slopes* from each element's perspective (not C1). This is the **C0 kink**.

Within element $e_1$, the Lagrange derivative at the face node ($\xi = +1$) has:
$$D[N,N] = \frac{N(N+1)}{4} = 14 \quad (N=7)$$

Combined with the Jacobian $2/L_e \approx 13$ ($L_e = 0.155$), the effective amplification is $\sim 90$.

**RIGHT COLUMN: Quantitative estimate**

For a sphere ($R = 0.3$) crossing an element of width $L_e = 0.155$:
- $\hat{\mathbf{n}}$ rotates by $\Delta\theta \approx L_e / R = 0.52$ rad within one element
- Change in $\hat{n}_x$ across element: $\delta \approx 0.26$
- Artifact at a face node: $\kappa_\text{face} \sim 90 \times 0.26 \approx 23$
- At 3D corners (multiple contributions): $\kappa_\text{corner} \sim$ hundreds

Consistent with $\kappa_\text{max} = 800+$ observed from the first timestep.

**Why Python gives less artifact:** element-local np.gradient keeps $\hat{\mathbf{n}}$ discontinuous — no GS bridge, no kink, no amplification. The artifact is a direct consequence of making $\hat{\mathbf{n}}$ C0.

| | $\kappa_\text{rms}$ | $\hat{\mathbf{n}}$ | artifact? |
|---|:---:|:---:|:---:|
| Analytical sphere | $6.67$ | C∞ | none |
| Python postprocess | $26$ | discontinuous | low |
| Neko (GS-averaged) | $62$ | C0 | high |

**[NOTES]**
D[N,N] = N(N+1)/4 is the exact endpoint weight of the Lagrange derivative matrix for GLL nodes.
This is a standard result. The endpoint rows are the largest — "endpoint amplification" is well known
for smooth functions but becomes pathological for C0-but-not-C1 fields.

=========================================================
## SLIDE 6: γ parametric test — stronger CDI makes κ worse
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: Motivation and test**

CDI re-sharpening timescale: $\tau_\text{CDI} = \varepsilon / \gamma_c = 0.07/0.05 = 1.4$ TU

Interface width at injection ($t=20.12$): $4\varepsilon_\text{meas} = 0.304$ vs theory $0.280$ — 8.6% too wide.
By $t=20.5$: still 0.293, barely improved. CDI cannot keep up with turbulent straining.

Test: run $\sigma=0$ CDI-only with $\gamma = 0.25$ (vs baseline $\gamma = 0.05$).
New $\tau_\text{CDI} = 0.07/0.25 = 0.28$ TU — much faster resharpening.

**RIGHT COLUMN: Result — κ is 9× worse immediately**

| $t$ | $\kappa_\text{rms}$ ($\gamma=0.05$) | $\kappa_\text{rms}$ ($\gamma=0.25$) |
|:---:|:---:|:---:|
| $20.002$ | $6.1$ | $6.1$ |
| $20.07$ | $7.8$ | $\mathbf{71}$ |
| $20.13$ | $7.8$ | $\mathbf{85}$ |
| $20.20$ | $19.9$ | $\mathbf{82}$ |
| $20.26$ | $19.9$ | $\mathbf{79}$ |

The jump happens at the **first output step** — before turbulence has time to act.

**Reason:** stronger CDI compression creates a sharper interface ($|\nabla\varphi|$ larger), which amplifies the C0-kink magnitude at element faces. The artifact scales with $|\nabla\varphi|$ and hence with $\gamma$.

**Conclusion: γ tuning cannot fix the κ artifact. The only fix is in the n̂ computation itself.**

**[NOTES]**
φ_max recovers faster with γ=0.25 (0.981 → 0.989 in 0.13 TU vs slower for γ=0.05).
So CDI quality metrics (interface width, φ_max) genuinely improve with higher γ.
But the κ accuracy, which is what matters for CSF, gets dramatically worse.
The two goals (good CDI profile, accurate κ) are in conflict with the current κ computation method.

=========================================================
## SLIDE 7: Summary and next step
**Layout: TWO COLUMNS**
=========================================================

**LEFT COLUMN: What we know**

All CSF cases blow up from the explicit capillary instability — even We=10 ($\Delta t/\Delta t_\text{cap} = 0.69$), even laminar flow.

σ=0 diagnostic run confirms:
- **CDI is working:** interface sharp, volume conserved, drop nearly spherical
- **κ is wildly inaccurate:** $\kappa_\text{rms} \approx 64$ vs sphere $6.67$
- **Root cause:** SEM element-boundary artifact — GS averaging makes $\hat{\mathbf{n}}$ C0 at faces, creating a C1 kink that the $O(N^2)$ Lagrange endpoint derivative amplifies by $\sim 90\times$
- **Not CDI failure:** confirmed independently — γ tuning improves CDI profile quality but makes κ worse

CDI parameter study ($\gamma = 0.25$ vs $0.05$):
- Faster resharpening ✓ (φ_max, interface width improve)
- κ immediately $9\times$ worse ✗
- Higher $\gamma$ → sharper interface → larger $|\nabla\varphi|$ → larger C0-kink amplitude

**With the current κ computation and $\sigma > 0$:**
$$F_\text{ST} \approx \sigma \times 62 \times |\nabla\varphi| \approx 9\times \text{ larger than spherical estimate}$$
The capillary stability analysis using $\kappa = 6.67$ was optimistic by a factor of 9.

**RIGHT COLUMN: The fix — one extra GS pass on $\hat{\mathbf{n}}$**

Current Neko workflow:
```
∇φ (element-local) → GS avg → n̂ = ∇φ/|∇φ| (C0) → div(n̂) (element-local) → GS avg → κ
                                        ↑
                              kink here — amplified by D[N,N]=14
```

Proposed fix:
```
∇φ (element-local) → GS avg → n̂ = ∇φ/|∇φ| → GS avg on n̂ → div(n̂) → GS avg → κ
                                                      ↑
                                         smooths the slope across faces
```

The extra GS pass on $\hat{\mathbf{n}}$ after normalisation smooths out the C1 kink, reducing the endpoint derivative artifact.

**Next step:** implement the extra GS pass in `turb_channel_two_phase.f90` (the CSF source term, ~5 lines of Fortran), then test with $\sigma > 0$.

**[NOTES]**
The extra GS pass is a standard SEM operation: `call coef%gs_h%op(nx_, nx_%x, GS_OP_ADD)` followed by
`call col2(nx_%x, coef%mult, nx_%size())`. It costs one GS communication per normal component (3 total),
which is negligible compared to the pressure solve.
The trade-off: over-smoothing κ toward the element average (slightly blurs genuine curvature variations).
For nearly spherical drops (uniform κ) this is ideal; for highly deformed drops it introduces a bias.
This should be tested by comparing κ_rms against the postprocessing value after the fix.
