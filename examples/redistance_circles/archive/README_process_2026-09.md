# redistance_circles: how we got here (archived process record, 2026-09-22 to 2026-09-29)

**This is an archive.** It is the case README as it stood on 2026-09-29, before it was rewritten
to describe the case on its own. It keeps the full process: the line-by-line audit against the
paper, the ε = H/N configuration and everything measured with it (the drift, node capture, the
apex failure, the N/4 arms, the eliminated hypotheses), the validation and rebuild, the reading
of the authors' public code, and the dated history. Section numbers below refer to this file.

The current description of the case is [`../README.md`](../README.md). The run-level record of the
2026-09-29 study of the authors' code is `../logs/saini_case/results.txt` (gitignored).

---

# redistance_circles — Saini et al. (2026) §4.4 in Neko

**This file is the record of the §4.4 reproduction.** It holds:
- the paper, equation by equation, against our code (§3);
- what the paper leaves unstated, and what we chose (§4);
- the Neko implementation of the paper's scheme and how it was validated (§5–§6);
- the results against Table 2 and Figs. 12–13 (§7);
- why they differ (§8);
- what is eliminated (§9) and what is open (§10).

Other documents summarise this page and point here. Read it before any work on the case.

## 1. Status

**Not reproduced. Every equation §4.4 prints is implemented as printed (§3), yet
$E_r$ is 1.4–3.4× above the paper in every cell.**

- $H{=}1/20$, $N{=}3$ fails at a cone apex, from the skewed IC only. The
  failure starts when the relaxation front collapses on the circle centre, at
  $\tau\approx1$ (§10).
- Started from the exact solution itself, the scheme drifts away from it. No
  parameter stops the drift.
- The drift ends when grid nodes beside the zero set are caught at $\psi=0$.
  $\varepsilon$ only sets its clock: at $\varepsilon\times3$ the same end state
  arrives three times later (§8.5).
- The paper does not state which $\varepsilon$ Eq. (46) uses in §4.4, nor whether
  any step holds the interface in place. Those two gaps (H8, H9′) are the likely
  explanations. **They are not pursued** (decision, 2026-09-28; §10).

What is settled, so a new session does not redo it:

| settled | how | where |
|---|---|---|
| The code implements Eqs. (7), (18), (22)–(36), (42)–(46), (83), (84) as printed | line-by-line audit against the rendered PDF, 2026-09-25 | §3 |
| The time integrator is not the gap | At §4.4's own settings, Saini's BDF2/EXT2 and BDF3/EXT3 in Neko are within 1.9% of SSP-RK3 in every holding Table 2 cell, at every $\tau$ to 24 (§7.4); BDF2 is within 2.2% over the whole Fig. 12 grid (§7.2). (With the weaker $N_{svv}=N/4$ at $H{=}1/5$, $N{=}7$, only RK3's split damping keeps the apex from failing, §7.5.) Step-map spectra: every integrator reproduces the operator's dynamics | §7.4 |
| The Neko code is the numpy replica's scheme | node by node to 1.3e-13 (RK3) and 4.7e-13 (BDF2, BDF3, exact IC) from the same $\psi_0$ | §6 |
| The metric is Eq. (84) as printed, $\int\psi_e$ | their Fig. 13 integrates to their $E_r$ | §7.3 |
| The gap is in our field, and it is a drift | from $\psi_e$, $E_r$ rises from 0 in every cell; the drift is the Eq. (31) SVV acting beside the zero set | §8 |
| $\varepsilon$ is a clock for the drift, not its size | from $\psi_e$ at $H{=}1/10$, $N{=}3$, $\varepsilon\times3$ at $\tau$ equals $\varepsilon$ at $\tau/3$ within 12% over $\tau=6$–72, down to which nodes are caught; the end states agree to 0.2% | §8.5 |
| §4.4 is not band-limited, and a band would not close the gap | the paper: §4.4 re-distances "in the entire domain" to $\tau=6$, and $E_r$ is "global" (p. 19). Their Fig. 13: 52–97% of $\int\lvert e\rvert$ lies beyond $2.5H$. A $2.5H$ mask leaves the drift as it is (H13, 2026-09-28) | §9 |
| The time-converged scheme fails at the apex with $N_{svv}=N/4$ ($H{=}1/5$, $N{=}7$) | BDF2 fails at every $\Delta\tau$ from $5\times10^{-4}$ to $1.25\times10^{-4}$; RK3+Lie holds at $\Delta\tau\ge2.5\times10^{-4}$ and fails at $1.25\times10^{-4}$ | §7.5 |
| The Eq. (44) solver has no implementation error | Independent audit, 2026-09-28, without the replica: the rendered paper, the SVV filter from Legendre modes, the implicit step against the unsymmetrised printed system, the gradient against exact gradients, Neko's own time integrator, and a manufactured solution | §6.1 |
| No shared misreading of a printed equation | A design from the rendered paper alone, written before the code was read, matches it term by term. A rebuild from Neko's own routines (`ax_helm` on a filtered space, `gmres`, `adv_no_dealias_t`, `opgrad`, `rhs_maker`) reproduces every logged $E_r$ and the $\tau=6$ fields to ≤5.4e-11 in every cell that holds, and the $H{=}1/20$, $N{=}3$ apex failure (2026-09-28) | §6.2 |

## 2. The problem, as printed (JCP pp. 18–19)

Eq. (44) **standalone**: no flow, no phase field, and no Neko time stepping. The
whole result is computed in the user `initialize` hook. The notation is theirs
transposed: their $\phi$ is our $\psi$ (`../../CDI_METHOD.md` §1).

$$\frac{\partial\psi}{\partial\tau} + \mathbf{w}\cdot\nabla\psi = \operatorname{sgn}(\psi) + S_{vv}(\psi), \qquad \mathbf{w} = \operatorname{sgn}(\psi)\,\frac{\nabla\psi}{|\nabla\psi|}, \qquad \operatorname{sgn}(\psi) = \tanh\!\Big(\frac{\psi}{2\varepsilon}\Big)$$

| | as printed | here |
|---|---|---|
| domain | $\Omega=[-2,2]^2$ | the same, one element thick in $z$, non-periodic in $x$ and $y$ |
| exact solution, Eq. (83) | signed distance to two circles, $r=1$, centres $(\pm a,0)$, $a=0.7$, negative inside; a corner branch where both radial projections land inside the other circle | `circle_distance` |
| initial condition, Eq. (83) | $\psi_0=((x-1)^2+(y-1)^2+0.1)\,\psi_e$ | `skewed_ic` |
| meshes | $H\in\{1/5,1/10,1/20\}$, $N=3$…8 | `box20/40/80.nmsh` |
| $\Delta\tau$, $\tau_{end}$ | $5\times10^{-4}$ fixed, to $\tau=6$ ("quasi-steady") | the same, 12000 steps |
| SVV | $N_{svv}=N/6$, $c_0=2$ "for all cases in this study" | the same |
| metric, Eq. (84) | $E_r=\int\lvert\psi-\psi_e\rvert\,d\Omega\,/\int\psi_e\,d\Omega$ | the same (`error_report`) |
| results | Table 2 (four cells), Fig. 12 ($E_r(6)$ over the grid), Fig. 13 ($\lvert e\rvert$ at $N{=}8$) | §7 |

Non-periodic boundaries are load-bearing. On $y=0$, $\psi_0$ is +3.03 at $x=-2$
against +0.63 at $x=+2$, so a periodic seam would be carried inward at unit speed.

## 3. The paper against the code, equation by equation

Transcribed from the rendered JCP pages on 2026-09-25 (page numbers are PDF
pages) and checked against the Fortran. The PDF is
`../../references/saini_tomboulides_2026.pdf` (gitignored, CC BY-NC-ND).
`paper-lookup` answers questions about it.

| paper | what it prints | code | agrees |
|---|---|---|---|
| Eq. (2)–(3), p. 3 | $\partial_tu+\mathbf c\cdot\nabla u=f+\mu\nabla\cdot\mathbf Q\nabla u$, with $\mu$ **outside** the divergence; the equation is pre-multiplied by $1/\mu$ | the SVV step | yes |
| Eq. (4), p. 3 | weak form, "periodic domains and domains with zero flux boundaries" | natural boundary conditions; no BC is imposed. The `.case` files declare zero-flux Neumann on `psi`, but the Eq. (44) solve bypasses Neko's scalar solver and never reads it | yes |
| Eq. (7), p. 4 | $\mathbf v^T[\mathbf D_\mu^{-1}(\mathbf M\dot u+\mathbf Cu-\mathbf Mf)+\mathbf S_{vv}u]=0 \Rightarrow \mathbf M\dot u=-\mathbf Cu+\mathbf Mf-\mathbf D_\mu\mathbf S_{vv}u$; $\mathbf C$ is Galerkin with GLL quadrature | $\mathbf D_\mu$ left-multiplies the assembled $\mathbf S_{vv}$; $\mathbf C$ is `conv1`, then $\mathbf B$, gather-scatter, $\mathbf B^{-1}$ | yes |
| Eq. (18), p. 5 | $(G_{lm})_{\hat k\hat k}=\rho_i\rho_j\rho_k\{\tilde G_{lm}J_e\}$ | Neko `G11`, `G22`; with the filter set to identity, `svv_local` equals Neko's stiffness matrix exactly (O2 audit); the off-diagonal $G_{lm}$ are checked to be zero | yes |
| Eq. (22), p. 6 | $\hat B_{ij}=L_{j-1}(r_i)$, "the basis matrix (or the Vandermonde matrix)" | Neko `elementwise_filter`, `nonBoyd` = plain Legendre | yes |
| Eq. (24), p. 6 | $\hat Q_{ii}=(k/N)^{N_{svv}}$, $k=i-1$; the exponent is $N_{svv}$ | `svv_init`: `trans(i)`, `nsvv = N/ratio` | yes |
| Eq. (25)–(28), p. 6–7 | $\tilde D_l=\hat B\hat Q\hat B^{-1}\hat D$ (derivative first, then filter); $\mathbf S^e_{vv}=\tilde D^TG\tilde D$, so $\hat Q$ enters once on each side | `svv_local`: $D$, $F$, $\nu G$, $F^T$, $D^T$ | yes |
| Eq. (29), (31), p. 7 | $\mu=c_0\lvert\mathbf c\rvert\mathbb H/N$ at the GLL points, $\mathbb H=2J_e^{1/d}$ | $\nu=c_0H/N$ (`svv_init`) times $\mathbf D_\mu=\lvert\operatorname{sgn}\psi\rvert=\lvert\mathbf w\rvert$. In 2D, $2J_e^{1/2}=H$ exactly | yes |
| Eq. (33), p. 7 | $\mathbf D^e_\mu\mathbf S^e_{vv}$ per element, then assembled | equal to left-multiplying the assembled operator, because $\lvert\operatorname{sgn}\psi\rvert$ is continuous and $\mathbb H$ is uniform | yes |
| Eq. (34)–(35), p. 7–8 | BDF$k$/EXT$k$, "BDF2/EXT2 … for all studies"; $\mathbf H=\tfrac{b_0}{\Delta t}\mathbf M+\mathbf D_\mu\mathbf S_{vv}$; the explicit terms $a_j(\mathbf Cu-\mathbf Mf)^{n+1-j}$ are extrapolated | `redistance_standalone` (BDF2/EXT2, BDF1/EXT1 on step 1) and `svv_step_eq31` ($\mathbf D_\mu$ from $\psi^n$), §5 | yes |
| p. 9; pp. 12–13 (before Algorithm 1) | the band is a pseudo-time budget: "iterated in pseudo time up to a limited distance"; $\Delta\tau_{tls}=H/(N{+}1)$, $N_{tls}=2.5H/\Delta\tau_{tls}$. Eqs. (44)–(47) contain no mask or clip | the coupled cases' `rd_niter = ceiling(band*H/dtau)`; not used here | yes |
| Eq. (36), p. 8 | $\epsilon=\xi H$, $\xi=\{1,1.5\}/N$, $H$ the largest element edge | $\varepsilon=H/N$ (`case.cdi.epsilon`) until 2026-09-29 | printed reading; **their code does not use it for Eq. (46)**, §3.1 |
| Eq. (42), p. 8 | $\mathbf n=\nabla\phi/\lvert\nabla\phi\rvert$ | `unit_normal`: gradient averaged at shared nodes | printed; how $\nabla\phi$ is formed is unstated, §4 |
| Eq. (44)–(46), p. 9 | as above | `rd_sgn`; §5 | yes |
| §3.4, p. 12 | "the SVV parameters are defined for each governing equation separately" | Eq. (44) has its own $c_0$ and $N_{svv}$ | yes |
| Eq. (83), p. 18 | the corner-branch condition $(a\mp x)/\sqrt{(a\mp x)^2+y^2}\ge a/r$ (both must hold); the branches $\max(-\sqrt{x^2+(y\pm\sqrt{r^2-a^2})^2})$ and $\min(\sqrt{(x\pm a)^2+y^2})-r$ | `circle_distance`, written without division. At the circle centres, where the printed form is 0/0, both branches give −1 | yes |
| Eq. (84), p. 19 | as above | `ensure_psie`, `error_report` | yes |
| §4.4, p. 18–19 | $\Omega$, $r$, $a$, $H$, $N$, $\Delta\tau$, $\tau$, $N_{svv}$, $c_0$ | the six committed `.case` files | yes |
| §4.4, p. 19 | "advanced to $\tau=6$ … which allows sufficient time for re-distancing in the entire domain"; Eq. (84) is "the global relative error". No band, no $N_{tls}$ | Eq. (44) over all of $\Omega$ to `tau_end`; $E_r$ over $\Omega$ | yes |

Two harmless oddities in the print:
- Eq. (30) reads $N/c_{sf}$ where Eq. (29) implies $N/c_0$.
- §3.4's automatic $\Delta\tau_{tls}=H/(N+1)$ and $N_{tls}=2.5H/\Delta\tau_{tls}$ are
  overridden in §4.4 by the fixed $5\times10^{-4}$ to $\tau=6$.

The paper mentions no dealiasing, over-integration or filtering anywhere. Their code
dealiases (§3.1).

### 3.1 Their code against ours (2026-09-29)

The authors' §4.4 case is public: `nandu90/nekLS_Examples`, branch `jcp`,
`intersectingCircles/` (commit `bceef31`). It runs on their Nek5000 fork,
`nandu90/Nek5000`, branch `nekLS` (commit `5e9b0ae`). Both were read, built and run here
(`logs/saini_case/`).

- **Their build reproduces `plot.py`** L83–85, the values behind Fig. 12 and Table 2, to 12–14
  digits in the seven cells run. That needs one fix: `useric` leaves `func` unset for the
  phase field. Left unset (0 here), the sign guard below freezes every node with $\psi<0$ and
  $E_r$ is 7.7.
- Their case against ours, line by line:

| their code | what it does | ours | same? |
|---|---|---|---|
| `lvlSet.f` `signls` | $\operatorname{sgn}\psi=\tanh(\psi/(2\cdot0.25))$: **a fixed $\varepsilon=0.25$** in Eq. (46). The $\varepsilon=H/N$ of `usrdat3` (`eps_cls`) enters only the phase field's Heaviside | `rd_sgn`, `svv_step_eq31`: $\varepsilon=H/N$ | **no** |
| `conv_tlsr` → `convect_new`; `SIZE` `lxd = 3 lx1/2` | $\mathbf C(\mathbf w)\psi$ **dealiased** on $\lfloor 3(N{+}1)/2\rfloor$ GL points; $\mathbf w$ re-formed from $\psi^n$ each step | `conv1`, GLL collocation | **no** |
| `constrainTLSR` | **sign guard**: a node whose $\psi^n$ disagrees in sign with the phase field ($=\operatorname{sgn}\psi_e$) keeps $\psi^n$ | none | **no** |
| `bdry_tlsr_fix` | flips an inward $\mathbf w$ at the walls | none | no, but it never fires here (bit-identical when removed) |
| `settime_cls` | BDF2 with Nek's EXT3 $[8/3,-7/3,2/3]$ | BDF2/EXT2 | no; ≤1e-4 relative |
| `cggo_cls` | PCG on the non-symmetric $\mathbf H$, `residualTol` 1e-8 | symmetrised CG, 1e-12 | no; below their dump's precision |
| `cls_normals` | $\mathbf n$ from the averaged gradient; $\mathbf n=0$ below $\lvert\nabla\psi\rvert=10^{-12}$ | averaged; floor $10^{-6}$ | yes (floor inert) |
| `svv.f` `setmu_svv`, `axhelm_svv`, `diffFilter1D`; `usrdat3` | $\mu=c_0\lvert\operatorname{sgn}\psi\,\mathbf n\rvert\,2J^{1/2}/N$ left-multiplying the element $\mathbf S^e$; $\hat Q=(k/N)^{N/6}$ (`svvcut` $=N/3$, halved); $c_0=2$ | `svv_init`, `svv_local`, `svv_step_eq31` | yes |
| `useric`, `getexact`, `lserrors.f` `ls_relerr` | Eq. (83) IC and $\psi_e$; $E_r=\int\lvert\psi-\psi_e\rvert/\int\psi_e$ | `skewed_ic`, `circle_distance`, `error_report` | yes |
| `.par`, `ls_init_maxiter` | $\Delta\tau=5\times10^{-4}$, 12000 steps; zero-Neumann walls | the same | yes |

- **What the three differences do** (`logs/saini_case/results.txt`):
  - **$\varepsilon$ is the gap.** With $\varepsilon=0.25$ and nothing else changed, our
    committed scheme gives $E_r(6)$ at 0.73–0.99× `plot.py` in all 18 Fig. 12 cells.
    $H{=}1/20$, $N{=}3$ holds. In the seven cells run in both codes, it equals their code
    with dealiasing off to 5 figures, and the fields to $3\times10^{-7}$ (the
    single-precision floor of their dump).
  - **Dealiasing** adds 12–37% to their $E_r$. Neko's dealiased advection equals their
    `convect_new` to 5 figures.
  - **The guard** acts only once dealiasing moves zero-set nodes: −19% at $H{=}1/20$,
    $N{=}3$, about 1e-5 at $N{=}8$.
- **With all three,** our scheme gives their $E_r$:
  - to 5 figures at $H{=}1/5$, $N{=}7$;
  - at $H{=}1/10$, $N{=}3$, 2.2% low from our own mesh, and exact once $\psi_0$ at the 18
    zero-set nodes is taken from their mesh. The field then agrees to $6\times10^{-8}$.
  - The genbox and genmeshbox coordinates differ at round-off, which flips the sign of
    $\psi_0\approx\pm10^{-15}$ at 9 of those 18 nodes. So at low $N$ their $E_r$ depends at
    the 2% level on the mesh generator's round-off.
- **The whole Fig. 12 grid with their configuration, from our own mesh** (`r2ge/`), as
  $E_r(6)$ ours / `plot.py`:

| | $N{=}3$ | $4$ | $5$ | $6$ | $7$ | $8$ |
|---|---|---|---|---|---|---|
| $H{=}1/5$ | 1.000 | fails (0.59) | 1.000 | 1.000 | 1.000 | fails (0.11) |
| $H{=}1/10$ | 0.978 | 0.995 | 1.000 | 0.999 | 1.000 | 1.000 |
| $H{=}1/20$ | 1.011 | 1.007 | 1.000 | 1.000 | 1.000 | 1.000 |

  - 12 cells equal theirs to 4–5 figures; the other holding cells are within 2.2%.
  - Both failures are the zero-set round-off. Seeded with their $\psi_0$ at the zero-set nodes,
    $H{=}1/5$ gives their value exactly at $N{=}4$ and $N{=}8$. With dealiasing but no guard,
    both hold (0.6–9% from their value).
  - The four Table 2 cells hold.

## 4. What the paper does not state: their code, and what we ship

| not stated | what their code does (2026-09-29, §3.1) | ours | tested |
|---|---|---|---|
| $\varepsilon$ in Eq. (46) for §4.4 | **fixed 0.25** (`signls`; $0.25L$ in the core code, with $L$ the domain size). The author's email gives $\varepsilon=H/N$, which is the phase field's. Asked him on 2026-09-29; open | 0.25 (D4) | yes: the whole gap (§3.1) |
| how $\nabla\phi$ in Eq. (42) is formed | averaged across elements (`cls_normals`; the author confirms) | the same | yes (as before) |
| the boundary condition in §4.4 | zero Neumann (the author confirms); plus a wall flip of inward $\mathbf w$ that never fires here | natural, zero flux | yes: the boundary is passive |
| the time level of $\mathbf D_\mu$ in $\mathbf H$ | $\psi^n$ (`sethlm_ls`) | the same | as before |
| the BDF startup and extrapolation | BDF1 on step 1; BDF2 with Nek's EXT3 | BDF2/EXT2 as printed | EXT2 = EXT3 within 5e-5 relative over the whole Fig. 12 grid |
| dealiasing (the paper mentions none) | on, $\lfloor 3(N{+}1)/2\rfloor$ points, "in all my cases" (the author) | on (D4) | yes: +12–37% in $E_r$ |
| any interface-preserving correction | none; only the sign guard `constrainTLSR` (named in the author's email, not in the paper) | the guard (D4) | yes: matters only with dealiasing, up to 19% at low $N$; makes low-$N$ cells round-off sensitive |
| how their code applies $\mathbf D_\mu$ | left-multiplying the element operator (`axhelm_svv`), $\lvert\mathbf c\rvert=\lvert\operatorname{sgn}\psi\,\mathbf n\rvert$ (the author confirms) | the same | as before |
| $c_0$, $N_{svv}$ | $c_0=2$; `svvcut` $=N/3$, halved internally to $N/6$ (the author) | the same | as before |

Whether §4.4 limits Eq. (44) or $E_r$ to Algorithm 1's $2.5H$ band is **stated**, so
it is not in this table: it does not (§3, p. 19; §9).

## 5. The implementation: Saini's scheme, in Neko

Neko's scalar solver is not used for Eq. (44): its advecting velocity is the
fluid's, an external field (`../../REDISTANCING.md` §9.0). Instead, the user file
runs its own pseudo-time loop built from Neko's operators. One step
$n\to n+1$ is Eqs. (34)–(35) in absolute form, with $\mathbf M=\mathbf B$ diagonal:

$$F^m=\operatorname{sgn}(\psi^m)-\mathbf B^{-1}\mathbf C(\mathbf w^m)\psi^m,\qquad \hat\psi=\tfrac{1}{b_0}\Big(\textstyle\sum_j\beta_j\psi^{n+1-j}+\Delta\tau\sum_j a_jF^{n+1-j}\Big),\qquad \big(\mathbf B+\tfrac{\Delta\tau}{b_0}\mathbf D_\mu(\psi^n)\mathbf S_{vv}\big)\psi^{n+1}=\mathbf B\hat\psi$$

The BDF2/EXT2 coefficients are $b_0=3/2$, $\beta=(2,-\tfrac12)$, $a=(2,-1)$.
Step 1 uses BDF1/EXT1.

**$\mathbf w$ is not a velocity field.**
- $\mathbf w^m$ is re-formed from $\psi^m$ at every level and applied only to $\psi^m$.
  The EXT history holds the products $F^m$, never $\mathbf w$.
- Why it matters: linearised at a distance function, $\mathbf w$'s own
  $\operatorname{sgn}'(\psi)\lvert\nabla\psi\rvert\,\delta$ cancels the source's
  $\operatorname{sgn}'(\psi)\,\delta$, leaving pure advection.
- A $\mathbf w$ decoupled from the $\psi$ it multiplies leaves growth at up to
  $1/2\varepsilon$ on the zero set. That is 15 at $H{=}1/10$, $N{=}3$.
- Measured in the replica: $\mathbf w$ frozen at the IC diverges, and a lag of 0.05 in
  $\tau$ is 6× worse.
- Extrapolating $\mathbf w$, subcycling it (OIFS), or routing it through a
  velocity-field solver would all break this.

**How $\mathbf w$ and $\mathbf C(\mathbf w)$ are formed.**
- The gradient is averaged at shared nodes (`unit_normal`, floor $10^{-6}$).
- $\mathbf w=\operatorname{sgn}(\psi)\,\mathbf n$ (`rd_sgn`).
- $\mathbf C(\mathbf w)\psi$ is Neko's `conv1` (element-local collocation), then $\mathbf B$,
  gather-scatter and $\mathbf B^{-1}$: the non-dealiased GLL Galerkin form.
- With GLL quadrature, $(\mathbf C\psi)_i$ sees only $\mathbf w_i$, so a node on the zero
  set is not moved.
- On this mesh it equals the identity $\operatorname{sgn}(\psi)\lvert\nabla\psi\rvert$ to
  1e-14 (O2 audit), since the shared-node GLL weights are equal.

**The implicit solve.**
- $\mathbf D_\mu=\lvert\operatorname{sgn}\psi^n\rvert$ vanishes on the zero set, and the
  operator is not symmetric.
- The substitution $\psi^{n+1}=\hat\psi+\mathbf D_\mu z$ turns it into
  $(\mathbf D_\mu\mathbf B+\tfrac{\Delta\tau}{b_0}\mathbf D_\mu\mathbf S\mathbf D_\mu)z=-\tfrac{\Delta\tau}{b_0}\mathbf D_\mu\mathbf S\hat\psi$,
  which is symmetric and consistent where $\mathbf D_\mu=0$.
- It is solved by CG to $10^{-12}$ relative (`svv_step_eq31`).
- It does not conserve mass. A distance function does not need to.
- **No node can change sign**: every term of $\psi_\tau$ vanishes as $\psi_i\to0$.

Routine by routine: `../../REDISTANCING.md` §9.

## 6. Validation (2026-09-25; `logs/sep25/i1_bdfk/results.txt`)

The reference is the numpy replica,
`../../../neko-multiphase/examples/saini_benchmarks/redistance_eq44/`. It is a few
hundred lines of numpy, with the whole operator inspectable.

| check | result |
|---|---|
| replica `bdfk.py` against the 09-23 `wterm.py bdf2` | bit-identical: field and 120 trace lines, 1.8392e-2 |
| BDF3/EXT3 order, scalar ODE | 3 from exact history (error ratio 8.1); 2 with the Nek ramp start |
| BDF3 at $\Delta\tau$ against $\Delta\tau/2$, $H{=}1/10$, $N{=}3$ | $E_r(1)$ differs by 0.013% |
| scratch build, default path | every logged line identical to the committed run (1.8150e-2) |
| Neko BDF2 / BDF3 / BDF2 from $\psi_e$ against the replica, same $\psi_0$ (Neko's round-off at the 18 zero-set nodes injected) | max node difference 3.9e-13 / 4.7e-13 / 4.1e-13 at $\tau=6$; all 120 logged $E_r$ equal to the printed figures |
| step-map spectra at $\psi_e$, $H{=}1/5$, $N{=}3$ (`stepspec.py`) | all integrators reproduce the operator's 38 growing modes; leading rate 0.6094 exactly (RK3+Lie 0.6089); worst rate error RK3+Lie 1.2e-2, BDF2 3.3e-4, BDF3 9e-7 |
| operator audit (O2, `logs/sep25/o2_audit/`) | averaged gradient equals the mass-weighted one to 5e-14; `conv1` Galerkin equals the identity to 5.6e-14; GPU equals CPU |

### 6.1 Independent audit (2026-09-28; `logs/audit/`: `PREREGISTERED.txt`, `results.txt`)

The checks above lean on the replica, which we wrote from the same reading of the paper.
This audit checks each ingredient against something the replica does not share: the
rendered paper, first principles, and Neko's own routines. It used a scratch build
generated from the committed file; with every audit key off, that build reproduces all 104
logged lines of `run_circles_h10_n3.log`. **No implementation error was found.**

| ingredient | check | result |
|---|---|---|
| the printed forms | Eqs. (7), (22)–(24), (27)–(29), (31), (33)–(35) re-read from the rendered pages | as coded. $N_{svv}$ is an exponent with no cutoff; the basis is plain Legendre; $\hat Q$ sits whole on each side; $\mathbf D_\mu$ left-multiplies the assembled operator; "second order extrapolation" |
| SVV filter and `svv_local` | `fh` applied to $L_k$ at the GLL points; the energy of fields with $\partial_ru=L_k$ against $\nu\sigma_k^2\sum G\,L_k^2$; $N=3$–8 | $\sigma_kL_k$ to 3e-16; energy ratio 1 to 1.3e-12; `fht`$=$`fh`$^T$; metric to 4e-14. Not $F^T$, not Boyd's basis |
| the implicit Eq. (31) step | the residual of the unsymmetrised printed system after each step; an independent damped Richardson solve at $\tau=0.5,3,6$; four Table 2 cells and $H{=}1/20$, $N{=}8$ | Richardson equals CG to a few ulps of $\lvert\psi\rvert$ (≤1.2e-14 absolute); zero-set nodes unchanged bit for bit; CG takes 9–24 iterations, never the 200 cap |
| gradient | `unit_normal` against exact gradients, $H\in\{1/5,1/10,1/20\}$, $N=3$–8 | round-off on degree-$N$ polynomials; rate $N$ in $H$ on a smooth field, with no face or wall defect. On $\psi_e$: rate about $N$ where it is analytic; it does not converge next to its cone points (the disc centres and the interface corners) or its kinks, as §9.9 describes |
| flat interface | tilted line, committed scheme, to $\tau=6$ | exact to 8e-13 (first time in Neko) |
| time stepping | Neko's own scalar solver integrating the same explicit operator, SVV off | node for node equal (2.7e-14 relative at 200 steps), once our loop uses Neko's extrapolation (below) |
| time order | Richardson in $\Delta\tau$, committed scheme, to $\tau=1$ | order 2, drifting to 1 at $H{=}1/5$. The drift is the $\mathbf D_\mu(\psi^n)$ lag: with $\mathbf D_\mu$ extrapolated it is 3.97–4.00 |
| whole solver | manufactured solution, $\varepsilon=0.2$ fixed | order 2 in $\Delta\tau$; spatial rate 3.2–3.4 at $N{=}3$ (above the time floor); the implicit SVV consistent to the same floor |

Three things were learned that are not bugs:
- **Neko's `time_order` 2 is BDF2 with a modified EXT3**, $[8/3,-7/3,2/3]$ from step 3
  (`src/time_schemes/time_scheme_controller.f90`). The paper prints second-order
  extrapolation, which is what this case does. The modified EXT3 changes $E_r(6)$ by 0.05%.
- **$\mathbf D_\mu$ taken from $\psi^n$ makes the scheme first order in $\Delta\tau$** (§4). With
  it extrapolated, $E_r(6)$ moves by at most 0.31% in the Table 2 cells, and
  $H{=}1/20$, $N{=}3$ still fails.
- **On a smooth field, the Eq. (31) SVV's consistency error converges at order 2 in $H$**,
  not 1 (§8.3).

### 6.2 Rebuild from the paper with Neko's own operators (2026-09-28; `logs/rebuild/`: `PREREGISTERED.txt`, `results.txt`)

§6.1 still read the committed code. This check asks whether the committed code and
the replica could share a misreading of the paper.

**The design came first.** It was written from the rendered pages (`paper-lookup`)
before `redistance_circles.f90` was opened. It agrees with the committed code in
every term of Eqs. (7), (24)–(35), (42)–(46), (83) and (84). The only differences
are in implementation.

**The rebuild uses Neko's routines wherever they fit.** Its core
(`rebuild_core.f90`) shares no routine with this case's file:

| piece | Neko routine |
|---|---|
| $\mathbf S_{vv}$, Eqs. (18), (22)–(28) | `ax_helm`, the stiffness kernel (CPU and CUDA), run on a second `space_t` whose derivative matrices are $F\hat D$, with $F$ = `elementwise_filter_t%fh` ("nonBoyd", $(k/N)^{N/6}$) |
| $\mathbf D_\mu$, Eqs. (29), (31) | $\frac{c_0}{N}\lvert\operatorname{sgn}\psi\rvert\,2J_e^{1/2}$ at the GLL points, from the in-plane Jacobian, left-multiplying the operator |
| $\mathbf H$ solve, Eqs. (34)–(35) | `gmres` on the unsymmetrised $\mathbf H$, in delta form, with the `jacobi` preconditioner |
| $\mathbf C(\mathbf w)\psi$ | `adv_no_dealias_t`, i.e. $-\mathbf B$`conv1` |
| $\nabla\psi$ | `opgrad`, gather-scatter, $\mathbf B^{-1}$ (mass-weighted average) |
| BDF/EXT | `bdf_`/`ext_time_scheme_t`, or `time_scheme_controller_t`; `rhs_maker_ext`/`_bdf`, `field_series` |
| Eqs. (83), (84) | written from the print, with Eq. (83)'s division |

**Operator checks.** A check build (the committed file, only read, plus the rebuild
core) compares each ingredient on the same fields: $\psi_e$, $\psi_0$ and the
committed $\tau=6$ field, at $H{=}1/10$, $N=3$ and 8. Its own relaxation reproduces
the committed logs.

| ingredient | result |
|---|---|
| Eq. (83) | bit-identical to `circle_distance` at every node of all 18 grid cells; the same exact zeros |
| Eq. (84) | bit-identical: 1.837871825798979e-2 and 3.416342603046530e-3 |
| explicit term $\operatorname{sgn}\psi-\mathbf B^{-1}\mathbf C(\mathbf w)\psi$ | ≤7.5e-16 relative |
| averaged gradient | ≤1.3e-13 relative (summation order; mass-weighted and arithmetic averages are equal on this mesh) |
| $\mathbf S_{vv}$ | absolute differences ≤1.9e-15, on outputs that are themselves small cancellation residuals. The excess is `ax_helm`'s $z$-direction terms, which are round-off on a $z$-invariant field and which `svv_local` omits. With them dropped, the $N{=}8$, $\tau=6$ difference falls from 1.9e-15 to 3.0e-19 |
| one implicit step | GMRES against the symmetrised CG: ≤1e-12 relative on all three fields. Zero-set nodes unchanged bit for bit by both |
| BDF/EXT coefficients | as printed; the controller gives Neko's modified EXT3 from step 3 |
| a tilted straight interface, to $\tau=6$ | exact: 4.0e-15 and 7.6e-15 |
| Neko's `cfl()` on $\mathbf w$, $H{=}1/20$, $N{=}8$ | 0.2822. The paper states "$CFL=0.28$ for the highest resolution case" (p. 19) |

**Full runs, from the skewed IC** (committed runs as reference; `cmp_fields.py`):

| cell | logged $E_r$ equal (5 figures) | max node difference, $\tau=6$ | $E_r(6)$, relative |
|---|---|---|---|
| $H{=}1/10$, $N{=}3$ | 26/26 | 2.4e-12 | 1.8e-11 |
| $H{=}1/5$, $N{=}7$ | 26/26 | 1.2e-11 | 2.0e-10 |
| $H{=}1/10$, $N{=}7$ | 26/26 | 5.4e-11 | 4.8e-11 |
| $H{=}1/5$, $N{=}8$ | 26/26 | 3.6e-13 | 1.2e-12 |
| $H{=}1/10$, $N{=}8$ | 26/26 | 6.1e-12 | 8.2e-11 |
| $H{=}1/20$, $N{=}8$ | 26/26 | 1.4e-12 | 3.9e-11 |
| $H{=}1/20$, $N{=}3$ | fails at the apex, identically to 5 figures through $\tau=4.5$ | | round-off takes over after (§10) |

The $N{=}8$ map split (§7.3) agrees within 5e-10 in all three cells.

With Neko's own `time_scheme_controller` (modified EXT3), the only change, $E_r(6)$
moves by −0.06% to +0.11% in the six holding cells, and $H{=}1/20$, $N{=}3$ still fails.

**What this settles.** A third implementation, sharing none of the hand-written
numerics, gives the committed numbers. So the gap to Table 2 is neither an
implementation error nor a shared misreading of any printed equation. It lies in
what §4.4 does not state (§4, §10).

**What it found on the way** (none of it a bug):
- **Their Fig. 1 confirms Eq. (24) as printed.** Digitised, every curve is
  $(k/N)^{N_{svv}}$ within 0.5%, with $\hat Q_0=0$ and no cutoff.
- **Figs. 11 and 13 carry no axes.** The L/R and T/B split of §7.3 rests on an
  inference: the circles lie side by side, and the flat sector of Fig. 11(a) is at
  the top right.
- **GMRES needs the `jacobi` preconditioner.** With the identity it stalls at
  $N{=}8$ (200 iterations, unconverged); with Jacobi it takes 9.8–29.6 per step on
  average over the full runs, and never fails to converge.
- **The CUDA `ax_helm` kernel is not bit-reproducible call to call at $N{=}8$**
  (6.5e-19). This is autotuning, not the plumbing.

## 7. Results

Numbers are $E_r$ by Eq. (84). "Theirs" is the authors' own values, `plot.py`
L83–85 in their public case (`nandu90/nekLS_Examples@jcp`, `intersectingCircles/`,
commit `bceef31`). Our earlier pixel reading of Fig. 12 was within 0.7% of them in every
cell, so no ratio changed at two figures (2026-09-29).

### 7.1 Table 2 — BDF2/EXT2, the paper's scheme

| cell | Saini Table 2 | ours, BDF2 | ratio |
|---|---|---|---|
| $H{=}1/10$, $N{=}3$ | 6.76e-3 | 1.8379e-2 | 2.7 |
| $H{=}1/5$, $N{=}7$ | 3.87e-3 | 8.2931e-3 | 2.1 |
| $H{=}1/10$, $N{=}7$ | 1.35e-3 | 3.7590e-3 | 2.8 |
| $H{=}1/20$, $N{=}3$ | 1.85e-3 | 0.168 | fails at the apex (§10) |

### 7.2 Fig. 12 — the grid, $E_r(6)$, theirs / ours (BDF2) / ratio

| BDF2 | $N{=}3$ | $4$ | $5$ | $6$ | $7$ | $8$ |
|---|---|---|---|---|---|---|
| $H{=}1/5$ | 1.562e-2 / 4.617e-2 / 3.0 | 1.030e-2 / 1.658e-2 / 1.6 | 8.131e-3 / 1.159e-2 / 1.4 | 4.125e-3 / 8.964e-3 / 2.2 | 3.874e-3 / 8.293e-3 / 2.1 | 1.680e-3 / 5.169e-3 / 3.1 |
| $H{=}1/10$ | 6.761e-3 / 1.838e-2 / 2.7 | 4.869e-3 / 1.018e-2 / 2.1 | 2.805e-3 / 7.578e-3 / 2.7 | 1.862e-3 / 4.276e-3 / 2.3 | 1.348e-3 / 3.759e-3 / 2.8 | 1.315e-3 / 3.416e-3 / 2.6 |
| $H{=}1/20$ | 1.858e-3 / 0.168 / fails | 2.017e-3 / 4.576e-3 / 2.3 | 1.923e-3 / 4.161e-3 / 2.2 | 9.810e-4 / 3.214e-3 / 3.3 | 5.815e-4 / 1.978e-3 / 3.4 | 5.116e-4 / 1.649e-3 / 3.2 |

- **1.4–3.4× above theirs in every cell**; $H{=}1/20$, $N{=}3$ fails.
- **The integrator makes no difference.** In every cell that holds, BDF2 is within
  −0.02% to +2.2% of the 2026-09-23 SSP-RK3 grid, which was 1.4–3.3× above theirs.
- **The curve shapes are ours, not theirs:**
  - Their $H{=}1/20$ curve is flat over $N=3$–5, then drops. Ours falls steadily
    from $N{=}4$.
  - Their $H{=}1/5$ curve drops most at $N{=}6$ and $N{=}8$. Ours drops most at
    $N=3\to4$.
- The figure is `evidence/fig12_error_decay.png`, with their values dashed.

### 7.3 Fig. 13 — $\lvert e\rvert$ at $N{=}8$, and the gap is in the field

- **Their maps integrate to their $E_r$.** Their Fig. 13, read back through its colour
  bar, integrates to 0.92–1.49× their Fig. 12 $E_r$ at $H{=}1/20$, where nothing
  saturates. So the metric is Eq. (84) as printed, and the difference is in the field.
- Their 2025 report (ANL/NSE-25/51) is not a target: its $E_r$ is about $10^3$× below
  the journal in every case.

The split below uses the "diamond", the corner branch of Eq. (83), with vertices
$(\pm0.7,0)$ and $(0,\pm0.714)$. The columns are $\int\lvert e\rvert$ inside and
outside it, and the left/right and top/bottom ratios of the outside part.
The skewed-IC rows are BDF2 (2026-09-26). The $\psi_e$ rows are RK3 (2026-09-24),
which BDF2 matches within 2% wherever both were run (§7.4).

| $N{=}8$ | $E_r(6)$ | diamond | outside | outside L/R, T/B |
|---|---|---|---|---|
| $H{=}1/5$, skewed IC | 5.169e-3 | 7.81e-3 | 9.86e-3 | 2.38, 0.58 |
| $H{=}1/5$, from $\psi_e$ | 3.620e-3 | 7.94e-3 | 4.44e-3 | 1.01, 1.00 |
| $H{=}1/5$, theirs | 1.68e-3 | ≥3.29e-3 (clipped) | 1.25e-3 | 1.0, 0.94 |
| $H{=}1/10$, skewed IC | 3.416e-3 | 7.44e-3 | 4.24e-3 | 2.14, 0.53 |
| $H{=}1/10$, from $\psi_e$ | 2.704e-3 | 7.41e-3 | 1.84e-3 | 1.00, 1.00 |
| $H{=}1/10$, theirs | 1.32e-3 | ≥3.45e-3 (clipped) | 2.1e-4 | 0.7, 1.21 |
| $H{=}1/20$, skewed IC | 1.649e-3 | 3.33e-3 | 2.30e-3 | 2.10, 0.57 |
| $H{=}1/20$, from $\psi_e$ | 1.264e-3 | 3.33e-3 | 9.95e-4 | 1.0, 1.0 |
| $H{=}1/20$, theirs | 5.1e-4 | 1.54e-3 | 5.2e-5 | 0.3, 2.77 |

What the comparison shows:
- Theirs keeps no memory of the IC's skew.
- Ours carries about 2.2× their error in the diamond, plus fans of rays outside it.
- The rays survive being pushed through their colour pipeline (`logs/fig13cmp.py`;
  BDF2: `logs/sep25/i1_bdfk/fig13cmp_bdf2.log`). Outside the diamond, at their
  colour resolution, ours carries 5.56e-3, 2.49e-3 and 1.43e-3 against their
  1.25e-3, 2.07e-4 and 5.2e-5 at $H{=}1/5$, $1/10$, $1/20$. The clipped diamond is
  1.8×, 1.7× and 2.1× theirs.
- The figure is `evidence/fig13_error_maps.png`.

### 7.4 The integrator, in Neko: RK3, BDF2, BDF3 from both ICs, to $\tau=24$

$E_r$ at $\tau=6$ and 24 (`logs/sep25/i1_bdfk/`; the $\tau=6$ runs are a bit-exact
prefix of the $\tau=24$ ones).

| cell | IC | SSP-RK3 + Lie, 6 / 24 | BDF2/EXT2, 6 / 24 | BDF3/EXT3, 6 / 24 | BDF2 against RK3, max over $\tau$ |
|---|---|---|---|---|---|
| $H{=}1/10$, $N{=}3$ | skewed | 1.8150e-2 / 3.4899e-2 | 1.8379e-2 / 3.4937e-2 | 1.8370e-2 / 3.4929e-2 | 1.6% |
| | $\psi_e$ | 1.5916e-2 / 3.5762e-2 | 1.5958e-2 / 3.5829e-2 | 1.5958e-2 / 3.5829e-2 | 0.38% |
| $H{=}1/5$, $N{=}7$ | skewed | 8.2351e-3 / 1.1689e-2 | 8.2931e-3 / 1.1756e-2 | 8.2946e-3 / 1.1756e-2 | 0.73% |
| | $\psi_e$ | 5.5281e-3 / 9.8969e-3 | 5.5689e-3 / 9.9447e-3 | 5.5689e-3 / 9.9447e-3 | 0.74% |
| $H{=}1/10$, $N{=}7$ | skewed | 3.7037e-3 / **NaN** | 3.7590e-3 / **NaN** | 3.7626e-3 / **NaN** | 1.6% to $\tau=15$ |
| | $\psi_e$ | 2.6145e-3 / **NaN** | 2.6568e-3 / **NaN** | 2.6569e-3 / **NaN** | 1.6% to $\tau=15$ |
| $H{=}1/20$, $N{=}3$ | skewed | 0.366 / **NaN** | 0.168 / **NaN** | 0.309 / **NaN** | fails; only "fails" is reproducible |
| | $\psi_e$ | 5.3491e-3 / 1.9529e-2 | 5.3739e-3 / 1.9672e-2 | 5.3739e-3 / 1.9672e-2 | 1.9% |

- **The integrator makes no difference.** Where the field holds, BDF2 and BDF3 are
  within 1.9% of RK3 at every reported $\tau$ to 24. From $\psi_e$, BDF2 and BDF3
  agree to every printed figure: both are converged in time.
- **The offset from RK3 is RK3's own error.** The 0.26–1.6% offset at $\tau=6$ is
  RK3's first-order Lie split: halving $\Delta\tau$ moves RK3 toward the BDF value,
  and the spectra show the same.
- **The map splits agree.** At $\tau=6$ and 24, in every cell that holds and from
  both ICs, the diamond is within 2.4%, the outside within 1.6%, and L/R and T/B
  are equal to 0.02 (`i1_bdfk/i1_table.py split`).
- **Nothing is steady at $\tau=24$.** From $\tau=6$ to 24, $E_r$ grows 1.4–3.6× in
  every cell that holds; the 3.6× is $H{=}1/20$, $N{=}3$ from $\psi_e$.
  - At $\tau=24$ the skewed and exact runs of a cell have nearly converged to each
    other ($H{=}1/10$, $N{=}3$: 3.49e-2 against 3.58e-2), so the IC memory decays
    and the drift remains.
- **$H{=}1/10$, $N{=}7$ blows up at $\tau\approx15.5$–16.5, under every integrator
  and from $\psi_e$ too.**
  - $E_r$ creeps from 2.6e-3 to 3.7e-3 by $\tau=15$ (from $\psi_e$), then jumps:
    0.11 at 15.5, 0.51 at 16, then NaN.
  - Beforehand the minimum $\lvert\nabla\psi\rvert$ grows exponentially from
    round-off, at about 1.1–1.6 per unit $\tau$: 8e-12 at $\tau=6$, 1.3e-6 at 10,
    3e-5 at 12, 8e-4 at 15. That is a symmetry-breaking mode, odd in $x$ and in $y$.
  - **But its amplitude does not set the time** (2026-09-26, `logs/sep26/results.txt`).
    The skewed run carries about 4e-3 of asymmetry from $\tau\le11$, about $10^{9}$ times
    the $\psi_e$ run's at $\tau=6$, yet starts only 0.2 earlier (onset 14.45 against
    14.65, from frames every 0.05). Something in the common state switches a local growth on at
    $\tau\approx14.3$–14.5, and it amplifies whatever asymmetry is there. Where: §10.
- **$H{=}1/20$, $N{=}3$ from the skewed IC** reaches NaN at $\tau\approx6.25$ under
  every integrator.
- All pre-registered predictions hold (`logs/sep25/i1_bdfk/PREREGISTERED.txt`):
  - **P1:** within 3% of RK3 wherever the field holds.
  - **P2:** $H{=}1/20$, $N{=}3$ fails from the skewed IC and holds from $\psi_e$.
  - **P3:** no cell within 1.4× of theirs.
  - The reopen criterion is not triggered.

### 7.5 The $N_{svv}=N/4$ arms: where the integrator does matter

The two committed $N/4$ arms are not §4.4 settings; the paper uses $N/6$ throughout
§4.4. They test the weaker filter the coupled cases inherit from Saini §4.5.
$E_r(6)$:

| $N_{svv}=N/4$ | SSP-RK3 + Lie (09-23) | BDF2/EXT2 (committed) | BDF3/EXT3 | BDF2 from $\psi_e$ |
|---|---|---|---|---|
| $H{=}1/10$, $N{=}3$ | 1.5113e-2 | 1.5256e-2 | — | — |
| $H{=}1/5$, $N{=}7$ | 9.169e-3, holds | **0.546, fails** | **0.529, fails** | 6.48e-3, holds |

- **The failure belongs to the scheme.** At $H{=}1/5$, $N{=}7$ with $N/4$ the scheme
  itself fails at the apex from the skewed IC: BDF2 and BDF3 agree with each other.
  - They part from RK3 at $\tau\approx1.25$, and the disc interiors fill in:
    min $\psi$ goes from −0.99 to −0.21, the failure mode of $H{=}1/20$, $N{=}3$.
- **RK3 held it** only because its first-order Lie split adds damping. Confirmed
  2026-09-28 by step size alone (`logs/sep26/results.txt`, $E_r(6)$):

  | $\Delta\tau$ | $5\times10^{-4}$ | $2.5\times10^{-4}$ | $1.25\times10^{-4}$ |
  |---|---|---|---|
  | RK3 + Lie | 9.169e-3, holds | 9.166e-3, holds | **0.484, fails** |
  | BDF2/EXT2 | 0.546, fails | 0.495, fails | 0.533, fails |

  - So the time-converged scheme fails, and RK3 holds only above a step between
    $1.25\times10^{-4}$ and $2.5\times10^{-4}$. More SVV does not help: BDF2 at
    $c_0=2.5$ also fails (0.492).
  - All the schemes agree to three figures until $\tau\approx1.25$ and part at
    1.3–1.5, just after the front collapses on the circle centres (§10).
- **What this says about the apex failure.** It is set by how much dissipation
  reaches the apex when the front collapses there. At $N/4$ the converged scheme
  lacks it, and neither $c_0=2.5$ nor a finer step supplies it. That fits the
  paper's remark that §4.4 "required stronger diffusion than the preceding cases"
  ($N/6$).
- **Consequence for the coupled cases.** They run $N/4$ in their redistancing SVV
  with RK3, so any apex they hold may be held by the split damping at their step.
  A smaller $\Delta\tau$ or BDF2 would remove it. That matters only where the seed
  is steeper than a distance function near a medial-axis point inside the band
  (§10); nobody has checked whether any coupled case has such a point.
- At $N/6$, the §4.4 setting, the integrator changes nothing (§7.4).
- Runs: `logs/sep25/i1_bdfk/runs/nsvv4_*`, `logs/sep26/q7_margin/`.

## 8. Where the gap comes from

These results were measured with RK3 in the replica and in Neko on 2026-09-24. They
carry over because the integrator changes nothing (§7.4).

### 8.1 The scheme does not hold the exact solution

- Started from $\psi_e$, $E_r$ rises from zero: 2.74e-3 at $\tau=0.5$, 1.59e-2 at 6,
  3.56e-2 at 24 ($H{=}1/10$, $N{=}3$).
- It ends 1.4–2.9× above their $E_r$ in all six cells tried.
- It is not steady at $\tau=24$. Saini call $\tau=6$ "quasi-steady".
- The skewed run's minimum at $\tau\approx1.5$–3 is where the decaying transient
  crosses this drift.

### 8.2 Three parts

Splitting $\int\lvert e\rvert$ by the diamond:
- **The diamond** does not depend on the IC.
  - A third of it is the valley kink lying on element faces: translating the geometry
    off the mesh lines removes 34%.
  - The rest follows where the interface corners settle, pinned toward a node or
    mesh line.
- **The arc drift outside it** is symmetric and grid-oriented, and present from
  $\psi_e$. A lone circle drifts too (8.9e-3 at $\tau=6$, $H{=}1/5$, $N{=}3$).
- **Memory of the skewed IC** is left- and bottom-heavy. It is the steep-side
  transient: ~12% of $\int\lvert e\rvert$ at $N{=}3$, about half the outside error at
  $N{=}8$.

### 8.3 The mechanism of the arc drift

- Every signed distance function is a steady state of Eq. (44). Nothing in the
  equation restores an interface's position; the paper's own introduction names
  this as the known weakness of PDE redistancing.
- $\mathbf D_\mu=\lvert\operatorname{sgn}\psi\rvert$ vanishes *on* the zero set, but
  not at the nodes beside it.
- There the SVV of a curved distance function is non-zero and moves both sides the
  same way. So the interpolated interface moves although no node changes sign.
- At $\psi_e$ that term is 18× the rest of the scheme on the arcs.
- **It is curvature, not dimension** (2026-09-25; replica, from $\psi_e$, $H{=}1/5$,
  $N{=}3$, $\tau=6$, `flat_vs_curved.py` in the replica directory).
  - A straight interface tilted 30° across the grid stays exact to round-off:
    $\int\lvert e\rvert$ 1.8e-13 at every $\tau$. Its distance function is linear, so the
    gradient is exact, $\lvert\nabla\psi\rvert=1$, and the SVV of it is zero ($\hat Q_0=0$). Every
    term vanishes, and it is an exact discrete steady state. Neko reproduces this to 8e-13
    (§6.1).
  - One circle drifts: $\int\lvert e\rvert$ 1.33e-2, 2.86e-2, 7.56e-2 at $\tau=0.5$, 2, 6.
  - Saini's two circles drift more: 2.82e-2, 4.08e-2, 1.29e-1.
  - The 1D analogue ($\psi_e=\lvert x\rvert-1$, skewed IC) converges and is
    quasi-steady (replica README §1; that run predates Eq. (31)).
  - From Saini's $r_f$ seed instead, the 1D run captures its walls at $N\ge7$
    (§10, "1D, reinit seed").
- **On a smooth field, the SVV's own error is second order in $H$** (2026-09-28,
  manufactured solution, $\varepsilon=0.2$, $\tau=0.5$, `logs/audit/results.txt` 4c).
  - The rms error is 6.7e-4, 1.6e-4, 4.1e-5 at $H=1/5$, $1/10$, $1/20$, $N{=}3$ (rate 2.0).
    It is 1.5–1.7 at $N{=}7$, and falls steeply with $N$ (2.0e-5 at $H{=}1/5$, $N{=}8$).
  - Pointwise the SVV forcing is $O(H)$, but it has zero mean per element, so its net
    effect on a smooth field is $O(H^2)$.
  - It was measured outside the §4.4 regime: a fixed $\varepsilon=0.2$, and a field that is
    not a distance function. So it neither explains nor bounds the §4.4 drift, which acts
    beside the zero set at $\varepsilon=H/N$ (§8.5).

### 8.4 No parameter stops the drift

$c_0$, $N_{svv}$ and $\varepsilon$ were varied one at a time from $\psi_e$. Each
changes the drift's rate, and none stops it.
- $\varepsilon$ is the strongest knob.
  - At $\varepsilon=H$ the $\tau=6$ values come to 0.81–1.51× theirs at $H{=}1/5$
    and $1/10$.
  - That reproduces their $H{=}1/5$ curve shape and $N{=}8$ map.
- But every cell is still rising at $\tau=6$, and the best-fit $\varepsilon$ changes
  with $N$: about $2H$ at $N{=}3$, $H$ at $N{=}7$.
- A fixed $\varepsilon=0.2$ (2026-09-25, `logs/sep25/n1_eps_fixed/`) met 3 of 4
  pre-registered predictions at $H{=}1/10$, $N{=}8$ ($E_r$ 0.91× theirs), but is
  not quasi-steady either.
  - At $H{=}1/20$, $N{=}8$ it meets both of its pre-registered bands (judged
    2026-09-26): $E_r(6)$ 4.77e-4, outside 8.08e-5.
  - Its diamond is 1.550e-3, against their 1.54e-3.
  - At fixed $\varepsilon$ the outside error converges at order 1.81 and 2.02 over
    $H=1/5\to1/10\to1/20$. Theirs is 2.59 and 1.99; ours at $\varepsilon=H/N$ is 1.23.
  - Only its outside L/R and T/B, 0.94 and 1.04, miss theirs (0.3 and 2.77).
- **No single fixed $\varepsilon$ gives their $H{=}1/20$ map** (2026-09-28,
  `logs/sep26/results.txt`, P15). At $N{=}8$, from the skewed IC:
  - $\varepsilon=0.25$ meets every pre-registered band except the $H{=}1/20$ outside
    error: 1.52e-4 against their 5.2e-5. Its L/R and T/B are 0.52 and 1.92. At
    $H{=}1/5$ it gives $E_r$ 1.40e-3, outside 1.10e-3, L/R and T/B 0.99 and 1.00.
  - $\varepsilon=0.30$ reproduces their $H{=}1/20$ asymmetry, 0.34 and 2.92, but with
    13× their outside error.
  - So the top-right-heavy signature has their shape only while our flat-sector lag
    is still large. The IC is flat there ($\psi_0=0.1\psi_e$ near $(1,1)$), and
    $\tanh(\psi/2\varepsilon)$ stays small.
  - Their flat sector relaxes faster at the same shape. A sign function normalised by
    $\lvert\nabla\psi\rvert$, or normalising $\psi_0$ first (Sussman–Fatemi,
    Appendix A), would do that. That is untested, and Eq. (46) prints
    $\tanh(\psi/2\varepsilon)$.
- So a matched $E_r(6)$ is no evidence for any $\varepsilon$. **The committed
  cases keep $H/N$.**
- §8.5 explains why $\varepsilon$ can only change when the drift shows, not where it
  ends.

### 8.5 How the drift ends: node capture, on an $\varepsilon$ clock

2026-09-26/28. Runs and analysis: `logs/sep26/` (`PREREGISTERED.txt`, `results.txt`),
and the replica's `capture_logs/`.

**Near the zero set every term is proportional to $\psi_i$.**
- $\mathbf D_\mu=\lvert\tanh(\psi/2\varepsilon)\rvert$ and the GLL-collocated
  $\mathbf C(\mathbf w)$ make every term of a node's rate proportional to $\psi_i$
  as $\psi_i\to0$. So a node beside the interface obeys

  $$\partial_\tau\psi_i\approx\psi_i\,\sigma_i,\qquad\sigma_i=\frac{(1-\lvert g_i\rvert)-\operatorname{sgn}(\psi_i)\,q_i}{2\varepsilon},\qquad q=\mathbf B^{-1}\mathbf S_{vv}\psi .$$

- With $\sigma_i<0$ the node is **caught**: $\psi_i\to0$, super-exponentially at
  the end, and the interface moves onto it.
- With $\sigma_i>0$ it is pushed off until tanh saturates.
- So the interface's position between nodes is not a neutral mode of the discrete
  scheme. It relaxes until enough nodes are caught.

**What the pattern of $q$ is.**
- The filter has $t_0=0$, so on a curved distance function each element keeps only
  the variation of its slope.
- The result has zero mean per element. It is diffusive at element-interior nodes
  and anti-diffusive, about $N(N+1)/2-1$ times stronger, on element-boundary lines.
- That is why the arc drift is grid-oriented, and why mean $\lvert dR\rvert$
  exceeds $\lvert$mean $dR\rvert$ in every run, by 1.2–12×.
- The mean is a small residual. At $H{=}1/10$, $N{=}3$ it turns outward once
  capture starts. At $H{=}1/5$, $N{=}7$ it stays inward to $\tau=24$.

**Measured.**
- **Captures accumulate while $E_r$ levels off.** At $H{=}1/10$, $N{=}3$ from
  $\psi_e$ the caught nodes number 4, 30, 32, 52 at $\tau=6$, 12, 18, 24.
  $E_r(48)/E_r(24)=1.04$, with a slow tail to 3.86e-2 at $\tau=72$.
- **Captures are how the drift ends, but which nodes go is not readable from
  $\psi_e$ alone.** The sign of $\sigma_i$ at $\psi_e$ predicts 93–96% of the caught
  nodes here, but only 33–88% on single circles at $H{=}1/5$.
- **The SVV on the nodes next to the interface holds them, it does not push them.**
  Zeroing $\mathbf D_\mu$ on just those nodes (their Eq. (44) rate kept) catches 289
  of 486 by $\tau=24$. The zero set thickens into a band, and $E_r$ settles at 0.132,
  3.7× the base. The push comes from the SVV on the next rings.
- **$\varepsilon$ is an exact clock.**
  - The nodes that decide the drift have $\lvert\psi\rvert\lesssim\varepsilon$,
    where $\tanh(\psi/2\varepsilon)\approx\psi/2\varepsilon$. So every term there
    scales as $1/\varepsilon$.
  - From $\psi_e$ at $H{=}1/10$, $N{=}3$ to $\tau=72$, $\varepsilon\times3$ at $\tau$
    equals $\varepsilon$ at $\tau/3$: outside within 12%, diamond within 6%, the
    corner shift and even the caught-node counts identical.
  - The end states agree: $E_r$ 3.567e-2 against 3.559e-2.
  - It holds at $H{=}1/5$, $N{=}7$ with $\varepsilon\times7$ too, within 6%.
- **The size of the end state goes roughly like h.** Circles end at mean
  $\lvert dR\rvert$ of 0.14–0.29 h (h = H/N). Halving H takes a circle's to 0.41–0.43×
  (strict ∝h would be 0.5×).
- **The clock has no clean H or curvature scaling.**
  - Circles of radius 0.5, 1 and 1.5 plateau at 1.9e-2, 1.75e-2 and 0.93e-2.
  - Halving H makes a circle's captures start sooner (τ≈12 against ≈36), but slows
    $H{=}1/20$, $N{=}3$ for Saini's geometry (outside ×6.7 from τ=6 to 24, against
    ×2.7 at $H{=}1/10$).

**What it means for §4.4.**
- A larger $\varepsilon$ does not reduce the drift, it postpones it. At a fixed
  $\varepsilon=0.2$ the drift reached at $\tau=6$ is what $\varepsilon=H/N$ reaches
  at $\tau=6\,(H/N)/0.2$: at $H{=}1/20$, $N{=}8$ that is $\tau\approx0.19$.
- So their "quasi-steady at $\tau=6$" fits a larger $\varepsilon$ with a field that
  would still drift later (H9′).
- Removing the drift itself needs something that stops capture. The Russo–Smereka
  anchor does (§10).

## 9. Eliminated — do not re-test

Details in the replica README §8–§10.

- **The time integrator:** §7.4 and the spectra, in Neko, with BDF2 and BDF3, from both ICs.
  Tested one at a time: Neko's own default extrapolation (modified EXT3) changes $E_r(6)$ by
  0.05%, and an extrapolated $\mathbf D_\mu$ by at most 0.31% (§6.1).
- **An implementation error in the Eq. (44) solver:** the independent audit of §6.1.
- **A misreading of a printed equation shared by the code and the replica:** the rebuild
  from the paper with Neko's own operators (§6.2) gives the committed numbers.
- **The time treatment of $\mathbf w$:** consistent treatments agree within 1.6%;
  frozen or stale ones fail (§5).
- **Dealiasing $\mathbf C(\mathbf w)$:** 3.3× worse, because it moves zero-set nodes
  at about 0.08.
- **$\mathbf w$ per element**, i.e. $\lvert\cdot\rvert$ before averaging: worse in every
  cell. It is 2.0× at $H{=}1/5$, $N{=}7$, about 9× at $H{=}1/20$, $N{=}3$, and an
  order-one failure at $H{=}1/10$, $N{=}3$ ($E_r$ 1.32; the disc interiors fill in).
- **The SVV inside the bilinear form** with $\lvert\mathbf c\rvert=1$ (the first
  version of this case): 14–920× off Table 2 across the four cells (221× at
  $H{=}1/10$, $N{=}3$), with spurious interfaces (2502 sign-change edges against 337
  there).
- **SUPG**, in each tried form.
- **A frozen $\operatorname{sgn}(\psi_0)$:** diverges.
- **Floors on $\mathbf D_\mu$:** do nothing.
- **$\mathbb H$, `cg_iters`, 2D against the 3D slab, the boundary** (from the skewed IC
  and $\psi_e$, which are steep at every wall; a flat reseed is another matter, §10).
- **The IC's sign convention:** Eq. (44) and the scheme are exactly odd in $\psi$;
  checked bit for bit.
- **Every $\varepsilon$ as a fix** (§8.4).
- **Saini's $2.5H$ band (H13)**, 2026-09-28 (`logs/band/`: `PREREGISTERED.txt`,
  `results.txt`). All three pre-registered falsifiers hold:
  - **The paper:** the band is a pseudo-time budget, $N_{tls}=2.5H/\Delta\tau_{tls}$
    (pp. 12–13). §4.4 runs to $\tau=6$ "in the entire domain", and $E_r$ is "global"
    (p. 19).
  - **Their Fig. 13:** 51.8%, 87.5% and 97.1% of their $\int\lvert e\rvert$ lies
    beyond $2.5H$ at $H=1/5$, $1/10$, $1/20$; 81–97% of that is in the diamond. So
    neither their field nor their metric was band-limited.
  - **The mechanism** (replica, $H{=}1/10$, $N{=}3$, from $\psi_e$, $\tau=24$): a mask
    $\lvert\psi_e\rvert\le2.5H$ on the whole update, $\psi$ frozen outside, leaves
    the drift as it is. Mean $\lvert dR\rvert$ is within 0.3% at $\tau=6$ and 3.6%
    at 24, the corner is identical, and the caught nodes are 4, 30, 44, 52 at
    $\tau=6$, 12, 18, 24 against 4, 30, 32, 52 unmasked, ending on the same set.
    $E_r$ falls 3.9× only because $\psi_e$ is held where the error is counted.
  - **The metric:** $E_r$ with only the numerator over the band lands at 0.63–1.14× in
    the four Table 2 cells, but only by discarding 47–99% of our error. Over the
    rest of the skewed-IC Fig. 12 grid it undershoots to 0.23×. It leaves our curve
    shapes and the $N{=}8$ L/R and T/B (2.1–2.5, 0.64–0.68) as they are. With both integrals over the band,
    $E_r$ is 11–300× theirs in the four Table 2 cells.
  - A pseudo-time cap gives $E_r\approx5.2$ at $\tau=2.5H$ ($H{=}1/10$, $N{=}3$). A
    frozen-outside mask from the skewed IC gives $E_r\ge9.79$.

## 10. Open

- **H8 and H9′ are not pursued** (decision, 2026-09-28). What is known stays here
  as a record; no further runs are planned for either.
  - **H8: an unstated interface-preserving step** (Sussman–Fatemi [19],
    Russo–Smereka [16], Saini & Bolotnov [20]). In the replica (2026-09-25,
    `preserve.py`, `preserve_logs/`), a Russo–Smereka subcell anchor stops the drift:
    $E_r$ 4.144e-3 at $H{=}1/10$, $N{=}3$ from the skewed IC (0.61× Table 2), steady
    to $\tau=24$. From the skewed IC at $N{=}8$ it fails at the apex (0.43 at
    $H{=}1/5$). Switched on at $\tau=1.5$ it holds the apex but keeps the IC memory
    (2.4× theirs; 2026-09-26/28, `logs/sep26/results.txt` P9). A step-to-step global volume constraint and a frozen
    $\operatorname{sgn}(\psi_0)$ are both worse than no anchor.
  - **H9′: a larger $\varepsilon$, with their field still drifting at $\tau=6$.**
    $\varepsilon$ is a clock for the drift (§8.5), so a larger one looks
    quasi-steady at 6. A fixed $\varepsilon=0.2$ gives their $N{=}8$ diamond and
    0.9–1.6× their outside error, but no single fixed $\varepsilon$ gives their
    $H{=}1/20$ L/R and T/B (§8.4).
- **$H{=}1/10$, $N{=}7$ blows up at $\tau\approx15$–16.5, from both ICs, under every
  integrator** (§7.4). Located 2026-09-26 from frames every 0.05 in $\tau$ written by
  the animation work, read only (`logs/sep26/q6_frames_analysis.log`):
  - **Where it starts.** Inside the left disc, at $(-1.10,+0.54)$ and its
    neighbours, where $\psi_e\approx-0.33$, on the element face $x=-1.1$.
    - It is the same site from both ICs. From $\psi_e$ it also appears at the point
      mirror $(+1.10,-0.54)$: the growing structure is odd in $x$ and in $y$.
    - It is not at the interface corners, and no node is caught there.
  - **What happens before it.** The fastest nodes are the exact zeros of $\psi_e$
    where the circle passes through element vertices, $(-1.3,\pm0.8)$ and
    $(-1.5,\pm0.6)$. The onset site lies inward from $(-1.3,0.8)$ along the normal.
  - **How it grows.** It e-folds in 0.13 $\tau$, then in 0.025. The disc interior
    fills in ($\psi\to0$ at the apex by $\tau\approx15.5$), then the domain corners.
  - **When.** The two ICs start 0.2 apart even though their asymmetry seeds differ
    by $10^{9}$ (§7.4). What switches the growth on is not found; the vertex zeros
    upstream are the lead.
- **The apex fails when the relaxation front collapses on the circle centre**
  (2026-09-26/28, `logs/sep26/results.txt`; fields every 0.05).
  - **The mechanism.** From the skewed IC, the region the characteristics have not
    reached yet keeps the IC's slope $s_0>1$ and rises at $s_0-1$. It reaches −1
    when the front arrives at $\tau\approx R=1$.
  - The vertex apex reads $\lvert\overline{\nabla\psi}\rvert\approx0.1$–0.15
    (`../../REDISTANCING.md` §9.9), so the apex lags while its surroundings rise.
  - **When.** Every failing run parts from its holding partner at $\tau$ 1.05–1.35:
    $H{=}1/20$, $N{=}3$ (left, steep apex); the $N{=}8$ anchor run; $N/4$ BDF2.
    Before $\tau=1$ each pair agrees to three or four figures.
  - **The necessary condition.** In the failing runs
    $\lvert\overline{\nabla\psi}\rvert$ at the apex node exceeds 1 right after the
    collapse: 1.60, 5.95 and 1.9. In two of the three holding runs it stays below 1.
  - But RK3 at $N/4$ also reaches 1.51, then recovers to 0.90 and holds. So
    $\lvert\overline{\nabla\psi}\rvert>1$ is necessary here, not sufficient.
  - Past $\tau\approx4.5$ the breakdown amplifies round-off. So $E_r(6)$ there is not
    a reproducible number under any integrator.
- **1D, reinit seed: wall capture** (2026-09-28, `logs/wall1d/results.txt`; the replica,
  1D, [−2, 2], natural BCs, $H{=}1/10$, the scheme of this case, from the seed
  $\psi_0=-r_f(\phi-\tfrac12)$, $r_f=0.1$).
  - **What happens.** At $N\le6$ it converges. At $N=7$ and 8, $\psi(\pm2)$ collapses when
    the front reaches the walls ($\tau\approx0.95$). The field then settles on the
    distance to $\{\pm1,\pm2\}$, with $E_r$ 0.250 against 4.5e-4 at $N{=}3$. $E_r$
    here uses $\int\lvert\psi_e\rvert$, because Eq. (84)'s $\int\psi_e$ is exactly 0 on
    this domain.
  - **The one deciding variable.** All 13 runs split on which way the wall node leans
    when the front arrives.
    - If it sits above its GLL neighbour (outflow), the wall gradient $\lvert g_b\rvert$
      peaks at 1.02–1.07 and the kink leaves the domain.
    - If it sits below (inflow), the one-sided, unaveraged endpoint gradient jumps to
      12–29. The node then falls at $\lvert g_b\rvert-1$.
    - Nothing on that node resists. $\mathbf D_\mu\approx1$, and the SVV term is
      $\le0.02$ in five of the seven failing runs. It is 0.80 and −0.09 in the two
      modified-seed runs, which fail anyway. $\psi_b$ reaches $\varepsilon$ in ~0.1 $\tau$,
      and capture (§8.5) holds it there.
  - **What sets the lean** is the arriving front's precursor at the wall.
    - In the $N$ sweep at otherwise base settings it grows at 40–44 per $\tau$ at
      $N=7$, 8 ($N=8$ with $\Delta\tau$ halved included).
    - It stays small or outflow-signed at $N\le6$.
    - The other $N=8$ variants shift the rate, from 23 (slope seed) to 69 (plateau
      seed), and all of them fail. It is tied to the front, not to time: on $[-3,3]$ the collapse waits for
    the front, at $\tau=1.95$.
  - **What it is not.**
    - It is not a linear wall mode: with the SVV, the linearised operator has none
      at any $N$.
    - It is not $\Delta\tau$: halving it changes nothing.
    - It is not the SVV amplitude: $c_0=4$ still fails. $N_{svv}=0.5$ at $N=8$ holds,
      so it is the filter's reach that matters.
    - It is not the wall value: the skewed IC starts at 9.1 at the left wall and
      converges, because it rises steeply toward both walls. A seed with slope 1 at
      the walls still fails.
  - **Periodic $[-2,2)$ at $N=8$ holds**, and so do the interior points where fronts
    meet, $x=0$ in every run. So the amplifier is the one-sided wall row.
  - **Against the apex failure above.** The sequence is the same: a front arrives
    where the characteristics end, $\lvert\nabla\psi\rvert$ overshoots 1, and the node is
    caught. But the 1D interior meeting points hold. The 2D apex, whose gradient is
    averaged, needs its own amplifier.
  - **The coupled cases are not exposed.** They are periodic, and their pseudo-time
    stops at $2.5H$ (`CDI_METHOD.md` §4.1c).
- **Questions for the authors** (drafted, not sent):
  - which $\varepsilon$ §4.4 used;
  - whether an interface-preserving correction is applied;
  - whether their $E_r(\tau)$ is flat at 6;
  - which boundary condition §4.4 uses (not stated; §5.3's zero flux is for the dam break);
  - the §4.4 `.usr`/`.par` files.

What this left open for the coupled cases (D2, D5) is in `../../NEXT_SESSION.md`.

## 11. Running it

```bash
source ../../setup-env-cuda.sh
genmeshbox -2 2 -2 2 -0.1 0 20 20 1 .false. .false. .true.  && mv box.nmsh box20.nmsh
genmeshbox -2 2 -2 2 -0.1 0 40 40 1 .false. .false. .true.  && mv box.nmsh box40.nmsh
genmeshbox -2 2 -2 2 -0.1 0 80 80 1 .false. .false. .true.  && mv box.nmsh box80.nmsh
./run.sh                      # the six committed cells, GPU-guarded
./run.sh circles_h10_n3       # or one cell
```

- The committed `.case` files are Table 2's four cells plus two $N_{svv}=N/4$ arms.
  They carry **Saini's own** $\varepsilon=H/N$.
- `logs/gen_cases.py <outdir> HN:1 <Hden>:<N> …` emits the rest of the Fig. 12
  grid. `logs/` is gitignored.
- **The drift control.** Add `"ic": "exact"` under `case.cdi` to start from
  $\psi_e$ instead of Eq. (83)'s skewed IC. The log's setup section then reads
  `initial condition : psi_e (DIAGNOSTIC)`. The default is `"skewed"`.
- Wall time per cell to $\tau=6$ on the RTX 3090: $H{=}1/5$ 0.5–2 min, $H{=}1/10$
  1–10 min, $H{=}1/20$ 4–50 min.

**Files.**
- `redistance_circles.f90` is the user file. Seventeen of its routines are shared
  byte-identically with the three coupled cases (`../../CLAUDE.md`).
- `logs/` (gitignored) holds the forensic runs by date:
  - `sep24/`: exact-IC and $\varepsilon=H$ grids;
  - `sep25/o2_audit/`: the operator audit;
  - `sep25/n1_eps_fixed/`: the fixed-$\varepsilon$ test;
  - `sep25/i1_bdfk/`: the integrator re-check, with `PREREGISTERED.txt`, `results.txt`,
    `make_bdfk.py`, `make_promote.py` (which produced the current `.f90`) and the runs;
  - `rk3_20260923/`: the SSP-RK3 runs and `.f90` the case used until 2026-09-26.
  - `sep26/`: the drift-mechanism and apex study (§7.5, §8.5, §10).
    - `PREREGISTERED.txt` and `results.txt`;
    - the Neko $N/4$ step-size runs (`q7_margin/`) and the fixed-$\varepsilon$ $N{=}8$
      runs (`n3_eps/`), made with `mkcases.py` and `chain.sh`;
    - `q6_frames_analysis.py`, the blow-up analysis.
    - The replica runs are in `capture_logs/` and `apex_logs/` in the replica
      directory.
  - `rebuild/`: the rebuild from the paper with Neko's own operators (§6.2), with
    `PREREGISTERED.txt` and `results.txt`. `rebuild_core.f90` is the solver,
    `rebuild.f90` its hooks, `make_check.py` the check build, `mkcases.py`,
    `runB.sh` and `cmp_fields.py` the runs and their comparison.
- The figure scripts are `logs/mkfigs.py`, `logs/figs.py` and `logs/fig13cmp.py`.

**Evidence.**
- `evidence/` holds the Fig. 11/12/13 recreations and the $E_r(\tau)$ traces of the
  Table 2 cells, from the BDF2 runs of 2026-09-26.
  - Fourteen of the 18 grid cells were run by the scratch build, which the committed
    build reproduces bit for bit.
- Fig. 12 adds their `plot.py` values, dashed. Fig. 13 uses their colour scale:
  0–7e-3, cool-to-warm, clipped.
- Their raster is not reproduced here (CC BY-NC-ND).
- **Animations** (2026-09-28) tell the §4.4 story to someone new to it, all with
  the scheme above (Eq. (31) SVV, BDF2/EXT2, $\varepsilon=H/N$). Scripts and run
  data: `logs/anim/` (gitignored); the Neko clips come from a scratch build that
  only adds a frame dump, and reproduce the committed logs line for line.
  - `anim1_1d.mp4`: 1D, the skewed IC relaxes onto $\lvert x\rvert-1$ ($H{=}1/10$,
    $N=3$ and 8; replica).
  - `anim2_flat_vs_curved.mp4`: from $\psi_e$, a tilted line stays exact while one
    circle and the two circles drift, to $\tau=24$ ($H{=}1/5$, $N{=}3$; replica).
  - `anim3_fig13_over_tau.mp4`: Fig. 13 over $\tau$ at $N{=}8$ on the three meshes,
    from both ICs, on their colour scale, with $E_r(\tau)$ against their value (Neko).
  - `anim3b_apex.mp4`: $H{=}1/20$, $N{=}3$, the steep-side apex failure from the skewed
    IC against the exact-start run (Neko).
  - `anim4a_1d_reinit.mp4`: 1D re-initialization in Saini's sense, the Eq. (47)
    seed then Eq. (44). At $N{=}8$ the seed's front creates a spurious zero set at
    the walls (replica; mechanism in §10, "1D, reinit seed").
  - `anim4b_blowup_h10_n7.mp4`: $H{=}1/10$, $N{=}7$ from $\tau=10$ to the blow-up,
    with the left–right and top–bottom asymmetric parts of $\psi$ (Neko).
  - `anim4c_anchor_n3.mp4`, `anim4c_anchor_n8_apex.mp4`: the Russo–Smereka anchor
    against none, ported to BDF2, at $H{=}1/10$, $N{=}3$ to $\tau=24$, and its apex
    failure at $H{=}1/5$, $N{=}8$ (replica).
  - `talk/`: the same stories for a presentation, one message and at most two
    panels per clip, in KTH style (Figtree, KTH palette, via `kthviz`). Each clip
    pauses at its key moments, and `talk/snapshots/` holds those exact frames as
    PNG stills. Made by `logs/anim/talk.py` from the same data; the detailed clips
    above are the reference. The one extra, `t4a2_not_fixes.png`, is a still from
    the wall-capture matrix (`logs/wall1d/`, §10). It shows that $N_{svv}=0.5$ and a
    periodic domain avoid the false interface at $N{=}8$. Neither is a fix: the first
    is not the printed $N/6$, and the second removes the walls.
- `visualize.ipynb` shows the $\tau=6$ error map at $H{=}1/10$, $N{=}3$ (holds) and
  $H{=}1/20$, $N{=}3$ (fails).

## History

- **2026-09-22 — the first form.** The SVV used `svv_step_imp`, with
  $\lvert\mathbf c\rvert=1$ inside the bilinear form. It acts on the zero set and
  created spurious interfaces: 14–920× off Table 2.
- **2026-09-23 — Eq. (31) as printed.** The left-multiplied form made the scheme
  sign-invariant (no node changes sign), at 2.1–2.7× Table 2. The case then used
  SSP-RK3 with a Lie-split backward-Euler SVV step. That scheme is 0.26–1.6% below the
  BDF2 one (§7.4).
- **2026-09-24 — the drift.** Found from $\psi_e$ (§8).
- **2026-09-25.**
  - Line-by-line audit against the rendered paper (§3).
  - Saini's BDF2/EXT2 integrator implemented and validated (§5–§6).
- **2026-09-26 — BDF2/EXT2 becomes the case's only scheme.**
  - The committed build reproduces the validated scratch build bit for bit in all
    four Table 2 cells.
  - `evidence/` and the notebook were rebuilt from BDF2 runs.
  - The RK3 runs are archived in `logs/rk3_20260923/` (gitignored), with the RK3
    `.f90`.
  - The $N_{svv}=N/4$ arm at $H{=}1/5$, $N{=}7$ now fails; RK3 held it (§7.5).
- **2026-09-26/28 — how the drift ends, and when the apex breaks** (`logs/sep26/`).
  - The drift ends in node capture, on an exact $\varepsilon$ clock (§8.5).
  - The apex failure starts at the front collapse, $\tau\approx1$ (§10).
  - The time-converged scheme fails at $N/4$; RK3 held it only at $\Delta\tau\ge2.5\times10^{-4}$ (§7.5).
  - The $H{=}1/10$, $N{=}7$ blow-up starts inside the left disc, on a clock common to
    both ICs (§7.4, §10).
  - Fixed $\varepsilon=0.2$ at $H{=}1/20$, $N{=}8$ judged: both bands met (§8.4).
- **2026-09-28 — Saini's $2.5H$ band eliminated (H13)** (`logs/band/`). The paper runs
  §4.4 over the whole domain, their Fig. 13 error lies mostly beyond $2.5H$, and a band
  mask leaves the drift as it is (§9).
- **2026-09-28 — rebuilt from the paper with Neko's own operators** (`logs/rebuild/`).
  The paper-only design matches the code term by term; the rebuild reproduces the
  committed runs. No shared misreading (§6.2).
- **2026-09-28 — H8 and H9′ not pursued** (decision). Their record is compressed into
  §10; no further runs are planned for either.
