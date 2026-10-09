# The method

This repo runs Neko's Conservative Diffuse Interface (CDI) equation for a phase field $\phi$, with
the compression term's interface normal taken from a **separately transported signed-distance
field** $\psi$ instead of from $\phi$'s own gradient. The $\psi$ machinery follows Saini &
Tomboulides (JCP 2026): transport with SVV (their Eq. 43), and re-initialization by the Eq. (47)
reseed followed by the Eq. (44) re-distancing equation (their Algorithm 1). The $\phi$ equation is
the one intended difference: Neko's fused CDI equation instead of their split conservative level
set (CLS). §6 lists their configuration beside ours, setting by setting.

**Attribution.** This is not "Saini's method". The transported-distance normal is, by their own
account, adopted from Al-Salami, Kamra & Hu (JCP 438, 2021) [their 29], and Eq. (44) is "the
traditional level-set (TLS) re-distancing equation" of Sussman, Smereka & Osher (JCP 114, 1994)
[their 15]. Cite Saini et al. for the equation forms and parameters we validate against.

**Cost.** CDI itself needs no re-initialization phase: the compression lives inside $\phi$'s
transport. The $\psi$ normal adds one transported scalar and, where it is maintained, one of
Saini's two periodic pseudo-time phases (Eq. 44, a fixed number of steps per event); his CLS
re-sharpening (Eq. 39) is not used.

The validated results transport $\psi$ from the exact signed distance and never redistance it
(§4). Real geometries have no analytic distance, so $\psi$ *built* by Eq. (44)
(`psi_init = "redistance"`) is the path that carries over. Why the $\phi$-gradient normal fails on
a spectral element method is the sibling repo's investigation (`../neko-multiphase/CDI_IN_SEM.md`).

## 1. Notation

$\phi$ is the phase field (Neko scalar `s`, 0 to 1, a tanh across the interface); $\psi$ is the
signed distance. **Saini et al.'s names are transposed**: their $\psi$ is the phase field and their
$\phi$ the distance. Their equation numbers are quoted as printed.

$$\xi \;=\; \frac{\varepsilon N}{H}, \qquad \varepsilon = \xi\,\frac{H}{N},$$

with $H$ the element edge and $N$ the polynomial order: **$\xi$ is the interface width in units of
$H/N$**, the resolution-independent way to state it. GLL nodes cluster at element ends, so three
spacings play three roles:

| quantity | at $N=7$, $H=1/50$ | in $H/N$ | role |
|---|---|---|---|
| $h_{\text{GLL,min}}$ (element ends) | 1.28e-3 | 0.449 | sets the time step: CFL, the compression limit, the SVV limit |
| $H/N$ (mean) | 2.86e-3 | 1.000 | what $\xi$ is defined on |
| $h_{\text{GLL,max}}$ (element centre) | 4.19e-3 | 1.465 | the worst-resolved point |

**Our $\xi$ is $N$ times Saini's.** Their Eq. (36) writes $\varepsilon=\xi H$ and they choose
$\xi=\{0.5,1,1.5\}/N$; our $\xi$ is that numerator, so our $\xi=1$ and their $\xi=1/N$ are the
same width. Never compare the bare symbols; quote $\varepsilon$ or $\varepsilon/(H/N)$.

## 2. The $\phi$ equation and how it is stepped

Neko's CDI, with $\Gamma(t)=\gamma\,u_{\max}(t)$:

$$\frac{\partial \phi}{\partial t} + \mathbf{u}\cdot\nabla\phi
= \nabla\cdot(\varepsilon\Gamma\nabla\phi) + \Gamma\,\nabla\cdot\big(-\phi(1-\phi)\,\mathbf{n}\big),
\qquad \mathbf n=\frac{\nabla\psi}{|\nabla\psi|}$$

- $u_{\max}$ is the flow's peak speed: constant on the slab and Zalesak, $|\cos(\pi t/8)|$ on
  Rider–Kothe. Both terms carry the same $\Gamma(t)$, so the equilibrium profile is a tanh of
  width $\varepsilon$ at every $\gamma$.
- **The diffusion is physical.** It balances the compression to hold that profile; it is not
  numerical stabilisation and is never replaced by SVV (§5).
- **$\gamma$ is a rate.** It scales both terms, so it sets how fast the interface relaxes, not to
  what. $\gamma=0$ removes the compression. On Zalesak $E_r$ still rises with $\gamma$, by 1.5–2×
  from 0.25 to 2 at $\xi\ge1.5$: a rigid rotation's exact solution contains no relaxation
  (`examples/zalesak_disk/README.md`).

Saini's CLS, for comparison only (their Eqs. 38–39, our symbols): transport
$\partial_t\phi+\mathbf v\cdot\nabla\phi=S_{vv}(\phi)$, plus a periodic pseudo-time re-sharpening
$\partial_\tau\phi+\nabla\cdot(\phi(1-\phi)\mathbf n)=\nabla\cdot(\varepsilon(\nabla\phi\cdot\mathbf n)\mathbf n)+S_{vv}(\phi)$.
His phase field reads the normal only in that re-sharpening, inside a band; ours reads it pointwise
everywhere, every step (§3).

**Stepping.** Neko's scalar scheme, BDF3/EXT3 (`time_order` 3 in every case; Saini uses BDF2/EXT2,
§6). The advection and the compression (a user source term) are explicit and extrapolated; the
diffusion is implicit, $\varepsilon\Gamma$ being the scalar's conductivity. So the step is limited
by the explicit terms:
- **The compression limit.** $C_{\text{comp}}=\Gamma\Delta t/h_{\text{GLL,min}}\le0.05$ is a hard
  error in the three coupled files, and `variable_timestep` is refused, because Neko would size
  $\Delta t$ from the advective CFL alone. The 0.05 is a safe value, not a measured stability
  boundary. The shipped cases run at up to 0.0485. This limit, not the advective CFL, sets
  Rider–Kothe's $\Delta t$: $8\times10^{-5}$ at $H=1/64$ and $\gamma=1$, against Saini's $4\times10^{-4}$.
- **The SVV limit.** With an explicit `svv_psi`, $\Delta t\,\rho\le0.4$ as well (BDF3/EXT3 carries a
  real negative eigenvalue explicitly up to 0.952). It binds first only at low $\gamma$ (Zalesak at
  $c_0=1$: below $\gamma\approx0.3$ at $N=7$, 0.35 at $N=5$).
- **An event that replaces $\psi$ restarts the history.** The coupled files replace $\psi$ in the
  `compute` hook, after Neko has updated its lags; without a restart BDF3 settles at
  $\psi_{old}+\tfrac{11}{6}(\psi_{new}-\psi_{old})$. They set `nadv = ndiff = 0`, as Saini's
  `ireset_ls` does, so the next step is BDF1/EXT1.

## 3. $\psi$ and the normal

**Transport** is Saini's Eq. (43): pure advection, with SVV as the only dissipation ($\psi$'s
conductivity is $10^{-16}$):

$$\frac{\partial \psi}{\partial t} + \mathbf{u}\cdot\nabla\psi = S_{vv}(\psi)$$

**The normal is read everywhere, every step.** Saini's scheme reads $\mathbf n$ only inside a
band-local relaxation (§2); the fused CDI term reads it at every node. Where $\psi$ is flat (outside
a built band) $|\nabla\psi|$ is round-off, about $10^{-19}$, and dividing by it gives a unit vector
in a random direction. Its divergence, about $1/h$, turns the compression term into a source
proportional to $\phi$ with a rate near $\pm2\times10^{3}$, which diverges. `unit_normal` therefore
floors $|\nabla\psi|$ at `case.cdi.grad_floor` $=10^{-6}$, so $\mathbf n\to0$ where there is no
gradient. **The floor is load-bearing**: at $10^{-30}$ a built $\psi$ diverged on Zalesak at $N=5$, $t\approx1.17$
(`examples/zalesak_disk/run_fx_none_n5.log`). For an analytic $\psi$, which is never flat, the floor
changes nothing: the run is bit-identical (measurements in `git show 80201cf13e6:CDI_METHOD.md`, §4.1).

**The initial $\psi$** (`case.cdi.psi_init`):
- `"exact"` (default): the analytic (periodic) signed distance. Every validated result uses it.
- `"redistance"`: Saini's Algorithm 1 line 2, Eq. (44) from the seed $r_f(\phi-\tfrac12)$. The
  analytic distance is then never evaluated. Its one call site sits in the `else` branch of
  `initial_conditions`.

Never pair `"exact"` with periodic redistancing: an event replaces a global distance with a field
that is a distance only near the interface.

## 4. Maintaining $\psi$

Saini maintains $\psi$ by **re-initialization**: every $\Delta t_{tls}$, reseed from the phase field
by Eq. (47) and relax by Eq. (44). The terms, what is known, the rules and the code are in
[`REDISTANCING.md`](REDISTANCING.md). What each case does:

| case | $\psi$ at $t=0$ | maintenance |
|---|---|---|
| `advecting_slab_1d` | exact | none, except the `psi_rd_*` testbed (built $\psi$; off, or 40 events with `seed` `"phi"` or `"psi"`) |
| `zalesak_disk` | exact (periodic) | none, except two variants (built $\psi$; timer 0.5, `seed` `"phi"` or `"psi"`) |
| `rider_kothe` | exact (periodic) | none in the shipped cases; re-initialization in Saini's configuration ran on scratch files only |
| `redistance_circles` | the problem's IC | Eq. (44) alone, standalone |

## 5. SVV

**SVV is a $\psi$-only knob.** $\phi$ carries the CDI equation; $\psi$ carries transport plus SVV.
There is no SVV on $\phi$: the coupled files stop at startup if a case sets `case.cdi.svv_phi`. SVV
on $\psi$ is always on: they stop if `normal = "psi"` has `svv_psi` off. Saini does the same, with
SVV on every scalar equation, his TLS transport included. Saini's names are transposed (§1), so his
§4.3 SVV sits on the phase field; we do not adopt that placement. SVV is not merged with $\phi$'s
diffusion either: that term defines the equilibrium profile, while SVV is a grid-scale filter. In
the sibling repo, SVV on $\phi$ only delayed the $\phi$-normal's growth, which lives in Legendre
modes 0–2, where SVV cannot act.

Two instances, with Saini's parameters:

| instance | $N_{svv}$ | $c_0$ | applied |
|---|---|---|---|
| `svv_psi`, $\psi$ transport (Eq. 43) | $N/2$ | 0.1, his TLS-transport value (§4.5, p. 21; §5 default). Zalesak still ships 1.0, his §4.3 value, which sits on his CLS transport (open, `NEXT_SESSION.md`) | explicit, through the source term |
| `svv_rd`, Eq. (44) | $N/4$ | 2.0, his §4.5 value | implicit, forced in `startup` (at $c_0=2$ the spectral radius would otherwise set $\Delta\tau$) |

His values per test, as printed (paper §4–§5; his public cases differ in places):
- **Transport:** $N/2$ (§4.1–§4.2 also try $N/4$), $c_0$ 0.1, except §4.3's 1.0 on his CLS transport
  and §4.1's $c_0$ sweep.
- **CLS re-sharpening:** $N/2$, 0.1 in §4.5; $N/4$, 1.0 by default from §5.1.
- **Eq. (44):** $N/6$, 2.0 in §4.4; $N/4$, 2.0 in §4.5, §5.2 and §5.3 ($N/6$ for the dam break at
  $N=9$); $N/6$, 1.0 as the §5.1 default.

**The viscosity** is $\nu=c_0|\mathbf c|\mathbb H/N$ (his Eq. 29) with $\mathbb H=2J_e^{1/d}$. We use
the element edge, which equals it for $d=2$ on these one-element-thick meshes; $d=3$ would fold in
the $z$ extent. Two deviations, both in the coupled files:
- **$|\mathbf c|$ is uniform.** Saini evaluates the local characteristic speed (Eq. 31). The
  coupled files use the flow's $u_{\max}$, taken once when the instance is built, for transport and
  for Eq. (44). On Rider–Kothe, Saini's local $|\mathbf u|$ on `svv_psi` alone gives $E_r(8)$
  0.0470 against 0.0452 (scratch user file, transport only).
- **Where $\nu$ sits.** Saini's Eqs. (7), (31), (33) left-multiply the assembled operator by a
  pointwise $\mathbf D_\mu$. `svv_local` applies $\nu$ inside the bilinear form. The two agree only
  for constant $\nu$. For Eq. (44) they do not: $|\mathbf c|=|\operatorname{sgn}\psi|$ vanishes on
  the zero set, and only the left-multiplied form keeps the SVV from moving the interface.
  `redistance_circles` uses the printed form (`svv_step_eq31`), and only it. The coupled files still
  use `svv_step_imp`; on Rider–Kothe it creates the events' spurious zero sets
  (`REDISTANCING.md` §3). Don't vary $\nu$ inside `svv_local` and call it Eq. (31).

**What SVV on $\psi$ buys.** On Zalesak, without it $E_r$ is 9–256× worse and $p$-refinement
reverses: a higher order gives a worse answer (`examples/zalesak_disk/README.md` §3).

## 6. Saini's configuration and ours

The $\psi$ machinery is meant to be Saini's. Sources: the paper (JCP 2026; page or equation), his
code (`nandu90/nekLS_Examples@jcp` with `nandu90/Nek5000@nekLS`; a local copy is in
`examples/rider_kothe/logs/d5/saini/`) and his email of 2026-09
(`references/saini_email_2026-09.md`, gitignored). Case-specific values (meshes, $\Delta t$, step
counts) are in each case README.

| setting | Saini | coupled files (committed) | `redistance_circles` |
|---|---|---|---|
| phase-field equation | split CLS: Eq. (38) transport with SVV, plus Eq. (39) re-sharpening every $\Delta t_{cls}$ (0.05 in §4.5); the phase field reads $\mathbf n$ only in Eq. (39) | fused CDI (§2), $\mathbf n$ read at every node every step. **The intended difference** | no phase field |
| physical time scheme | BDF2, extrapolated with Nek's EXT3 coefficients. The paper prints BDF3/EXT3 for §4.5 at $H=1/128$ (p. 21), but his public `circVortex/cv.par` uses BDF2 at every resolution | BDF3/EXT3 in every case. Neko's `time_order` 2 is BDF2 with a modified EXT3 | — |
| dealiased transport | yes, "for all my cases" (email) | yes, `dealias: true` in every `.case` | — |
| $\Delta t$ | advective CFL (§4.5: $\approx0.4$ at $N=6$) | CDI's compression limit (§2) | — |
| `svv_psi` | $N/2$, $c_0=0.1$, local $|\mathbf u|$, left-multiplied, implicit in the Helmholtz operator | $N/2$, 0.1 (Zalesak and six slab cases 1.0); explicit, through the source term; uniform $u_{\max}$ inside the bilinear form (§5) | — |
| $\psi$ at $t=0$ | built by Eq. (44) from the Eq. (47) seed (Algorithm 1, line 2) | `psi_init`: `"exact"` (validated runs) or `"redistance"` | the problem's IC, Eq. (83) |
| re-initialization | every $\Delta t_{tls}\in[0.01,0.5]$ (0.5 in §4.5, "including at the initial step"), always Eq. (47) then Eq. (44) | off in the validated cases; on, on the timer, in Zalesak's two variants and three slab testbed cases. `seed = "phi"` is his operation, `"psi"` our in-place variant | Eq. (44) alone |
| reseed, Eq. (47) | $r_f(\phi-\tfrac12)$, $r_f=0.1$; his code scales it by $L$, the smallest domain extent | $r_f(\phi-\tfrac12)$, $r_f=0.1$ | — |
| sign width in Eq. (46) | **0.25**, fixed, applied to $\psi/L$ (his `signls`; not printed; intentional and recommended, email) | $\varepsilon$, the phase-field width | 0.25 |
| pseudo-time extent | $2.5H$ printed (§3.4); **$25H$** in his `circVortex` | `band`$\times H$, band 2.5 | $\tau=6$ |
| $\Delta\tau$ | $H/(N{+}1)$ (§3.4): 150 steps over $25H$ at $N=5$, each capped at a Nek CFL of 1 (no subcycling) | $0.1\,h_{\text{GLL,min}}$ (`cfl` 0.1; `dtau` overrides) | $5\times10^{-4}$ |
| pseudo-time scheme | BDF2/EXT2, Eqs. (34)–(35), SVV implicit in the Helmholtz operator; his code extrapolates with Nek's EXT3 coefficients | SSP-RK3, then one implicit SVV step (Lie split) | BDF2/EXT2 |
| Eq. (44) SVV | printed Eq. (31): $\mathbf D_\mu=c_0|\operatorname{sgn}\psi\,\mathbf n|\mathbb H/N$ left of the assembled operator, $\mathbf n=0$ where $|\nabla\psi|<10^{-12}$; §4.5 $N/4$, $c_0=2$; §4.4 $N/6$, $c_0=2$ | `svv_step_imp`: uniform $\nu=c_0u_{\max}H/N$ inside the bilinear form; $N/4$, $c_0=2$ | `svv_step_eq31`: the printed form, $N/6$, $c_0=2$; $\mathbf D_\mu$ omits $|\mathbf n|$, which differs only where $|\nabla\psi|<10^{-12}$ |
| $\mathbf C(\mathbf w)\psi$ | dealiased on $\lfloor3(N+1)/2\rfloor$ Gauss points (`convect_new`; email) | not dealiased: `rd_rhs` evaluates $\operatorname{sgn}(\psi)\lvert\nabla\psi\rvert$ on the GLL points, equal to the GLL Galerkin $\mathbf B^{-1}\mathbf C(\mathbf w)\psi$ to $10^{-14}$ | dealiased (`adv_dealias_t`) |
| sign guard | `constrainTLSR`: a node whose $\psi^n$ disagrees in sign with $\phi-\tfrac12$ keeps $\psi^n$ (email; `lvlSet.f`) | none | the same, against the sign of $\psi_0$ |
| convergence test | none; a fixed step count | none; a fixed `rd_niter` | none |
| history after an event | restarted (`ireset_ls`) | restarted (§2) | — |
| $E_r$ on Rider–Kothe | his `ls_relerr` divides by the area *outside* the disk, about 13× ours; recompute on his dumps with `examples/norms.py` | Eqs. (79)–(81) as printed | Eq. (84), $\int\psi_e$ |

`redistance_circles` uses his configuration and reproduces his §4.4 values in 16 of the 18 Fig. 12
cells; two low-$N$ cells fail on our mesh through zero-set round-off (its README §5). The coupled files'
events path is not his configuration; adopting it on Rider–Kothe is `NEXT_SESSION.md` item 1.

## 7. Non-goals

- Whether CDI works with the $\phi$-gradient normal at all. It does not on a spectral element
  method; that is the sibling repo's investigation (`../neko-multiphase/CDI_IN_SEM.md`).
- Jain's (2022) algebraic normal: tried there, and it does not reach both accuracy and boundedness
  (`CDI_IN_SEM.md` §6.2).
- The strong- vs weak-form question for the compression term
  (`../neko-multiphase/examples/cdi_face_artifact/`). This repo inherits the strong form.
- Parameter sweeps beyond what each case README documents. The full $\xi$/$\Gamma$/SVV ladder is
  in `../neko-multiphase/examples/saini_benchmarks/`.
