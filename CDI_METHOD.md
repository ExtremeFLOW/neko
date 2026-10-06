# The method this repo showcases

This repo demonstrates Neko's Conservative Diffuse Interface (CDI) method for
two-phase interface capture, once its compression term is given an interface
normal it can actually use: computed from a **separately transported
signed-distance field**, $\psi$, instead of straight from the sharp phase
field's own gradient. The method also supports *periodic redistancing* of
$\psi$ (§3–§4). Note up front that **the validated results predate it**: in
all three, $\psi$ is seeded once from the exact analytic signed distance and
thereafter only advected. That is a validation baseline, not a
recommendation — **real geometries have no analytic distance**, so
`psi_init = "redistance"` (§4.1), which builds $\psi$ from $\phi$ by Saini's
Eq. (44) and never evaluates an analytic distance, is the path that carries
over.

The normal fix — and every upstream number quoted below — was found and
measured in the sibling investigation repo, `neko-multiphase` (see
`../neko-multiphase/CDI_IN_SEM.md` §0 for the full derivation). This repo
exists to show the working result cleanly, not to re-run the investigation.

**A naming note, deliberately:** this repo is not called anything with
"Saini" in it, and `CDI_METHOD.md` is not "Saini's method." Neko's CDI
equation (below) is structurally different from the method in Saini et al.
(2026) — theirs splits transport and re-sharpening into two separate
equations with SVV in each; CDI fuses compression into the transport equation
itself, every step, no SVV needed there. What *is* shared is narrower: the
idea of computing the interface normal from a separately transported,
periodically redistanced signed-distance field. Saini et al. themselves
attribute that idea to earlier work — their own text (§3.1) reads: *"In this
work we adopt a similar approach as described by Salami et al. [29]"*, and
the redistancing PDE itself (their Eq. 44) is introduced as *"the traditional
level-set (TLS) re-distancing equation [15]"* — i.e. Saini's own paper treats
it as a classical technique (almost certainly Sussman-style
reinitialization), not something they originated. Their [15] is Sussman,
Smereka & Osher, J. Comput. Phys. 114 (1994) 146–159, and their [29] is
Al-Salami, Kamra & Hu, J. Comput. Phys. 438 (2021) 110376. Cite the technique as
"a transported-signed-distance-field normal with periodic redistancing," and
cite Saini et al. (2026) only for the specific equation forms and parameters
this repo's implementation was validated against — not as the technique's
origin.

**Treat status claims as provisional:** check them against the case files, logs
and notebook outputs (§7 and each case's `README.md`).

## 1. Symbol convention

$\phi$ = the phase field (Neko scalar `s`, 0–1, tanh profile across the
interface). $\psi$ = the signed-distance field. **This is the opposite of
Saini et al.'s own notation** (they call the phase field $\psi$ and the
distance $\phi$) — their equation numbers read transposed against this
project's code and this doc. Jain (2022)'s algebraic surrogate for a
distance-like field, $\varepsilon\ln\!\big(\tfrac{\phi+\delta}{1-\phi+\delta}\big)$,
is called `jain` here and never $\psi$, since it's a pointwise transform of
$\phi$, not a transported field — and it does not work well enough to use
(see §9).

### $\xi$: the interface width, and a second collision with Saini's notation

Everywhere in this repo,

$$\xi \;=\; \frac{\varepsilon N}{H}, \qquad\text{equivalently}\qquad \varepsilon = \xi\,\frac{H}{N},$$

with $H$ the element edge length and $N$ the polynomial order. In words: **$\xi$
is the interface width $\varepsilon$ measured in units of $H/N$** — the nominal
node spacing, element edge divided by polynomial order. It is the one
resolution-independent way to state an interface width, which is why every
sweep here varies $\xi$ rather than $\varepsilon$.

GLL nodes are not uniformly spaced — they cluster at element ends — so there are
**three** spacings with three distinct roles, and none of them is "the GLL
spacing":

| quantity | at $N=7$ | in units of $H/N$ | role |
|---|---|---|---|
| $h_{\text{GLL,min}}$ (element ends) | 1.28e-3 | 0.449 | sets the **timestep** — CFL and the SVV explicit limit |
| $H/N$ (mean) | 2.86e-3 | 1.000 | what **$\xi$ is defined on** |
| $h_{\text{GLL,max}}$ (element centre) | 4.19e-3 | 1.465 | the **worst-resolved point** — sets how thin $\varepsilon$ can usefully go |

All three appear within a few lines of the same startup report. Saini et al.
write $\epsilon = \xi H$ and never call $H/N$ a GLL spacing at all; they use GLL
only for the quadrature.

**How thick is a given $\xi$, concretely.** For $\phi = \tfrac12(1+\tanh(d/2\varepsilon))$
the 5%–95% band has width $4\,\text{artanh}(0.9)\,\varepsilon = 5.89\varepsilon$. At
$H=1/50$, $N=7$:

| $\xi$ | 5–95% band | in $H/N$ | in **elements** | nodes across the band at $h_{\text{GLL,max}}$ |
|---|---|---|---|---|
| 0.5 | 8.4e-3 | 2.9 | 0.42 | **2.0** |
| 1.0 | 1.7e-2 | 5.9 | 0.84 | 4.0 |
| 1.5 | 2.5e-2 | 8.8 | 1.26 | 6.0 |
| 2.0 | 3.4e-2 | 11.8 | 1.68 | 8.0 |
| 2.8 | 4.7e-2 | 16.5 | **2.36** | 11.3 |

This is the clearest statement of what $\xi$ means physically, and it explains
the stability boundary measured in `examples/zalesak_disk`: at $\xi=0.5$ the
interface spans **two nodes** at the worst-resolved point and the run diverges;
at $\xi=1$ it is four and violations are small but real; from $\xi=1.5$ up it is
six or more and the field is bounded to round-off. It also shows that the
$\xi=2.8$ reference configuration smears the interface across **more than two
whole elements** — see §7 and that case's README for what that costs.

**$\xi$ here is $N$ times Saini et al.'s $\xi$, and the numbers are a trap.**
Their Eq. (36) defines $\epsilon = \xi H$ with $H$ the *maximum* element edge,
so their $\xi$ is a fraction of an element; they then choose
$\xi = \{0.5, 1, 1.5\}/N$. This repo's $\xi$ is the numerator of that
expression. So:

The values coincide with their **numerator**, so "$\xi = 1.5$" here and
"$\xi = 1.5/N$" in their paper denote the **same physical interface width** —
but read as bare numbers, theirs is $N$ times smaller. Quote $\varepsilon$ or
$\varepsilon/(H/N)$ explicitly in any comparison, never the symbol alone. This is
the second notational collision with that paper, alongside the $\phi$/$\psi$
transposition above.

**Why we keep our form rather than adopting theirs.** Saini never use their bare
$\xi$ as a controlled variable — they always write it as $c/N$ and vary $c$
(their Eq. 36 discussion; $c = \{0.5, 1, 1.5\}$ for Zalesak). Their effective
knob *is* our $\xi$. Ours also reads directly as a resolution criterion — "the
interface spans $\xi$ nominal node spacings" — independent of $N$, which is what
makes this repo's $h$- and $p$-refinement rows commensurable at a glance.

### Where our $\xi$ values sit against theirs

At $H = 1/50$, $N = 7$ — a setup both this repo and they ran:

| our $\xi$ | $=$ their $\xi_S$ | $\varepsilon$ | in Saini et al. | measured here |
|---|---|---|---|---|
| 0.5 | $0.5/N = 0.0714$ | 1.43e-3 | sharpest studied | **diverges** at $\gamma \ge 1$; $2.5\times10^{-2}$ otherwise |
| 1.0 | $1/N = 0.1429$ | 2.86e-3 | **robust**, their default | bounded only to $4\times10^{-6}$–$6\times10^{-5}$ |
| 1.5 | $1.5/N = 0.2143$ | 4.29e-3 | **robust**, upper | **exactly bounded**, $2\times10^{-10}$ |
| 2.0 | 0.2857 | 5.71e-3 | beyond their range | exactly bounded, $4\times10^{-10}$ |
| 2.8 | 0.4000 | 8.00e-3 | beyond their range | exactly bounded, $7\times10^{-10}$ |

The mapping is confirmed independently by their constant-thickness control
$\epsilon_c = 1/150$, which is $\xi_S = 1/N$ at $H=1/50$, $N=3$ — our $\xi = 1.0$
at $N=3$, i.e. $6.667\times10^{-3}$ exactly.

Two points worth carrying:

- **Our sharpest exactly-bounded setting, $\xi = 1.5$, is exactly the top of
  their robust range.** Independent agreement, not a tuned result.
- **Where we differ, their own text explains it.** They do not claim
  boundedness — "the high order spatial discretization and EXTk/BDFk temporal
  discretization does not guarantee boundedness of $\psi \in [0,1]$", with
  overshoots "within 2%". So $\xi = 1$ failing to be *exactly* bounded here is
  consistent with their account, and our $6\times10^{-5}$ is ~300× tighter than
  their 2%. Their sharpest setting diverging for us at $\gamma \ge 1$ is a real
  difference, plausibly because their scheme splits transport and re-sharpening
  with SVV on both, while Neko's CDI fuses compression into transport (§2).

## 2. The equation actually run here, vs. the split scheme it's compared against

Neko's CDI (fused, this repo), with $\Gamma(t)=\gamma\,u_{\max}(t)$:
$$
\frac{\partial \phi}{\partial t} + \mathbf{u}\cdot\nabla\phi
= \nabla\cdot(\varepsilon\Gamma\nabla\phi)
+ \Gamma\,\nabla\cdot\!\big(-\phi(1-\phi)\,\mathbf{n}\big)
$$

$u_{\max}$ is the flow's peak speed. It is constant in the slab and Zalesak, and
$|\cos(\pi t/8)|$ in Rider–Kothe. Both halves carry the same $\Gamma(t)$, so the
equilibrium width stays $\varepsilon$; $\gamma$ is dimensionless (§6).

Saini's CLS (split, for comparison only — not implemented here):
$$
\text{transport:}\quad \frac{\partial \phi}{\partial t} + \mathbf{v}\cdot\nabla\phi = S_{vv}(\phi)
\qquad\text{(their Eq. 38, our symbols)}
$$
$$
\text{re-sharpening:}\quad \frac{\partial \phi}{\partial \tau} + \nabla\cdot\!\big(\phi(1-\phi)\mathbf{n}\big)
= \nabla\cdot\!\big(\varepsilon(\nabla\phi\cdot\mathbf{n})\,\mathbf{n}\big)
\qquad\text{(their Eq. 39, pseudo-time, periodic)}
$$

The $\varepsilon\Gamma\nabla\phi$ diffusion term in Neko's CDI equation is
**physical, not a stabilization add-on**: it's what balances the compression
term to hold a $\tanh$ profile of half-width $\varepsilon$ at equilibrium. It
is not a free knob for numerical hygiene (see §5).

**How it is stepped.** Neko's scalar scheme is BDF3/EXT3 (`time_order` 3 in every case).
- **Explicit (extrapolated):** the advection by the prescribed velocity, and the compression
  term, which the user file adds as a source.
- **Implicit:** the diffusion. $\varepsilon\Gamma$ is the scalar's conductivity, inside the
  Helmholtz solve.
- **A time-varying conductivity has to reach the solver.** The solve reads `s_lambda_tot`, which
  Neko copies from `s_lambda` only at initialisation unless a turbulence model is set. So the
  user file's `material_properties` fills both. Rider–Kothe filled only `s_lambda` until
  2026-10-02, which froze its diffusion (§4.1d).

So the step is limited by the explicit terms, not the diffusion:
- **The compression guard.** $C_{\text{comp}}=\Gamma\Delta t/h_{\text{GLL,min}}\le0.05$ is a
  hard error in all three coupled files. `variable_timestep` is refused, because Neko sizes
  $\Delta t$ from the advective CFL alone and would step straight past it.
  - The 0.05 is called "empirically measured" in the sibling repo's READMEs, but no
    measurement is cited there or here. So it is a safe value, not a known stability boundary.
  - The shipped cases run at up to 0.0485; Rider–Kothe at 0.044–0.045.
- **The SVV limit.** With an explicit `svv_psi`, $\Delta t\,\rho\le0.4$ also applies. BDF3/EXT3
  carries a real negative eigenvalue explicitly up to 0.952; the rest of the budget is left to
  advection and compression. It binds first only at low $\gamma$ (Zalesak, $\gamma<0.27$ at
  $c_0=1$).
- **What implicit diffusion buys.** An explicit diffusion would need roughly
  $\Delta t\lesssim h_{\text{GLL,min}}^2/(\varepsilon\Gamma)$, which scales as $1/(\xi\gamma)$.
  The compression limit scales as $1/\gamma$. So implicitness matters only at large $\xi$. On
  Rider–Kothe at $\xi=1$, $h_{\text{GLL,min}}^2/(\varepsilon\Gamma)$ is about $13\Delta t$.
- **What would buy larger steps:** measuring where the explicit compression actually goes
  unstable, not changing how the diffusion is treated.

The interface normal is where the two schemes actually meet:
$$
\mathbf{n} = \frac{\nabla\psi}{|\nabla\psi|}
$$
with $\psi$ built and maintained as described in §3 — not
$\nabla\phi/|\nabla\phi|$, which is the baseline that fails
(`../neko-multiphase/CDI_IN_SEM.md` §0, §9).

## 3. The $\psi$ field: transport, redistancing, re-seeding

Three distinct operations, each with a distinct name — see §4 for why they're
not interchangeable, and why they're not free of a real HPC cost either.

**Transport** (pure advection, no compression, no physical diffusion):
$$
\frac{\partial \psi}{\partial t} + \mathbf{u}\cdot\nabla\psi = S_{vv}(\psi)
$$
$\psi$'s material diffusivity is effectively zero ($\lambda \sim 10^{-16}$);
SVV, if enabled, is its *only* dissipative term (§5).

**Redistancing** (pseudo-time PDE relaxation, restores $|\nabla\psi|\approx 1$):
$$
\frac{\partial \psi}{\partial \tau} + \mathbf{w}\cdot\nabla\psi = \operatorname{sgn}(\psi) + S_{vv}(\psi),
\qquad \mathbf{w} = \operatorname{sgn}(\psi)\,\mathbf{n}
$$
which is algebraically
$$
\frac{\partial \psi}{\partial \tau} = \operatorname{sgn}(\psi)\big(1 - |\nabla\psi|\big) + S_{vv}(\psi),
\qquad \operatorname{sgn}(\psi) = \tanh\!\left(\frac{\psi}{2\varepsilon}\right)
$$
**The $\varepsilon$ in that sign function is the same symbol as $\phi$'s
interface width, and that is not obviously right.** Saini define it once, at
Eq. (36), as $\epsilon = \xi H$ — a *phase-field smearing width* — and Eq. (46)
reuses it. But §4.4 solves Eq. (44) with no phase field present at all, and there
$\varepsilon$ is the only thing setting the collapse rate
$(1-\lvert\nabla\psi\rvert)/2\varepsilon$. Their code in fact uses a fixed 0.25 there,
which the paper does not state; `redistance_circles` ships it (§4.4), and whether the
coupled cases follow is D5 (`NEXT_SESSION.md`).
Integrated with **SSP-RK3** in pseudo-time. Writing $L(\psi)$ for the
right-hand side above, the three stages are

$$\psi^{(1)} = \psi^{n} + \Delta\tau\,L(\psi^{n})$$
$$\psi^{(2)} = \tfrac{3}{4}\,\psi^{n} + \tfrac{1}{4}\Big(\psi^{(1)} + \Delta\tau\,L(\psi^{(1)})\Big)$$
$$\psi^{n+1} = \tfrac{1}{3}\,\psi^{n} + \tfrac{2}{3}\Big(\psi^{(2)} + \Delta\tau\,L(\psi^{(2)})\Big)$$

Each stage is one explicit Euler step, convexly averaged back toward $\psi^n$.
That is exactly what **strong-stability-preserving** means: the scheme is built
only from Euler steps and convex combinations of them, so whatever bound Euler
respects at a given $\Delta\tau$ — monotonicity, a maximum principle, no new
extrema — the full third-order scheme respects at the *same* $\Delta\tau$ (the
SSP coefficient of this tableau is 1). It buys third-order accuracy without
buying the spurious oscillations a general third-order scheme would allow near
the kinks, which this equation has by construction along the medial axis.

Explicit Euler is not an option: a central, continuous-Galerkin advection
operator has an essentially imaginary spectrum, and Euler's stability region
touches the imaginary axis only at the origin, so it is unstable there at every
$\Delta\tau$. SSP-RK3's region contains an interval of the imaginary axis
($|\Delta\tau\,\lambda| \le \sqrt3$), which is what makes the hyperbolic stage
integrable at all.

Neko has this scheme built in — `rk_order = 3` in
`src/time_schemes/runge_kutta_scheme.f90`, in Butcher form ($A_{21}=1$,
$A_{31}=A_{32}=\tfrac14$, $b=(\tfrac16,\tfrac16,\tfrac23)$, which is the tableau
the three stages above expand to) — but that class is wired only into the
compressible fluid solver, and the redistancing relaxation is not a Neko-stepped
equation: it runs to convergence inside a user hook. So each user file writes the
three stages out directly, in the Shu–Osher form above. **The integrator is not
what fails on Saini §4.4.** Their own BDF2/EXT2, and BDF3/EXT3, run in Neko
within 1.9% of SSP-RK3 in the Table 2 cells to $\tau=24$, and BDF2 within 2.2% over the
Fig. 12 grid (`examples/redistance_circles/archive/README_process_2026-09.md` §7.2, §7.4).

Pseudo-timestep
$\Delta\tau$ is set from $h_{\text{GLL,min}}$ (not the paper's own $H/(N+1)$
formula, which is a pseudo-CFL of 2.75 on `advecting_slab_1d`'s mesh —
$H=1/10$, $N=10$ — and is unstable under repeated events, §4.2);
iteration count set to cover a band of $2.5H$ (element edges) from the
interface, matching the "banded, not global" scope the technique is meant to
run at.

**Re-seeding** (only happens as part of a redistancing event, and only if
`seed="phi"`): discards the transported $\psi$ and rebuilds the pseudo-time
initial condition from the phase field,
$\psi_0 = r_f(\phi - 0.5)$, $r_f = 0.1$.

## 4. Redistancing is not reinitialization — and it is not free

> **See [`REDISTANCING.md`](REDISTANCING.md)** for when to use which: what the
> normal actually requires, why the validated results need neither, why a
> $\psi$ *built* by Eq. (44) is the path that carries over, and what the
> ablation measured. Its §9 walks through the Fortran that solves Eq. (44).

**What Saini's own words are**, because this repo's split does not map onto them
cleanly. Their Eq. (44) is *"the traditional level-set (TLS) re-distancing
equation"* — "re-distancing" names the **PDE**. The *operation* that uses it is
called **re-initialization** (Algorithm 1, line 19's margin note), and it is
always both steps together: line 20 reads *"Re-initialize $\phi^{n+1}$ using
Eq. (44), with the initial condition from $\psi^{n+1}$, Eq. (47)"* — reseed from
the phase field, **then** relax. Line 2 initialises the same way.

**Saini have no in-place mode.** The transported $\phi$ (their Eq. 43, line 5 of
the loop) is *discarded* at every reconstruction; it only has to survive one
$\Delta t_{tls}$ interval, never the whole run. So `seed="psi"` is this repo's
invention with no counterpart in their method, and `seed="phi"` is their
Algorithm 1 as written.

That distinction is worth a number. On `advecting_slab_1d`, 40 events, everything else held
fixed, worst violation over every frame (`examples/advecting_slab_1d/README.md`):

| | $E_r$ | $E_v$ | worst violation |
|---|---|---|---|
| no redistancing at all | 0.00006 | $1.7\times10^{-12}$ | $1.2\times10^{-10}$ |
| `seed="psi"` — relax in place | **0.00006** | $-8.1\times10^{-10}$ | $1.2\times10^{-10}$ |
| `seed="phi"` — Eq. (47) reseed, then relax | **0.00006** | $-3.1\times10^{-7}$ | $1.4\times10^{-10}$ |

**On an isometry neither half costs anything.** Until 2026-10-02 the reseed row read 0.00023, "4×
in $E_r$ and 3600× in boundedness". That was the event's stale-history bug (§4.1d), not the
reseed. Where $\phi$'s interface is itself wrong, the reseed does cost: on Zalesak and
Rider–Kothe (§4.1c–d).

**But `seed="psi"` is a diagnostic, not a method, and should not be read as a
recommendation.** Saini state the reseed's purpose three times, and it is not
noise reduction — it is *registration*: the TLS field is rebuilt from the CLS
field *"ensuring congruence in the interface location, as given by the zero and
0.5 isocontour of the respective fields"* (§3.1; also §1 and the conclusions).
Relaxing in place **cannot** do that by construction, because Eq. (44)'s
$\operatorname{sgn}(\psi)$ pins $\psi$'s *own* zero contour: whatever drift has
accumulated between $\psi = 0$ and $\phi = 0.5$ is preserved, not corrected. On an
isometry both fields stay exact, no drift exists, and the option looks free —
which is precisely why the table above does not generalise. The value of
`seed="psi"` is that it isolates which half of Saini's operation costs what; it
is not a candidate for a case with strain, where re-registration is the whole
point.

With that said, these get conflated in the literature generally, and this repo
draws its own line explicitly:

- **Reinitialization** = discarding the current $\psi$ and rebuilding it from
  another source (the re-seed above, or an exact analytic distance at
  $t=0$). A hard reset.
- **Redistancing** = the pseudo-time PDE relaxation (§3) that pushes an
  *existing* field toward $|\nabla\psi|=1$. A soft correction.

You can redistance **without** reinitializing — iterate the currently
transported $\psi$ in place (`seed="psi"`) — or reinitialize **and then**
redistance (`seed="phi"`, Saini's full Eq. 47 procedure). The in-place option
exists specifically to break a measured feedback loop: re-seed from $\phi$ →
redistance sets $\mathbf{n}$ → $\mathbf{n}$ perturbs $\phi$ via the
compression term → the next event re-seeds from a $\phi$ that's already been
perturbed by the last one. The "$\sim$1.34× gain per event" once measured in the 1D testbed
was the stale-history bug of §4.1d. With the order restart, $\lVert dn\rVert$ there stays at
$10^{-7}$–$10^{-5}$. The loop is real where $\phi$'s interface fragments (Zalesak, §4.1c).

**Why this matters beyond correctness.** One of CDI's original selling
points over classical level-set/CLS methods is that it needs **no
reinitialization phase at all**: the compression term lives inside the same
fused transport equation as $\phi$, so there's nothing to periodically
re-solve. Saini's own CLS scheme, by contrast, runs *two* separate periodic
pseudo-time phases per Algorithm 1 — TLS redistancing for their distance
field and CLS re-sharpening for their phase field — each iterating on the
order of 100 pseudo-steps to convergence, each pseudo-step needing its own
gradient evaluation (halo exchange across ranks) and often a global
reduction for stability. That's a burst of many small, synchronization-heavy
operations layered onto the regular timestep — a poor shape for GPUs and
large-scale MPI, which want large, infrequent, well-overlapped work, not many
tiny serialized steps. This is a genuine, recognized weakness of classical
level-set reinitialization in HPC contexts, not an incidental detail of one
paper.

Once this repo adds $\psi$-transport and redistancing to fix the normal, it
does give some of that advantage back — but only partially, not completely.
The $\phi$ equation itself stays fully fused; Saini's Eq. 39 (CLS
re-sharpening) was never reintroduced. Only *one* of their two periodic
phases (redistancing) came back, not both. How much of a regression that
actually is in practice depends entirely on how often redistancing needs to
fire — which `REDISTANCING.md` §6 and §8 take up — and is why `zalesak_disk`'s ablation
(see its `README.md`) is framed as "does adding redistancing to an already-
working configuration help, hurt, or do nothing," not as a foregone
conclusion that we need it.


### Which of these each case actually uses

Stated plainly, because the distinction only matters if it is applied:

| Case | $\psi$ initialization | $\psi$ maintenance |
|---|---|---|
| `advecting_slab_1d` (all but `psi_rd_*`) | exact analytic slab distance, once | **none** — advection plus SVV on the $\psi$ equation, 20 flow-throughs |
| `advecting_slab_1d` `psi_rd_*` | built by Eq. (44) | off, or 40 events (`seed="phi"`/`"psi"`); the redistancing testbed (§4.2) |
| `zalesak_disk` primary | exact *periodic* (nine-image) distance, once | **none** — advection, plus SVV on the $\psi$ equation |
| `zalesak_disk` `_redistance_phi` | built by Eq. (44) | redistancing on, timer $\Delta t_{tls}=0.5$, `seed="phi"` — **reinitializes** from $\phi$ ($\psi_0 = r_f(\phi-0.5)$), then relaxes |
| `zalesak_disk` `_redistance_psi` | built by Eq. (44) | redistancing on, timer $\Delta t_{tls}=0.5$, `seed="psi"` — relaxes **in place**, no reinitialization |
| `rider_kothe` | exact periodic distance, once | **none** — advection plus SVV on the $\psi$ equation |

So the three validated results — the 1D slab, the Zalesak reference
configuration and Rider–Kothe — use neither redistancing nor reinitialization.
Only the ablation variants do, and they exist to test whether adding it helps,
not because it was needed.

### How ours differs from Saini et al.'s, and why the ablation is confounded

Our implementation follows their Eqs. (44)–(47) — the pseudo-time TLS equation,
the $\psi_0 = r_f(\phi-0.5)$ seed with $r_f=0.1$, the smoothed
$\text{sgn} = \tanh(\psi/2\varepsilon)$, a $2.5H$ band, SVV on the relaxation.
The first difference matters:

| | Saini et al. | this repo |
|---|---|---|
| $\psi$ at $t=0$ | **built** by solving Eq. (44) from the $\phi$ seed (their Algorithm 1, line 2) | `psi_init`: `"exact"`, the analytic periodic distance (the validated runs), or `"redistance"`, Algorithm 1 line 2 |
| re-init seed | Eq. (47), $r_f = 0.1$ | same |
| band | $2.5H$ | $2.5H$ |
| $\Delta\tau_{tls}$ | $H/(N{+}1)$ | $0.1\,h_{\text{GLL,min}}$ |
| $N_{tls}$ | $2.5(N{+}1)$ — 10/15/20 at $N=3/5/7$ | ~213 |
| trigger | timer, $\Delta t_{tls}$ | timer, $\Delta t_{tls}$ (`redistance.dt_tls`) |

They never use an analytic distance: Algorithm 1 line 2 *constructs* $\psi$ with
the same solve that later maintains it, so their field is **band-consistent from
$t=0$**. An analytic $\psi$ is a global exact distance, and each redistancing
event would replace it discontinuously with a field that is distance-like only
within $2.5H$.

**So `zalesak_disk`'s recorded redistancing ablation is confounded** — its
strongly negative result condemned *the pairing it ran*, an analytic-distance
initial condition with band-limited redistancing, not redistancing as such. The
shipped variants now build $\psi$ by Eq. (44) (`psi_init = "redistance"`); the
coherent configuration is examined in §4.1.

### 4.1 Algorithm 1 done coherently — and the bug it exposed

$H=1/50$, $N\in\{3,5,7\}$, $\xi=1$, $\gamma=1$, $\text{svv}_\psi$ on, $t=20$.
Three arms differing **only** in $\psi$ ($\phi$'s initial condition is identical;
$E_r(0)$ agrees to $10^{-10}$):

| | $\psi$ at $t=0$ | periodic redistancing |
|---|---|---|
| **A** | exact analytic periodic distance | off |
| **B** | built by Eq. (44) from the $\phi$ seed | off |
| **C** | built by Eq. (44) | on, $\Delta t_{tls}=0.5$ |

**The built-$\psi$ runs diverged at $N=5$ and $7$ until a defect in this repo's
code was fixed.** The defect is invisible while $\psi$ is analytic.

**The bug.** `unit_normal` floors $|\nabla\psi|$ before dividing it out:

```fortran
call field_sqrt(w, n)                  ! w = |grad psi|
call field_cpwmax2(w, grad_floor, n)   ! w = max(w, floor)
call field_invcol2(g1, w, n)           ! n = grad psi / w
```

The floor was $10^{-30}$. A numerically **flat** field's gradient is not zero but
round-off, $\sim10^{-19}$ — so the floor never engaged, and wherever $\psi$ was
flat the code divided round-off by round-off and produced a **unit vector in a
random direction**. Its own comment stated the intent correctly ("guards 0/0
where the field is flat") and the value simply did not implement it. For a global
analytic $\psi$, which is never flat, this is unreachable dead code. Algorithm 1
line 2 bands $\psi$ by design, leaving 67.6% of the domain flat, and the dead
code became live.

**Why it is fatal and not negligible.** Expanding the compression term,

$$\gamma\nabla\cdot\!\big(-\phi(1-\phi)\mathbf{n}\big)
= -\gamma(1-2\phi)\,\mathbf{n}\cdot\nabla\phi \;-\; \gamma\,\phi(1-\phi)\,\nabla\cdot\mathbf{n}$$

Outside the band $\phi(1-\phi)\sim10^{-16}$, so the *forcing* is negligible — the
reason this was initially dismissed. But for a random unit field
$\nabla\cdot\mathbf{n}\sim 1/h\sim10^{3}$, so the second term is a **source
proportional to $\phi$ with rate $\sim\gamma u_{\max}\nabla\cdot\mathbf{n}
\approx\pm2\times10^{3}$**. That is an exponential instability, not a forcing: a
$10^{-16}$ seed reaches $\mathcal{O}(1)$, and the onset time depends on $N$
through how $\phi$'s diffusion $\varepsilon\gamma u_{\max}$ damps grid scales.

**Why Saini never hit it.** Their CLS transport, Eq. (38), is pure advection plus
SVV — **no compression term and no normal at all**. $\mathbf{n}$ appears only in
Eq. (39), the periodic re-initialization, and their own text says the band
gradient must be smooth "to compute smooth normal vector, *for the CLS
re-initialization equation, Eq. (39)*". So their normal is never evaluated
outside a band-localised relaxation, and a band-limited $\psi$ is exactly
sufficient. Neko's CDI fuses compression into transport and therefore reads
$\mathbf{n}$ pointwise over the whole domain every step — it needs a guard their
scheme does not. The guard existed and was set to a value that never fired.

**The fix, and its verification.** `case.cdi.grad_floor`, default $10^{-6}$:
where there is no gradient there is no direction, so $\mathbf{n}\to0$ and the
compression term vanishes — which is the answer the mathematics has there.
Measured at $N=5$, $t=1.5$ (arm B dies at $t=1.17$ without it):

| floor | redistancing band | outcome | $\phi$ range |
|---|---|---|---|
| $10^{-30}$ | $2.5H$ | **diverges** $t{=}1.166$ | $[-3.25,\ 1.85]$ |
| $\mathbf{10^{-6}}$ | $2.5H$ | survives | $[-3.4\times10^{-9},\ 1.000]$ |
| $10^{-30}$ | $40H$ (global) | survives | $[-1.9\times10^{-10},\ 1.000]$ |
| $10^{-6}$ | $40H$ | survives | identical to the row above |

The last two rows being **bit-identical** is the confirmation: building $\psi$
globally removes the plateau, so the floor has nothing left to do. Two
independent fixes, one mechanism.

At $N=7$, the worst order, the floor holds: bounded to $1.7\times10^{-11}$ at
$t=1.5$ with Saini's $2.5H$ band. The global band also survives but is worse on
both counts that matter ($\phi_{\max}=0.9936$ rather than $1.000$; mass drift
$7\times10^{-7}$ rather than $1.5\times10^{-9}$) and costs 240 pseudo-steps
instead of 15 — so **keep Saini's $2.5H$ band and fix the floor**.

**The floor is a safe default, not an opt-in.** On the analytic-$\psi$
configuration that produced every validated result in this document, raising it
from $10^{-30}$ to $10^{-6}$ leaves the run **bit-identical** — `bnd` at
$N=5$, $t=1.4944$ reads `-0.2594E-09  0.1000E+01  0.1332E-07` either way. It
cannot perturb what already worked, because that field is never flat.

A banded $\psi$ and a fused pointwise-normal scheme are therefore compatible: a
fused scheme reads $\mathbf{n}$ where a split scheme never does, so it must handle
$\nabla\psi=0$ explicitly, and the floor does that.

### 4.1b The redistancing band must cover the interface band: $N \ge 3.68\,\xi$

The compression term reads $\mathbf{n}$ wherever $\phi(1-\phi) > 10^{-4}$, an
interface half-width of $9.21\varepsilon = 9.21\,\xi H/N$. Redistancing only
rebuilds $\psi$ within `redistance.band` $\times H$. If the first exceeds the
second, part of the interface sits in the clamped plateau, where `grad_floor`
correctly sets $\mathbf{n} \to 0$ — but $\phi(1-\phi)$ there is **not**
negligible, so **compression is switched off over part of the interface** and
$\phi$ smears. With Saini's band of $2.5H$:

$$9.21\,\xi H/N \;\le\; 2.5H \quad\Longrightarrow\quad N \;\ge\; 3.68\,\xi$$

| $\xi$ | 0.5 | 1.0 | 1.5 | 2.0 | 2.8 |
|---|---|---|---|---|---|
| smallest usable $N$ at band $2.5H$ | 2 | **4** | 6 | 8 | 11 |

**It is a real coverage statement, and it is *not* what kills arm C.** With the
event's history restart (§4.1d), arm C at band $2.5H$ diverges at $N=3$ ($t=11.2$,
criterion violated), completes at $N=5$ (satisfied) but 17× worse than without
events, and diverges at $N=7$ (satisfied; §4.1c).

Raising the band at $N=3$ does improve the *build* — the band mean moves
1.259 → 1.027, converging by $\approx10H$ — so the coverage effect is real. A
band-$4.0H$ arm C run at $N=3$ also diverged, but it had the history bug and has
not been repeated. So treat $N \ge 3.68\,\xi$ as a **quality** consideration for
the built field, not a stability requirement, and see §4.1c for what actually
fails.

Saini's §4.5 runs $N \in \{4,5,6\}$ at $\xi = 1/N$, all of which satisfy it, so
the question does not arise for them either way.

### 4.1c What actually fails on Zalesak: the reseed copies $\phi$'s fragments

Arm C re-run on 2026-10-02 with the event's order restart (§4.1d). The cases are the
originals with only paths changed; predictions and verdicts are in
`examples/zalesak_disk/logs/armC_2026-10-02/PREREGISTERED.txt` (gitignored, local). Arm B runs no events and so
never had the bug.

| | $N=3$ | $N=5$ | $N=7$ |
|---|---|---|---|
| **B** — built $\psi$, no events | $t=20$ ✓, $E_r$ 0.090 | $t=20$ ✓, $E_r$ 0.0336 | stopped at $t=9.4$, clean |
| **C**, order restart | **diverged $t=11.2$** | completes, **$E_r$ 0.574**, violation $7.9\times10^{-2}$ | **diverged $t=4.5$** |
| $\lVert dn\rVert$ per event, C restarted | 32 → 77 | 57 → ~95, flat from event ~20 | 166 → 575 |

**The restart was not the whole story here**, unlike in 1D, where it made the reseed free (§4).
Measured on frames written right after an event, at $N=5$ in the compression band
(`examples/zalesak_disk/logs/armC_2026-10-02/rd_quality_C_B_n5.txt`, gitignored, local):
- **The relaxation converges.** $|\nabla\psi|$ has mean 1.01–1.17, and $\psi$ is within
  0.25–0.5$\varepsilon$ of the signed distance to $\phi$'s 0.5 contour, its target.
- **The target fragments.** $\phi$'s 0.5 contour is in 2 pieces at $t=1$ and 8 at $t=20$; arm B
  keeps 1.
- **Each reseed turns every spurious piece into a zero set of $\psi$**, which the compression
  then maintains. $\psi$'s normal goes from 13° to 73° off the exact interface; arm B stays near 3°.

That is the $\phi$ → reseed → $\psi$ → $\mathbf n$ → $\phi$ loop of §4, measured. A 1D slab cannot
fragment, which is why the same operation is free there. Where the first spurious piece comes
from (the slot's corners are the obvious suspect) is not measured; `NEXT_SESSION.md`.

**Walls are a separate risk, and none of these cases has one.** In the 1D replica of
`redistance_circles`, the reseed front from a flat $r_f(\phi-0.5)$ seed captures a
natural-BC wall at $N\ge7$ when it arrives with the wall node below its neighbour. The
periodic domain holds (`examples/redistance_circles/archive/README_process_2026-09.md` §10, "1D, reinit seed").
All three coupled cases are periodic. Their pseudo-time also stops at
$\text{band}\times H = 2.5H$, so a front can reach a wall only if an interface lies
within $2.5H$ of it. A future walled case would need that check.

### 4.1d An event must restart the time history, a varying diffusion must reach the solver, and a reseed is only as good as $\phi$ (2026-10-01/05)

On Rider–Kothe, the settings of Saini's circVortex case diverged in every events arm, including
his complete configuration. Two bugs were found and fixed. Tables are in
[`examples/rider_kothe/README.md`](examples/rider_kothe/README.md); the record is in
`examples/rider_kothe/logs/d5/PREREGISTERED.txt` (gitignored, local).

**1. A redistancing event that replaces $\psi$ must restart the time history.**
- **The bug.** The coupled files run the event in the `compute()` hook, after Neko's
  `slag%update()`. So the next BDF3 steps combined the rebuilt $\psi$ with the pre-event $\psi$
  still in the lags, and settled at $\psi_{old}+\tfrac{11}{6}(\psi_{new}-\psi_{old})$. That is
  an O($\Delta\psi$) error, not O($\Delta t$).
- **The fix.** All three coupled files now set `fluid%ext_bdf%nadv = ndiff = 0` after an event,
  as Saini's fork does (`ireset_ls`), so the next step is BDF1/EXT1.
- **Checks.** With redistancing off the output is bit-identical. On a Rider–Kothe events run it is
  bit-identical to the scratch key that was validated against copying the lags.
- **Effect.** In 1D it removed the whole "cost of the reseed" (§4). On Zalesak it did not
  (§4.1c).

**2. A time-varying diffusion must reach the solver.**
- **The bug.** `rider_kothe.f90` sets the diffusion $\varepsilon\gamma u_{\max}(t)$ in
  `material_properties`, but Neko's scalar solve reads `s_lambda_tot`. Neko copies that from
  `s_lambda` only at initialisation unless a turbulence model is set (Neko v1.1.0,
  `src/scalar/scalar_scheme.f90`). So the diffusion stayed at $\varepsilon\gamma$ while the
  compression followed $|\cos(\pi t/8)|$.
- **What it did.** The equilibrium width grew as $\varepsilon/|\cos(\pi t/8)|$. That predicts the
  dissolved filament's core to 0.02 at $t=2.88$ and 3.36. At $t=3.84$, where the compression is
  nearly off, $\phi$ lags behind that equilibrium (core 0.31–0.41 against 0.20–0.27).
- **The fix.** `material_properties` fills both fields.
- **Effect.** Transport-only $E_r(8)$ went from 0.226 to 0.0446 (0.0452 after item 4's velocity fix); Saini's, in our normalisation,
  is 0.041.
- **Scope.** `zalesak_disk` and `advecting_slab_1d` have a constant $u_{\max}$ and are unaffected.

**3. With both fixed, Saini's full configuration completes on Rider–Kothe** ($E_r$ 0.0407, worst
violation $8.9\times10^{-4}$, $\lVert dn\rVert$ 86–159 per event, rising 1.9× over 16 events). But its rebuilt $\psi$ is a
worse normal source than the transported one:
- **The solve does not reach a distance in the compression band.** Its sign-function width is
  0.25, so the pseudo-velocity there is below about 0.1 and $\psi$ stays close to the steep seed.
  Right after each event $|\nabla\psi|$ has mean 1.66, and $\psi$ is $2.6\varepsilon$ off its
  target on average.
- **Its normal is 4–8° off the exact interface; the transported $\psi$'s is 0.6–0.8°.**
- **Where the filament is thinner than $2\varepsilon$** (5–6% of its length at maximum stretch),
  $\phi$ has no 0.5 contour, and the rebuilt $\psi$ is wrong across it at 76–79% of points.
- **$\phi$'s normal error rises to 28° at $t=8$**; transport alone ends at 2°.

**A reseed takes $\phi$'s interface as the truth.** It is only safe while $\phi$'s interface is
better than the transported $\psi$'s. On a straining flow with an under-resolved tail, and
around corners, it is not.

**4. A prescribed velocity belongs to the explicit terms' time level (2026-10-05).** Neko's scalar
step applies advection and the compression source to $s^n$ and extrapolates them to $t_{n+1}$
(`src/scalar/scalar_pnpn.f90`), Saini's Eq. (34) with $\mathbf C\,u^{n+1-j}$. So the velocity
must be $\mathbf u(t_n)$. `rider_kothe.f90` and `zalesak_disk.f90` set it in `compute`, after the
step, which gave $\mathbf u(t_n)$ from step 2 on but $u=0$ on step 1. Both now set it in
`preprocess` at `time%tlag(1)` $=t_n$, as Saini's circVortex `userchk` does (it runs once before
step 1). Rider–Kothe's implicit diffusion now reads $u_{\max}(t_{n+1})$. Every Rider–Kothe
$E_r(8)$ moved up by 0.5–3.0% (0.0446 → 0.0452 at $H=1/64$, $\xi=1$); boundedness and filament
area did not change; Zalesak moved 0.1%. Records: `examples/rider_kothe/logs/vel_fix_2026-10-05/`
(gitignored, local).

### 4.2 Pseudo-timestep: fine steps, and why a one-shot test misleads

Saini's $\Delta\tau_{tls} = H/(N{+}1)$ against this repo's
$\Delta\tau = 0.1\,h_{\text{GLL,min}}$. Both cover the same total pseudo-time
$2.5H$, so both target the same fixed point; they differ only in how many steps
they take to get there ($2.5(N{+}1)$ against $\lceil 2.5H/\Delta\tau\rceil$).

**On a single build, Saini's looks better** — measured on Zalesak, $\xi=1$:

| $N$ | ours: $N_{tls}$, band mean, max, wall | Saini: $N_{tls}$, band mean, max, wall |
|---|---|---|
| 3 | 91, 1.259, 2.664, 4 s | **10**, 1.243, 2.380, 3 s |
| 5 | 213, 0.995, 2.617, 13 s | **15**, 0.988, **1.698**, 6 s |
| 7 | 390, 1.029, 4.403, 36 s | **20**, **0.987**, **2.254**, 11 s |

**It is wrong even for the build, and unstable under repeated events.** Measured
on `advecting_slab_1d` (four minutes a run), holding the build $\Delta\tau$ and
the event $\Delta\tau$ as separate controlled variables:

| $\psi$ at $t=0$ | build $\Delta\tau$ | periodic redist. | $E_r(t{=}20)$ | worst violation |
|---|---|---|---|---|
| exact analytic | — | off | **0.00006** | $1.2\times10^{-10}$ |
| built by Eq. (44) | `cfl = 0.1` | off | **0.00006** | $1.2\times10^{-10}$ |
| built by Eq. (44) | `cfl = 0.1` | on, 40 events | **0.00006** | $1.4\times10^{-10}$ |
| built by Eq. (44) | $H/(N{+}1)$ | off | 0.06458 | $2.6\times10^{-3}$ |
| built by Eq. (44) | $H/(N{+}1)$ | on, 40 events | **1.283** | destroyed |

(Every frame, with the event's history restart, re-measured 2026-10-05. Before the restart, the
40-event `cfl = 0.1` row read 0.00023; that was the history bug, §4.1d.)

The coarse $\Delta\tau$ does not merely cost accuracy — it produces a *wrong
build*: $\lvert\nabla\psi\rvert$ spans 0.067–3.165 where the fine build gives
0.999–1.001. Each subsequent event then overshoots the fixed point and the
overshoot compounds, beyond $\lvert\nabla\psi\rvert\sim10^{30}$ by event 40.
$H/(N{+}1)$ is a pseudo-CFL of 2.75 on the slab's mesh ($H=1/10$, $N=10$); on
Zalesak it is 0.905/1.419/1.949 at $N=3/5/7$.

**With a fine build, giving up the analytic distance costs nothing at all** in
1D: $E_r$ and the worst violation are identical to the analytic reference. That
is the result the whole `psi_init` exercise exists to establish.

**Use `redistance.cfl = 0.1`** (the default) and do not set `redistance.dtau`.

**Saini's $\Delta\tau$ is stable in his own configuration for a reason that costs him
elsewhere.** His sign-function width of 0.25 keeps the pseudo-velocity
$|\operatorname{sgn}\psi|$ below about 0.1 near the interface, so $H/(N{+}1)$ is a small
pseudo-CFL there. On Rider–Kothe his full configuration ran 16 events cleanly. The same width is
why his solve does not reach a distance in the compression band (§4.3).

### 4.3 There is no convergence criterion on the pseudo-time solve

`redistance()` runs a fixed `rd_niter` steps with no residual check and no
stopping test. Since Eq. (44)'s residual is $\operatorname{sgn}(\psi)(1-|\nabla\psi|)$,
the band mean already logged *is* a convergence measure, and read that way the
build converges to 1.2–1.3% at $N=5,7$ but only 24% at $N=3$ — the relaxation
does not reach $|\nabla\psi|=1$ at low order, in either pseudo-timestep setting.
That is not what caused the divergences (those were at $N=5,7$, where the
residual is small), but nothing in the code would have told us either way.

**A fixed step count can leave the solve unconverged where it matters.** With Saini's Rider–Kothe
settings (sign-function width 0.25, $25H$, $\Delta\tau=H/(N{+}1)$, 150 steps), right after each
event the compression band still has $|\nabla\psi|$ with mean 1.66 and 5th–95th percentiles of
1.1–2.2. $\psi$ is $2.6\varepsilon$ off its target there on average (§4.1d). The cause: the
pseudo-velocity $|\operatorname{sgn}\psi|\approx\psi/2\varepsilon_{sgn}$ is below about 0.1 in
the band, so $\psi$ there barely moves from the seed. With the sign-function width equal to the
phase field's (Zalesak arm C) the same check converges to $0.25$–$0.5\varepsilon$ (§4.1c). A
residual on the compression band, not the build band, is the measure to use.

A caution on the band **minimum**, after `REDISTANCING.md` §4(d): it is
$7.5\times10^{-14}$ at $t=0$ for the *exact analytic* $\psi$ at $N=3$ — lower
than any built field — in the configuration that runs ten rotations cleanly.
Gather-scatter averaging cancels the gradient across kinks. **The band minimum
diagnoses nothing.**

### 4.4 Saini §4.4: reproduced — and *where the SVV acts* is what matters

**The record of this benchmark is
[`examples/redistance_circles/README.md`](examples/redistance_circles/README.md)**:
the configuration and where each setting comes from (§2), the paper against the
code (§3), the results (§4) and the limitations (§5). This section keeps only the
lessons that belong to the method.

**What the case is.** Saini's own rigour test for Eq. (44), the one equation this
repo adopted from them, run **standalone**. The skewed $\psi_0$ of Eq. (83) relaxes
toward the distance around two intersecting circles on $[-2,2]^2$, with
$\Delta\tau=5\times10^{-4}$ to $\tau=6$, $N_{svv}=N/6$, $c_0=2$ and their BDF2/EXT2.
It uses the authors' configuration, taken from their public case. Three of its
settings are not in the paper:
- the **sign-function** $\varepsilon=0.25$, fixed, in Eq. (46) (their `signls`). The
  phase-field width $\varepsilon=\xi H/N$ of Eq. (36) is a different quantity, and
  this case has no phase field;
- $\mathbf C(\mathbf w)\psi$ **dealiased** on $\lfloor3(N+1)/2\rfloor$ Gauss points;
- their **sign guard**: a node whose $\psi^n$ disagrees in sign with $\psi_0$ keeps
  $\psi^n$.

**Status (2026-09-29): reproduced.** 12 of the 18 Fig. 12 cells are within 0.1% of
their values, four more within 2.2%. Two cells ($H{=}1/5$, $N{=}4$ and 8) fail on our mesh through
zero-set round-off; started from their zero-set $\psi_0$, our code gives their
values exactly (README §4–§5).

**Lessons for the method:**

- **Where the SVV acts decides everything.** Every signed distance function is a
  steady state of Eq. (44), so the interface position is held by nothing but
  $\operatorname{sgn}(0)=0$.
  - Saini's Eqs. (2), (7), (31) and (33) put $\mathbf D_\mu=\operatorname{diag}(c_0|\mathbf c|\mathbb H/N)$
    **left** of the assembled $\mathbf S_{vv}$. For Eq. (44), $|\mathbf c|=|\operatorname{sgn}\psi|$
    vanishes on the zero set, and the scheme becomes sign-invariant: no node can
    change sign.
  - $\nu$ inside the bilinear form with $|\mathbf c|=1$, the coupled cases'
    `svv_step_imp`, acts on the zero set. Here it was 14–920× off Table 2, with spurious
    interfaces (at the earlier sign-function $\varepsilon=H/N$, 2026-09-23). Moving the coupled cases to the printed form is a separate decision
    (`NEXT_SESSION.md`, D2).
- **$\mathbf w$ is a function of $\psi$, not a velocity field.** A $\mathbf w$
  decoupled from the $\psi$ it multiplies leaves the source's
  $\operatorname{sgn}'(\psi)\delta$ uncancelled on the zero set: growth at up to
  $1/2\varepsilon$.
  - Frozen $\mathbf w$ diverges.
  - Every consistent time treatment gives the same result: SSP-RK3 with a Lie
    split, and BDF2/EXT2 and BDF3/EXT3 with each level's own product extrapolated,
    agree within about 2% in Neko at §4.4's settings.
- **Dealiasing $\mathbf C(\mathbf w)$ lets the zero set move.** The non-dealiased
  GLL Galerkin form equals $\operatorname{sgn}(\psi)|\nabla\psi|$ to 1e-14 and leaves
  zero-set nodes alone. Dealiasing evaluates $\mathbf w$ between nodes, so they move
  and can change sign. The sign guard then freezes them by the sign of $\psi_0$,
  which is $\pm10^{-15}$ there, so a low-$N$ result depends on round-off (README §5).
  The coupled cases keep the non-dealiased form until D5.
- **Eq. (44) alone drifts.** Every signed distance function is a steady state, so
  nothing restores an interface's position, the weakness the paper's introduction
  names. The SVV beside the zero set moves it slowly, and dealiasing makes even a
  straight interface drift (README §5). That matters for any coupled use of
  redistancing (§4.1c).
- What was learnt at the earlier sign-function $\varepsilon=H/N$ (the drift ending in
  node capture, $\varepsilon$ as its clock, a step size holding the apex at
  $N_{svv}=N/4$, a Russo–Smereka anchor) is in
  `examples/redistance_circles/archive/README_process_2026-09.md` §7.5, §8, §10.

The open items (D2, D5) are in `NEXT_SESSION.md`.

## 5. SVV is a $\psi$-only knob

**There is no SVV on $\phi$, and it is never a variable in any study.** The split is
fixed: **$\phi$ carries the CDI equation** — transport, compression, and the
physical diffusion $\varepsilon\gamma u_{\max}$ — and **$\psi$ carries transport
plus SVV**. `svv_phi` no longer exists in the code; the coupled files stop at startup if a
case sets `case.cdi.svv_phi`. **SVV on $\psi$ is always on** (user rule, 2026-10-02): the three coupled files stop
at startup if `normal = "psi"` has `svv_psi` off. Saini does the same: his TLS transport,
Eq. (43), carries $S_{vv}$, as do all his scalar equations (JCP p.3, p.9). The only SVV knob that
ever varies is `case.cdi.svv_psi.c0`.

Note that Saini et al.'s naming is transposed (§1), so their §4.3 SVV sits on the
*phase* field, i.e. on our $\phi$. **We do not adopt or test that
placement.** Their scheme is structurally different anyway (§2): one equation with
SVV holding the interface together, against our two with a compression term doing
that job.

Both SVV and $\phi$'s diffusion are diffusive operators, but they are not
interchangeable and this repo does not merge them:

- $\phi$'s $\varepsilon\gamma\nabla\phi$ term (§2) is **physical** — it's the
  other half of the compression/diffusion balance that defines the
  equilibrium interface profile. Removing or replacing it changes what shape
  the interface relaxes to, not just how noisy the solution is.
- SVV is a **numerical, grid-scale-targeted** high-pass spectral filter
  ($S_{vv} = \widetilde{D}^{\mathsf T} G \widetilde{D}$, per-equation
  $c_0$/`nsvv_ratio` knobs — see `references/saini_2026_test_cases.md` §1.4
  for the operator in full).

Measured evidence for keeping them separate rather than merging: SVV on the
$\phi$ equation delays but does not cure CDI's grid-scale growth problem — it
wins the first Zalesak rotation and has lost the advantage by the fifth, and
never helps the 1D slab at all, because the growth lives in Legendre modes
0–2, which SVV cannot touch by construction (mode 0 passes untouched, mode
$N$ passes at full strength). Once the normal comes from $\psi$ instead of
$\phi$, no stabilization is needed on the $\phi$ equation at all
(`../neko-multiphase/CDI_IN_SEM.md` §0, §6.1). The $\psi$ equation is the
opposite case: it has *no* physical diffusion of its own (pure advection), so
if transported $\psi$ develops grid-scale ringing (e.g. near medial-axis
kinks), SVV is the textbook tool for exactly that equation. This is the one
place in the whole method SVV earns its place, and it is always on there.

### Two instances, two settings — both Saini's own

The SVV parameters are **not** the same on every equation, and the differences
are taken from the paper rather than tuned here. $\phi$ has no instance (§5 above):

| instance | $N_{svv}$ | $c_0$ | source | applied |
|---|---|---|---|---|
| `svv_psi` ($\psi$ transport, Eq. 43) | $N/2$ | **0.1**; Zalesak still 1.0 | 0.1 is his TLS-transport value: §4.5, *"The stabilization parameters for the advection equation, of both CLS and TLS fields, ... are $N_{svv}=N/2$; $c_0=0.1$"* (p.21), and his §5 default (p.24). Zalesak's 1.0 is his §4.3 value, which sits on the **CLS** transport (our $\phi$), so it was borrowed from the wrong equation; whether to move Zalesak is open | explicit, via the source term |
| `svv_rd` ($\psi$ redistancing, Eq. 44) | $N/4$ | $2.0$ | their §4.5 TLS setting, because the Eq. (44) fixed point has kinks at medial axes; only their standalone §4.4 test goes stronger, to $N/6$ | implicit, inside each pseudo-step |

So the redistancing equation gets a **broader and twice-stronger** filter than the
transport equation. That is deliberate on their part and reproduced here.
`svv_rd` is forced implicit in `startup` regardless of the case file: at $c_0=2$
the operator's spectral radius would otherwise set the pseudo-timestep instead of
the CFL.

The viscosity is $\nu = c_0|\mathbf{c}|\mathbb{H}/N$ (their Eq. 29) with
$\mathbb{H}=2J_e^{1/d}$. This repo substitutes the element edge, which is
**exact** when $d$ is the true spatial dimension — $2J_e^{1/2}=H$ for the
Zalesak element, $0.02$ either way. Reading $d=3$ on these one-element-thick
meshes would fold in a meaningless $z$ extent and inflate $\mathbb{H}$ by 1.71×.
The one genuine deviation is $|\mathbf{c}|$ — and it is two deviations, not
one.

- **Pointwise vs uniform.** Saini evaluate $|\mathbf c|$ pointwise, from the
  *local* characteristic speed (Eq. 31). The coupled files use a uniform bound,
  the flow's $u_{\max}$, for transport **and** for redistancing (whose
  characteristic speed is $|\mathbf w|\le1$, so on Zalesak their redistancing SVV
  is $\pi/\sqrt2\approx2.2\times$ Saini's). `redistance_circles` uses 1, with
  the pointwise $|\operatorname{sgn}\psi|$ applied as $\mathbf D_\mu$.
- **Where $\mu$ sits.** They left-multiply the assembled operator,
  $\mathbf D_\mu\mathbf S_{vv}$ (Eqs. 7, 33, 35), whereas `svv_local` applies
  `nu` *inside* the bilinear form. The two agree only while $\mu$ is constant.

For transport the difference is a mild over-diffusion. **For the redistancing
equation it is not mild.** There $|\mathbf c|=|\operatorname{sgn}\psi|$ vanishes
on the zero set, and in the left-multiplied form that is what stops the SVV
dragging the interface. At the earlier sign-function $\varepsilon=H/N$ (2026-09-23) it
turned Saini §4.4 from 14–920× off into 2.1–2.7× for three of four cells; with the
authors' full configuration the case reproduces theirs (§4.4). `redistance_circles` therefore uses the printed form
(`svv_step_eq31`, its only SVV step). The coupled cases still use
`svv_step_imp`: there the same instance also carries $\psi$ *transport*, where
$|\mathbf c|$ is the flow speed, so moving them over is a separate decision
(`NEXT_SESSION.md`, D2).

### What SVV on $\psi$ actually buys — measured

Saini et al. argue for SVV largely from theory. Their only on/off comparison is §4.1's Fig. 2, for a
generic advected scalar. They never show the level-set fields without it. This is that result,
on their own §4.3 Zalesak setup, with SVV on $\psi$ at their §4.3 (CLS) strength:
$H=1/50$, $N \in \{3,5,7\}$, $\xi = 1$ (their $\xi = 1/N$), $N_{svv}=N/2$,
$c_0=1.0$, $\gamma=1$, ten rotations. Each on/off pair uses an **identical**
$\Delta t$, so `svv_psi.c0` is the only thing that differs.

| $N$ | $E_r(t{=}20)$ on | off | ratio | worst violation on | off |
|---|---|---|---|---|---|
| 3 | 0.1012 | 0.9316 | 9.2× | $2.6\times10^{-8}$ | $1.6\times10^{-3}$ |
| 5 | 0.0208 | 0.9747 | 46.8× | $6.4\times10^{-6}$ | $1.4\times10^{-2}$ |
| 7 | **0.0040** | **1.0266** | **256×** | $9.5\times10^{-6}$ | $3.1\times10^{-2}$ |

**SVV is what makes $p$-refinement work at all.** With it, $E_r$ falls **25×**
from $N=3$ to $N=7$ — reproducing the convergence Saini's Fig. 7 reports.
Without it, $E_r$ *rises*: 0.93 → 0.97 → 1.03. Raising the polynomial order makes
the answer **worse**, so the entire reason to use a high-order method is
forfeited. At $N=7$, $E_r > 1$ means the error exceeds the exact solution's own
integral — the field has effectively dissolved.

The mechanism is the one the $\phi$-normal analysis already established: higher
$N$ gives grid-scale noise more modes to grow in. On the $\psi$ equation SVV
damps exactly those modes, so the benefit **widens with $N$** (9× → 47× → 256×)
rather than being a fixed factor. That is why a single on/off pair at one
resolution understates the case.

Three qualifications, because "SVV fixes everything" would be wrong:

- **It is not needed for mass conservation.** $|E_v|$ is $10^{-9}$–$10^{-7}$
  *either way* — the CDI form is conservative regardless. SVV buys one or two
  orders there, on a quantity already negligible.
- **Its boundedness benefit is $\xi$-dependent.** At $\xi=2.8$ boundedness is
  exactly zero with or without SVV (see `examples/zalesak_disk/README.md`): the
  compression flux's $\phi(1-\phi)$ factor carries it alone. At $\xi=1$ — a
  sharper interface, and Saini's own setting — SVV is worth 3–4 orders of
  magnitude. So the older claim that "SVV buys accuracy, not boundedness" holds
  only at thick interfaces.
- **Without SVV nothing crashes.** All three off-runs completed ten rotations and
  produced plausible-looking output that is up to 256× wrong. That is a worse
  failure mode than divergence, and the reason this comparison had to be run
  rather than assumed.

### The mechanism, measured on $\psi$ itself

The above is $\phi$'s error. The cause is visible directly in $\psi$. Measured in
the interface band $\phi(1-\phi) > 10^{-4}$ — the only place the compression term
reads the normal — at $t=0$ and after ten rotations:

| $N$ | mean $\lvert\nabla\psi\rvert$, SVV on | SVV off | band min, SVV on |
|---|---|---|---|
| 3 | 0.999 → 0.797 | 0.999 → **26.0** | 0.104 → 0.046 |
| 5 | 1.000 → 0.940 | 1.000 → **51.6** | 0.143 → 0.016 |
| 7 | 0.998 → 0.956 | 0.998 → **69.5** | 0.171 → 0.0047 |

**Without SVV, $\lvert\nabla\psi\rvert$ does not drift — it explodes**, to a band
mean of 26–70 and a maximum of 3515 at $N=7$, and it gets worse with $N$. The
normal $\mathbf{n} = \nabla\psi/\lvert\nabla\psi\rvert$ is then read off a field
whose gradient is grid-scale noise, which is exactly why $E_r$ saturates near 1.
This is the mechanism, not an inference from $\phi$'s error.

**A caveat on the isometry argument (§3).** That argument — rigid rotation
preserves $\lvert\nabla\psi\rvert$ because $\nabla\mathbf{u}$ is antisymmetric —
covers **pure advection only**. $\psi$'s equation carries $S_{vv}$ on the
right-hand side, so it does not apply to the SVV-on runs, and indeed the band
mean drifts down by 4% ($N=7$) to 20% ($N=3$). Two reasons that drift is benign:
the compression term uses the **unit** normal, so a magnitude change alters no
direction; and the drift *shrinks* with $N$ while the band **minimum** degrades
*most* at $N=7$ — the order with the best $E_r$. So $\lvert\nabla\psi\rvert$
conditioning is not what drives the remaining error.

Figures: `examples/zalesak_disk/evidence/saini_fig{5,6,7}_*.png`.

### At $c_0=0.1$: the slab and Rider–Kothe (2026-10-02)

- **`advecting_slab_1d`.** $E_r$ is unchanged for $\xi\le1.5$, every $\gamma$ and $N=6$–12.
  - **Boundedness, at $\xi=1$:** the worst violation over every frame goes from
    $1.5\times10^{-3}$ to $1.2\times10^{-10}$. The slab's tables had sampled every 10th frame,
    which aliased with the element spacing and hid the SVV-off violations.
  - **At large $\xi$ it costs shape:** 0.00001 → 0.00005 at $\xi=2$ and 0.00003 → 0.00042 at
    $\xi=2.8$. The mechanism is not measured.
- **`rider_kothe`, $H=1/64$:** $E_r$ 0.0498 → 0.0446 ($\xi=1$), 0.0708 → 0.0631 ($\xi=1.5$,
  $\gamma=1$) and 0.0688 → 0.0615 ($\xi=1.5$, $\gamma=2$), at unchanged boundedness (both
  columns before the 2026-10-05 velocity fix, §4.1d).
- **What this supersedes.** The 2026-09-10 Rider–Kothe finding "SVV does not help the shape
  error" ran with the frozen diffusion (§4.1d).

## 6. $\gamma$ is a rate — and the rate itself costs accuracy

$\gamma$ multiplies **both** terms of the right-hand side — the diffusivity is
$\varepsilon\gamma u_{\max}$ and the compression term carries $\gamma u_{\max}$
— so it sets how fast the interface relaxes toward equilibrium, not which
equilibrium. $\varepsilon$ alone fixes the equilibrium tanh width. The trap: on a
test with no strain, lowering $\gamma$ does not improve the method, it *removes*
it. At $\gamma=0$ there is no compression at all.

**But "a rate" does not mean "free".** Accuracy is flat in $\gamma$ above 0.25
on the 1D slab, but **not on Zalesak**: over the $(\xi,\gamma)$
map, $E_r$ rises monotonically with $\gamma$ — about $2\times$ from 0.25 to 2 at
every $\xi \ge 1.5$. The reason is exactly the "rate" framing taken seriously: a
rigid rotation's exact solution contains **no relaxation at all**, so once the
interface has reached its equilibrium profile, further relaxation is pure
perturbation and more of it is worse. This was confirmed against a fixed-$\Delta t$
control row ($\Delta t$ moves $E_r$ by $\le 5\%$ over a $6.8\times$ change, while
$\gamma$ moves it $2.06\times$), so it is not the $\Delta t \propto 1/\gamma$
coupling of the sweep.

Boundedness does improve with $\gamma$, but **the sign inverts with $\xi$**:
better at $\xi \ge 2$ ($4.3\times10^{-10} \to 8.8\times10^{-12}$ at $\xi=2$),
an order of magnitude *worse* at $\xi=1$, and fatal at $\xi=0.5$, where
$\gamma \ge 1$ diverges outright. The ceiling on $\gamma$ is the compression CFL,
$\gamma u_{\max}\Delta t/h_{\text{GLL,min}} \le 0.05$, enforced as a hard error —
and with `svv_psi` explicit there is a second ceiling, $\Delta t\,\rho \le 0.4$,
which binds *first* below $\gamma \approx 0.27$.

Full tables and the two heatmaps: `examples/zalesak_disk/README.md`.

**Recommended settings live in the root [`README.md`](README.md)**, which is the
single source for them. They are regime-dependent — strain inverts the $\xi$
recommendation — so don't restate them in this file.


## 7. What's actually validated, and what's narrower than it first looked

Each row states exactly what was tested — not what the headline number might
imply.

| Case | Configuration tested | Result | Evidence |
|---|---|---|---|
| 1D slab | $\psi$-normal, $\xi=\varepsilon N/H=1.0$ | $E_r$: **≈1.2** ($\phi$-normal) → **0.00006** | `examples/advecting_slab_1d/README.md`, `evidence/` |
| 1D slab | $\psi$-normal, $\xi=0.5$ (Saini's sharpest setting) | $\phi$-normal **diverges** ($t\approx6.8$) → **0.0004** | same |
| 2D Zalesak | $\psi$-normal, $\xi=2.8$, $\text{svv}_\psi.c_0=1.0$, **redistancing OFF** | $E_r=0.022$, zero boundedness violations, 10 full rotations | `examples/zalesak_disk/README.md`, `evidence/` |
| 2D Zalesak | same, matched-SVV control at $\xi=0.5$ | **diverges** at $t=1.74$ (step 34869) | isolates $\xi$ as the cause, holding SVV fixed |
| 2D Zalesak | full $(\xi,\gamma)$ grid, $\xi\in\{0.5,1,1.5,2,2.8\}\times\gamma\in\{0.25,0.5,1,2\}$, one rotation, $N=7$ | **gap closed.** $\xi\ge1.5$ bounded to round-off at every $\gamma$; boundary between $\xi=1$ and $\xi=0.5$; $\gamma$ **not** neutral ($E_r$ ~2× over the row) | `examples/zalesak_disk/evidence/zalesak_xi_gamma_svv_on.png`; 20 runs, one rotation only |
| 2D Zalesak | the two redistancing variants, $\xi=2.8$, $N=5$: `seed="phi"` and in-place `seed="psi"` | recorded: `seed="phi"` completes at $E_r$ 0.7087, worst violation $1.06\times10^{-1}$; `seed="psi"` $E_r$ 1.667 at $t=3$, stopped. **Not citable**: from the earlier configuration (analytic $\psi$, $\lvert\nabla\psi\rvert$ trigger, before the `grad_floor` fix, with the history bug). The shipped `.case` files now use `psi_init = "redistance"` and the 0.5 timer, and have not been re-run (`NEXT_SESSION.md`); $\xi=2.8$ at $N=5$ violates $N\ge3.68\xi$ (§4.1b) | `examples/zalesak_disk/README.md` |
| 2D Zalesak | arm C: built $\psi$ + periodic reseed, $\xi=1$, $N\in\{3,5,7\}$, with the history restart (2026-10-02) | $N=3$ diverges $t=11.2$; $N=5$ completes at $E_r$ 0.574 (arm B 0.034); $N=7$ diverges $t=4.5$ | §4.1c: the relaxation converges, $\phi$'s contour fragments and the reseed sustains the fragments |
| 2D Zalesak | SVV on vs off, $\xi=1$ (Saini's $\xi=1/N$), $N\in\{3,5,7\}$, ten rotations, identical $\Delta t$ per pair | SVV-on $E_r$ **0.101 → 0.021 → 0.0040**; SVV-off **0.932 → 0.975 → 1.027** | single-variable; the gap grows 9× → 47× → 256× with $N$ (§5) |
| Rider–Kothe | re-run 2026-10-02 with the diffusion fix and 2026-10-05/06 with the velocity fix (§4.1d): $\psi$-normal, `svv_psi` $c_0=0.1$, redistancing **OFF**; $\xi\in\{1,1.5,2\}$ × $\gamma$ at $H=1/64$; $h$-series to $H=1/128$ at $\xi=1$ | $E_r(t{=}8)$ 0.0452 → 0.0104 under $h$-refinement (rate ≈2.1); $\xi$ trades shape (0.0452 at $\xi=1$, 0.0710 at $\xi=2$) against boundedness ($2.5\times10^{-3}$ → 0); band $\lvert\nabla\psi\rvert$ to 15 at $t=4$, back to 1.008 at $t=8$, normal within 0.6–1.3° of exact; SVV on $\psi$ worth 10–11% in $E_r$. Saini's full periodic-reseed configuration (scratch user file) completes, $E_r$ 0.0407, but its $\psi$ normal is 4–8° off | `examples/rider_kothe/README.md`; `evidence/` predates the fix |
| Saini §4.4 circles | Eq. (44) standalone with the authors' configuration from their public case: sign-function $\varepsilon=0.25$ in Eq. (46), dealiased $\mathbf C(\mathbf w)$, their sign guard; their BDF2/EXT2, $c_0{=}2$, $N_{svv}{=}N/6$, $\tau{=}6$, over the 18-cell Fig. 12 grid | **reproduced** (2026-09-29). 12 of 18 cells within 0.1% of their values, four more within 2.2%; the four Table 2 cells are all within 2.2%. Two cells ($H{=}1/5$, $N{=}4$ and 8) fail on our mesh through zero-set round-off, and match exactly when started from their zero-set $\psi_0$ | `examples/redistance_circles/README.md` §4–§5, `evidence/` |

The $\phi$-normal runs predate the `grad_floor` fix (§4.1) and were not re-run;
read them for direction only.

The $(\xi,\gamma)$ grid ran **one** rotation, not ten: it maps the operating
envelope (the stability boundary, $\gamma$ not neutral) but does not measure
accumulated degradation.

Error norms ($E_r$, $E_v$, $E_s$, boundedness) are Saini et al.'s Eqs.
79–81 — see `references/saini_2026_test_cases.md` §2 for definitions.

## 8. Evidence local to this repo

Each case's `evidence/` folder holds figures and animations produced **by this
repo's own runs**, generated by that case's `visualize.ipynb` and regenerable by
rerunning it — except `redistance_circles`, whose figures and animations come
from scripts under its gitignored `logs/saini_case/` (local only).

**`advecting_slab_1d/`** — `slab_1d_methods.mp4` (one panel per method:
$\gamma=0$, $\phi$-normal, $\psi$-normal, with the exact slab overlaid);
`slab_1d_cdi_off_vs_phi.mp4` (its first two panels only);
`slab_1d_psi_field.mp4`; `slab_1d_grad_psi.png`.

**`zalesak_disk/`** — `zalesak_methods.mp4` (four panels, each differing from the
$\xi=2.8,\ \gamma=1$ reference in exactly one setting, colour scale running past
$[0,1]$ so a boundedness violation shows as saturation);
`zalesak_xi_gamma_svv_on.png` (the $(\xi,\gamma)$ map, §6);
`saini_fig5_contours.png`, `saini_fig6_slot_zoom.png`, `saini_fig7_norms.png` —
Saini et al. §4.3 recreated, with Fig. 7 carrying the SVV-off family their paper
does not show (§5).

**`rider_kothe/`** — `rider_kothe_methods.mp4`; `rider_grad_psi.png`
($|\nabla\psi|$ drifting and returning); `rider_xi_and_h.png`.

**`redistance_circles/`** — `fig11_ic_and_exact.png`, `fig12_error_decay.png`,
`fig13_error_maps.png`, `fig13_vs_saini.png` (Saini Figs. 11–13 recreated with the
authors' configuration, §4.4; Fig. 12 with their values dashed, Fig. 13 on their
colour scale and against their own field); `er_vs_tau.png` ($E_r(\tau)$ to 24 in
the Table 2 cells, against their code); and the `anim*.mp4` clips.
The list and captions are in that case's README §7.

Material carried over from the upstream investigation has been **removed** where
it showed a method this repo does not implement — Jain's algebraic normal, the
$\Gamma$ and SVV sweep points, the `psi_init` ladder, the 1D redistancing
overlay. Those live in `../neko-multiphase/` and belong there.

## 9. Explicit non-goals

This repo does not re-litigate whether CDI works without $\psi$ at all — it
doesn't (see `../neko-multiphase/CDI_IN_SEM.md`), and that's the sibling
repo's investigation to own, not this one's to repeat. Specifically out of
scope here:

- Jain (2022)'s algebraic normal — tried, doesn't thread the
  accuracy/boundedness needle (`../neko-multiphase/CDI_IN_SEM.md` §6.2).
- The strong-form vs. weak-form divergence question for the compression
  term — an open, separate investigation
  (`../neko-multiphase/examples/cdi_face_artifact/`). This repo inherits the
  strong-form discretization as already validated and does not take a
  position on it.
- Parameter sweeps beyond what each case's `README.md` documents, and any
  ablation set beyond `zalesak_disk`'s one sanctioned redistancing ablation.
  Anyone wanting the full $\xi$/$\Gamma$/SVV sweep ladder should look at
  `../neko-multiphase/examples/saini_benchmarks/` directly rather than
  expecting it reproduced here.
