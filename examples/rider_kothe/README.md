# rider_kothe — Saini et al. (2026) §4.5 with Neko's CDI

A disk is stretched into a thin spiral filament by $t=4$, then the flow reverses and it unwinds back
to a disk at $t=8$. It is the only case with real strain.

**Status.**
- **Transport only** (the shipped `.case` files, no re-initialization) is validated on the current
  build. The tables are from the re-runs of 2026-10-05/06, after the velocity fix. Don't quote
  Rider–Kothe numbers from before 2026-10-02. It works here only because the flow is time-reversible
  on a periodic box: $\psi$ starts as the exact distance and is carried back to it. It is not an
  option for a general flow.
- **Re-initialization in Saini's configuration** (§3.4) ran on a scratch user file and is not in
  `rider_kothe.f90`. Adopting it is `../../NEXT_SESSION.md` item 1.

$\xi=\varepsilon N/H$ (`../../CDI_METHOD.md` §1). Paths under `logs/` are gitignored and local.

## 1. The problem

Saini §4.5. A disk of radius 0.15 at $(0.5,0.75)$ in $[0,1]^2$ is advected by

$$u = \sin^2(\pi x)\sin(2\pi y)\cos(\pi t/T),\qquad v = -\sin(2\pi x)\sin^2(\pi y)\cos(\pi t/T),\qquad T=8$$

(his Eqs. 85–86; divergence-free, non-penetrating, time-reversible). The flow returns the disk to its
start, so the error is quoted at $t=8$ only. Boundedness, mass and $|\nabla\psi|$ need no reference
and are valid throughout.

What his figures measure:
- **Fig. 16:** $\phi$ (his CLS) and $\psi$ (his TLS) along $y=0.75$ at $t=8$, $H=1/128$, $N=4$–6.
- **Fig. 17:** $E_r$ at $t=8$ under $h$- and $p$-refinement. **Fig. 18:** the $L^1$ norm, $|E_v|$
  and $E_s$, against other methods.
- **Figs. 19–22:** mean thickness $l_{avg}/l_0$, the $L_\infty$ of his CLS, the 0.5 contours at $t=8$,
  and $|E_v|$ against $t$.
- **Table 3:** pairs at equal GLL count, $N=3$ against $N=7$.

He notes (p. 21) that re-distancing "slows down the convergence rate of the coupled algorithm".

## 2. Configuration

**Shipped cases.** All are transport-only: `svv_psi` $c_0=0.1$, $N_{svv}=N/2$, $\psi_0$ the exact
periodic distance, BDF3/EXT3.

| `.case` | $\xi$ | $\gamma$ | $H$ | $\Delta t$ | note |
|---|---|---|---|---|---|
| `rider_kothe` | 1.5 | 0.5 | 1/64 | $1.6\times10^{-4}$ | bounded to round-off |
| `rider_kothe_xi10` | 1.0 | 1.0 | 1/64 | $8\times10^{-5}$ | Saini's width; the baseline below |
| `rider_x15_g10`, `rider_x15_g20`, `rider_x20_g10` | 1.5, 1.5, 2.0 | 1, 2, 1 | 1/64 | $\le8.3\times10^{-5}$ | the $\xi$ cross |
| `rider_h96`, `rider_h128`, `rider_h192` | 1.0 | 1.0 | 1/96, 1/128, 1/192 | $5.5$, $4.1$, $2.8\times10^{-5}$ | $h$-series; `h192` has never run |

$N=5$ throughout. $\Delta t$ is set by CDI's compression limit,
$\gamma u_{\max}\Delta t/h_{\text{GLL,min}}\le0.05$ (`../../CDI_METHOD.md` §2), not by the advective
CFL.

**Both halves of the CDI balance follow the flow**, at rate $\gamma u_{\max}(t)$ with
$u_{\max}(t)=|\cos(\pi t/8)|$:
- the compression is a source term;
- the diffusion $\varepsilon\gamma u_{\max}$ is the implicit conductivity. `material_properties`
  fills both `s_lambda` and `s_lambda_tot`; that and the velocity's time level are Neko gotchas in
  `../../CLAUDE.md`.
- The `preprocess` hook prescribes $\mathbf u(t_n)$ before each step, as his `userchk` does.

**Saini §4.5 against ours.** These are the case-specific rows. Rows common to all cases are in
`../../CDI_METHOD.md` §6.

| | Saini §4.5 (paper) | his `circVortex`, $H=1/64$, $N=5$ | ours |
|---|---|---|---|
| $\xi$ | his $1/N$ (our 1) | same | 1 (baseline, $h$-series); 1.5–2 (cross) |
| mesh, order | $H\in\{1/32,1/64,1/128\}$, $N\in\{4,5,6\}$ | | $H\in\{1/64,1/96,1/128\}$, $N=5$; his configuration also $H\in\{1/64,1/128\}$, $N=5$ |
| $\Delta t$ | $\{8,4,2\}\times10^{-4}$ (CFL $\approx0.4$ at $N=6$) | $4\times10^{-4}$ | $8\times10^{-5}$ at $H=1/64$, $4.1\times10^{-5}$ at $1/128$ |
| time scheme | BDF3/EXT3 at $H=1/128$, else BDF2/EXT2 | BDF2 (his shipped $H=1/128$ case too), extrapolated with EXT3 coefficients | BDF3/EXT3; his configuration: BDF2 with Neko's modified EXT3 |
| phase-field re-sharpening, Eq. (39) | every $\Delta t_{cls}=0.05$, $\Delta\tau_{cls}=0.1H/(N{+}1)$, $N_{cls}=\varepsilon/\Delta\tau_{cls}$ | same, also at $t=0$; 12 steps | none: the CDI compression acts every step |
| re-initialization interval $\Delta t_{tls}$ | 0.5, "including at the initial step" | same | none in the shipped cases; his configuration: 0.5, and a $t=0$ build |
| $\Delta\tau_{tls}$, steps | $H/(N{+}1)$, $N_{tls}=2.5H/\Delta\tau_{tls}$ (§3.4) | $H/(N{+}1)$, 150 steps | his configuration: $H/(N{+}1)$, 150 steps (151 at $H=1/128$) |
| pseudo-time extent | $2.5H$ | $25H$ | his configuration: $25H$ |
| Eq. (44) SVV | $N/4$, $c_0=2$ | same, $\lvert\mathbf c\rvert=\lvert\operatorname{sgn}\psi\rvert$ | his configuration: the same, printed form |
| pseudo-CFL | $\approx0.24$ at $N=6$ | 0.20 at $N=5$ for Eq. (39); Eq. (44) at $\lvert\mathbf c\rvert=1$ would be 10× | — |

## 3. Results

$E_r$, $E_s$ and $E_v$ are Saini's Eqs. (79)–(81) as `../norms.py` computes them, at $t=8$.
- **Over every output frame** (every 0.08, or 0.04 at $H=1/96$ and 1/128): the worst violation, the
  mass drift and the areas.
- **Sources, transport only:** `logs/d5/gamma_sweep.py`, collected in `logs/table_velfix_2026-10-06.txt`.
  $E_s$ comes from `logs/d5_2026-10-06/tables/full_*.txt` and `logs/d6_2026-10-07/tables/summary_*.txt`.
- **Sources, re-initialization (§3.4–§3.6):** `logs/s45_2026-10-12/measure.sh` on every frame, and
  `logs/s45_2026-10-12/tables/`. Its `PREREGISTERED.txt` holds every prediction and verdict.
- **$|E_v|$** is reported at $t=8$ and as the maximum over the run. Below about $10^{-5}$ it reflects
  our solver tolerance (an absolute residual of $10^{-6}$), not the method.
- **End times:** `rider_kothe_xi10` and the $H=1/64$ re-initialization runs at $\Delta t=8\times10^{-5}$
  end one step past the reversal, at $t=8.00008$, because Neko's summed time is a round-off below 8 after 100000 steps. The
  other runs end within $4\times10^{-5}$ of 8.

### 3.1 Transport only, $\xi=1$

$|\nabla\psi|$ drifts by an order of magnitude and comes back. Band $\phi(1-\phi)>10^{-4}$, from
`run_rider_kothe_xi10.log`:

| $t$ | band $\lvert\nabla\psi\rvert$ min | mean | max |
|---|---|---|---|
| 0.08 | 0.805 | 1.026 | 1.28 |
| 2.0 | 0.0494 | 5.48 | 10.7 |
| **4.0** (max stretch) | **0.0394** | **7.69** | **15.2** |
| 6.0 | 0.0480 | 5.42 | 10.0 |
| **8.0** | **0.778** | **1.008** | **1.19** |

**Its direction never leaves.** The band-weighted mean of $\psi$'s normal error against the exact
interface is 0.62–1.45° from $t=0.24$ on, and less before
(`logs/d6_2026-10-07/tables/rdq_all_xi10.txt`; the maximum is at $t=7.76$).
- That mean includes the 5–6% of the filament thinner than $2\varepsilon$.
- There, 1–10% of the points have a normal more than 25° off at $t=3.5$–4.5.

**Saini's local $|\mathbf c|$ on `svv_psi`** (scratch file, `logs/d6_2026-10-07/tables/`).
- **The difference.** Ours takes $|\mathbf c|=u_{\max}$ once, at step 1. His is the local
  $|\mathbf u(\mathbf x,t)|$, left of the operator.
- **The cost.** It gives $E_r(8)$ 0.0470 against 0.0452 and $E_s$ 0.0247 against 0.0237, at the same
  worst violation.
- **The normal is no worse:** 0.62–0.70° against 0.64–0.81° at $t=2$–5, and 0.83° against 1.25° at
  $t=8$.

### 3.2 $\xi$ trades shape against boundedness

$H=1/64$, $N=5$:

| $\xi$ | $\gamma$ | $E_r(t{=}8)$ | worst violation | area($\phi>0.5$)/$A_0$ at $t=4$ | $\phi$ peak at $t=4$ |
|---|---|---|---|---|---|
| **1.0** | 1.0 | **0.0452** | $2.6\times10^{-3}$ | 0.978 | 0.996 |
| 1.5 | 0.5 | 0.0632 | $\sim10^{-30}$ | 0.932 | 0.976 |
| 1.5 | 1.0 | 0.0637 | $\sim10^{-70}$ | 0.931 | 0.977 |
| 1.5 | 2.0 | 0.0618 | **0** | 0.930 | 0.977 |
| 2.0 | 1.0 | 0.0710 | **0** | 0.872 | 0.940 |

The mass drift $|m/m_0-1|$ is at most $4.9\times10^{-6}$ in every row.

- **Why.** At CDI equilibrium a filament of thickness $d$ holds a core of $1-e^{-d/2\varepsilon}$.
  The exact filament at $t=4$ is $7.3\,H/N$ thick at its median: $7.3\varepsilon$ at $\xi=1$, and
  $3.7\varepsilon$ at $\xi=2$. So a larger $\varepsilon$ keeps less of it above 0.5 (0.978 → 0.872),
  and a smaller one is less well bounded.
- **Without strain the trade differs.** On the slab, $\xi=1.5$ is both the most accurate and exactly
  bounded. That is why the slab and Zalesak, both isometries, don't settle the recommendation.

### 3.3 $h$-refinement

$\xi=1$, $\gamma=1$, $N=5$:

| $H$ | $\varepsilon$ | $E_r(t{=}8)$ | worst violation | area($\phi>0.5$)/$A_0$ at $t=4$ | $\phi$ peak at $t=4$ |
|---|---|---|---|---|---|
| 1/64 | $3.13\times10^{-3}$ | 0.0452 | $2.6\times10^{-3}$ | 0.978 | 0.996 |
| 1/96 | $2.08\times10^{-3}$ | 0.0215 | $3.0\times10^{-3}$ | 0.990 | 1.000 |
| 1/128 | $1.56\times10^{-3}$ | **0.0104** | $2.4\times10^{-3}$ | 0.998 | 1.000 |

- **$E_r$ falls 4.3× per doubling**, an observed rate of about 2.1. At fixed $\xi$ the filament gets
  thicker in units of $\varepsilon$ as $H$ shrinks.
- **The worst violation does not shrink.** It stays at $2$–$3\times10^{-3}$. For exact boundedness,
  use $\xi\ge1.5$.

### 3.4 Re-initialization in Saini's configuration

**The configuration.** The scratch user file is `logs/s45_2026-10-12/rider_kothe_s45.f90`; the rows are
`../../CDI_METHOD.md` §6.
- $\psi$ is built at $t=0$ by Eq. (44) from the Eq. (47) seed, then reseeded and relaxed every 0.5.
- Eq. (44):
  - sign width 0.25 and extent $25H$;
  - $\Delta\tau=H/(N{+}1)$, 150 steps (151 in the $H=1/128$ cases, whose $\Delta\tau$ is written a hair
    below $H/6$);
  - BDF2 with the printed Eq. (31) SVV ($N/4$, $c_0=2$);
  - dealiased $\mathbf C(\mathbf w)\psi$ and his sign guard.
- `svv_psi`: $N/2$, $c_0=0.1$, local $|\mathbf u|$ left of the operator.
- BDF2 in physical time (`time_order` 2, which extrapolates with Neko's modified EXT3, as his code does).
- **Still not his**, besides the intended CDI $\phi$ and CDI $\Delta t$:
  - `svv_psi` is explicit; his is implicit.
  - Eq. (44) extrapolates with EXT2 rather than his EXT3 coefficients. With EXT3 (`s6_ext3`) the same
    run gives $\psi$'s normal to 0.05° through $t=6$. At $t=8$ its errors are at the edge of the noise
    band below: $E_r$ 0.0408, $E_s$ 0.0188, but a worst violation of $2.0\times10^{-3}$.

**Saini's three metrics at $t=8$**, $\xi=1$, $\gamma=1$, $N=5$. His values are his code's dumps
rescored with `../norms.py` (§3.6).

| $H$ | run (source) | $E_r$ | $E_s$ | $\lvert E_v\rvert$ at $t=8$ / max | worst violation |
|---|---|---|---|---|---|
| 1/64 | ours (`s1_h64`) | 0.0379 | 0.0182 | $4.3\times10^{-6}$ / $8.1\times10^{-6}$ | $9.2\times10^{-4}$ |
| 1/64 | his code (`d5/saini/cv64n5`) | 0.0410 | **0.0121** | $1.5\times10^{-6}$ / $2.0\times10^{-6}$ | $8.6\times10^{-2}$ |
| 1/128 | ours (`s2_h128`) | **0.00852** | $4.60\times10^{-3}$ | $3.5\times10^{-7}$ / $1.8\times10^{-6}$ | $9.3\times10^{-4}$ |
| 1/128 | his code (`d5/saini/cv128n5`) | 0.0119 | $\mathbf{1.63\times10^{-3}}$ | $1.9\times10^{-8}$ / $3.5\times10^{-8}$ | $4.8\times10^{-2}$ |

Sources: `logs/s45_2026-10-12/tables/full_*.txt` and `summary_all_h64.txt`.
- **The noise band.** Two repeats of the $H=1/64$ run change nothing physical: $\Delta t=7.99\times10^{-5}$
  (`n1_dt`) and $r_f=0.1001$ (`n2_rf`). With the run itself they give $E_r$ 0.0379–0.0407 and $E_s$
  0.0182–0.0192, about ±4% and ±3%. Smaller differences are not results.
- **$E_r$** is 7.6% below his at $H=1/64$, at the edge of the band, and 28% lower at $H=1/128$.
- **$E_s$** is 1.5× his at $H=1/64$ and 2.8× at $H=1/128$. His falls 7.4× per doubling, ours 4.0× (at
  $\gamma=0.5$, 3.5×; why: §3.8).
- **$|E_v|$** is at our solver tolerance (an absolute residual of $10^{-6}$). His solves converge to
  $10^{-8}$.
- **Boundedness.** Ours stays within about $10^{-3}$ of $[0,1]$: $9.2\times10^{-4}$–$2.0\times10^{-3}$
  over the repeats. His phase field leaves it by 5–14% at every resolution (§3.6).

**Where the shape is lost.** These are $E_s$ against the marker-advected exact interface at every
frame (`logs/s45_2026-10-12/es_t.py`, tables `es_t_*`).
- **From the start.** Ours is already 1.3–3.0× his up to $t=2$ at $H=1/64$, including before the
  first event, and 1.5–3.2× his throughout at $H=1/128$. Part of the gap is the phase field itself.
- **Mid-run, at some events.** $E_s$ steps up at some events, by +30% across $t=3$ at $H=1/64$.
- **Late, between events.** From $t\approx7$ it climbs between events, about +50% per half
  interval at $t=7$–8, while $\psi$ is rebuilt from an already degraded $\phi$.

**$\psi$ at the event instant**, $H=1/64$. This is the $25H$ field as the run uses it, written by the
probe of §3.5 (`tables/event_psi_s1_h64.txt`).
- **Under strain the rebuild is worse than its seed.** At $t=2$–4, $\psi$'s normal is 4.3–8.6° off the
  exact interface, while $\phi$'s own is 0.6–4.1°. Against $\phi$'s 0.5 contour, $\psi$'s normal is
  4.3–7.7° off.
- **$\phi$ follows within a frame.** At $t=3.5$, $\phi$'s normal is 2.3° at the event and 6.3° at the
  next frame, 0.02 later.
- **Late in the run the cause reverses.** From $t\approx7$, $\phi$'s own normal is the worse
  (8.9–24°), and $\psi$ tracks it.
- **Refinement shrinks the event error.** At $H=1/128$, $\psi$'s normal on the frames after events is
  1.3–2.6° at $t=2$–5, and 7.0° at $t=8$.
- **$\psi$ is not a distance in the band**, as in his code (1.64 after his build at $H=1/64$). Band
  $|\nabla\psi|$ after each event is 1.64–1.66 at $H=1/64$ and 2.54–2.69 at $H=1/128$.
- **The thin tail.** At $t=3.5$–4.5 the exact filament is thinner than $2\varepsilon$ over 5–6% of its
  length. There $\phi$ has no 0.5 contour for a reseed to start from.
  - At $t=4$, $H=1/128$, our 0.5 contour stops $5.7H$ ($\gamma=0.5$) and $7.8H$ ($\gamma=1$) short of
    the exact tail tip. His stops $0.8H$ short (`logs/s45_2026-10-12/tables/tail_tip_t4_h128.txt`,
    `evidence/rider_saini45_tail_h128.png`).
  - **The tip itself is a small part of $E_s$.** Split by Lagrangian origin (each misplaced node traced back
    to $t=0$ along the exact flow; `logs/s45_2026-10-12/es_origin.py`, `tables/es_origin_*.txt`), the
    tip's origin arc (74–99° about the disk centre) holds 2–9% of the gap to his code at $t=4$ and none
    net at $t=8$.
  - **The gap starts on the top of the disk.** At $t=1$–2, 64–76% of it comes from material starting at
    60–115°, the arm that becomes the outer filament, and it is mostly missing area there. At $t=8$
    it sits at 20–130°, with both signs.
  - **There $\phi$ and $\psi$ disagree most** (`logs/s45_2026-10-12/align_t.py`, `tables/align_*.txt`).
    At the events of $t=1$–2 the angle between $\nabla\phi$ and $\psi$'s normal is 2.6–5.3° on that arm,
    against his 1.1–1.9°, and 0.5–0.8° elsewhere.
  - **Saini's Eq. (39) diffusion form recovers the early shape but not the run** (appendix, $H=1/128$).
    His diffusion is projected onto $\mathbf n$ (`../../CDI_METHOD.md` §2); ours is isotropic and so
    also diffuses along the interface where the two normals disagree. Projected, our $E_s$ equals his
    code's to $t\approx1$ and the tail shortfall falls from $5.7H$ to $3.1$–$3.3H$. But variations
    along $\psi$'s contours are then undamped: $\phi$ leaves $[0,1]$ by up to $3\times10^{-2}$, or
    diverges, and $E_r$, $E_s$ at $t=8$ are worse than isotropic. SVV on $\phi$ at his phase-field
    transport strength changes nothing measurable.

### 3.5 The knobs on Saini's configuration

**$\gamma$, the CDI rate.** It is dimensionless: the applied rate is $\Gamma(t)=\gamma\,u_{\max}(t)$, with
$u_{\max}$ the flow's peak speed at $t=0$ times $|\cos(\pi t/8)|$. It scales both halves of the $\phi$
equation, and so sets how fast $\phi$ relaxes to the profile that $\psi$'s normal asks for. It is our
parameter; Saini has none. $\Delta t$ is kept at $8\times10^{-5}$ at $H=1/64$ ($\gamma=2$ needs
$4\times10^{-5}$) and at $4.1\times10^{-5}$ at $H=1/128$.

| $H$ | $\gamma$ (source) | $E_r$ | $E_s$ | $\lvert E_v\rvert$ at $t=8$ / max | worst violation |
|---|---|---|---|---|---|
| 1/64 | 0.25 (`s7_g025`) | 0.0177 | 0.0080 | $7.2\times10^{-6}$ / $9.8\times10^{-6}$ | $2.9\times10^{-3}$ |
| 1/64 | 0.5 (`s7_g05`) | **0.0268** | **0.0114** | $4.6\times10^{-6}$ / $5.7\times10^{-6}$ | $9.7\times10^{-4}$ |
| 1/64 | 1 (`s1_h64`) | 0.0379 | 0.0182 | $4.3\times10^{-6}$ / $8.1\times10^{-6}$ | $9.2\times10^{-4}$ |
| 1/64 | 2 (`s7_g20`) | 0.0475 | 0.0222 | $3.2\times10^{-6}$ / $6.0\times10^{-6}$ | $1.7\times10^{-2}$ |
| 1/128 | 0.5 (`s8_g05_h128`) | **0.00575** | $\mathbf{3.25\times10^{-3}}$ | $1.6\times10^{-6}$ / $2.5\times10^{-6}$ | $1.0\times10^{-3}$ |
| 1/128 | 1 (`s2_h128`) | 0.00852 | $4.60\times10^{-3}$ | $3.5\times10^{-7}$ / $1.8\times10^{-6}$ | $9.3\times10^{-4}$ |

- **$\gamma=0.5$ beats $\gamma=1$ in $E_r$ and $E_s$ at both resolutions:** −29% and −38% at $H=1/64$,
  −33% and −29% at $H=1/128$. That is 7–12× the noise band measured at $H=1/64$, at an unchanged worst
  violation.
- **Against his code at $\gamma=0.5$:** its $E_s$ at $H=1/64$ is below his (0.0114 against 0.0121). At
  $H=1/128$ it is still 2.0× his.
- **A weaker pull hurts less.** With a lower $\gamma$, $\phi$ takes in the degraded post-event normal
  more slowly.
- **$\gamma=2$ is worse on every measure but $|E_v|$.** Its worst violation is 18× larger.
- **Boundedness at $\gamma<1$.** Exact boundedness is not required, but it is nice to have. $\gamma<1$
  means $\Gamma<\max|\mathbf u|$, below the sufficient condition usually quoted for CDI.
  - In this repo boundedness against $\gamma$ depends on $\xi$. The slab needs $\gamma\ge1.5$ at $\xi=1$
    for exact boundedness, and on Zalesak the trend changes sign with $\xi$ (their READMEs).
  - Here, $\gamma=0.5$ stays within about $10^{-3}$ at $N=5$, $H\ge1/64$. It does not on the $H=1/32$
    mesh or at $N=6$ (§3.7), but no $\gamma=1$ run exists at those cells to blame $\gamma$; the
    $\gamma=1$ run at $H=1/32$, $N=6$ diverges too.
- **$\gamma=0.25$ trades boundedness for shape.** $E_r$ and $E_s$ fall further, by another 34% and
  30%. But the worst violation triples, and the maximum $|E_v|$ reaches $10^{-5}$. $\phi$ also sheds
  more short pieces: up to 9 at $t=3$–6, against 6–7 at $\gamma=0.5$–1.

**The extent is his $25H$, and a $25H$ solve is not converged under strain.** The probe
(`logs/s45_2026-10-12/probe/`) repeats the $H=1/64$ run of §3.4 with every Eq. (44) solve carried on
to $400H$ and logged every $25H$, with the $25H$ field restored afterwards. Its output is
byte-identical to that run's.
(`tables/probe_s1p_h64.txt`)
- **$|\nabla\psi|$** in the band goes from 1.66 at $25H$ to 1.22 at $50H$ and 1.09 at $75H$, at every
  event.
- **The band normal barely moves on the near-circle:** 0.14° and 0.19° from the $25H$ one at
  $t\le0.5$.
- **Under strain it does move:** 3.1–5.8° by $50H$ and 4.9–8.8° by $75H$, at $t=2$–5.
- **The steady state** is, by construction, the distance to $\phi$'s 0.5 contour, not the exact
  interface. Whether the turn moves the normal towards it, or towards the exact interface, was not
  measured.
- **From $75$–$100H$** the far field goes unstable: from $75H$ at the events from $t=5$, from $100H$
  before. The likely cause, not tested, is EXT2 at $\Delta\tau=H/6$ once $|\operatorname{sgn}\psi|$
  grows.
- **Why his code needs $25H$ rather than the printed $2.5H$.** At width 0.25 the seed's sign function
  is $\operatorname{sgn}(\psi_0)\approx0.2(\phi-\tfrac12)\le0.1$, so characteristics move at about 0.1,
  not 1. With the printed $2.5H$, $\psi$ stays the seed: its normal is 37° off by $t=0.5$, and an
  earlier full run diverged at $t=7.8$ (`logs/d5/run_rk_a.log`, before the velocity fix).

**From our earlier scratch settings to his.** The earlier settings were our $\Delta\tau$, uniform
$|\mathbf c|$ on `svv_psi` and BDF3. Changing them one at a time (appendix):
- local $|\mathbf c|$ gives +11% in $E_r$;
- then his $\Delta\tau$ gives −16%, just below the band;
- then BDF2 gives +1%.

The net is −6%, at the band's edge. The two large steps are not explained by what they change.

Re-initialization intervals other than 0.5 are in the appendix: 1.0 is mildly better, and 0.25
diverges at $t=5.69$.

### 3.6 Saini's own code on this problem

His fork (`nandu90/Nek5000@nekLS` 5e9b0ae; a local copy is in `logs/d5/saini/`) reproduces his
paper over his whole grid (below), and his Table 3's $H=1/128$, $N=3$ value (1.966e-3 against 1.96e-3).

**His $E_r$ is not on our scale.** `ls_relerr` divides by $\int\psi_e$ with $\psi_e=1$ *outside* the
disk: 0.929, against our disk area of 0.0707. With `norms.E_r` on his $t=8$ dumps he gets **0.0410**
($H=1/64$, $N=5$) and **0.0258** ($H=1/128$, $N=3$).
- **His filament retention matches ours**: his CLS keeps 0.978 of its $\phi>0.5$ area at $t=4$.
- **His low $E_r$ comes from his phase-field re-sharpening, Eq. (39).** Without it
  (`userParam02 = 0`) he gets 0.208.
- **His CLS is not bounded**: $[-0.09,1.15]$ with re-sharpening, $[-0.21,1.16]$ without
  ($H=1/64$, $N=5$).
- **Fed the same $\phi$, our Eq. (44) build and his TLS agree to 0.10° in the band-weighted mean.**
  - The check is one event at $t=0.5$, against his single-precision dump; the 95th percentile is 26°.
  - It is the one direct check of our build against his.
- **His Fig. 16b** (his TLS along $y=0.75$ at $t=8$, $H=1/128$) plateaus at 0.074.
  - That is the seed's $\pm r_f/2$ grown at width 0.25 over $25H$:
    $\tfrac12\operatorname{asinh}(\sinh(0.1)\,e^{2\cdot25H})=0.0738$. The printed $2.5H$ would give
    0.052.
  - Its slope near 1 is at $t=8$. Right after his build his code's band mean is 2.32 ($H=1/128$,
    $N=3$), and at $t=8$ it is 0.94.

**His code over his whole grid.** These are his shipped case at each $H$ and $N$ of his Figs. 17–18,
changing only $H$, `lx1` and his paper's $\Delta t$ ($\{8,4,2\}\times10^{-4}$). They are rescored in our
normalisation at $\varepsilon=H/N$ (`logs/d5/saini/mk_saini_grid.sh`, `logs/s45_2026-10-12/tables/full_saini_cv*.txt`;
$H=1/64$, $N=5$ in `logs/d5_2026-10-06/tables/full_saini_cv64n5.txt`).

| $H$ | $N$ | $E_r$ | $E_s$ | $\lvert E_v\rvert$ at $t=8$ | worst violation |
|---|---|---|---|---|---|
| 1/32 | 4, 5, 6 | 0.272, 0.182, 0.186 | 0.096, 0.062, 0.050 | $3.3\times10^{-4}$, $3.7\times10^{-5}$, $4.3\times10^{-5}$ | 0.08–0.14 |
| 1/64 | 4, 5, 6 | 0.0550, 0.0410, 0.0317 | 0.0156, 0.0121, 0.0075 | $1.0\times10^{-5}$, $1.5\times10^{-6}$, $9.4\times10^{-7}$ | 0.06–0.09 |
| 1/128 | 4, 5, 6 | 0.0159, 0.0119, 0.0099 | 0.0024, 0.0016, 0.0010 | $3\times10^{-9}$–$2\times10^{-8}$ | 0.05 |

- **They reproduce his Fig. 17a.** His logged $E_r$, in his normalisation, matches the plotted values
  within 3% in every cell but $H=1/32$, $N=6$ (5%).
- **His $E_s$ shares his $E_r$'s denominator.** Ours over his plotted values is 12.4–13.3× for $E_r$
  and 12.0–14.4× for $E_s$, in every cell. That is the outside-to-inside area ratio of 13.15, so his
  printed $E_s$ is not comparable to ours either.
- **His Fig. 18a plots $E_{L1}$, not $E_r$.** He doesn't define it; it is about 0.91–0.95× his $E_r$.
- **His Table 3 pairs equal GLL counts:** $(1/64,N{=}3)$ with $(1/32,N{=}7)$, and $(1/128,N{=}3)$ with
  $(1/64,N{=}7)$.

### 3.7 His Figs. 17–18 grid

Ours is Saini's configuration at $\gamma=0.5$, against his code at the same cell, both in our
normalisation.
- $\Delta t$ per cell is the one that gives a compression CFL of about 0.0436 at $\gamma=1$, kept at
  $\gamma=0.5$: from $2.35\times10^{-4}$ at $H=1/32$, $N=4$ down to $2.88\times10^{-5}$ at $H=1/128$,
  $N=6$.
- $\varepsilon=H/N$.
- $\Delta\tau=H/(N{+}1)$, with $25(N{+}1)$ steps (151 in the $H=1/128$, $N=5$ cell). Sources: `logs/s45_2026-10-12/tables/full_gd*_g05.txt`
(and `full_s7_g05`, `full_s8_g05_h128` for $N=5$), his cells as §3.6. Figure:
`logs/s45_2026-10-12/figs/s45_grid.png` (local).

| $H$ | $N$ | $E_r$ ours / his | $E_s$ ours / his | $\lvert E_v\rvert$ at $t=8$ / max, ours | worst violation, ours / his |
|---|---|---|---|---|---|
| 1/32 | 4 | 0.249 / 0.272 | 0.091 / 0.096 | $2.0\times10^{-6}$ / $5.3\times10^{-6}$ | 0.39 / 0.084 |
| 1/32 | 5 | 0.146 / 0.182 | 0.061 / 0.062 | $5.9\times10^{-6}$ / $8.1\times10^{-6}$ | 0.036 / 0.14 |
| 1/32 | 6 | diverged at $t=7.28$ / 0.186 | — / 0.050 | | above $10^{-2}$ from $t=5.3$ / 0.079 |
| 1/64 | 4 | 0.0392 / 0.0550 | 0.0202 / 0.0156 | $2.2\times10^{-6}$ / $3.5\times10^{-6}$ | $6.4\times10^{-4}$ / 0.059 |
| 1/64 | 5 | 0.0268 / 0.0410 | 0.0114 / 0.0121 | $4.6\times10^{-6}$ / $5.7\times10^{-6}$ | $9.7\times10^{-4}$ / 0.086 |
| 1/64 | 6 | 0.0336 / 0.0317 | 0.0137 / 0.0075 | $3.2\times10^{-6}$ / $8.3\times10^{-6}$ | 0.071 / 0.056 |
| 1/128 | 4 | 0.0078 / 0.0159 | 0.0038 / 0.0024 | $3.3\times10^{-7}$ / $1.1\times10^{-6}$ | $8.8\times10^{-4}$ / 0.049 |
| 1/128 | 5 | 0.0057 / 0.0119 | 0.0033 / 0.0016 | $1.6\times10^{-6}$ / $2.5\times10^{-6}$ | $1.0\times10^{-3}$ / 0.048 |
| 1/128 | 6 | 0.0065 / 0.0099 | 0.0032 / 0.0010 | $2.8\times10^{-6}$ / $3.4\times10^{-6}$ | $6.1\times10^{-3}$ / 0.052 |

- **$h$-refinement works wherever the runs complete.** $E_r$ and $E_s$ fall with $H$ at $N=4$ and 5
  from $H=1/32$, and at $N=6$ from $H=1/64$ (it diverges at 1/32).
- **Our $E_r$ is below his** except at $H=1/64$, $N=6$ and the diverged cell. At $H=1/128$ it is about
  half his at $N=4$–5.
- **Our $E_s$ is level with or below his at $H=1/32$** and at $H=1/64$, $N=5$, and above it elsewhere:
  1.6×, 2.0× and 3.2× at $H=1/128$, $N=4$–6. His falls faster under refinement.
- **$p$-refinement stops at $N=5$.** At $N=6$, $E_r$ is worse than at $N=5$ at $H=1/64$ and 1/128, $E_s$
  does not improve, and boundedness is lost: $7.1\times10^{-2}$ at $H=1/64$, $6.1\times10^{-3}$ at
  1/128, and divergence at 1/32. His $E_s$ falls monotonically in $N$, and so does his $E_r$ at
  $H=1/64$ and 1/128.
- **The coarsest mesh is unbounded at every $N$**, late in the run: 0.39 at $N=4$ (from $t\approx6.6$),
  0.036 at $N=5$ (from $t=7.6$).
- **No single change fixes $N=6$** (at $H=1/32$, one change each from the diverged cell):
  - EXT3 in pseudo-time diverges at $t=7.04$;
  - $\gamma=1$ diverges at $t=6.08$;
  - half the $\Delta t$ diverges at $t=6.16$;
  - the Eq. (44) SVV at $N/6$ diverges later, at $t=7.84$.

  The explicit terms are far inside their stability limits at $N=4$, 5 and 6: `svv_psi` $\Delta t\rho$
  0.011–0.013 against 0.4, compression CFL 0.022, advective CFL 0.043.

  In every run the violation passes $10^{-3}$ at $t=1.28$–1.36, long before the divergence. In the two
  frames located (late in the run, at $H=1/64$ and 1/32), the worst values sit within about 1.5$H$ of
  the interface, where $|\nabla\psi|$ is ordinary. The mechanism is open (`../../NEXT_SESSION.md`).

### 3.8 Why our $E_s$ converges more slowly

At $t=8$, $N=5$, $\gamma=0.5$, his $E_s$ falls 7.4× from $H=1/64$ to $1/128$ and ours 3.5× (§3.7). This
section takes that rate apart. It uses the outputs of §3.4–§3.7 plus repeats of his code.
- Scripts and tables are in `logs/s45_2026-10-12/`: `drift_split.py`, `persist.py` with
  `persist_stats.py`, `mass_disp.py` and `hb_predict.py`.
- `rate_tables.py` prints the window rates, the drift and the retention (→ `tables/rate_tables.txt`).
- Its `PREREGISTERED.txt` holds the predictions and verdicts.
- "His" is his code at the same cell, rescored as in §3.6.
- "In $E_s$ units" means a misplaced area divided by $2P_eA_e$, as in $E_s$.

**His rate is real.** Repeats of his runs that change nothing physical (`logs/d5/saini/mk_s45_rate.sh`;
`tables/full_saini_cv64n5_{d04,n1dt}.txt`, `full_saini_cv128n5_n1dt.txt`):
- at $H=1/64$, the dump interval alone, which perturbs his run at round-off from the first step;
- at both $H$, a 0.08% smaller $\Delta t$.

With the original runs:
- $E_s(8)$ is 0.0121–0.0138 at $H=1/64$ and $1.63$–$1.71\times10^{-3}$ at $1/128$; $E_r$ is 0.0410–0.0438
  and 0.0119.
- So his rate is 7.1–8.4× per doubling.
- At $H=1/64$ his runs agree to four digits through $t=4$, then part. His band there is ±7%, against our
  ±3% (§3.4).

**A single frame is a noisy rate.** Between frames 0.04–0.08 apart, $E_s$ changes by 7–13% (rms of the log
ratio over $t=6.5$–8) in both codes. The rates below are ratios of means over time windows. His $1/32$ run
has frames only every 0.5.

| coarse → fine | $t$ 0–2 | 2–4 | 4–6 | 6–7.5 | 7.5–8 |
|---|---|---|---|---|---|
| ours, $\gamma=0.5$, $1/32\to1/64$ | 3.30 | 4.17 | 4.16 | 4.64 | 5.35 |
| his, $1/32\to1/64$ | 3.65 | 4.26 | 4.35 | 5.25 | 4.50 |
| ours, $\gamma=0.5$, $1/64\to1/128$ | 3.25 | 3.62 | 4.12 | 4.33 | 4.30 |
| ours, $\gamma=1$, $1/64\to1/128$ | 3.16 | 3.12 | 3.85 | 4.64 | 4.81 |
| his, $1/64\to1/128$ | 3.97 | 5.40 | 5.81 | 5.53 | 6.30 |

At the coarse doubling the two codes agree within about 20%. From the coarse doubling to the fine one, his
rate rises in every window after $t=2$ and ours does not.

**Per interval, $\phi$ does the same thing in both codes.** Each reseed puts $\psi$'s zero set on $\phi$'s
contour, and transport then carries it nearly materially. So the area between $\phi$'s 0.5 contour and
$\psi$'s zero set at the last frame before the next event, $D_{drift}$, is $\phi$'s non-material motion
over that interval (`tables/drift_*.txt`).
- **Its size matches.** Ours at $\gamma=0.5$ is within a factor of about 1.6 of his at every interval to
  $t=7.5$ (2.2 at $\gamma=1$). At $H=1/128$ it is $(2.1$–$6.4)\times10^{-4}$ in ours and
  $(1.7$–$4.6)\times10^{-4}$ in his, in $E_s$ units.
- **It falls at the same rate:** per-interval ratios 3.3–6.1× (median 4.4) for ours at $\gamma=0.5$, and
  3.5–5.4× (median 4.5) for his.
- **The exception is our last interval** ($t=7.5$–8): 2.2× against his 5.5×.

**The difference is how much error carries over.** Retention is
$r_k=(E_{s,k}-D_{drift,k})/E_{s,k-1}$ at successive pre-event frames: the share of the earlier error still
present. Means over the events at $t=1$–7.5:

| | $H=1/64$ | $H=1/128$ |
|---|---|---|
| ours, $\gamma=0.25,\ 0.5,\ 1,\ 2$ | 0.60, 0.66, 0.70, 0.75 | –, 0.71, 0.75, – |
| his | 0.45 | 0.23 |

- **This accounts for the rates.** With the drift falling about 4.5× in both codes,
  $E_s\approx D_{drift}/(1-r)$ gives about 3.8× for ours and 6.3× for his, close to the window rates.
- **His error is not noise.** The misplaced area, binned by where its material started at $t=0$, has the
  same pattern in repeats of either code (correlation 0.92–1.00). It also stays in place across an event
  (`tables/persist_stats.txt`).
- **Ours locks in more with refinement.** Across an event, the pattern's correlation at $t=2.5$–6 goes from
  0.73 to 0.90 in ours, and from 0.75 to 0.78 in his.
- **Our early error lasts.** The flow is reversible, so the frame at $8-t$ has the geometry of the frame at
  $t$. At $H=1/128$ the $t\le2$ pattern correlates +0.70 with its mirror in ours and +0.27 in his ($H=1/64$:
  +0.39 and +0.32).

**Ruled out.**
- **$\psi$'s quality.** The drift per interval is the same in both codes. As a diagnostic only, CDI fed a
  $\psi$ transported from the exact distance also converges at 4.4–4.7× in every window.
- **A curvature–strain drift.** A profile that the flow compresses and that CDI relaxes fast sheds
  $-(\pi^2/3)\kappa S_n\varepsilon^2$ of area per unit length and time. Summed unsigned along the exact
  interface (`tables/hb_predict.txt`):
  - it matches the first interval at $H=1/64$: $7.5\times10^{-4}$ in $E_s$ units, measured
    $7.5\times10^{-4}$ at $\gamma=0.5$;
  - it gives $(0.3$–$2.5)\times10^{-4}$ per interval at $t=2$–6. The measured drift there is
    $(13.5$–$19.4)\times10^{-4}$ at $\gamma=0.5$, and about the same at every $\gamma$ from 0.25 to 2.
- **His error as a deviation tied to the current strain.** Neither code's pattern changes sign between $t$
  and $8-t$.
- **His error as profile shape**, a distorted profile whose mass is still where pure transport puts it.
  - At $t\ge7.4$ the contour error per material wedge follows the mass moved non-materially into that wedge,
    in both codes. At $t=8$ the regression slope is 0.7–0.8 and the correlation 0.81–0.94
    (`tables/mass_*.txt`).
  - His moved mass at $t=8$ is 0.74× ours at $H=1/128$.
  - The test is blind at $t\approx2.5$–6. The carried initial profile is then sharper than the node spacing,
    and its sampling sets the moved mass to within about 20% in every run.
  - It is also weak early: at $t\approx0.5$ only about 20 nodes are misplaced per frame.
- **The levers in our $\phi$ equation.** At $H=1/128$, with means over the events at $t=1$–6:
  - retention is 0.67 with isotropic diffusion and 0.67 with SVV on $\phi$;
  - it is 0.56 with the projected Eq. (39) diffusion, projected everywhere or only inside the band;
  - his code gives 0.14;
  - $\gamma=0.25$ gives 0.56 at $H=1/64$.

**Open.** Per interval, $\phi$ moves about as far from $\psi$'s zero set in both codes. At $t=8$ both
codes' contour errors are mostly mass moved non-materially. But his error carries over less: retention
0.45 and 0.23, against about 0.7 in ours. If his mass error also accumulates less, his non-material flux
partly reverses from one interval to the next while ours keeps its sign. This has not been measured
(`../../NEXT_SESSION.md` item 2).

## 4. Limitations and open items

- **Shape.**
  - In Saini's configuration at $\gamma=1$, our $E_s$ is 1.5× his at $H=1/64$ and 2.8× at $H=1/128$; at
    $\gamma=0.5$, 0.94× and 2.0× (§3.4–§3.5).
  - His $E_s$ falls 7.1–8.4× per doubling over his repeats.
  - Per interval, $\phi$ moves as much in both codes. Ours keeps about 70% of the earlier error at every
    $H$; his keeps 45%, then 23% (§3.8). What lets his $\phi$ shed it is open.
  - At $H=1/128$ our thin tail stops 5.7–7.8$H$ short of the tip at $t=4$; his stops $0.8H$ short (§3.4).
  - The gap starts on the top of the disk by $t=1$. His projected Eq. (39) diffusion closes it early but
    costs boundedness and the $t=8$ result (§3.4).
- **$p$-refinement.** In Saini's configuration $N=6$ loses boundedness and is no better than $N=5$;
  at $H=1/32$ it diverges. The cause is open (§3.7).
- **The coarsest mesh**, $H=1/32$, is unbounded late in the run at every $N$ (§3.7).
- **Re-initialization interval.** At $\Delta t_{tls}=0.25$, inside his range of 0.01–0.5, the run
  diverges at $t=5.69$ after 22 events (appendix).
- **$\psi$ under re-initialization.** Under strain the rebuilt $\psi$'s normal is worse than that of
  the $\phi$ it is seeded from (§3.4).
- **Not in the code.** Saini's configuration runs only on the scratch file. The committed events path
  is not his (appendix, first row). Porting it is `../../NEXT_SESSION.md` item 1.
- **Boundedness of the shipped cases.** Transport at $\xi=1$ is bounded only to $2$–$3\times10^{-3}$,
  at every $H$.
- **Housekeeping.** `rider_h192` has never run. `rider_kothe_xi10` ends at $t=8.00008$.

## 5. Running

```bash
genmeshbox 0 1 0 1 -0.1 0 64 64 1 .true. .true. .true.   # box.nmsh, H = 1/64
./run.sh                    # rider_kothe and rider_kothe_xi10
./run.sh rider_h128         # or any case by name
```

`rider_kothe_xi10` is 100000 steps to $t=8$, so launch chains detached
(`setsid nohup ./run.sh > chain.log 2>&1 < /dev/null &`). The finer meshes are `box96.nmsh`,
`box128.nmsh` and `box192.nmsh`.

## 6. Evidence

Every animation also ships as a GIF.

**Saini's configuration against his code, over his grid** (`logs/s45_2026-10-12/figs_saini45.py`; §3.4, §3.7).
Ours is CDI with his $\psi$ configuration at $\gamma=0.5$. His is his own code on the same cell. Error
norms are in his normalisation: ours divided by 13.147.
- **`rider_saini45_map.png`:** where our method does well. Ours and his values for $E_r$, $E_s$, $|E_v|$
  and the worst violation, coloured by ours over his, per cell of $H\in\{1/32,1/64,1/128\}\times N\in\{4,5,6\}$.
  - $E_r$ is 0.48–0.92× his in seven of eight cells.
  - The worst violation is 0.01–0.02× his at $H\ge1/64$, $N\le5$.
  - $E_s$ is 0.94–3.2× his.
  - $|E_v|$ is at our solver tolerance.
  - $H=1/32$, $N=6$ diverged.
- **`rider_saini45_snapshots.png`:** $\phi$ and $\psi$ at $t=4$ and $t=8$, $H=1/64$, $N=5$, ours at $\gamma=0.5$,
  against the exact interface (dashed); the layout of `rider_snapshots.png` (`logs/s45_2026-10-12/snapshots_s45.py`).
  - $\phi$ is bounded to $3\times10^{-12}$ at $t=4$; at $t=8$ it spans $[-3.2\times10^{-4},\,1+5.7\times10^{-4}]$,
    at the walls and inside the disk.
  - At $t=8$ $\psi$'s zero contour is ragged across the top-left of the disk.
- **`rider_saini45_tail_h128.png`:** the open problem. $\phi=0.5$ at $H=1/128$, $N=5$ against the exact
  interface, for ours at $\gamma=0.5$ and 1 and his code.
  - At $t=4$ (the whole spiral, and the tail zoomed), his contour reaches the tail tip and ours stops
    5.7–7.8$H$ short.
  - At $t=8$ ours is ragged across the top-left of the disk.
- **`rider_saini45_h128.mp4`:** $H=1/128$, $N=5$, every 0.04 in every column: ours at $\gamma=0.5$, ours at
  $\gamma=1$, and his code. Top row $\phi$, bottom row $\psi$ with its zero contour; dashed is the exact
  interface.
  - His code was rerun with a dump every 0.04 (`logs/d5/saini/cv128n5_d04`), otherwise unchanged. Its
    $t=8$ $E_r$ and $E_s$ equal the original run's to every logged digit.
  - Made by `logs/anim/anim_rk.py`, which reads his dumps with `saini=<dir>`.
  - The positions agree at every frame. Ours loses the tail tip at $t\approx4$ where his keeps it.
  - His $\phi$ leaves $[0,1]$ over much of the domain: $[-0.048,\,1.010]$ at $t=4$.
- **`rider_saini45_fig17.png`:** his Fig. 17, $E_r$, $|E_v|$ and $E_s$ at $t=8$ against $N$, one curve
  per $H$. Ours solid, his dotted, his $N=3$ at $H=1/128$ included.
- **`rider_saini45_fig18.png`:** his Fig. 18, the same three measures against $1/H$, one curve per $N$; x marks ours
  diverged ($H=1/32$, $N=6$). Our $|E_v|$ stays near $10^{-7}$ as $H$ falls, his drops to $10^{-9}$–$10^{-10}$.
- **`rider_saini45_fig19.png`:** his Fig. 19, the mean interface thickness $l_{avg}/l_0$ against $t$, per $H$ and
  $N$. At $H=1/128$ ours stays within about 5% of 1; his falls to about 0.8 and rises to about 1.3 after $t=6$.
- **`rider_saini45_fig20.png`:** his Fig. 20 as the worst violation against $t$, per $H$ and $N$, on a
  log scale. His figure plots $\max|\psi|$, overshoot only.
  - At $N=4$–5, $H\ge1/64$, ours has a per-frame median of $10^{-11}$–$10^{-8}$ and a maximum of about
    $10^{-3}$.
  - At $N=6$, 40% of frames exceed $10^{-3}$ at $H=1/64$, and 19% at $H=1/128$. At $H=1/32$ the
    violation grows late in the run.
  - His sits at $10^{-2}$–$10^{-1}$ throughout.
  - Source: `logs/s45_2026-10-12/tables/full_*.txt`, column 3.
- **`rider_saini45_fig22.png`:** his Fig. 22, $|E_v|$ against $t$, in our normalisation for both. Ours climbs to
  $10^{-6}$–$10^{-5}$ by $t\approx4$ and then holds; at $H=1/128$ his stays near $10^{-8}$.
- **`rider_saini45_fig21.png`:** his Fig. 21, the 0.5 contour at $t=8$ for $H=1/32$, 1/64 and 1/128
  against the exact circle, per $N$. Ours is on top, his below; ours at $H=1/64$, $N=6$ is the ragged
  one.

**From `visualize.ipynb`** (`rider_kothe_xi10`, `rider_h128`):
- **`rider_snapshots.png`:** $\phi$ and $\psi$ at $t=4$ and $t=8$, $H=1/64$, against the
  marker-advected exact interface (dashed).
  - $\phi<0$ shows yellow: up to $9\times10^{-4}$ at $t=4$, and the run's worst, $2.6\times10^{-3}$,
    at $t=3.2$.
  - The pale streak in $\psi$ at $t=4$ is small and negative, not a zero set.
- **`rider_grad_psi.png`:** band $|\nabla\psi|$ against $t$ at $H=1/64$ and 1/128 (mean, and the
  5–95% range shaded).

## Appendix: other configurations

These runs are not in Saini's configuration. Scratch user files: `logs/d5_2026-10-06/rider_kothe_0c.f90`,
`logs/d6_2026-10-07/rider_kothe_d6.f90`, `logs/s45_2026-10-12/rider_kothe_s45.f90` and
`logs/s45_2026-10-12/cdi_forms/rider_kothe_forms.f90` (the $\phi$-diffusion rows). Columns:
$E_r$, $E_s$, $|E_v|$ at $t=8$, worst violation, $\phi$ pieces at $t=3.04$, 4, 5.04 and 6, and
$\psi$'s normal on those frames' events at $t=2$–5 / $t=8$. Rows without $|E_v|$ predate its
measurement. $H=1/64$, $\xi=1$, $\gamma=1$ unless stated.

| run (source) | how it differs from Saini | numbers | note |
|---|---|---|---|
| $\Delta t_{tls}=1.0$ (`s45_2026-10-12/s3_tls10`) | 8 events instead of 16 | 0.0358, 0.0163, $8.7\times10^{-7}$, $9.3\times10^{-4}$, 1–2, 4.4–7.8° / 21° | $E_s$ −10% against 0.5, outside the ±3% noise band |
| $\Delta t_{tls}=0.25$ (`s45_2026-10-12/s3_tls025`) | 32 events | diverged at $t=5.69$ after 22 events | a $\phi$ undershoot within $0.01$ of the exact interface at (0.26, 0.155) grows from $-0.11$ at $t=5.52$ to $-2.1$ at 5.68 |
| our $\Delta\tau$, uniform $\lvert\mathbf c\rvert$, BDF3 (`d6_2026-10-07/runA_s025_25H_dg`) | $\Delta\tau=0.1h_{\text{GLL,min}}$ (2129 steps), uniform $\lvert\mathbf c\rvert$ on `svv_psi`, BDF3 | 0.0402, 0.0173, $3.4\times10^{-6}$, $7.1\times10^{-4}$, 2–3, 4.3–8.2° / 29° | |
| the same, local $\lvert\mathbf c\rvert$ (`s45_2026-10-12/r1_runB`) | our $\Delta\tau$, BDF3 | 0.0446, 0.0210, $3.3\times10^{-8}$, $9.9\times10^{-4}$, 3–4, 4.3–8.5° / 30° | outside the noise band |
| the same, his $\Delta\tau$ (`s45_2026-10-12/r2_runB_dtau`) | BDF3 | 0.0374, 0.0178, $1.1\times10^{-8}$, $9.9\times10^{-4}$, 3, 4.3–8.6° / 25° | just below the band |
| his $\Delta\tau$, uniform $\lvert\mathbf c\rvert$, BDF3 (`d5/rk_eLf`) | as the first of these three, with his $\Delta\tau$ | 0.0407, 0.0179, $3.6\times10^{-6}$, $8.9\times10^{-4}$, 3, 4.3–8.2° / 29° | before the velocity fix |
| committed events path (`events_2026-10-05/rk_events`) | SSP-RK3, uniform-$\lvert\mathbf c\rvert$ SVV (`svv_step_imp`), width $\varepsilon$, $2.5H$, $\Delta\tau=0.1h_{\text{GLL,min}}$, no dealiasing or guard | 0.835, 0.340, $1.0\times10^{-5}$, 0.41, 9–60, 12–47° / 86° | its SVV step makes $\psi$'s spurious zero sets; each reseed copies $\phi$'s fragments (98 pieces at $t=8$) |
| run 1 (`d5_2026-10-06/rider_bdf2_eq31`) | BDF2 and the printed Eq. (31), but width $\varepsilon$, $2.5H$, no dealiasing or guard | 0.0960, 0.0367, $1.4\times10^{-5}$, $6.4\times10^{-3}$, 2–7, 5.7–12° / 54° | $\phi$'s thin tail breaks from $t\approx2.5$; the reseeds copy the breaks |
| run 2 (`d5_2026-10-06/rider_bdf2_eq31_s025_25H`) | Saini's settings without dealiasing and guard | 0.0342, 0.0147, $3.4\times10^{-6}$, $6.2\times10^{-4}$, 2–3, 1.0–2.6° / 29° | the best $\phi$ of all runs; $\psi$'s band $\lvert\nabla\psi\rvert$ 1.66 |
| run 2 at $H=1/128$ (`d5_2026-10-06/rider_h128_bdf2_eq31_s025_25H`) | the same | 0.0085, 0.0045, $1.2\times10^{-7}$, $1.3\times10^{-3}$, 1, 0.5–1.1° / 7.3° | band $\lvert\nabla\psi\rvert$ 2.7 |
| width 0.25 at the printed $2.5H$ | extent | normal 37° off by $t=0.5$; 9.5° and 19 crossed nodes after event 2 | the measured runs stop at $t\approx1$; an earlier full run (before the velocity fix) diverged at $t=7.8$ |
| after event 2 ($t=1.00016$): committed path | as above | 3 $\psi$ pieces, 228 crossed nodes, zero set up to $3.2\varepsilon$ from $\phi$'s, 6.3° | |
| after event 2: SSP-RK3 + printed Eq. (31) | | 1, 8, $0.41\varepsilon$, 4.3° | the printed SVV form cuts the crossings from 228 to 8 |
| after event 2: run 1 | | 1, 7, $0.41\varepsilon$, 4.3° | BDF2 equals RK3 here |
| after event 2: run 1 + dealiasing and guard | width $\varepsilon$ | 4, 131, $5.2\varepsilon$, 5.9° | |
| $\phi$ diffusion projected, $\nabla\cdot(\varepsilon\Gamma(\nabla\phi\cdot\mathbf n)\mathbf n)$ (`s45_2026-10-12/cdi_forms/cdi_proj_h128`) | his Eq. (39) form in our CDI, an implicit split sub-step; $H=1/128$, $\gamma=0.5$ | diverged at $t=7.27$; violation $>10^{-3}$ from $t=1.36$ | tail $3.3H$ short at $t=4$ ($5.7H$ isotropic) |
| the same with SVV on $\phi$ (`cdi_forms/cdi_proj_svv_h128`) | plus `svv_phi` $N/2$, $c_0=0.1$, local $\lvert\mathbf u\rvert$, explicit | diverged at $t=7.43$ | tail $3.1H$ |
| projected inside the band, isotropic outside (`cdi_forms/cdi_proj_band_h128`) | $\mathbf D=\mathbf I-w(\mathbf I-\mathbf n\mathbf n)$, $w=\min(1,\max(0,4\phi(1-\phi)/0.01))$ | 0.00783, $3.69\times10^{-3}$, $1.4\times10^{-11}$, $3.0\times10^{-2}$ | isotropic: 0.00575, $3.25\times10^{-3}$, $1.6\times10^{-6}$, $1.0\times10^{-3}$; tail $3.1H$; $\lvert E_v\rvert$ is below Neko's $10^{-6}$ tolerance because $\phi$'s own solve has no diffusion |
| SVV on $\phi$, isotropic diffusion (`cdi_forms/cdi_svv_h128`) | `svv_phi` as above | 0.00589, $3.19\times10^{-3}$, $1.7\times10^{-6}$, $9.4\times10^{-4}$ | within noise of isotropic; tail $5.7H$ |
| SVV off, $\xi=1$, $\gamma=1$ (`lambda_fix_nosvv_2026-10-02/`) | no `svv_psi` | 0.0498 (on: 0.0446), $2.5\times10^{-3}$ both | before the velocity fix |
| SVV off, $\xi=1.5$, $\gamma=1$ and 2 | no `svv_psi` | 0.0708 and 0.0688 (on: 0.0631, 0.0615) | before the velocity fix |
