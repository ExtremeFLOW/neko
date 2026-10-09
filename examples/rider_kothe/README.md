# rider_kothe — Saini et al. (2026) §4.5 with Neko's CDI

A disk is stretched into a thin spiral filament by $t=4$, then the flow reverses and it unwinds back
to a disk at $t=8$. It is the only case with real strain.

**Status.**
- **Transport only** (the shipped `.case` files, no re-initialization) is validated on the current
  build. The tables are from the re-runs of 2026-10-05/06, after the velocity fix. Don't quote
  Rider–Kothe numbers from before 2026-10-02.
- **Re-initialization in Saini's configuration** ran on scratch user files only, and is not in
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
| mesh, order | $H\in\{1/32,1/64,1/128\}$, $N\in\{4,5,6\}$ | | $H\in\{1/64,1/96,1/128\}$, $N=5$ |
| $\Delta t$ | $\{8,4,2\}\times10^{-4}$ (CFL $\approx0.4$ at $N=6$) | $4\times10^{-4}$ | $8\times10^{-5}$ at $H=1/64$, $4.1\times10^{-5}$ at $1/128$ |
| time scheme | BDF3/EXT3 at $H=1/128$, else BDF2/EXT2 | BDF2 (his shipped $H=1/128$ case too), extrapolated with EXT3 coefficients | BDF3/EXT3 |
| phase-field re-sharpening, Eq. (39) | every $\Delta t_{cls}=0.05$, $\Delta\tau_{cls}=0.1H/(N{+}1)$, $N_{cls}=\varepsilon/\Delta\tau_{cls}$ | same, also at $t=0$; 12 steps | none: the CDI compression acts every step |
| re-initialization interval $\Delta t_{tls}$ | 0.5, "including at the initial step" | same | none in the shipped cases; 0.5 on the scratch files |
| $\Delta\tau_{tls}$, steps | $H/(N{+}1)$, $N_{tls}=2.5H/\Delta\tau_{tls}$ (§3.4) | $H/(N{+}1)$, 150 steps | scratch: $0.1\,h_{\text{GLL,min}}$ (2129 steps) or $H/(N{+}1)$ |
| pseudo-time extent | $2.5H$ | $25H$ | scratch: $25H$ |
| Eq. (44) SVV | $N/4$, $c_0=2$ | same, $\lvert\mathbf c\rvert=\lvert\operatorname{sgn}\psi\rvert$ | scratch: the same, printed form |
| pseudo-CFL | $\approx0.24$ at $N=6$ | 0.20 at $N=5$ for Eq. (39); Eq. (44) at $\lvert\mathbf c\rvert=1$ would be 10× | — |

## 3. Results

$E_r$, $E_s$ and $E_v$ are Saini's Eqs. (79)–(81) as `../norms.py` computes them, at $t=8$.
- **Over every output frame** (every 0.08, or 0.04 at $H=1/96$ and 1/128): the worst violation, the
  mass drift and the areas.
- **Sources:** `logs/d5/gamma_sweep.py`, collected in `logs/table_velfix_2026-10-06.txt`. $E_s$ comes
  from `logs/d5_2026-10-06/tables/full_*.txt` and `logs/d6_2026-10-07/tables/summary_*.txt`.
- **End times:** `rider_kothe_xi10` ends one step past the reversal, at $t=8.00008$, because Neko's
  summed time is a round-off below 8 after 100000 steps. The other runs end within $4\times10^{-5}$
  of 8.

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

**The scratch settings.**
- Sign width 0.25 and extent $25H$; BDF2/EXT2 with the printed Eq. (31) SVV ($N/4$, $c_0=2$);
  dealiased $\mathbf C(\mathbf w)\psi$ and his sign guard.
- Reseed every 0.5, and a $t=0$ build.
- **Still not his:**
  - the uniform $|\mathbf c|$ on `svv_psi`;
  - BDF3 in physical time;
  - EXT2 rather than his EXT3 coefficients in pseudo-time;
  - the CDI $\Delta t$.
- The angle columns give $\psi$'s normal against the exact interface on the frames right after the
  events at $t=2$–5, $\phi(1-\phi)$-weighted over the compression band.

| configuration (source) | $E_r(8)$ | $E_s(8)$ | worst violation | $\phi$ pieces, $t=3$–6 | angle, $t=2$–5 | at $t=8$ |
|---|---|---|---|---|---|---|
| transport only (`rider_kothe_xi10`) | 0.0452 | 0.0237 | $2.6\times10^{-3}$ | 1–2 | 0.6–0.8° | 1.3° |
| **Saini's settings, our $\Delta\tau$** (`d6_2026-10-07/runA_s025_25H_dg`) | **0.0402** | 0.0173 | $7.1\times10^{-4}$ | 2–3 | 4.3–8.2° | 29° |
| the same, his $\Delta\tau=H/(N{+}1)$, build before the velocity fix (`d5/rk_eLf`) | 0.0407 | 0.0179 | $8.9\times10^{-4}$ | 3 | 4.3–8.2° | 29° |
| his own code (`d5/saini/cv64n5`), our normalisation | 0.0410 | 0.0121 | $8.6\times10^{-2}$ | | | |

There is no run of his configuration at $H=1/128$ yet. The records are `logs/d5/PREREGISTERED.txt`
and `logs/d6_2026-10-07/tables/`.

- **Against transport** it is lower in $E_r$, $E_s$ and the worst violation.
- **Against his code** its $E_r$ is level. His $E_s$ is lower, and his worst violation far higher:
  his phase field is unbounded.
- **$\psi$ is not a distance in the band.** Its band-mean $|\nabla\psi|$ is 1.65–1.66 after each
  event (1.61 after the last), against 1.64 in his code right after its build (2.32 at $H=1/128$, $N=3$).
- **Its normal is 4–8° off**, against 0.6–0.8° for transport. Two causes:
  1. **The reseed takes $\phi$'s 0.5 contour as the interface.** $\phi$'s own normal is 3.0–5.5° off at
     $t=3$–5 (0.9–1.0° under transport). Once $\psi$ follows $\phi$, nothing pulls $\phi$ back,
     whereas a transported $\psi$ is independent of $\phi$ and keeps correcting it.
  2. **Dealiasing moves nodes across zero, and the guard freezes them where they land.** It does not
     put them back.
     - After the second event ($t=1.00016$) $\psi$ is one piece but 18 nodes stay crossed: its zero set is
       up to $0.41\varepsilon$ from $\phi$'s contour, and its normal 1.30° off.
     - Without the pair: 0 nodes, $0.29\varepsilon$, 0.74°.
     - Over the run, without the pair, $E_r(8)$ is 0.0342 and the normal 1.0–2.6° off. So the pair
       costs about 18%, and which of the two carries it is not measured.
- **The thin tail.** At $t=3.5$–4.5 the exact filament is thinner than $2\varepsilon$ over 5–6% of
  its length. There even transport-only $\phi$ has a core of only 0.33–0.35, so there is no 0.5
  contour to reseed from. In `rk_eLf`, $\psi$'s normal across that tail is wrong at 76–79% of points.
- **Each event at $t=1$–5 raises $\psi$'s normal error**, by ×1.4–6.5, and transport brings it back
  to 0.23–0.86 of that by the next event.
  - The first event lowers it: 1.07° → 0.49°.
  - From $t\approx6$ it grows between events, to 16.1–16.5° at $t=7.52$ and 29° at $t=8$ (transport
    1.4° and 1.3°). His TLS shows the same rise, 12–17° at $t=7.5$–8.
  - Each event resets the band $|\nabla\psi|$ to 1.65–1.66 (1.61 at the last).

### 3.5 The two knobs: pseudo-time and frequency

- **Extent matters.** At width 0.25 the seed's sign function is
  $\operatorname{sgn}(\psi_0)\approx0.2(\phi-\tfrac12)\le0.1$, so its characteristics move at about
  0.1, not the 1 that §3.4 of the paper assumes.
  - With the printed $2.5H$, $\psi$ stays the seed: its normal is 37° off by $t=0.5$, and $\phi$ is
    57% thicker before the first event. An earlier full run at these settings
    (`logs/d5/run_rk_a.log`, before the velocity fix) diverged at $t=7.8$.
  - His code's $25H$, ten times the time at a tenth of the speed, gives roughly the paper's reach.
  - In 1D the slope at $\psi=0$ then follows the logistic law
    $\dot y=\tfrac{1}{2\varepsilon}y(1-y)$ from the seed's 19.1 to 2.79, whatever the SVV.
  - In 1D the reach in units of $H$ shrinks with refinement (2.1$H$ → 1.5$H$ from $H=1/40$ to
    1/80).
  - 0.25 slows every near-interface process by $0.25/(H/N)$: 80 here at $H=1/64$, $N=5$.
    That is 3.75–40 on the §4.4 grid.
  - Source for the 1D facts: `git show 80201cf13e6:examples/redistance_circles/eps_1d/README.md`, §4–5.
- **Step size barely matters.** 2129 steps ($0.1\,h_{\text{GLL,min}}$) and 150 ($H/(N{+}1)$) over
  $25H$ give 0.0402 and 0.0407, though the two runs also differ in build and file. $H/(N{+}1)$ is
  safe at width 0.25 because $|\mathbf w|$ is small near the interface (in 1D within $4\times10^{-4}$
  of a fine-$\Delta\tau$ run). At width $\varepsilon$ it is not: on the slab it is a pseudo-CFL of
  2.75, and 40 events destroy the run (`../advecting_slab_1d/README.md` §3).
- **Frequency has never been varied.** $\Delta t_{tls}=0.5$ in every production run. The per-event pattern of
  §3.4 is what makes it the next thing to try (`../../NEXT_SESSION.md`).

### 3.6 Saini's own code on this problem

His fork (`nandu90/Nek5000@nekLS` 5e9b0ae; a local copy is in `logs/d5/saini/`) reproduces his
paper:
- as shipped ($H=1/128$, $N=3$): $E_r$ 1.966e-3, against Table 3's 1.96e-3;
- at our resolution ($H=1/64$, $N=5$, $\Delta t=4\times10^{-4}$): 3.13e-3, against Fig. 17's ≈3.3e-3.

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

## 4. Limitations and open items

- **Boundedness.** Transport at $\xi=1$ is bounded only to $2$–$3\times10^{-3}$, at every $H$.
- **Re-initialization.** The committed events path is not Saini's (appendix, first row), and his
  configuration is not in the code yet; that is the next step (`../../NEXT_SESSION.md`).
- **$\psi$ under re-initialization.** The rebuilt $\psi$'s normal is worse than the transported one's
  (§3.4).
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

**From `visualize.ipynb`** (`rider_kothe_xi10`, `rider_h128`):
- **`rider_snapshots.png`:** $\phi$ and $\psi$ at $t=4$ and $t=8$, $H=1/64$, against the
  marker-advected exact interface (dashed).
  - $\phi<0$ shows yellow: up to $9\times10^{-4}$ at $t=4$, and the run's worst, $2.6\times10^{-3}$,
    at $t=3.2$.
  - The pale streak in $\psi$ at $t=4$ is small and negative, not a zero set.
- **`rider_grad_psi.png`:** band $|\nabla\psi|$ against $t$ at $H=1/64$ and 1/128 (mean, and the
  5–95% range shaded).

**From `logs/anim/anim_rk.py` and `logs/figs/figs_rk.py`** (local). In these legends, "run 2" is
Saini's re-initialization settings without dealiasing and guard (appendix), "Run A" is §3.4's
"Saini's settings, our $\Delta\tau$", and "run 1" is BDF2 with the printed SVV at width $\varepsilon$
(appendix).
- **`rider_transport.mp4`:** transport only, $H=1/64$ and 1/128 side by side, every 0.08. Yellow
  marks anything below 0; the patch at $H=1/128$ from $t\approx6.7$ is below $10^{-6}$.
- **`rider_convergence.png`:** $E_r(8)$ against $H$ (with run 2 beside it), and against $\xi$ with
  each point's worst violation.
- **`rider_redistancing_normals.png`:** $\psi$'s normal against the exact interface on ten frames
  next to events, and $E_r(8)$, for transport only and four events runs.
- **`rider_redistancing_during.png`:** $\psi$'s normal at every frame, and the band-mean
  $|\nabla\psi|$ before and after each event, for transport only, run 2 and Run A.
- **`rider_redistancing_contours.png`:** $\phi=0.5$ against the exact interface at $t=4$ (whole domain
  and the tail, whose tip no run reaches) and at $t=8$, for transport only, run 2 and Run A.
- **`rider_redistancing_dealias.mp4`:** transport only, run 2 and Run A. At $t=4$ Run A's $\psi=0$
  wiggles along the outer arm; at $t=8$ both rebuilt zero contours are ragged across the top of the
  disk.
- **`rider_redistancing_events.mp4`:** transport only, the committed events path, run 1 and run 2.
  The top row is $\phi$; the bottom row is $\psi$, with its zero contour solid. At the tail tip the
  rebuilt zero contour stops where $\phi$'s does, while transport's follows the exact interface.

## Appendix: other configurations

These runs are not in Saini's configuration. Scratch user files: `logs/d5_2026-10-06/rider_kothe_0c.f90`
and `logs/d6_2026-10-07/rider_kothe_d6.f90`. Columns: $E_r$, $E_s$, worst violation, $\phi$ pieces
at $t=3$–6, and $\psi$'s normal at $t=2$–5 / $t=8$.

| run (source) | how it differs from Saini | numbers | note |
|---|---|---|---|
| committed events path (`events_2026-10-05/rk_events`) | SSP-RK3, uniform-$\lvert\mathbf c\rvert$ SVV (`svv_step_imp`), width $\varepsilon$, $2.5H$, $\Delta\tau=0.1h_{\text{GLL,min}}$, no dealiasing or guard | 0.835, 0.340, 0.41, 9–60, 12–47° / 86° | its SVV step makes $\psi$'s spurious zero sets; each reseed copies $\phi$'s fragments (98 pieces at $t=8$) |
| run 1 (`d5_2026-10-06/rider_bdf2_eq31`) | BDF2 and the printed Eq. (31), but width $\varepsilon$, $2.5H$, no dealiasing or guard | 0.0960, 0.0367, $6.4\times10^{-3}$, 2–7, 5.7–12° / 54° | $\phi$'s thin tail breaks from $t\approx2.5$; the reseeds copy the breaks |
| run 2 (`d5_2026-10-06/rider_bdf2_eq31_s025_25H`) | Saini's settings without dealiasing and guard | 0.0342, 0.0147, $6.2\times10^{-4}$, 2–3, 1.0–2.6° / 29° | the best $\phi$ of all runs; $\psi$'s band $\lvert\nabla\psi\rvert$ 1.66 |
| run 2 at $H=1/128$ (`d5_2026-10-06/rider_h128_bdf2_eq31_s025_25H`) | the same | 0.0085, 0.0045, $1.3\times10^{-3}$, 1, 0.5–1.1° / 7.3° | transport: 0.0104; band $\lvert\nabla\psi\rvert$ 2.7 |
| width 0.25 at the printed $2.5H$ | extent | normal 37° off by $t=0.5$; 9.5° and 19 crossed nodes after event 2 | the measured runs stop at $t\approx1$; an earlier full run (before the velocity fix) diverged at $t=7.8$ |
| after event 2 ($t=1.00016$): committed path | as above | 3 $\psi$ pieces, 228 crossed nodes, zero set up to $3.2\varepsilon$ from $\phi$'s, 6.3° | |
| after event 2: SSP-RK3 + printed Eq. (31) | | 1, 8, $0.41\varepsilon$, 4.3° | the printed SVV form cuts the crossings from 228 to 8 |
| after event 2: run 1 | | 1, 7, $0.41\varepsilon$, 4.3° | BDF2 equals RK3 here |
| after event 2: run 1 + dealiasing and guard | width $\varepsilon$ | 4, 131, $5.2\varepsilon$, 5.9° | |
| SVV off, $\xi=1$, $\gamma=1$ (`lambda_fix_nosvv_2026-10-02/`) | no `svv_psi` | 0.0498 (on: 0.0446), $2.5\times10^{-3}$ both | before the velocity fix |
| SVV off, $\xi=1.5$, $\gamma=1$ and 2 | no `svv_psi` | 0.0708 and 0.0688 (on: 0.0631, 0.0615) | before the velocity fix |
