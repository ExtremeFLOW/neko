# rider_kothe

Vortex-in-a-box: a disk is stretched into a thin spiral filament by $t=4$, then
the flow reverses and it unwinds back to a disk at $t=8$.

**Status: complete with redistancing off.**
- **The tables** are from the re-runs of 2026-10-02.
- **The redistancing comparison** (`rk_Af`, `rk_eLf`) ran on 2026-10-01 on the scratch user
  file `logs/d5/rider_kothe_d5.f90`, which already filled `s_lambda_tot`; the committed
  `rider_kothe.f90` got that fix on 2026-10-02.
- **The measurements of the reseeded $\psi$** are from 2026-10-05.

It is the only case with genuine strain, which makes it the only one that can
settle whether redistancing is needed and what $(\xi, \gamma)$ to use.

Two bugs affected the earlier Rider–Kothe numbers:
- **Every number before 2026-10-02, apart from the 2026-10-01 scratch-file runs, ran with the CDI diffusion frozen at its $t=0$ value**,
  while the compression followed $u_{\max}(t)=|\cos(\pi t/8)|$. That dissolved the filament
  (see "Configuration").
- **Every redistancing run before then also broke $\psi$'s time history** at each event
  (`../../CDI_METHOD.md` §4.1d).

Both are fixed. The old results are not quoted here.

$\xi = \varepsilon N/H$ throughout — see [`CDI_METHOD.md`](../../CDI_METHOD.md) §1.

Paths under `logs/` are gitignored and local to this workstation.

## How $\psi$ is maintained here: transport and SVV only — and here that is a *result*

$\psi$ is seeded once from the exact periodic signed distance, then advected with SVV
(`svv_psi` $c_0=0.1$, $N_{svv}=N/2$, Saini's TLS-transport value; every `.case` carries it).
There is no redistancing.

In the other two cases that choice is justified by an identity: they are
isometries, so $|\nabla\psi|$ cannot drift
($D|\nabla\psi|/Dt = -|\nabla\psi|(\mathbf{n}\cdot\nabla\mathbf{u}\cdot\mathbf{n})$,
and $\nabla\mathbf{u}$ is antisymmetric or zero). **Here it does drift, hugely** —
and then comes back. `rider_kothe_xi10`, band $\phi(1-\phi)>10^{-4}$, from the run log's own
band report (`run_rider_kothe_xi10.log`, lines `|grad psi| t=`):

| $t$ | band $\|\nabla\psi\|$ min | mean | max |
|---|---|---|---|
| 0.08 | 0.805 | 1.026 | 1.28 |
| 2.0 | 0.0255 | 5.48 | 10.7 |
| **4.0** (max stretch) | **0.0391** | **7.69** | **15.2** |
| 6.0 | 0.0496 | 5.42 | 10.0 |
| **8.0** | **0.778** | **1.008** | **1.19** |

The strain integrates back to zero over the cycle, so $|\nabla\psi|$ returns to 1 (band mean
1.008; the SVV keeps it from returning exactly). **Its direction never leaves.** At every frame
measured ($t=2$–8), $\psi$'s normal is within 0.6–1.3° of the exact interface in the compression band
(`logs/d5/ev/rd_quality_eLf_Af.txt`, run `rk_Af`, the same configuration). That includes the
5–6% of the filament that is thinner than $2\varepsilon$ at $t=4$.

**Redistancing does worse here** (next-but-one section): its rebuilt $\psi$ points 4–8° off the
exact interface, and is wrong across that thin tail. See
[`../../REDISTANCING.md`](../../REDISTANCING.md).

## Results

Error is quoted **only at $t=8$**. The flow returns the disk to its start, so the
exact solution there is the initial condition; at intermediate times the shape is
a spiral with no closed form. Boundedness, mass and $|\nabla\psi|$ need no
reference and are valid throughout. Every row has `svv_psi` $c_0=0.1$, $N/2$, and no
redistancing. The numbers come from `logs/d5/gamma_sweep.py` (gitignored, local), over every output
frame (spacing 0.08; 0.04 at $H=1/96$ and 1/128), and are collected in `logs/table_lambda_svv_2026-10-05.txt`.

### $\xi$ trades shape against boundedness

$H=1/64$, $N=5$:

| $\xi$ | $\gamma$ | $E_r(t{=}8)$ | worst violation | area($\phi>0.5$)/$A_0$ @ $t{=}4$ | $\phi$ peak @ $t{=}4$ |
|---|---|---|---|---|---|
| **1.0** | 1.0 | **0.0446** | $2.5\times10^{-3}$ | 0.978 | 0.996 |
| 1.5 | 0.5 | 0.0620 | $\sim10^{-30}$ | 0.933 | 0.976 |
| 1.5 | 1.0 | 0.0631 | $\sim10^{-69}$ | 0.931 | 0.977 |
| 1.5 | 2.0 | 0.0615 | **0** | 0.930 | 0.977 |
| 2.0 | 1.0 | 0.0705 | **0** | 0.873 | 0.940 |

The maximum mass drift $|m/m_0-1|$ is at most $4.8\times10^{-6}$ in every row.

$E_r$ rises with $\xi$, and boundedness improves with it. At $\xi=1.5$, $\gamma$ hardly matters
(0.0615–0.0631 over $\gamma=0.5$–2).

**Why.** A filament of thickness $d$ holds, at CDI equilibrium, a core of
$1-e^{-d/2\varepsilon}$. That falls below 0.5 at $d<2\ln2\,\varepsilon$.
- The exact filament at $t=4$ is $7.3\,H/N$ thick at its median (marker-advected interface,
  `logs/d5/PREREGISTERED.txt`, step 4). That is $7.3\varepsilon$ at $\xi=1$ and $3.7\varepsilon$
  at $\xi=2$.
- So a larger $\varepsilon$ holds less of the filament above 0.5 (0.978 → 0.873), and less of
  the shape survives the return.
- A smaller $\varepsilon$ is less well bounded.

This is the opposite of the strain-free conclusion, where larger $\xi$ won on both counts. That
is why a recommendation tuned on `advecting_slab_1d` or `zalesak_disk` does not transfer: both
are isometries, and nothing thins.

### $h$-refinement improves everything at once

At $\xi=1.0$, $\gamma=1.0$, $N=5$:

| mesh | $\varepsilon$ | $E_r(t{=}8)$ | worst violation | area($\phi>0.5$)/$A_0$ @ $t{=}4$ | $\phi$ peak @ $t{=}4$ |
|---|---|---|---|---|---|
| 1/64 | $3.13\times10^{-3}$ | 0.0446 | $2.5\times10^{-3}$ | 0.978 | 0.996 |
| 1/96 | $2.08\times10^{-3}$ | 0.0211 | $3.0\times10^{-3}$ | 0.991 | 1.000 |
| 1/128 | $1.56\times10^{-3}$ | **0.0101** | $2.4\times10^{-3}$ | 0.998 | 1.000 |

- **Accuracy:** $E_r$ falls 4.4× over one doubling of the mesh, an observed rate of about 2.1. At
  fixed $\xi$, refining shrinks $\varepsilon$ and the mesh together, so the filament is
  $2\times$ more $\varepsilon$-widths thick at $H=1/128$.
- **Boundedness does not improve:** the worst violation at $\xi=1$ stays at
  $2$–$3\times10^{-3}$ at every $H$. For exact boundedness, use $\xi\ge1.5$.

### SVV on $\psi$

Three single-variable pairs, the fixed diffusion and $H=1/64$ in all. The SVV-off $E_r$ are
recomputed from `logs/lambda_fix_nosvv_2026-10-02/` with `norms.E_r`; SVV off is not a
configuration of this method (`../../CLAUDE.md`):

| $\xi$, $\gamma$ | $E_r$ SVV off | $E_r$ `svv_psi` $c_0=0.1$ | worst violation (off / on) |
|---|---|---|---|
| 1.0, 1.0 | 0.0498 | **0.0446** | $2.5\times10^{-3}$ / $2.5\times10^{-3}$ |
| 1.5, 1.0 | 0.0708 | **0.0631** | $\sim0$ / $\sim0$ |
| 1.5, 2.0 | 0.0688 | **0.0615** | 0 / 0 |

SVV on the $\psi$ transport is worth 10–11% in $E_r$ here at no cost in boundedness.

**The trap worth knowing:** an explicit SVV instance is initialised lazily by its own
source-term hook. So a `.case` whose `psi` scalar lacks `"source_terms": [{"type": "user"}]`
never applies the term, **while the header still prints it as on**.
- Every `.case` here carries the entry.
- `rider_kothe.f90` stops after step 1 if an explicit SVV term was requested but never fired.
- At startup it stops if `normal = "psi"` has `svv_psi` off.

## Redistancing with Saini's settings

**With both fixes, the periodic reseed completes; it does not beat transport alone in what
matters, the normal.** The validated configuration stays redistancing off.

The redistancing settings of Saini's own circVortex case (`nandu90/nekLS_Examples@jcp`) were
run against transport alone. Both runs use $H=1/64$, $N=5$, $\xi=1$, $\gamma=1$ and `svv_psi`
$c_0=0.1$, $N/2$, the fixed diffusion, and run to $t=8$.

| run | configuration | $E_r(t{=}8)$ | outcome |
|---|---|---|---|
| `rk_Af` | transport only, exact $\psi_0$ (= `rider_kothe_xi10`, bit-identical) | 0.0446 | completes |
| `rk_eLf` | his full configuration: built $\psi_0$; every 0.5 an Eq. (47) reseed and an Eq. (44) solve over $25H$ with $\Delta\tau=H/(N{+}1)$, sign-function $\varepsilon=0.25$, Eq. (31) SVV, BDF2/EXT2, dealiased $\mathbf C(\mathbf w)$, his sign guard | 0.0407 | completes, worst violation $8.9\times10^{-4}$, $\lVert dn\rVert$ 86–159 per event |

`rk_eLf` ran on the scratch user file `logs/d5/rider_kothe_d5.f90` (gitignored, local).
- Of its keys, only the history restart is in `rider_kothe.f90`.
- The sign-function width, BDF2, dealiasing and the guard are scratch-only (decision D2/D5 in
  `../../NEXT_SESSION.md`).
- The committed events path (SSP-RK3, sign-function $\varepsilon$ equal to the phase field's,
  $2.5H$) has not been run here since the fixes.

**How good is the reseeded $\psi$?** All of the following was measured on frames written right
after an event (`logs/d5/rd_quality.py`, `logs/d5/filament_core.py`; tables in `logs/d5/ev/`; gitignored, local):
- **How often.** Every 0.5, so 16 events, as Saini re-distances his TLS. He also re-sharpens his
  phase field every 0.05; this repo has no such step.
- **The solve does not reach a distance where it matters.** In the compression band, right after
  each event at $t=2$–5 (at $t=6$ and 8 the mean $|\nabla\psi|$ is 1.409 and 1.484,
  `logs/d5/ev/rd_quality_eLf_Af.txt`):
  - $|\nabla\psi|$ has mean 1.60–1.67, with 5th and 95th percentiles of 1.1 and 2.2;
  - $|\psi-d_\phi|$, where $d_\phi$ is the distance to $\phi$'s 0.5 contour (the solve's own
    target), has mean $2.6\varepsilon$ and 95th percentile $5.6\varepsilon$.

  The likely reason: with the sign-function width 0.25, the pseudo-velocity
  $|\operatorname{sgn}\psi|\approx2|\psi|$ is below about 0.1 inside the band. So $\psi$ there
  stays close to the steep seed $r_f(\phi-\tfrac12)$, whose gradient is $r_f/4\varepsilon=8$.
- **So the normal is $\phi$'s, and worse than transport's.** The rebuilt $\psi$ points 4–8° off
  the exact interface at $t=2$–5. The transported $\psi$ in `rk_Af` points 0.6–0.8° off.
- **The thin tail loses its contour.**
  - At $t=3.5$–4.5, the exact filament is thinner than $2\varepsilon$ over 5–6% of its length.
    There, even transport-only $\phi$ has a core of only 0.33–0.35, so no 0.5 contour.
  - A reseed then has nothing to build from. In `rk_eLf`, $\psi$'s normal across that tail is
    wrong at 76–79% of points, and at about 20% of those in the 2–4$\varepsilon$ band.
  - $\phi$'s 0.5 contour is in 3–4 pieces at $t=3$–6, against 1–2 for transport alone.
- **The end of the run degrades.** $\phi$'s own normal is 17° off the exact interface at
  $t=7.52$ and 28° at $t=8$; transport alone ends at 2.0°. Saini's TLS shows the same late rise,
  12–17° at $t=7.5$–8.

### Saini's own code on the same problem

His circVortex in his fork (`nandu90/Nek5000@nekLS` 5e9b0ae) reproduces his paper:
- as shipped ($H=1/128$, $N=3$): $E_r$ 1.966e-3, against Table 3's 1.96e-3;
- at our resolution ($H=1/64$, $N=5$, $\Delta t=4\times10^{-4}$): $E_r$ 3.13e-3, against Fig. 17's ≈3.3e-3.

**His $E_r$ is not on our scale.** His `ls_relerr` divides by $\int\psi_e$ with $\psi_e=1$
*outside* the disk (his CLS is 1 outside), about 0.929 against our disk area 0.0707. Recomputed
with our `norms.E_r` on his $t=8$ dumps, he gets **0.0258** ($H=1/128$, $N=3$) and **0.0410**
($H=1/64$, $N=5$).

At $H=1/64$, $N=5$ we get 0.0446 by transport alone and 0.0407 with his events, so the methods
are level.
- **Filament retention matches.** His CLS keeps 0.978 of its $\phi>0.5$ area at $t=4$, as ours
  now does. Its core follows the same equilibrium curve $1-e^{-d/2\varepsilon}$.
- **His low $E_r$ comes from his phase-field re-sharpening.** Without it (`userParam02 = 0`), his
  $E_r(8)$ is 0.208.
- **His CLS is not bounded:** with re-sharpening its range is $[-0.09, 1.15]$; without it,
  $[-0.21, 1.16]$.
- Fed the same $\phi$, our Eq. (44) build and his TLS agree to 0.10° in the band.

All of this is in the gitignored `logs/d5/`:
- runs and `PREREGISTERED.txt`, with every prediction and verdict;
- the scratch user file;
- `saini/`: his fork, his runs and the comparison scripts;
- `ev/`: the measurement tables.

## Configuration

Mesh `genmeshbox 0 1 0 1 -0.1 0 64 64 1 .true. .true. .true.` ($H=1/64$),
$N=5$, so $H/N = 3.125\times10^{-3}$ and $\varepsilon = \xi\,H/N$. Velocity

$$
u = \sin^2(\pi x)\sin(2\pi y)\cos(\pi t/T), \qquad
v = -\sin(2\pi x)\sin^2(\pi y)\cos(\pi t/T), \qquad T = 8
$$

which is $\mathbf{u} = \nabla\times\psi_s\hat{z}$ for
$\psi_s = \frac{1}{\pi}\sin^2(\pi x)\sin^2(\pi y)$ — verified divergence-free to
$4\times10^{-10}$, non-penetrating on all walls, and time-reversible over $[0,T]$
to $6\times10^{-17}$.

**The CDI rate follows the flow.** Both halves of the CDI balance scale with
$\gamma u_{\max}(t)$, $u_{\max}(t)=|\cos(\pi t/T)|$:
- the compression, an explicit source term;
- the diffusion $\varepsilon\gamma u_{\max}$, implicit, in the scalar's Helmholtz solve.

**Neko's solve reads `s_lambda_tot`**, which it copies from `s_lambda` only at initialisation
unless a turbulence model is set, so `material_properties` fills both every step
(`../../CDI_METHOD.md` §4.1d).

$\Delta t$ is set by the explicit compression term's guard,
$\gamma u_{\max}\Delta t/h_{\text{GLL,min}} \le 0.05$ (`../../CDI_METHOD.md` §2).

`.case` files are named for their parameters: `rider_x15_g20` is $\xi=1.5$,
$\gamma=2.0$; `rider_h128` is the $H=1/128$ refinement.

## Evidence

**All files in `evidence/` predate the 2026-10-02 fixes.** They were made on 2026-09-10 with the
frozen diffusion, and show the dissolving filament it produced. They are to be regenerated from
`visualize.ipynb`, in the kthviz style (`../../NEXT_SESSION.md`).
- Without SVV: `rider_kothe_methods.mp4`, `rider_grad_psi.png`, `rider_xi_and_h.png`.
- With SVV on $\psi$ ($H=1/128$, $N=6$): `rider_svv_snapshots.png`, `rider_svv_methods.mp4`,
  `rider_svv_c01_filmstrip.png`.

Animations of the fixed runs are in the gitignored `logs/anim/`, each as MP4 and GIF:
- `rk_frozen_vs_fixed_diffusion`
- `rk_reinit_events_phi_psi`
- `rk_xi_series`
- `rk_h_series`

## Reference implementation (read, don't copy)

`../../../neko-multiphase/examples/saini_benchmarks/rider_kothe_saini/rider_kothe_sdf.f90`
