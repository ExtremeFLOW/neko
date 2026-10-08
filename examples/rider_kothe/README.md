# rider_kothe

Vortex-in-a-box: a disk is stretched into a thin spiral filament by $t=4$, then
the flow reverses and it unwinds back to a disk at $t=8$.

**Status: complete with redistancing off.**
- **The tables** are from the re-runs of 2026-10-05/06, after the velocity fix below.
- **The redistancing runs** are the committed events path (2026-10-05/06) and two sets on scratch
  user files, D5 (2026-10-06/07) and D6 (2026-10-07/08), all on the velocity-fix build, plus
  `rk_eLf` (2026-10-01).

It is the only case with genuine strain, which makes it the only one that can
settle whether redistancing is needed and what $(\xi, \gamma)$ to use.

**Three defects were fixed. Results from before them are not quoted, except `rk_eLf`, which is
labelled where it appears** (`../../CDI_METHOD.md` §4.1d):
- until 2026-10-02 the CDI diffusion was frozen at its $t=0$ value while the compression followed
  $u_{\max}(t)=|\cos(\pi t/8)|$, which dissolved the filament (see "Configuration"). The 2026-10-01
  scratch-file runs already had the fix;
- until 2026-10-02 every redistancing event left the old $\psi$ in the BDF lags;
- until 2026-10-05 the velocity was prescribed after each step, so step 1 ran with $u=0$. Fixing it
  moved every $E_r(8)$ up by 0.5–3.0% (`logs/vel_fix_2026-10-05/PREREGISTERED.txt`). The SVV on/off
  pairs and `rk_eLf` ran before this fix, and each place that quotes them says so.

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
| 2.0 | 0.0494 | 5.48 | 10.7 |
| **4.0** (max stretch) | **0.0394** | **7.69** | **15.2** |
| 6.0 | 0.0480 | 5.42 | 10.0 |
| **8.0** | **0.778** | **1.008** | **1.19** |

The strain integrates back to zero over the cycle, so $|\nabla\psi|$ returns to 1 (band mean
1.008; the SVV keeps it from returning exactly). **Its direction never leaves.** At every frame
measured ($t=1$–8), $\psi$'s normal is within 0.6–1.3° of the exact interface in the compression band
(`logs/vel_fix_2026-10-05/rdq_output_v2.txt`, this run). That includes the
5–6% of the filament that is thinner than $2\varepsilon$ at $t=4$.

**Redistancing gives a worse normal here from $t\approx2$ on, in every configuration tried**
("Redistancing" below). The best of them, Saini's sign width and extent with the printed Eq. (31)
SVV, does give a better $\phi$ shape than transport at $H=1/64$ and $1/128$. But its $\psi$ is not a
distance, and its normal is 0.9–2.1° off at maximum stretch ($t=4$) where transport's is
0.3–0.8°. See
[`../../REDISTANCING.md`](../../REDISTANCING.md).

## Results

Error is quoted **only at $t=8$**. The flow returns the disk to its start, so the
exact solution there is the initial condition; at intermediate times the shape is
a spiral with no closed form. Boundedness, mass and $|\nabla\psi|$ need no
reference and are valid throughout. Every row has `svv_psi` $c_0=0.1$, $N/2$, and no
redistancing. The numbers come from `logs/d5/gamma_sweep.py` (gitignored, local), over every output
frame (spacing 0.08; 0.04 at $H=1/96$ and 1/128), and are collected in `logs/table_velfix_2026-10-06.txt`.
At $H=1/64$ the run ends one step past the reversal, at $t=8.00008$: Neko's accumulated time is
a round-off below 8 after 100000 steps of $8\times10^{-5}$.

### $\xi$ trades shape against boundedness

$H=1/64$, $N=5$:

| $\xi$ | $\gamma$ | $E_r(t{=}8)$ | worst violation | area($\phi>0.5$)/$A_0$ @ $t{=}4$ | $\phi$ peak @ $t{=}4$ |
|---|---|---|---|---|---|
| **1.0** | 1.0 | **0.0452** | $2.6\times10^{-3}$ | 0.978 | 0.996 |
| 1.5 | 0.5 | 0.0632 | $\sim10^{-30}$ | 0.932 | 0.976 |
| 1.5 | 1.0 | 0.0637 | $\sim10^{-70}$ | 0.931 | 0.977 |
| 1.5 | 2.0 | 0.0618 | **0** | 0.930 | 0.977 |
| 2.0 | 1.0 | 0.0710 | **0** | 0.872 | 0.940 |

The maximum mass drift $|m/m_0-1|$ is at most $4.9\times10^{-6}$ in every row.

$E_r$ rises with $\xi$, and boundedness improves with it. At $\xi=1.5$, $\gamma$ hardly matters
(0.0618–0.0637 over $\gamma=0.5$–2).

**Why.** A filament of thickness $d$ holds, at CDI equilibrium, a core of
$1-e^{-d/2\varepsilon}$. That falls below 0.5 at $d<2\ln2\,\varepsilon$.
- The exact filament at $t=4$ is $7.3\,H/N$ thick at its median (marker-advected interface,
  `logs/d5/PREREGISTERED.txt`, step 4). That is $7.3\varepsilon$ at $\xi=1$ and $3.7\varepsilon$
  at $\xi=2$.
- So a larger $\varepsilon$ holds less of the filament above 0.5 (0.978 → 0.872), and less of
  the shape survives the return.
- A smaller $\varepsilon$ is less well bounded.

This is the opposite of the strain-free conclusion, where larger $\xi$ won on both counts. That
is why a recommendation tuned on `advecting_slab_1d` or `zalesak_disk` does not transfer: both
are isometries, and nothing thins.

### $h$-refinement improves everything at once

At $\xi=1.0$, $\gamma=1.0$, $N=5$:

| mesh | $\varepsilon$ | $E_r(t{=}8)$ | worst violation | area($\phi>0.5$)/$A_0$ @ $t{=}4$ | $\phi$ peak @ $t{=}4$ |
|---|---|---|---|---|---|
| 1/64 | $3.13\times10^{-3}$ | 0.0452 | $2.6\times10^{-3}$ | 0.978 | 0.996 |
| 1/96 | $2.08\times10^{-3}$ | 0.0215 | $3.0\times10^{-3}$ | 0.990 | 1.000 |
| 1/128 | $1.56\times10^{-3}$ | **0.0104** | $2.4\times10^{-3}$ | 0.998 | 1.000 |

- **Accuracy:** $E_r$ falls 4.3× over one doubling of the mesh, an observed rate of about 2.1. At
  fixed $\xi$, refining shrinks $\varepsilon$ and the mesh together, so the filament is
  $2\times$ more $\varepsilon$-widths thick at $H=1/128$.
- **Boundedness does not improve:** the worst violation at $\xi=1$ stays at
  $2$–$3\times10^{-3}$ at every $H$. For exact boundedness, use $\xi\ge1.5$.

### SVV on $\psi$

Three single-variable pairs, the fixed diffusion and $H=1/64$ in all, both columns on the build
before the 2026-10-05 velocity fix (hence 0.0446, not 0.0452). The SVV-off $E_r$ are
recomputed from `logs/lambda_fix_nosvv_2026-10-02/` with `norms.E_r`; SVV off is not a
configuration of this method (`../../CLAUDE.md`):

| $\xi$, $\gamma$ | $E_r$ SVV off | $E_r$ `svv_psi` $c_0=0.1$ | worst violation (off / on) |
|---|---|---|---|
| 1.0, 1.0 | 0.0498 | **0.0446** | $2.5\times10^{-3}$ / $2.5\times10^{-3}$ |
| 1.5, 1.0 | 0.0708 | **0.0631** | $\sim0$ / $\sim0$ |
| 1.5, 2.0 | 0.0688 | **0.0615** | 0 / 0 |

SVV on the $\psi$ transport is worth 10–11% in $E_r$ here at no cost in boundedness.

**Saini's local $|\mathbf c|$** (D6, 2026-10-08, velocity-fix build). Our `svv_psi` takes
$|\mathbf c|=u_{\max}$ once, at step 1, so here it never follows $|\cos(\pi t/8)|$. His is the
local $|\mathbf u(\mathbf x,t)|$, left of the assembled operator. On the scratch file that is a
pointwise multiply of the explicit SVV term by $|\mathbf u(\mathbf x,t_n)|/u_{\max}$; it changes
nothing else (`logs/d6_2026-10-07/tables/`):

| `svv_psi` $\lvert\mathbf c\rvert$ | $E_r(t{=}8)$ | $E_s(t{=}8)$ | worst violation | $\psi$ normal vs exact, $t=2$–5 | at $t=8$ |
|---|---|---|---|---|---|
| uniform $u_{\max}$ (`rider_kothe_xi10`) | 0.0452 | 0.0237 | $2.6\times10^{-3}$ | 0.64–0.81° | 1.25° |
| local $\lvert\mathbf u\rvert$ | 0.0470 | 0.0247 | $2.6\times10^{-3}$ | 0.62–0.70° | 0.83° |

The local form is the weaker filter. It costs $\phi$'s shape 4%, leaves the worst violation and
the band $|\nabla\psi|$ (within 0.7% at every frame) unchanged, and makes $\psi$'s normal no worse
(0.2–0.4° better at $t=4$ and $t=7.5$–8). The committed files keep the uniform $u_{\max}$.

**The trap worth knowing:** an explicit SVV instance is initialised lazily by its own
source-term hook. So a `.case` whose `psi` scalar lacks `"source_terms": [{"type": "user"}]`
never applies the term, **while the header still prints it as on**.
- Every `.case` here carries the entry.
- `rider_kothe.f90` stops after step 1 if an explicit SVV term was requested but never fired.
- At startup it stops if `normal = "psi"` has `svv_psi` off.

## Redistancing

**Redistancing is off in every shipped `.case`, and that stays the validated configuration.**
`rider_kothe.f90` carries one events path, and it fails here. Every other configuration below ran
on a scratch user file and is not in the committed code (`../../NEXT_SESSION.md`).

### How an event works

Every $\Delta t_{tls}=0.5$ the `compute` hook, after the scalar step, replaces $\psi$. The same
solve builds $\psi_0$ in the `initialize` hook (`psi_init = "redistance"`, Saini's Algorithm 1
line 2), without step 3:
1. **Reseed**, Eq. (47): $\psi\leftarrow r_f(\phi-\tfrac12)$, $r_f=0.1$.
2. **Relax** by Eq. (44), $\partial_\tau\psi=\operatorname{sgn}(\psi)-\mathbf C(\mathbf w)\psi$,
   with $\mathbf w=\operatorname{sgn}(\psi)\,\nabla\psi/|\nabla\psi|$ re-formed from each level's own
   $\psi$ and $\operatorname{sgn}(\psi)=\tanh(\psi/2\varepsilon_s)$ (Eq. 46). It runs for
   `band`$\times H$ of pseudo-time at $\Delta\tau=0.1\,h_{\text{GLL,min}}$, with the Eq. (44) SVV
   ($c_0=2$, $N_{svv}=N/4$, implicit).
3. **Restart** the scalar's time history (`nadv = ndiff = 0`), so the next step is BDF1/EXT1.

The configurations differ only in step 2, and in one transport setting:

| setting | committed `rider_kothe.f90` | scratch key | Saini's `circVortex` |
|---|---|---|---|
| sign width $\varepsilon_s$, Eq. (46) | $\varepsilon$, the phase field's | `sgn_eps`: 0.25 | 0.25 |
| extent | $2.5H$ (213 steps) | `band` 25: $25H$ (2129 steps) | $25H$ (the paper prints $2.5H$) |
| $\Delta\tau$ | $0.1\,h_{\text{GLL,min}}$ | unchanged | $H/(N{+}1)$ (150 steps) |
| pseudo-time scheme | SSP-RK3, SVV Lie-split after each step | `scheme = "bdf2"`: BDF2/EXT2 (Eqs. 34–35), SVV in the implicit solve, the lag holding $F^{n-1}$ | BDF2/EXT2 |
| Eq. (44) SVV | `svv_step_imp`: $\nu=c_0u_{\max}H/N$ inside the bilinear form, nonzero on the zero set | `svv_form = "eq31"`: $\mathbf D_\mu=\lvert\operatorname{sgn}\psi^n\rvert$ left of the assembled operator (`svv_step_eq31`), zero on the zero set | the Eq. (31) form |
| $\mathbf C(\mathbf w)\psi$ | on the GLL points, equal to $\operatorname{sgn}(\psi)\lvert\nabla\psi\rvert$ (`rd_rhs`) | `dealias`: on $\lfloor3(N{+}1)/2\rfloor$ Gauss points (`adv_dealias_t`) | dealiased |
| sign guard | none | `guard`: after each implicit solve, a node whose $\psi^n$ disagrees in sign with $\phi-\tfrac12$ keeps $\psi^n$ | `constrainTLSR`, the same |
| `svv_psi`'s $\lvert\mathbf c\rvert$ (transport) | $u_{\max}$, taken at step 1 | `svv_psi_local`: $\lvert\mathbf u(\mathbf x,t_n)\rvert$, left of the operator | local $\lvert\mathbf u\rvert$ |
| phase-field re-sharpening, Eq. (39) | none: the CDI compression acts every step | none | every 0.05 |

The keys sit under `case.cdi.redistance`, except `svv_psi_local` under `case.cdi`; `dealias` and
`guard` need `scheme = "bdf2"`. Two scratch files carry them:
- `logs/d5_2026-10-06/rider_kothe_0c.f90` (D5, 2026-10-06/07): `sgn_eps`, `svv_form`, `scheme`;
- `logs/d6_2026-10-07/rider_kothe_d6.f90` (D6, 2026-10-07/08): the same, plus `dealias`, `guard`
  and `svv_psi_local`.

With every key at its default, both reproduce the committed file's events run byte for byte in
frames f00000–f00007 (to $t=0.56$: the build and the first event; against
`logs/events_2026-10-05/output_rk_events`). The D6 file was also checked against the transport-only
baseline in the same frames and, with `bdf2` and `eq31` set, against D5's own BDF2 run
(`bdf2_e2post`, all 14 frames to $t=1.00016$).

### Results

All at $N=5$, $\xi=1$, $\gamma=1$, `svv_psi` $c_0=0.1$, $N/2$. The $H=1/64$ runs use
$\Delta t=8\times10^{-5}$ and end at $t=8.00008$ (16 events; the 16th fires on that last step, so the
$E_r$ frame is post-reseed); the $H=1/128$ runs use $\Delta t=4.1\times10^{-5}$. Every events run
builds $\psi_0$ by Eq. (44); transport only starts from the exact distance. The angle is $\psi$'s
normal against the exact interface on the frames right after the events at $t=2$–5,
$\phi(1-\phi)$-weighted over the compression band. For run 2 and Run A, the two runs measured at
every frame, those frames sit at or near each interval's peak (item 7 below). Cases are in
`logs/`.

| run | settings beyond the committed path | $E_r(8)$ | $E_s(8)$ | worst violation | $\phi$ pieces, $t=3$–6 | angle, $t=2$–5 | at $t=8$ |
|---|---|---|---|---|---|---|---|
| transport only (`rider_kothe_xi10`) | no events, exact $\psi_0$ | 0.0452 | 0.0237 | $2.6\times10^{-3}$ | 1–2 | 0.6–0.8° | 1.3° |
| committed events path (`events_2026-10-05/rk_events`) | — | 0.835 | 0.340 | 0.41 | 9–60 | 12–47° | 86° |
| run 1 (`d5_2026-10-06/rider_bdf2_eq31`) | `bdf2`, `eq31` | 0.0960 | 0.0367 | $6.4\times10^{-3}$ | 2–7 | 5.7–12° | 54° |
| run 2 (`d5_2026-10-06/rider_bdf2_eq31_s025_25H`) | run 1 + width 0.25, $25H$ | **0.0342** | **0.0147** | $6.2\times10^{-4}$ | 2–3 | 1.0–2.6° | 29° |
| Run A (`d6_2026-10-07/runA_s025_25H_dg`) | run 2 + `dealias`, `guard`: Saini's re-distancing settings but his $\Delta\tau$ | 0.0402 | 0.0173 | $7.1\times10^{-4}$ | 2–3 | 4.3–8.2° | 29° |
| `rk_eLf` (`d5/rk_eLf`, 2026-10-01, before the velocity fix) | the same with his $\Delta\tau=H/(N{+}1)$, on its own scratch file | 0.0407 | 0.0179 | $8.9\times10^{-4}$ | 3 | 4.3–8.2° | 29° |
| $H=1/128$, transport only (`rider_h128`) | no events | 0.0104 | 0.0053 | $2.4\times10^{-3}$ | 1 | 0.3° | 0.3° |
| $H=1/128$, run 2's settings (`d5_2026-10-06/rider_h128_bdf2_eq31_s025_25H`) | as run 2 | **0.0085** | **0.0045** | $1.3\times10^{-3}$ | 1 | 0.5–1.1° | 7.3° |

Saini's own code gives 0.0410 at $H=1/64$, $N=5$ in our normalisation (below). The records are
in `logs/d5/PREREGISTERED.txt` (every prediction and verdict), `logs/d5_2026-10-06/tables/` and
`logs/d6_2026-10-07/tables/`, measured by `logs/d5/rd_quality.py` and
`logs/d5_2026-10-06/{full_measures,frag_first,phi_normal_exact}.py`.

**One event**, measured one step after event 2 ($t=1.00016$); $\phi$ is one piece in every row:

| settings | $\psi$ pieces | nodes with $\operatorname{sgn}\psi\ne\operatorname{sgn}(\phi-\tfrac12)$ | furthest $\psi=0$ from $\phi=0.5$ | angle |
|---|---|---|---|---|
| committed (SSP-RK3, `svv_step_imp`, width $\varepsilon$, $2.5H$) | 3 | 228 | $3.2\varepsilon$ | 6.3° |
| SSP-RK3 + `eq31` | 1 | 8 | $0.41\varepsilon$ | 4.3° |
| run 1 (`bdf2`, `eq31`) | 1 | 7 | $0.41\varepsilon$ | 4.3° |
| run 1 + `dealias`, `guard` | 4 | 131 | $5.2\varepsilon$ | 5.9° |
| run 2 (+ width 0.25, $25H$) | 1 | 0 | $0.29\varepsilon$ | 0.74° |
| run 2 + `dealias`, `guard` | 1 | 18 | $0.41\varepsilon$ | 1.30° |

BDF2 against RK3 (row 3) and dealiasing with the guard (rows 4 and 6) change the band means of
$|\nabla\psi|$ and $|\psi-d_\phi|$ by at most 2.5%. `eq31` makes $|\psi-d_\phi|$ 21% smaller (row 2).
Width 0.25 with $25H$ (row 5) raises them to 1.65 and $9.3\times10^{-3}$, against 1.11 and
$1.3\times10^{-3}$: the seed's slow relaxation (item 3 below).

### What each setting does

1. **The committed path's SVV makes the spurious zero sets.**
   - Within an event an SSP-RK3 stage of Eq. (44) cannot change a node's sign unless
     $|\nabla\psi|>1+2\varepsilon/\Delta\tau\approx35$; the log's maximum is 10.8. The reseed's
     zero set is $\phi$'s contour. So `svv_step_imp`, nonzero on the zero set, makes the crossings.
   - $\psi$ breaks first: two extra loops in the thin tail at $t=1.04$, with $\phi$ still one piece.
     $\phi$'s first extra piece forms on one of them at $t=1.12$.
   - Every later reseed turns $\phi$'s pieces into zero sets, and the compression maintains them:
     the Zalesak arm C mechanism (`../../CDI_METHOD.md` §4.1c).
2. **The printed Eq. (31) SVV removes them** (run 1), and BDF2/EXT2 gives the same event as RK3.
   $\psi$ then converges onto $\phi$'s contour: its band-mean distance from it is 0.19–0.22$\varepsilon$
   at $t=1$–5. From
   $t\approx2.5$ $\phi$'s thin tail breaks first, and every reseed copies the breaks.
3. **Saini's width 0.25 with $25H$ gives the best $\phi$ of all runs** (run 2). $E_r$, $E_s$ and the
   worst violation are below transport's at both resolutions, and the $h$-rate is 2.0 against
   transport's 2.1. $|E_v|$ is not below at $H=1/64$ ($3.4\times10^{-6}$ against $1.0\times10^{-8}$),
   but transport's own reaches $1.8\times10^{-6}$ at $t=4$. But:
   - **That $\psi$ is not a distance.** At width 0.25 the pseudo-velocity
     $|\operatorname{sgn}\psi|\approx2|\psi|$ is below about 0.1 in the band. So $\psi$ relaxes only
     logistically from the steep seed, whose gradient is $r_f/4\varepsilon=8$
     (`../redistance_circles/eps_1d/README.md` §4). Its band $|\nabla\psi|$ is 1.66 (2.7 at $H=1/128$),
     and it sits $2.6\varepsilon$ ($7$–$9\varepsilon$) from the distance to $\phi$'s contour. It is a
     smoothed copy of $\phi$'s interface, refreshed every 0.5. Saini's own code gives the same right
     after its $t=0$ build: 1.64, and 2.32 at $H=1/128$, $N=3$.
   - **Its normal follows $\phi$, and $\phi$ drifts.** Against $\phi$'s 0.5 contour it is 0.9–1.5° off
     at $t=2$–5, closer than transport's (1.2–2.7°). But $\phi$'s own normal is 0.5–3.1° off the exact
     interface at $t=2$–6 and 27° at $t=8$, against 0.7–1.0° and 2.0° for transport-only $\phi$
     (`logs/d5_2026-10-06/tables/phin_*.txt`). Once $\psi$ follows $\phi$, nothing pulls $\phi$
     back. Transport's $\psi$ is independent of $\phi$ and keeps correcting it.
   - The width 0.25 alone, at $2.5H$, keeps the zero set but leaves $\psi$ the seed. Its normal is 37°
     off by $t=0.5$, and $\phi$'s interface is 57% thicker before the first event. It was not run
     further.
4. **Saini's dealiasing with his guard is the whole gap between run 2 and his settings, and it
   costs** (Run A, D6).
   - Run A reproduces `rk_eLf`: $E_r$ 0.0402 against 0.0407, $E_s$ within 4%, the normal within
     0.1° at every post-event frame to $t=7$, the same $\phi$ pieces to $t=5$. It does so despite
     `rk_eLf`'s 14× larger $\Delta\tau$, older build and separate file, which together change little
     here.
   - Dealiasing evaluates $\mathbf w$ between the nodes, so it moves nodes across zero. The guard
     then freezes a node that has crossed where it lands; it does not put it back. At width
     $\varepsilon$ that makes three $\psi$ islands, 1–5$\varepsilon$ off $\phi$'s contour, in one
     event. At run 2's settings it makes no island, but 18 nodes stay crossed and the normal is 0.6°
     further from both $\phi$'s contour and the exact interface.
   - Over 16 events $\psi$'s normal is 2.4–4.2× further off at $t=2$–5 and $\phi$'s own 2–2.6×
     (3.0–5.5° against 1.4–2.3° at $t=3$–5). $\psi$ has up to 11 pieces ($t=4$), and $E_r$ and
     $E_s$ are about 18% higher.
   - Run A is still below transport in $E_r$, $E_s$ and the worst violation, and level with Saini's
     code.
   - The two ran as a pair, so which of them carries the cost is not measured.
5. **The thin tail has no contour to rebuild from.** At $t=3.5$–4.5 the exact filament is thinner
   than $2\varepsilon$ over 5–6% of its length. There even transport-only $\phi$ has a core of only
   0.33–0.35, so no 0.5 contour. In `rk_eLf`, $\psi$'s normal across that tail is wrong at 76–79%
   of points (`logs/d5/ev/`, 2026-10-05).
6. **The runs with width 0.25 and $25H$ degrade at the end.** $\psi$'s normal is 12.7–16.5° off at
   $t=7.52$ and 29° at $t=8$ at $H=1/64$, against transport's 1.4° and 1.3°. Saini's TLS shows the
   same rise, 12–17° at $t=7.5$–8.
7. **Each event from $t=1$ to 5 adds to $\psi$'s normal error, and transport removes part of it
   before the next** (every output frame, run 2 and Run A; `logs/d6_2026-10-07/tables/rdq_all_*.txt`,
   `evidence/rider_redistancing_during.png`).
   - The first event lowers the error: 1.06° → 0.45° in run 2, 1.07° → 0.49° in Run A. Before it,
     transport alone had raised it from 0.02° to 1.04° (transport only) and from 0.46° to 1.06° (run 2).
   - The events at $t=1$ to 5 raise it, ×1.3–2.5 in run 2 and ×1.4–6.5 in Run A. By the next event
     the transport has lowered it to 0.23–0.86 of the post-event value, except run 2 at $t=5$ (0.99).
     Run A at $t=2$: 0.66° before, 4.3° after, 0.97° at $t=2.48$. The event at $t=5.5$ is in between
     (×1.08, ×1.33).
   - From $t\approx6$ the events change it by ×0.8–1.05, and it grows between them: the end-of-run
     rise of item 6.
   - Transport alone changes by at most 8% across the same pairs of frames.
   - The troughs rise with $\phi$'s own normal error (item 3: 0.5–2.3° at $t=2$–5), which the reseed
     imports; the frames are 0.08 apart, so the event and the transport after it are not separated.
     Run 2's error still rises for 0.16 after its event at $t=4$.
   - Each event also resets the band $|\nabla\psi|$ to 1.61–1.67. Between events the strain moves
     it, as it moves transport's to 7.7 at $t=4.08$.

### Saini's own code on the same problem

His circVortex in his fork (`nandu90/Nek5000@nekLS` 5e9b0ae, gitignored copy in `logs/d5/saini/`)
reproduces his paper:
- as shipped ($H=1/128$, $N=3$): $E_r$ 1.966e-3, against Table 3's 1.96e-3;
- at our resolution ($H=1/64$, $N=5$, $\Delta t=4\times10^{-4}$): $E_r$ 3.13e-3, against Fig. 17's ≈3.3e-3.

**His $E_r$ is not on our scale.** His `ls_relerr` divides by $\int\psi_e$ with $\psi_e=1$
*outside* the disk (his CLS is 1 outside), about 0.929 against our disk area 0.0707. Recomputed
with our `norms.E_r` on his $t=8$ dumps, he gets **0.0258** ($H=1/128$, $N=3$) and **0.0410**
($H=1/64$, $N=5$), level with transport alone (0.0452) and with Run A (0.0402).
- **Filament retention matches.** His CLS keeps 0.978 of its $\phi>0.5$ area at $t=4$, as ours does.
  Its core follows the same equilibrium curve $1-e^{-d/2\varepsilon}$.
- **His low $E_r$ comes from his phase-field re-sharpening.** Without it (`userParam02 = 0`), his
  $E_r(8)$ is 0.208.
- **His CLS is not bounded:** with re-sharpening its range is $[-0.09, 1.15]$; without it,
  $[-0.21, 1.16]$.
- Fed the same $\phi$, our Eq. (44) build and his TLS agree to 0.10° in the band.
- **His Fig. 16b** (p. 23) is his TLS, our $\psi$, along $y=0.75$ at $t=8$, $H=1/128$, $N=4$–6.
  Its slope near 1 is therefore not to be compared with the 2.32 above, which is right after his
  $t=0$ build; his code's band mean at $t=8$ is 0.94 (`cv128`, $N=3$). A band mean is not the slope
  along one line, and the paper does not say whether the frame is post-event.
  - Its plateau of 0.074 is the seed's $\pm r_f/2$ relaxed at width 0.25 where $|\nabla\psi|=0$:
    $\tfrac12\operatorname{asinh}(\sinh(0.1)\,e^{2\cdot25H})=0.0738$ at $H=1/128$. Our $H=1/128$
    events run has exactly that; the printed $2.5H$ would give 0.052.
  - A plateau survives transport, so it does not show that the frame is post-event.

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

**Time levels.** The scalar step applies the explicit terms (advection, compression) to $s^n$
and extrapolates them to $t_{n+1}$, as in Saini's Eq. (34). So the `preprocess` hook prescribes
$\mathbf u(t_n)$ and $u_{\max}(t_n)$ before each step, including $\mathbf u(0)$ before step 1,
as his circVortex `userchk` does. The implicit diffusion reads $u_{\max}$ at $t_{n+1}$. Until
2026-10-05 the velocity was set in `compute`, after the step, so step 1 ran with $u=0$.

$\Delta t$ is set by the explicit compression term's guard,
$\gamma u_{\max}\Delta t/h_{\text{GLL,min}} \le 0.05$ (`../../CDI_METHOD.md` §2).

`.case` files are named for their parameters: `rider_x15_g20` is $\xi=1.5$,
$\gamma=2.0$; `rider_h128` is the $H=1/128$ refinement.

## Evidence

Made on 2026-10-08 (the events animation on 2026-10-07), in the kthviz style, from this case's
runs: the shipped transport-only cases and the scratch-file runs of "Redistancing". Every animation
also ships as a GIF.

**The validated configuration, transport only** (`visualize.ipynb`, from `rider_kothe_xi10` and
`rider_h128`):
- **`rider_snapshots.png`** shows $\phi$ and $\psi$ at maximum stretch ($t=4$) and at the return
  ($t=8$), $H=1/64$, against the marker-advected exact interface (dashed).
  - $\phi$ values below 0 show in yellow. The specks at $t=4$ are an undershoot of up to
    $9\times10^{-4}$ at that time (the run's worst, $2.6\times10^{-3}$, is at $t=3.2$).
  - The pale streak in $\psi$ at $t=4$ (left edge) is small but negative, at most $-0.023$. It is not
    a zero set: $\psi=0$ is drawn solid and follows the filament.
- **`rider_grad_psi.png`** plots band $|\nabla\psi|$ against $t$ at $H=1/64$ and $1/128$: the mean,
  and its 5–95% range shaded. It rises to about 8 at the reversal and returns to 1.00–1.01.
- **`rider_transport.mp4`** shows the whole run at $H=1/64$ and $1/128$ side by side, every 0.08
  (`logs/anim/anim_rk.py`). Anything below 0 shows yellow, however small; each panel's box gives
  $\phi$'s range. The yellow patch at $H=1/128$ from $t\approx6.7$ on is an undershoot below
  $10^{-6}$ (`logs/d5_2026-10-06/tables/full_h128.txt`).
- **`rider_convergence.png`** plots $E_r(8)$ against $H$, with run 2's redistancing beside it, and
  against $\xi$, with each point's worst violation (`logs/figs/figs_rk.py`, from
  `logs/table_velfix_2026-10-06.txt` and the D5 tables).

**Redistancing** (`logs/figs/figs_rk.py`, from the tables of "Redistancing"):
- **`rider_redistancing_normals.png`** plots $\psi$'s normal against the exact interface on ten frames
  next to events ($t=0.56$, 1.04, 2, 3.04, 4, 5.04, 7.04, 7.52 just after; 6 and 8 just before), and
  $E_r(8)$, for transport only and the four events runs. For run 2 and Run A the post-event frames
  at $t=2$–5 sit at or near each interval's peak (item 7).
- **`rider_redistancing_during.png`** plots, for transport only, run 2 and Run A, $\psi$'s normal at
  every output frame (every 0.08) and the band mean $|\nabla\psi|$ from the run logs, with each
  event's before and after values (`logs/d6_2026-10-07/tables/band_log_*.txt`).
- **`rider_redistancing_contours.png`** draws $\phi=0.5$ against the exact interface for transport
  only, run 2 and Run A: at $t=4$ over the whole domain and at the thin tail, whose tip no run
  reaches, and at $t=8$.
- **`rider_redistancing_dealias.mp4`** shows transport only, run 2 and Run A, so the effect of
  dealiasing with the guard: at $t=4$ Run A's $\psi=0$ wiggles along the outer arm, and at $t=8$ both
  rebuilt zero contours are ragged across the top of the disk (`logs/anim/anim_rk.py`).
- **`rider_redistancing_events.mp4`** (2026-10-07, `logs/anim/anim_rk.py`) shows transport only, the
  committed events path, run 1 and run 2. The top row is $\phi$, the bottom row $\psi$ with its zero
  contour solid. At the tail tip the rebuilt $\psi$'s zero contour stops where $\phi$'s does, while
  transport's follows the exact interface.

## Reference implementation (read, don't copy)

`../../../neko-multiphase/examples/saini_benchmarks/rider_kothe_saini/rider_kothe_sdf.f90`
