# advecting_slab_1d

1D linear advection of a $\tanh$ slab through a periodic domain, twenty
flow-throughs, under Neko's CDI compression term. The case exists to isolate one
change — where the compression term's interface normal comes from — with
everything else held fixed.

**Status: built and validated.** Both $\psi$-normal targets reproduce.
- **SVV on $\psi$:** since 2026-10-02 every $\psi$-normal case carries `svv_psi`. The 17 sweep
  cases use $c_0=0.1$, $N/2$, and the tables below are from those runs. The six SVV and
  redistancing variants (`psi_svvon`, `psi_rd_*`) use $c_0=1$.
- **The events cases** (`psi_rd_on*`) were re-run on 2026-10-05 with the event's history restart.

$\xi = \varepsilon N/H$ throughout — see [`CDI_METHOD.md`](../../CDI_METHOD.md) §1.

Paths under `logs/` are gitignored and local to this workstation.

## How $\psi$ is maintained here: it does not need to be, and provably so

The reference runs seed $\psi$ once from the exact analytic signed distance
$\psi_0 = a - |x - x_c|$ and then advect it, with SVV, for twenty flow-throughs. No
redistancing. That rests on two facts below.

**SVV on $\psi$ is always on** (`case.cdi.svv_psi`; $c_0=0.1$, $N/2$, Saini's TLS-transport value,
in the sweep cases): the user file stops at startup if `normal = "psi"` has it off
(`../../CLAUDE.md`).
**Redistancing** (`case.cdi.redistance`, with its `psi_init` companion) is supported and off by
default. It was ported from `../zalesak_disk/` so that this case, four minutes a run against
hours there, can serve as the cheap regression test for the machinery the other cases depend
on. It is not needed here; see below for what it measures.

1. **The normal ignores $|\nabla\psi|$.** For strictly increasing $f$,
   $\nabla f(\psi)/|\nabla f(\psi)| = \nabla\psi/|\nabla\psi|$, so only $\psi$'s
   *level sets* matter. In 1D this collapses to
   $\mathbf{n} = \operatorname{sign}(\partial_x\psi)$ — magnitude is irrelevant
   entirely.
2. **There is no drift to correct.**
   $D|\nabla\psi|/Dt = -|\nabla\psi|(\mathbf{n}\cdot\nabla\mathbf{u}\cdot\mathbf{n})$,
   and $\mathbf{u}$ is uniform, so $\nabla\mathbf{u} = 0$ and the right-hand side
   is **exactly zero**. Redistancing would have nothing to undo.

### Measured, not just argued

`psi_xi10`, over every frame of the run, band $\phi(1-\phi)>10^{-4}$
(`scripts/slab_table.py`, 2026-10-05):

| quantity, interface band | value |
|---|---|
| $\|\nabla\psi\|$ range | 0.65 – 1.48 |
| $\|\nabla\psi\|$ band *mean* | stays in $[0.994, 1.009]$ |
| $\max\|\psi - \psi_{\text{exact}}\|$ in band / outside | 0.00067 / **0.0056** |
| $\max\|\phi - \phi_{\text{exact}}\|$ | 0.00084 |

$\psi$'s error is **8× larger outside the band than inside**, i.e. concentrated where the flux
$-\phi(1-\phi)\mathbf{n}$ vanishes and cannot reach $\phi$. Without SVV (the archived run of
2026-09) $|\nabla\psi|$ ranged over 0.001–3.6, nearly four decades, at the same $E_r$: the
invariance in (1) is not theoretical.

It also shows why the $|\nabla\psi|$-tolerance criterion that the code used to offer
would have measured nothing: here the band *mean* never leaves $[0.80, 1.25]$, while the
*minimum* would have fired at nearly every check and redistanced a correct solution. See
[`../../REDISTANCING.md`](../../REDISTANCING.md).

The same check, run on `psi_xi10` (then without SVV) on 2026-10-01 with the normal-error tool of the
redistancing study (`../rider_kothe/README.md`), finds no wrong sign at any band node,
with band $\phi(1-\phi)>10^{-4}$, at $t=0,1,\dots,20$. $E_r$ stays at 0.00006 throughout.
Redistancing is not needed here.

### What SVV buys here: nothing in $E_r$, seven orders in boundedness

$\xi=1$, $\gamma=1$, $\psi$-normal, twenty flow-throughs, every frame
(`logs/svv_off_2026-10-02/slab_rows.txt`). The SVV-off run is the archived one of 2026-09; it
is not a configuration of this method:

| `svv_psi` | $E_r$ | worst violation | violating nodes | $\max\lvert\psi-\psi_{ex}\rvert$ in band | band $\lvert\nabla\psi\rvert$ |
|---|---|---|---|---|---|
| off | 0.00006 | $1.5\times10^{-3}$ | 3872 | 0.0031 | 0.001 – 3.6 |
| $c_0=0.1$ (`psi_xi10`) | 0.00006 | $\mathbf{1.2\times10^{-10}}$ | 484 | 0.00067 | 0.65 – 1.48 |
| $c_0=1$ (`psi_svvon`) | 0.00006 | $1.2\times10^{-10}$ | **363** | **0.00027** | **0.975 – 1.02** |

$E_r$ is **identical to five figures**: this case is already exact to $10^{-4}$, and SVV cannot
improve that. What SVV buys is boundedness: seven orders of magnitude. $c_0=1$ also holds
$|\nabla\psi|$ within 2.5% of 1.
- Compare `../zalesak_disk`, where SVV is worth 9–256× in $E_r$ itself. Uniform advection
  generates little grid-scale noise for SVV to damp.
- **At large $\xi$ it costs shape** (see the sweep below): $E_r$ 0.00001 → 0.00005 at
  $\xi=2$, 0.00003 → 0.00042 at $\xi=2.8$. The mechanism is not measured.

### What a built $\psi$ costs here: nothing, if the build is resolved

This is the cheap testbed for the machinery the other cases depend on. All runs
$\xi=1$, $\gamma=1$, SVV on, twenty flow-throughs; the only difference is how
$\psi$ is produced. **None of the "built" rows evaluates an analytic distance
anywhere** — that is the point, since real geometries have none.

Every frame, `svv_psi` $c_0=1$ in these cases. The events rows are from the shipped outputs
re-run on 2026-10-05 with the event's history restart (`logs/svv_off_2026-10-02/slab_rows.txt`):

| $\psi$ at $t=0$ | build $\Delta\tau$ | periodic redist. | $E_r(t{=}20)$ | worst violation |
|---|---|---|---|---|
| exact analytic | — | off | **0.00006** | $1.2\times10^{-10}$ |
| built by Eq. (44) | `cfl = 0.1` | off | **0.00006** | $1.2\times10^{-10}$ |
| built by Eq. (44) | `cfl = 0.1` | on, 40 events | **0.00006** | $1.4\times10^{-10}$ |
| built by Eq. (44) | $H/(N{+}1)$ | off | 0.06458 | $2.6\times10^{-3}$ |
| built by Eq. (44) | $H/(N{+}1)$ | on, 40 events | **1.283** | destroyed |

**Giving up the analytic distance costs nothing.** A $\psi$ built from $\phi$
alone matches the closed-form distance to five significant figures in $E_r$ and
in boundedness. That is the result this repo needs, because the geometries it is
aimed at have no closed form.

**The pseudo-timestep is the only thing that matters, and it must be fine.**
At Saini's $\Delta\tau_{tls} = H/(N{+}1)$ — a pseudo-CFL of 2.755 here — the
*build itself* is already wrong ($|\nabla\psi|$ spans 0.067–3.165 instead of
0.999–1.001), and 40 repeated events then drive it beyond $10^{30}$. At
`cfl = 0.1` the build is exact and each event is nearly a no-op
($\lVert dn\rVert = 9\times10^{-8}$). See `../../CDI_METHOD.md` §4.2.

**Periodic reconstruction is not needed here, and with the event's history restart it costs
nothing either.** Uniform advection is an isometry, so $|\nabla\psi|$ is provably conserved and
there is nothing to correct. Splitting Saini's operation into its two halves, 40 events either
way, everything else fixed:

| | $E_r$ | $E_v$ | worst violation |
|---|---|---|---|
| no redistancing at all | 0.00006 | $1.7\times10^{-12}$ | $1.2\times10^{-10}$ |
| `seed="psi"`: relax in place | **0.00006** | $-8.1\times10^{-10}$ | $1.2\times10^{-10}$ |
| `seed="phi"`: Eq. (47) reseed, then relax (Saini Alg. 1) | **0.00006** | $-3.1\times10^{-7}$ | $1.4\times10^{-10}$ |

On an isometry neither half costs anything. Under strain and with corners they do (Zalesak arm
C and Rider–Kothe, `../../CDI_METHOD.md` §4.1c–d).

## Results

All at $t = 20$ (twenty flow-throughs), against the exact slab evaluated at
**the frame's own time** — Neko writes its last frame one step past `end_time`,
and over a half-interface-width that is worth a factor of ~5 in $E_r$.
"worst violation" is $\max(\phi - 1, -\phi, 0)$ over every frame of the run (spacing 0.05), not
just the last. The $\psi$ rows have `svv_psi` $c_0=0.1$ and come from
`logs/svv_off_2026-10-02/slab_table.txt`. The $\phi$-normal rows are the shipped outputs, measured
the same way on 2026-10-05.

| `.case` | $\xi$ | $\gamma$ | normal | $E_r$ | worst violation |
|---|---|---|---|---|---|
| `gamma0_xi10` | 1.0 | 0 | — | 0.01486 | $1.7\times10^{-2}$ |
| `gamma025_xi10` | 1.0 | 0.25 | $\psi$ | 0.00007 | $6.0\times10^{-6}$ |
| `psi_xi10` | 1.0 | 1 | $\psi$ | **0.00006** | $1.2\times10^{-10}$ |
| `psi_xi15` | 1.5 | 1 | $\psi$ | 0.00001 | **0** (exactly) |
| `psi_xi28` | 2.8 | 1 | $\psi$ | 0.00042 | **0** (exactly) |
| `psi_xi05` | 0.5 | 1 | $\psi$ | **0.00044** | $1.3\times10^{-2}$ |
| `phi_xi10` | 1.0 | 1 | $\phi$ | **≈1.2** | $3.1\times10^{-2}$ |
| `phi_xi05` | 0.5 | 1 | $\phi$ | **diverged** at $t=6.80$ | — |

The $\phi$-normal runs predate the `grad_floor` fix (`CDI_METHOD.md` §4.1) and
were not re-run; read them for direction only.

At both $\xi$ the $\psi$-normal is far better: at $\xi=1.0$ the $\phi$-normal is
destroyed, and at $\xi=0.5$ it diverges.
- Swapping $\nabla\phi$ for $\nabla\psi$ in the normal was the only difference between each
  pair.
- Since 2026-10-02 the $\psi$ cases also carry `svv_psi` and the ψ scalar's user source term.
  The $\phi$-normal runs were not changed.

**Read the $\phi$-normal column as one significant figure.** That failure is an
instability, not a converged wrong answer, so its endpoint amplifies round-off:
changing only the scalar solver tolerance from `1e-10` to `1e-11` moves it
between 1.11 and 1.32 (and moves the *upstream* reference implementation between
1.23 and 1.27), with $\max|\Delta\phi| \approx 0.8$ pointwise. At $\xi=0.5$ it
diverged at every tolerance tried, at $t \approx 5$–$7$. The $\psi$-normal side
is the opposite: before `svv_psi` was added it reproduced upstream's stored output to
$1.8\times10^{-14}$.

Note also that $\phi$-normal at $\xi=1.0$ *ends* inside $[0,1]$ — it violates
transiently and then smears into a bounded blob. Endpoint counts alone mislead,
which is why the table reports the worst violation over the whole run.

## Parameter sweep: where boundedness turns on

$\phi \in [0,1]$ is worth more than a small edge in $E_r$, and both knobs move
it the same way. Neither endpoint alone shows this — it took the sweep.

**$\xi$ at $\gamma=1$** ($\psi$-normal; $\varepsilon = \xi \cdot 0.01$):

| $\xi$ | 0.5 | 0.6 | 0.8 | 1.0 | **1.5** | 2.0 | 2.8 |
|---|---|---|---|---|---|---|---|
| $E_r$ | 0.00044 | 0.00046 | 0.00023 | 0.00006 | **0.00001** | 0.00005 | 0.00042 |
| worst violation | 1.3e-2 | 6.5e-3 | 3.3e-4 | 1.2e-10 | **0** | 0 | 0 |

**Exact boundedness switches on between $\xi=1.0$ and $\xi=1.5$** — not one node
outside $[0,1]$ at any frame from 1.5 upward — and $\xi=1.5$ is also the clear $E_r$ optimum:
40× better than 2.8. With SVV on $\psi$, accuracy turns over sharply above 1.5. Without it,
2.0 and 2.8 were 0.00001 and 0.00003. So $\xi \approx 1.5$ is the sweet spot for this case: the
sharpest interface that is still exactly bounded. Reading only the endpoints ($\xi$ = 0.5 and
2.8) would have missed both the threshold and the turnover.

$\gamma$ multiplies **both** halves of the CDI right-hand side — the diffusivity
is $\varepsilon\gamma u_{\max}$ and the compression term carries
$\gamma u_{\max}$ — so it is a rate, setting how fast the interface relaxes to
equilibrium, not what equilibrium it relaxes to. The trap is that on a test with
**no strain**, turning $\gamma$ down does not improve the method, it removes it:
at $\gamma=0$ there is no compression *and* no balancing diffusion, and the slab
merely advects. That row is the null the rest of the table has to beat, and it
does — by 250× in $E_r$ and $10^8$ in boundedness.

**$\gamma$ at $\xi=1.0$** ($\psi$-normal):

| $\gamma$ | 0 | 0.25 | 0.5 | 1.0 | 1.5 |
|---|---|---|---|---|---|
| $E_r$ | 0.01486 | 0.00007 | 0.00006 | 0.00006 | 0.00006 |
| worst violation | 1.7e-2 | 6.0e-6 | 1.8e-5 | 1.2e-10 | **0** |

$E_r$ is **flat** for $\gamma \ge 0.25$ while boundedness improves by five orders of magnitude,
to exactly bounded at $\gamma=1.5$. So raising $\gamma$ is close to free accuracy-wise and buys
boundedness. The
ceiling is the compression CFL, not accuracy: $C_{\text{comp}} = 0.0303$ at
$\gamma=1$ against the hard 0.05 guard, so $\gamma \le 1.65$ at this $\Delta t$
($\gamma=1.5$ runs at $C_{\text{comp}} = 0.0455$).

**$\gamma$ at $\xi=1.5$**, the boundedness sweet spot. $\Delta t$ headroom is
$0.05/C_{\text{comp}}$ — how much larger the timestep could be before the
compression CFL guard fires:

| $\gamma$ | **0.25** | 0.5 | 1.0 | 1.5 |
|---|---|---|---|---|
| $E_r$ | **0.00001** | 0.00001 | 0.00001 | 0.00001 |
| worst violation | **0** | 0 | 0 | 0 |
| $\Delta t$ headroom | **6.6×** | 3.3× | 1.65× | 1.1× |

**The two knobs trade against each other, and $\xi$ is the better buy.** At $\xi=1.0$ exact
boundedness needs $\gamma=1.5$, which leaves only 1.1× of timestep headroom. At $\xi=1.5$ it
holds at every $\gamma\ge0.25$, with the same $E_r$ and up to **6.6×** the timestep. Without SVV
on $\psi$, $\gamma=0.25$ at $\xi=1.5$ violated by $1.1\times10^{-6}$. So raising $\xi$ lets you
*lower* $\gamma$ and come out strictly ahead on all three of accuracy, boundedness and cost.

At $\xi=1.5$, $E_r$ is flat at 0.00001 for every $\gamma \ge 0.25$: there,
$\gamma$ buys boundedness only, and nothing else.

The mechanism behind both is the compression flux's own $\phi(1-\phi)$ factor,
which vanishes at each bound: it is what actively pushes $\phi$ back inside, so
weakening it ($\gamma \to 0$) removes the enforcement along with the sharpening,
and resolving it better (larger $\xi$) lets it act on more nodes.

### The setting this picks out

$$
\xi \approx 1.5, \qquad \gamma \approx 0.25\text{–}0.5
$$

Exactly bounded, best measured $E_r$, and 3.3–6.6× the timestep headroom. This is
the repo's recommended starting point — the root README's
[Recommended settings](../../README.md#recommended-settings) carries it and states what it does and does not establish beyond this
case (it is measured with no strain, and on a normal that carries only a sign).
`../rider_kothe` runs at it.

**The shipped `.case` set still keeps $\xi$ = 0.5 and 1.0 as the validated
comparison points**, and deliberately: those are what the $\phi$-vs-$\psi$
targets are quoted at, and 0.5 is Saini's sharpest setting, the hard case on
purpose. Changing the primary to the recommended point would make the headline
numbers no longer comparable with upstream's. The sweep points ship alongside
them so the tables above can be regenerated.

## p-refinement isolates the operator effect

At **fixed** $\xi=1.0$ ($\varepsilon = H/N$, so the interface always spans the
same number of points), only the polynomial order changes. And 1D has **no
geometry to under-resolve** — two flat interfaces, no curvature, no corners — so
this isolates the operator:

| $N$ | ψ $E_r$ | ψ worst viol | φ $E_r$ | φ worst viol |
|---|---|---|---|---|
| 6 | 0.00027 | **0** | 0.762 | $8.7\times10^{-3}$ |
| 8 | 0.00017 | **0** | 1.211 | $2.7\times10^{-2}$ |
| 10 | 0.00006 | $1.2\times10^{-10}$ | 1.214 | $3.1\times10^{-2}$ |
| 12 | **0.00001** | $4.0\times10^{-6}$ | 1.189 | $3.4\times10^{-2}$ |

(Every frame, re-measured 2026-10-05: `logs/svv_off_2026-10-02/slab_table.txt` for $\psi$,
`slab_rows_phi.txt` for $\phi$.)

**ψ converges cleanly — 27× from $N=6$ to 12.** φ degrades from $N=6$ and then
**saturates**, destroyed, with boundedness violations growing. (Read the φ column to one significant figure: that endpoint is an
instability, as established above.)

With no geometry in play, the only thing p-refinement changes is the discrete
operator — so the φ-normal's failure to converge is an **operator** property, not
geometric under-resolution. Raising the order gives the grid-scale instability
more modes to grow in; the φ-normal feeds it, the ψ-normal does not.
`../zalesak_disk` reproduces the same contrast in 2D at $N=5\to7$.

## Configuration

Mesh `genmeshbox 0 1 0 0.1 0 0.1 10 1 1 .true. .true. .true.`, $N=10$, so
$H/N = 0.01$, so $\varepsilon = \xi \cdot 0.01$ across
the whole sweep ($\xi$ from 0.5 to 2.8). $\Delta t = 10^{-4}$, `end_time` 20,
uniform $u = 1$ with `case.fluid.freeze = true`.

The `.case` files are named for their $\xi$ or $\gamma$: `psi_xi15` is
$\xi=1.5$, `gamma15_xi10` is $\gamma=1.5$ at $\xi=1.0$.

**CPU only.** Neko's CUDA kernels cap at `lx = 10` (`opr_dudxyz`, `opr_cfl`,
`opr_conv1`), i.e. `polynomial_order` $\le 9$, and this case is settled at
$N=10$. The `.f90` itself is backend-agnostic; it is the polynomial order that
rules the GPU out. `../zalesak_disk` at $N=5$ runs on either.

## Evidence

All produced from this repo's own runs by `visualize.ipynb`, which ships executed so the
numbers are visible without rerunning it. **They predate 2026-10-02**: the $\psi$-normal panels
are the SVV-off runs now archived in `logs/svv_off_2026-10-02/`, and the style is not yet kthviz.
They are to be regenerated (`../../NEXT_SESSION.md`).

- `evidence/slab_1d_methods.mp4` — 402 frames over the full run, one panel per
  method ($\gamma=0$, $\phi$-normal, $\psi$-normal without SVV, $\psi$-normal
  with SVV, all at $\xi=1.0$), with the exact slab dashed behind and $E_r$ and
  the range of $\phi$ annotated per frame.
- `evidence/slab_1d_cdi_off_vs_phi.mp4` — the same animation with only the
  first two panels, $\gamma=0$ (CDI off) and $\phi$-normal.
- `evidence/slab_1d_grad_psi.png` — $|\nabla\psi|$ over the interface band
  against $t$, **for both SVV states**. Without SVV the mean is pinned at 1 while
  the min–max envelope spans 0 to 3.6; with SVV the envelope collapses onto 1.
  The run is correct either way, which is the point: the envelope is not a
  failure signal.
- `evidence/slab_1d_redistancing.png` — the redistancing result in three
  panels: $\psi$ at $t=20$ for the five constructions, $E_r$ against $t$ (the
  fine-build runs sitting on the analytic reference, the coarse-build ones
  departing), and band $|\nabla\psi|$ against $t$, where the coarse-$\Delta\tau$
  run climbs to $10^{30}$ while `cfl = 0.1` sits on 1.
- `evidence/slab_1d_psi_field.mp4` — $\phi$, $\psi$ and $\psi$'s error, band
  shaded, **with SVV off and on overlaid**. At field scale $\psi$ looks clean;
  magnified $\times500$ its error is grid-scale oscillation spiking at its two
  kinks, both *outside* the band. That is how $|\nabla\psi|$ swings decades
  while $\phi$ stays exact — and the SVV-on curve shows what removing that
  oscillation does and does not buy.

## Reference implementation (read, don't copy)

`../../../neko-multiphase/examples/saini_benchmarks/distance_normal/advect_1d_sdf.f90`
and its `.case` files.
