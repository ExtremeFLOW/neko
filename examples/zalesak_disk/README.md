# zalesak_disk

The classic slotted-disk solid-body-rotation benchmark — ten full rotations —
run with Neko's CDI compression term and the $\psi$-normal. It is this repo's
most developed 2D case and the one carrying a deliberate ablation.

**Status: validated.** Every upstream number reproduces.

**On $\xi=2.8$, $\text{svv}_\psi=1.0$ — "the reference configuration".** That
pairing appears throughout this README as a fixed point to measure against,
because it is *the configuration this case reproduces* — the first run anywhere
to complete ten rotations. It is **not** a recommendation, and this repo does not
have a settled recommended configuration: the investigation is still open.
Nothing measured makes $\xi=2.8$ preferable. Over ten rotations $\xi=1.5$
matches it on $E_r$ (0.02205 vs 0.02204) and on exact boundedness, costs the
same, and gives a visibly sharper interface ($\phi$ reaching 1.000 rather than
0.997, an interface 4% of the disk radius rather than 7.5%). Over one rotation at
$N=7$ the $(\xi,\gamma)$ map below puts $\xi=1.5,\ \gamma=0.25$ ahead of it by
$2.1\times$ in $E_r$ (0.00368 vs 0.00769, both from that map). Under strain, lower $\xi$ is better still — see
`../rider_kothe`. Treat $\xi=2.8,\ \gamma=1$ as the **reproduction target and
common comparison point**, nothing more.

$\xi = \varepsilon N/H$ throughout ($H$ element edge, $N$ polynomial order), so
$\varepsilon = \xi H/N$ — canonical definition and the Saini comparison in
[`CDI_METHOD.md`](../../CDI_METHOD.md) §1.

## Results at a glance

Every number below is measured in this repo; the section named has the detail.

| | result | section |
|---|---|---|
| $\psi$-normal vs $\phi$-normal, one key changed | $E_r$ **0.0220 vs 1.0738**; **0 vs ~55 000** violating nodes | [φ-vs-ψ](#the-phi-vs-psi-comparison-measured-here) |
| SVV on the $\psi$ equation | $E_r$ **0.0220 vs 1.3497** at $\xi{=}2.8$; **256×** at $\xi{=}1$, $N{=}7$, and the gap grows with $N$ | [SVV](#svv-on-the-ψ-equation-is-load-bearing) |
| $p$- and $h$-refinement | $\psi$ improves and stays exactly bounded; $\phi$ **degrades** | [Resolution](#resolution-ψ-converges-φ-does-not) |
| $(\xi,\gamma)$ map, $N=7$, one rotation | $\xi\ge1.5$ bounded to round-off at every $\gamma$; boundary sharp between $\xi=1$ and $0.5$ | [the map](#the-xi-gamma-map-svv-on) |
| $\gamma$ | **not** neutral — $E_r$ rises ~2× from $\gamma=0.25$ to 2, confirmed at fixed $\Delta t$ | [γ](#gamma-is-not-flat--and-that-contradicts-what-we-predicted) |
| interface thickness | $\xi=2.8$ spans **2.36 elements** — too thick, 2.1× worse $E_r$ for nothing | [too thick](#the-interface-here-has-been-too-thick-and-the-map-says-so) |
| adding redistancing | hurts badly. Done *coherently* (Saini Alg. 1, arm C, $\xi=1$) with the event's history restart: $N=3$ and 7 diverge, and $N=5$ completes 17× worse than without events, because $\phi$'s contour fragments and every reseed sustains the fragments (`CDI_METHOD.md` §4.1c) | [ablation](#the-redistancing-ablation) |

**There is no recommended configuration yet.** $\xi=2.8,\ \gamma=1$ is the
reproduction target and the common comparison point, nothing more; the
investigation is open.

## Configuration

Mesh `genmeshbox 0 1 0 1 -0.1 0 50 50 1 .true. .true. .true.`, $N=5$, so
$H/N = 0.004$. $\Delta t = 5\times10^{-5}$,
`end_time` 20. Rotation $u = \pi(0.5-y)$, $v = \pi(x-0.5)$ — one turn per
$t=2$, so ten rotations. $C_{\text{comp}} = 0.0473$ against the hard 0.05 guard.

`case.cdi.report_every` (default 1000 steps) logs $\min\phi$, $\max\phi$ and
the mass drift to stdout. It is how boundedness is measured over a *whole* run
rather than at the frames — see the $(\xi,\gamma)$ map above for what that
distinction turned out to be worth.

`zalesak.case` and `zalesak_phi_normal.case` write every $t=0.05$ (401 frames)
because they feed the animation; the ablation variants write every $t=0.2$,
since their deliverable is a number. Frames are ~25 MB, so that is a deliberate
storage trade rather than an oversight.

Runs on either backend. CPU and GPU were checked against each other here and
agree to 1.1e-9 in $\phi$ and 1.3e-13 in $\psi$ after 4000 steps, with $E_r$
identical to six decimals; the GPU is ~8.8× faster (~75 min per run).

## How $\psi$ is maintained here

For the primary and the $\phi$-normal baseline: **it isn't**. $\psi$ is seeded
once from the exact *periodic* signed distance and thereafter only advected,
with SVV on its own equation. No reinitialization, no redistancing.

The two ablation variants are the exception, and they differ in a way worth
keeping straight (`CDI_METHOD.md` §4):

| `.case` | $\psi$ maintenance |
|---|---|
| `zalesak` (primary) | none — transport + SVV only |
| `zalesak_phi_normal` | none (and $\psi$ is not read at all — the normal comes from $\phi$) |
| `zalesak_redistance_phi` | redistancing on, `seed="phi"`: **reinitializes** from $\phi$ ($\psi_0 = r_f(\phi-0.5)$), then relaxes in pseudo-time — Saini Eq. 47 |
| `zalesak_redistance_psi` | redistancing on, `seed="psi"`: relaxes **in place**, no reinitialization |

Redistancing, when it fires, is SSP-RK3 in pseudo-time with
$\Delta\tau = \text{cfl}\cdot h_{\text{GLL,min}} = 2.35\times10^{-4}$ over 213
pseudo-steps, covering a band of $2.5H$ from the interface.

**Transport alone keeps the normals (measured 2026-10-01).** The error is the
$\phi(1-\phi)$-weighted mean angle between $\psi$'s normal and the normal of the run's
own $\phi=0.5$ contour, over the compression band, every 0.2 to $t=20$:
- at $\xi=1$, $N=5$ (`s43_n5_svvon`), it stays between 1.01° and 1.26° over ten rotations;
- at the primary's $\xi=2.8$ it drifts slowly, 1.16° → 2.90°, while $E_r$ stays 0.022 and
  $\phi$ stays bounded.

Saini's own Zalesak case does not redistance at all: the paper calls the re-distancing
equation "inactive for the Zalesak case" (§4.5), and his public case never calls it. The
method and the tool are in `../rider_kothe/README.md`, the one case where redistancing
was then tried with his settings.

## Primary result

$\psi$-normal, $\xi = \varepsilon N/H = 2.8$
($\varepsilon = 0.0112$), $\text{svv}_\psi.c_0 = 1.0$, redistancing **off**,
ten rotations to $t=20$:

| | measured | upstream |
|---|---|---|
| $E_r$ | **0.02204** | 0.022 |
| $\|E_v\|$ | **4.28e-08** | 4.3e-08 |
| $E_s$ | **0.00460** | 0.0046 |
| $\phi$ range | **[0.00000, 0.99668]** | [0.000, 0.997] |
| nodes outside $[0,1]$ | **0, at every one of 101 frames** | zero |

$E_r$ grows slowly and monotonically — 0.0133 after one rotation, 0.0187 at
five, 0.0220 at ten — with mass conserved to $4\times10^{-8}$ and **not one
boundedness violation anywhere in the run**. That last point is the one worth
emphasising: the interface stays inside $[0,1]$ exactly, for ten rotations.

## The $\phi$-vs-$\psi$ comparison, measured here

`zalesak_phi_normal.case` is identical to `zalesak.case` except
`case.cdi.normal`, so this is a genuine single-variable result:

| $t=20$ | $E_r$ | $E_s$ | worst violation | nodes outside $[0,1]$ |
|---|---|---|---|---|
| $\psi$-normal | **0.0220** | **0.0046** | **0** | **0** |
| $\phi$-normal | 1.0738 | 0.2962 | $5.5\times10^{-3}$ | 54 980 |

**49× in $E_r$, and exact boundedness against ~55 000 violating nodes**, from
changing one key.

The number usually quoted for this benchmark, "$E_r$ 3.13 → 0.022", is a larger
factor but a **three-variable** one: it compares upstream's `sdf_phi` against
`sdf_psi_xi28`, which differ in $\varepsilon$ (0.002 vs 0.0112), the normal,
*and* $\text{svv}_\psi$ (0 vs 1.0). That it went from broken to working is a
fair statement about the method as a whole; it is not an isolation of the
normal. The table above is.

## SVV on the ψ equation is load-bearing

The reference configuration sets $\xi=2.8$ **and** $\text{svv}_\psi = 1.0$
together, and
upstream flagged the cell that separates them as never tested. Run here, over ten
rotations:

| $\text{svv}_\psi$ | $E_r$ | $E_s$ | worst violation |
|---|---|---|---|
| 1.0 | **0.0220** | 0.0046 | **0** |
| 0.0 | **1.3497** | 0.3491 | **0** |

**61× in $E_r$ from SVV alone** — without it the ψ-normal is no better than the
φ-normal. But note the last column: at **this** $\xi$, boundedness is exactly zero
either way. The two knobs do different jobs:

- **SVV** buys *accuracy*, by damping grid-scale noise in $\psi$ that would
  otherwise corrupt the normal.
- **Boundedness** comes from the compression flux's own $\phi(1-\phi)$ factor,
  which vanishes at both bounds — which is why $\xi$ and $\gamma$ move it and SVV
  does not.

With upstream's $\xi=0.5$ results this closes the square: **both** $\xi=2.8$ and
SVV are necessary; neither alone suffices.

### …and the benefit grows with polynomial order

The pair above is one point at one resolution. Swept over $N$ at $\xi=1$
(Saini et al.'s own §4.3 setting), $H=1/50$, $\gamma=1$, ten rotations, with an
identical $\Delta t$ in each on/off pair:

| $N$ | $E_r$ SVV on | SVV off | ratio | worst violation on | off |
|---|---|---|---|---|---|
| 3 | 0.1012 | 0.9316 | 9.2× | $2.6\times10^{-8}$ | $1.6\times10^{-3}$ |
| 5 | 0.0208 | 0.9747 | 46.8× | $6.4\times10^{-6}$ | $1.4\times10^{-2}$ |
| 7 | **0.0040** | **1.0266** | **256×** | $9.5\times10^{-6}$ | $3.1\times10^{-2}$ |

**SVV is what makes $p$-refinement work.** With it $E_r$ falls 25× from $N=3$ to
$7$; without it $E_r$ *rises* (0.93 → 0.97 → 1.03), so raising the order makes
the answer worse and the point of a high-order method is lost.

**This also qualifies the "zero either way" line above.** That was measured at
$\xi=2.8$, where the compression flux's $\phi(1-\phi)$ factor carries boundedness
by itself. At $\xi=1$ — a sharper interface — SVV is worth 3–4 orders of
magnitude in boundedness. The claim is $\xi$-dependent, not general.

And **nothing crashes without SVV**: all three off-runs completed ten rotations
and produced plausible-looking output that is up to 256× wrong.

Full argument and the volume-error caveat: [`CDI_METHOD.md`](../../CDI_METHOD.md)
§5. Figures: `evidence/saini_fig{5,6,7}_*.png`.

## Resolution: ψ converges, φ does not

One rotation ($t=2$), $\xi=2.8$, matched settings — each ψ/φ pair differs only in
`case.cdi.normal`:

| normal | $N$ | mesh | dofs | $E_r$ | worst violation |
|---|---|---|---|---|---|
| ψ | 5 | 1/50 | 540 k | 0.01333 | **0** |
| ψ | **7** | 1/50 | 1.28 M | **0.00769** | **0** |
| ψ | 5 | **1/100** | 2.16 M | **0.00705** | **0** |
| φ | 5 | 1/50 | 540 k | 0.29033 | $3.7\times10^{-5}$ |
| φ | **7** | 1/50 | 1.28 M | **0.38371** | $1.1\times10^{-2}$ |

**Refining the polynomial order makes ψ better and φ worse.** ψ improves 1.7×
and stays *exactly* bounded; φ degrades and its boundedness violation grows 290×.
$h$-refinement improves ψ similarly (1.9×), also exactly bounded.

`../advecting_slab_1d` shows the same under p-refinement with **no geometry to
under-resolve** (ψ converges 27× from $N=6$ to 12; φ saturates at
$E_r \approx 1.2$), identifying this as an **operator** effect rather than
geometric under-resolution: raising the order gives the grid-scale instability
more modes to grow in, and the φ-normal feeds it while the ψ-normal does not.

## $\xi=1.5$ is as good as $\xi=2.8$ here

Ten rotations, $N=5$: $E_r$ 0.02205 at $\xi=1.5$ against 0.02204 at $\xi=2.8$,
both exactly bounded. On a rigid rotation the thicker interface buys nothing, so
$\xi \approx 1.5$ is preferable — same accuracy, sharper, cheaper. Under
**strain** the picture changes completely; see `../rider_kothe`.

## The $(\xi, \gamma)$ map, SVV on

The $\xi=1.5$ vs $\xi=2.8$ comparison above is two points. This is the plane:
$\xi \in \{0.5, 1, 1.5, 2, 2.8\}$ against $\gamma \in \{0.25, 0.5, 1, 2\}$, one
rotation ($t=2$) at $N=7$ on the same $50^2$ mesh, $\text{svv}_\psi.c_0 = 1.0$
throughout. Twenty runs, ~9.9 h on the GPU. Every cell is `zalesak_p7.case`
with only $\varepsilon$, $\gamma$, $\Delta t$ and the output cadence changed.

![xi-gamma map, SVV on](evidence/zalesak_xi_gamma_svv_on.png)

**The good region is broad.** All twelve cells with $\xi \ge 1.5$ are bounded to
$\le 4.3\times10^{-10}$ — round-off — at every one of their ~400 logged steps,
across the full $8\times$ range of $\gamma$. Nothing in that block is fragile.

**The boundary is sharp, and it sits between $\xi=1$ and $\xi=0.5$.**

| $\xi$ | worst violation | outcome |
|---|---|---|
| $\ge 1.5$ | $\le 4.3\times10^{-10}$ | exactly bounded, all $\gamma$ |
| 1.0 | $3.9\times10^{-6}$ – $6.1\times10^{-5}$ | small but *real*; four orders worse |
| 0.5 | $2.5\times10^{-2}$, or divergence | fails |

$\xi=0.5$ survives one rotation only at $\gamma \le 0.5$, and even then with
$\phi \in [-0.025, 1.025]$. At $\gamma \ge 1$ the scalar solver diverges outright
($t=0.89$ and $t=0.34$). This locates upstream's "$\xi=0.5$ diverged at $t=1.74$"
as one point on a boundary that $\gamma$ also moves.

### The interface here has been too thick, and the map says so

$\xi$ is an abstraction; this is what it means physically. For
$\phi = \tfrac12(1+\tanh(d/2\varepsilon))$ the 5%–95% band is $5.89\varepsilon$
wide. GLL nodes cluster at element ends, so the **worst**-resolved point is the
element *centre*, where the spacing is $1.465\,H/N$ ($N=7$) — that, not
$h_{\rm GLL,min}$, is what limits how thin $\varepsilon$ can usefully be.

| $\xi$ | 5–95% band | in $H/N$ | in **elements** | nodes across it at the worst point | measured here |
|---|---|---|---|---|---|
| 0.5 | 8.4e-3 | 2.9 | 0.42 | **2.0** | diverges at $\gamma\ge1$ |
| 1.0 | 1.7e-2 | 5.9 | 0.84 | 4.0 | violations $4$–$60\times10^{-6}$ |
| 1.5 | 2.5e-2 | 8.8 | 1.26 | 6.0 | exactly bounded |
| 2.0 | 3.4e-2 | 11.8 | 1.68 | 8.0 | exactly bounded |
| 2.8 | 4.7e-2 | 16.5 | **2.36** | 11.3 | exactly bounded |

**This explains the boundary rather than just locating it.** At $\xi=0.5$ the
interface spans **two nodes** where the grid is coarsest — the tanh profile
cannot be represented there, and the run diverges. At $\xi=1$ it is four:
representable, but violations are small and real. From $\xi=1.5$ up it is six or
more and $\phi$ stays inside $[0,1]$ to round-off.

**And it shows the $\xi=2.8$ reference configuration is genuinely too thick.**
It smears the interface across **more than two whole spectral elements** —
$1.9\times$ Saini et al.'s thickest setting and $5.6\times$ their sharpest (they
study $\xi = 0.5$–$1.5$ in this repo's units; `CDI_METHOD.md` §1). It buys
nothing for it: $\xi=1.5$ is $2.1\times$ better in $E_r$ at identical round-off
boundedness, and $\xi=1$ is $3.2\times$ better still if a $6\times10^{-6}$
violation is acceptable. Before this sweep, Zalesak had only ever been run here
at $\xi=2.8$, so the comparison had never been available.

### $\gamma$ is **not** flat — and that contradicts what we predicted

An earlier plan argued the $\gamma$ axis would come out nearly flat on
$E_r$: rigid rotation is an isometry, $\nabla\mathbf{u}$ is antisymmetric, so
$|\nabla\psi|$ cannot drift and $\gamma$'s re-sharpening job is not exercised.
**The measurement disagrees.** At every $\xi \ge 1.5$, $E_r$ rises monotonically
with $\gamma$, by about $2\times$ across the row:

| $\xi$ | $\gamma$=0.25 | 0.5 | 1 | 2 |
|---|---|---|---|---|
| 1.0 | 0.00267 | **0.00240** | 0.00256 | 0.00308 |
| 1.5 | **0.00368** | 0.00432 | 0.00572 | 0.00719 |
| 2.0 | 0.00472 | 0.00589 | 0.00736 | 0.00829 |
| 2.8 | 0.00542 | 0.00666 | 0.00769 | 0.00806 |

The premise was right and the conclusion did not follow. $\gamma$ is not only a
re-sharpening rate: it multiplies the compression flux **and** the physical
diffusion together — $\lambda = \varepsilon\gamma u_{\max}$
(`zalesak_disk.f90:1238`) against $\gamma u_{\max}\,\nabla\!\cdot\!(\phi(1-\phi)\mathbf{n})$.
The equilibrium tanh width is therefore set by $\varepsilon$ alone and is
$\gamma$-independent, while $\gamma$ sets the *rate* of the whole CDI relaxation
relative to advection. The exact solution of a rigid rotation contains no
relaxation at all, so once the interface has reached its equilibrium profile,
further relaxation is pure perturbation and $E_r$ grows with $\gamma$. That is a
genuine $\gamma$ effect on an isometry — just not the one §5 was reasoning about.

**The obvious objection, and the control that answers it.** $\gamma$ and
$\Delta t$ are coupled by construction in the map above
($\Delta t \propto 1/\gamma$, to hold $C_{\text{comp}}$ fixed), so its $\gamma$
axis also varies $\Delta t$ by $8\times$ and does **not** isolate $\gamma$. So the
$\xi=1.5$ row was re-run at a **fixed** $\Delta t = 1.2990\times10^{-5}$ — the
$\gamma=2$ value, the largest legal for every $\gamma$ in the row. Three new
runs; the $\gamma=2$ point needs none, being the same cell already.

| $\gamma$ | 0.25 | 0.5 | 1.0 | 2.0 | span |
|---|---|---|---|---|---|
| $E_r$, **fixed** $\Delta t$ | **0.00349** | 0.00438 | 0.00572 | 0.00719 | **2.06×** |
| $E_r$, swept $\Delta t$ | 0.00368 | 0.00432 | 0.00572 | 0.00719 | 1.95× |
| $\Delta t$ change | 6.8× smaller | 4× | 2× | — | |
| $E_r$ change | −5.3% | +1.4% | +0.1% | — | |
| worst violation | 2.95e-10 | 2.20e-10 | 1.99e-10 | 2.12e-10 | |

**$\Delta t$ moves $E_r$ by at most 5% across a $6.8\times$ change — and not
monotonically — while $\gamma$ moves it $2.06\times$.** At fixed $\Delta t$ the
$\gamma$ effect is slightly *stronger* than in the swept map, which is the right
direction: removing the $\Delta t$ variation removes a small opposing
contribution. Boundedness is untouched at round-off throughout, so this does not
disturb the boundedness half of the map either. **The $\gamma$ result is real.**

$\gamma$ also moves boundedness, **in opposite directions depending on $\xi$**:
at $\xi \ge 2$ higher $\gamma$ is better bounded ($4.3\times10^{-10} \to
8.8\times10^{-12}$ at $\xi=2$), at $\xi=1$ it is an order of magnitude worse
($3.9\times10^{-6} \to 6.1\times10^{-5}$), and at $\xi=0.5$ it is fatal. So "$\gamma$
matters for boundedness" holds, but not with a single sign.

### What this changes

**$\xi=1.5$, $\gamma=0.25$ is the best exactly-bounded cell on the grid:**
$E_r = 0.00368$ against 0.00769 at the $\xi=2.8,\ \gamma=1$ reference — **$2.1\times$
better shape error at the same round-off boundedness.** The existing "prefer
$\xi \approx 1.5$" recommendation survives and now carries a $\gamma$ to go with
it: prefer the low end of $\gamma$ too, on a strain-free flow.

$\xi=1.0$ gives the best $E_r$ anywhere on the grid (0.00240 at $\gamma=0.5$,
$3.2\times$ better than the $\xi=2.8,\ \gamma=1$ reference) but is **not** exactly bounded
($5.9\times10^{-6}$). That is the same shape-vs-boundedness trade
`../rider_kothe` finds under strain, appearing here without any strain at all,
and it is a concrete instance of the open question of whether a violation that
small is worth a third off the shape error. This repo does not answer it.

### Two smaller results the step-resolved log made visible

**"Exactly bounded" was a frame-resolution claim.** Boundedness here is read
from the solver's own stdout every 200 steps (`case.cdi.report_every`), not from
frames — at 62 MB a frame, writing often enough to catch a transient would cost
hundreds of GB and would still only sample it. Doing so shows a
$7.4\times10^{-10}$ undershoot **at step 1**, decaying to $2\times10^{-22}$ by
step 200 and never recurring. Frames are written every ~770 steps and never see
it. It is numerically negligible — seven orders below the $\phi$-normal's
$5.5\times10^{-3}$ — but at $\xi \ge 2$ it is the *entire* content of the "worst
violation" column, so the honest statement is that those cells are bounded to
round-off with a startup transient, not that no node ever left $[0,1]$.

**The $t=0$ self-check is $\varepsilon$-dependent, not a fixed floor.** It reads
$3.5\times10^{-17}$ at $\xi=0.5$, $3.7\times10^{-14}$ at $\xi=1.5$ and
$2.6\times10^{-9}$ at $\xi=2.8$ — monotone and steep in $\varepsilon$. So it is a
property of representing the tanh initial condition on the GLL grid rather than a
defect in the norm, which is more than an earlier plan could say when it
recorded the number as unexplained. It remains four or more orders below any
measured $E_r$; still not worth chasing.

### $\Delta t$ has two limits when SVV is on

Worth recording because the sweep's first attempt died on it. With
$\text{svv}_\psi$ explicit, Neko's guard requires
$\Delta t\,\rho(B^{-1}S_{vv}) \le 0.4$, and $\rho = 4089.607$ here — fixed across
the grid, since it depends only on $c_0$, $u_{\max}$, the mesh, $N$ and
`nsvv_ratio`. So

$$\Delta t = \min\!\left(\frac{0.9\cdot 0.05\,h_{\text{GLL,min}}}{\gamma u_{\max}},\ \frac{0.9\cdot 0.4}{\rho}\right)$$

and **below $\gamma \approx 0.27$ the SVV limit binds first**, not the compression
CFL. The $\gamma=0.25$ row runs at $\Delta t = 8.80\times10^{-5}$
($C_{\text{comp}} = 0.038$, $\Delta t\rho = 0.36$) rather than the
$1.04\times10^{-4}$ the compression limit alone would permit. Relaxing $\Delta t$
was preferred to setting `implicit = true` on that row: the latter would change
the operator between rows and destroy $\gamma$ as a single-variable axis.

### What Mirjalili, Ivey & Mani (2020) predicts — and where it stops

Worth reading next to this map, because it is the closest thing in the
literature to a theory of it: *A conservative diffuse interface method for
two-phase flows with provable boundedness properties*, JCP 401 (2020) 109006.
Their Eq. (1) is **structurally Neko's CDI equation** — one equation, diffusion
against a $\phi(1-\phi)$ compression flux, no split — with their $\gamma$ a
*velocity*, i.e. $\gamma_{\text{theirs}} = \gamma\,u_{\max}$ here. They prove
that a bounded initial condition stays bounded provided

$$\frac{\varepsilon}{\Delta x} \ \ge\ \frac{1}{2}\!\left(\frac{u}{\gamma_{\text{theirs}}} + 1\right)
\qquad\text{and}\qquad \Delta t \le \frac{\Delta x^2}{2\gamma_{\text{theirs}}\varepsilon},$$

which in this repo's variables is $\xi \ge \tfrac12(1/\gamma + 1)$ — **a
$\gamma$-dependent lower bound on $\xi$, decreasing in $\gamma$.** That is
exactly the shape of the boundary this map measures, which is why it is worth
putting the two side by side.

**It does not transfer, and the numbers say so.** Their proof is for
second-order central differences on a staggered grid with explicit Euler/RK, and
— decisively — with the compression normal taken from $\nabla\phi/|\nabla\phi|$.
This case is spectral ($N=7$), BDF3, SVV on $\psi$, and takes its normal from a
*separately transported* $\psi$, which is the entire point of the repo. There is
not even an unambiguous $\Delta x$ to put in their formula: $H/N$ (the mean node
spacing, which $\xi$ is defined on) and $h_{\rm GLL,min}$ (the minimum, which the
timestep guards use) differ by $2.23\times$ here. The table below uses $H/N$.

| $\gamma$ | 0.25 | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|
| $\xi$ required, $\tfrac12(1/\gamma+1)$ | 2.50 | 1.50 | 1.00 | 0.75 |
| lowest $\xi$ measured clean ($<10^{-9}$) | 1.5 | 1.5 | 1.5 | 1.5 |

**The criterion predicts the threshold sweeps $3.3\times$ across this $\gamma$
range; measured, it sits flat at $\xi \approx 1.5$.** Their second condition is
satisfied everywhere on our grid (tightest cell at 0.56 of the limit), so it is
the $\varepsilon/\Delta x$ condition that is being tested, and it is far more
$\gamma$-sensitive than what we see. The likely reason is the one thing the two
schemes most differ in: their threshold is derived from how
$\nabla\phi/|\nabla\phi|$ behaves at the grid scale, and our normal does not come
from $\phi$ at all.

Its *directional* prediction — more $\gamma$ relaxes the $\xi$ requirement —
splits by regime rather than holding or failing outright. Above the boundary
($\xi \ge 2$) more $\gamma$ **is** better bounded ($4.3\times10^{-10} \to
8.8\times10^{-12}$ at $\xi=2$), as it implies. At $\xi=1$, right at the edge, more
$\gamma$ is an order of magnitude **worse** ($3.9\times10^{-6} \to
6.1\times10^{-5}$), the opposite.

**Read it for the form, not the number.** What it legitimately supplies is the
expectation that boundedness is governed by a threshold in the
$(\varepsilon/\Delta x,\ \gamma/u)$ plane rather than by $\varepsilon$ alone — the
advective term does not scale with $\gamma$ while diffusion and compression both
do, which is why the ratio $u/\gamma$ appears — and their explicit three-way
trade between boundedness, interface resolution and timestep stiffness. It is
not a design rule for this solver, and $\xi \ge \tfrac12(1/\gamma+1)$ should not
be quoted as one.

## The redistancing ablation

Per `CLAUDE.md` this case carries the repo's one sanctioned variant set. The
validated result has redistancing **off**, so the question was whether *adding*
it to an already-working configuration helps, hurts, or does nothing. It had
never been tried at $\xi=2.8$ — at $\xi=0.5$ upstream both seed modes crashed,
but that is a thinner-interface regime.

**Both seed modes hurt badly, and neither crashed — which is worse than
crashing, because a run that finishes invites being believed.**

**But be careful what this is evidence for.** All three runs start from the
*exact global* periodic distance, while redistancing produces a field that is
distance-like only within $2.5H$ of the interface. Those two choices are
incompatible, and each event replaces one with the other discontinuously. Saini
et al. never hit this: their Algorithm 1 line 2 *builds* $\psi$ by solving the
redistancing equation from a $\phi$ seed, even where an analytic distance
exists, so their field is band-consistent from $t=0$. Upstream showed that same
solve works exactly as advertised when used that way. So this measures the
incoherent pairing we ran, not redistancing as such — see
[`../../REDISTANCING.md`](../../REDISTANCING.md) §4.

**That gap has since been closed, and the ablation is doubly confounded.**
`case.cdi.psi_init = "redistance"` now implements Algorithm 1 line 2, and the
coherent version runs. Getting it to run exposed a defect that also affected
these very numbers: `unit_normal` floored $|\nabla\psi|$ at $10^{-30}$, below the
round-off gradient of a flat field, so wherever $\psi$ was flat the compression
term received a *random* unit normal. With `grad_floor = 1e-6` — now the default,
and bit-identical on analytic-$\psi$ runs — the built-$\psi$ configuration is
bounded to $10^{-11}$. **So the table below should not be cited for anything:**
it pairs an analytic initial condition with banded redistancing *and* runs with
the broken floor. See `../../CDI_METHOD.md` §4.1. The two redistancing variants also ran with the
events' broken time history (§4.1d), and have not been re-run since it was fixed. The coherent arm C,
re-run with the fix, is in §4.1c.

| | $E_r$ | $E_s$ | worst violation | events | wall time |
|---|---|---|---|---|---|
| primary (off) | **0.0220** | **0.0046** | **0** | — | 74 min |
| `seed="phi"` (reinit + relax) | 0.7087 | 0.1969 | $9.1\times10^{-2}$ | 122 | 87 min |
| `seed="psi"` (in place) | 1.667 *at $t=3$* | — | — | 452 *by $t=3$* | ~5× projected |

- **`seed="phi"`** ran to $t=20$ and landed **32× worse in $E_r$** than the
  primary — nearly as bad as abandoning the $\psi$-normal altogether (0.71 vs
  1.07). It destroyed boundedness outright: from *exact* to a worst violation of
  $9.1\times10^{-2}$, an order of magnitude worse than even the $\phi$-normal
  failure mode. The mechanism is visible in the event log: the interface band
  grows from 113 000 nodes at the first event to 496 000 at the last — 92% of
  the domain — and $|\nabla\psi|$ *after* each relaxation sits at mean 0.92 with
  max ~10, so the pseudo-time solve is not achieving its own goal. Cost: 18%
  more wall time.
- **`seed="psi"`** was stopped at $t=3$ (16.6%, 62 min) after **452 events** —
  one every ~146 steps, i.e. firing at nearly every check because the relaxation
  never restores $|\nabla\psi|$ to tolerance before the next trigger. At $t=3$
  it was already at $E_r = 1.667$ against the primary's 0.0147 — **113× worse**,
  and 21× worse than `seed="phi"` at the same instant. Projected wall time was
  ~5× the primary. It was stopped deliberately rather than run out; the thrashing
  rate and the $t=3$ error are the finding, and five more hours would only have
  refined a number that was already conclusive.

So the in-place variant, which exists specifically to break the
re-seed → normal → $\phi$ → re-seed feedback loop, is the *worse* of the two
here. Whether that is redistancing itself or a defect in the pseudo-time solve
at these settings remains unestablished — upstream said the same at $\xi=0.5$,
and this repo cannot separate them either.

**The practical conclusion for this repo:** the reference run's redistancing-off
configuration is the right choice *for this case*, and nothing measured here
argues for adding a periodic phase back. This case still cannot settle whether
redistancing is ever *needed* — rigid rotation preserves $|\nabla\psi|$ exactly,
so the drift it exists to counter is absent by construction, and
`../rider_kothe` remains the only case that can decide that. What it now does
settle is narrower and firmer: **Saini's redistancing, as written, cannot be
used with this scheme at all**, needed or not, because the banded $\psi$ it
produces has no gradient to give over most of the domain (§4.1).

## Evidence

`evidence/zalesak_methods.mp4` — the disk over ten rotations, **four panels,
each differing from the reference configuration in exactly one setting**:
$\phi$-normal; the reference run; SVV off; $\xi=1.5$. Annotated per frame with $E_r$ (legitimate
here — rigid rotation has a known exact solution at every instant) and the range
of $\phi$. The colour scale runs past $[0,1]$, so cyan marks $\phi<0$ and green
$\phi>1$.

It makes three things visible at once that the tables state separately: the
$\phi$-normal panel is drowned in cyan (tens of thousands of nodes below zero);
the SVV-off panel has all but dissolved ($\phi$ peaks at 0.36 by $t=19$); and
$\xi=1.5$ is **visibly sharper than $\xi=2.8$ at identical $E_r$** — 0.0220 on screen in
both (0.02205 and 0.02204), but $\phi$ reaches 1.000 rather than 0.997. That is the case for
preferring $\xi \approx 1.5$ here, in one frame.

`evidence/zalesak_xi_gamma_svv_on.png` — the two $(\xi, \gamma)$ heatmaps,
$E_r$ and worst boundedness violation, with diverged cells hatched and labelled
by time of death. Independent colour scales: $E_r$ spans only $3.5\times$ and
reads better linear, boundedness spans nine decades and needs a log scale.

`evidence/saini_fig5_contours.png`, `saini_fig6_slot_zoom.png`,
`saini_fig7_norms.png` — Saini et al. §4.3 recreated at $\xi=1$,
$N\in\{3,5,7\}$, with Fig. 7 carrying the SVV-off curve family their paper does
not show (`../../CDI_METHOD.md` §5).

All five are generated by `visualize.ipynb`, which ships executed. Material
carried over from the upstream investigation has been removed now that this repo
produces its own.

## Reference implementation (read, don't copy)

`../../../neko-multiphase/examples/saini_benchmarks/zalesak_saini/zalesak_sdf.f90`
and its `sdf_psi_xi28.case`.
