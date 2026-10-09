# Maintaining $\psi$: re-initialization and Eq. (44)

Where $\psi$ is maintained, it is maintained as in Saini's Algorithm 1. Every $\Delta t_{tls}$ the
transported $\psi$ is discarded, reseeded from $\phi$ by Eq. (47), and relaxed by Eq. (44). The
validated results transport $\psi$ from the exact distance and need none of this
([`CDI_METHOD.md`](CDI_METHOD.md) §4). A $\psi$ *built* by Eq. (44) is what carries over, because
real geometries have no analytic distance.

This file covers the terms (§1), when $\psi$ needs maintaining (§2), what is known (§3), the rules
(§4) and how the code solves Eq. (44) (§5). Saini's configuration against ours, setting by setting,
is `CDI_METHOD.md` §6.

## 1. Terms

**Saini's usage.**
- Eq. (44) is "the traditional level-set (TLS) re-distancing equation": the word names the PDE.
- The operation is **re-initialization** (Algorithm 1, lines 19–20): reseed from the phase field by
  Eq. (47), *then* relax by Eq. (44). Line 2 builds $\psi$ at $t=0$ the same way.
- There is no in-place mode. The transported $\psi$ is discarded at every event, so it only has to
  survive one $\Delta t_{tls}$.
- The reseed is **registration**, "ensuring congruence in the interface location, as given by the
  zero and 0.5 isocontour of the respective fields".

**This repo.**
- `seed = "phi"` is Saini's re-initialization.
- `seed = "psi"` relaxes the transported $\psi$ in place, with no reseed. It is our own diagnostic,
  not a method: Eq. (44) pins $\psi$'s own zero set, so it cannot re-register $\psi$ to $\phi$.
- "Redistancing" on its own means the Eq. (44) relaxation.

Say whose sense you mean when it matters.

## 2. When $\psi$ needs maintaining

**The normal needs a direction, not a distance.** For any increasing $f$,
$\nabla f(\psi)/|\nabla f(\psi)|=\nabla\psi/|\nabla\psi|$. So $|\nabla\psi|=1$ is not a
requirement. What is required is weaker:
- the zero set sits on the interface;
- $\nabla\psi\neq0$ near it.

**Transport keeps the zero set exact in exact arithmetic.** Level sets of
$\partial_t\psi+\mathbf u\cdot\nabla\psi=0$ are material, so maintaining $\psi$ is a numerical
correction only.

**Only normal strain changes $|\nabla\psi|$.** Along a particle path
$D|\nabla\psi|/Dt=-|\nabla\psi|\,(\mathbf n\cdot\nabla\mathbf u\cdot\mathbf n)$. Rigid rotation and
uniform advection preserve $|\nabla\psi|$ exactly. So **the slab and Zalesak cannot test whether
maintenance is needed.** Rider–Kothe can, with normal strain rates up to $\pm2.5$. Transported
alone there, $\psi$ survives:
- its band-mean $|\nabla\psi|$ peaks at 7.7 (max 15) at maximum stretch and returns to 1.008 at
  $t=8$;
- its normal stays 0.62–1.45° off the exact interface from $t=0.24$ on, and less before
  (`examples/rider_kothe/README.md` §3).

**Errors sit where they cannot act.** The compression flux vanishes outside the interface band. In
1D, $\max|\psi-\psi_{\text{exact}}|$ is 0.00067 inside the band and 0.0056 outside
(`examples/advecting_slab_1d/README.md` §3.4).

**So there is no $|\nabla\psi|$ trigger.** The normal does not depend on $|\nabla\psi|$; the only
cadence is a timer, as in Saini's code.

## 3. What is known

**§4.4, Eq. (44) alone.** Saini's configuration reproduces his values
([`examples/redistance_circles/README.md`](examples/redistance_circles/README.md)).

**Rider–Kothe, Saini's re-initialization settings** (`examples/rider_kothe/README.md` §3).
- **The settings.** Sign width 0.25, $25H$, BDF2/EXT2, the Eq. (31) SVV, dealiasing and the sign
  guard. They ran on scratch user files and are not in `rider_kothe.f90`.
- **The result.** $E_r(8)$ is 0.0402 with our $\Delta\tau$. Transport only gives 0.0452, and his own
  code 0.0410 in our normalisation. A run with his $\Delta\tau$, on an older build, gave 0.0407.
- **Not yet his:** local $|\mathbf c|$ on `svv_psi`, BDF2 in physical time, and any run at
  $H=1/128$.
- **$\psi$ is not a distance in the band.** Its band-mean $|\nabla\psi|$ is 1.65–1.66 after each event
  (1.61 after the last), against 1.64 in his own code after its build. At width 0.25 the pseudo-velocity
  $|\operatorname{sgn}\psi|$ is below about 0.1 in the band, so $\psi$ there relaxes only slowly from
  the steep seed.
- **Its normal is 4–8° off** the exact interface at $t=2$–5, against 0.6–0.8° for transport. Two
  causes:
  1. **The reseed takes $\phi$'s 0.5 contour as the interface.** Once $\psi$ follows $\phi$, nothing
     corrects $\phi$'s own drift. $\phi$'s normal is 3.0–5.5° off at $t=3$–5, against 0.9–1.0°
     under transport.
  2. **Dealiasing moves nodes across zero, and the guard freezes them where they land.** Without the
     pair, the same run's normal is 1.0–2.6° off and $E_r(8)$ is 0.0342, so the pair costs about
     18%. Which of the two carries the cost is not measured.
- **The thin tail.** Where the filament is thinner than $2\varepsilon$ (5–6% of its length at
  maximum stretch), $\phi$ has no 0.5 contour, and no reseed can rebuild $\psi$ there.
- **The end of the run.** $\psi$'s normal rises to 16.1–16.5° at $t=7.52$ and 29° at $t=8$, against
  1.4° and 1.3° for transport. Saini's own code shows the same rise (12–17° at $t=7.5$–8).

**The two knobs to check when the rebuilt $\psi$ is not good enough.**
- **Pseudo-time extent matters.**
  - At width 0.25 the printed $2.5H$ leaves $\psi$ as the seed: its normal is 37° off by $t=0.5$,
    and $\phi$ is 57% thicker before the first event.
  - His code's $25H$ is what the results above use.
- **The step size barely matters.**
  - Over $25H$, 2129 steps of $0.1\,h_{\text{GLL,min}}$ and 150 of $H/(N{+}1)$ give 0.0402 and 0.0407.
    The second run also differs in build and file.
  - $H/(N{+}1)$ is safe only with width 0.25, where the pseudo-velocity near the interface is small.
    At width $\varepsilon$ it is a pseudo-CFL of 2.75 on the slab, and 40 events destroy the run.
- **Re-initialization frequency has never been varied.** It is $\Delta t_{tls}=0.5$ in every run.
  - At $t=1$–5 each event raises $\psi$'s normal error (×1.4–6.5 with dealiasing and the guard), and
    transport brings it back to 0.23–0.86 of that before the next event.
  - The first event lowers the error. From $t\approx6$ it grows between events.

**Zalesak.** Saini does not re-initialize on §4.3. Our only events runs on it used the committed
path.

**The committed coupled path is not Saini's, and it fails.**
- It uses SSP-RK3, the uniform-$|\mathbf c|$ SVV step, sign width $\varepsilon$, $2.5H$, and no
  dealiasing or guard.
- On Rider–Kothe, $E_r(8)$ is 0.835. On Zalesak (arm C), the run diverges at $N=3$ and 7 and is 17×
  worse than without events at $N=5$.
- The mechanism: its SVV step acts on the zero set and makes spurious zero sets, and every reseed
  then copies $\phi$'s fragments.
- Each case README lists these runs in its appendix.

**A reseed takes $\phi$'s interface as the truth.** It is safe only while $\phi$'s interface is
better than the transported $\psi$'s.

## 4. Rules

- **$\mathbf w=\operatorname{sgn}(\psi)\nabla\psi/|\nabla\psi|$ is a function of $\psi$, not a
  velocity.**
  - Re-form it from each level's own $\psi$ and apply it only to that $\psi$.
  - BDF/EXT lags hold the products $\operatorname{sgn}(\psi^m)-\mathbf C(\mathbf w^m)\psi^m$, never
    $\mathbf w$.
  - A frozen, lagged or extrapolated $\mathbf w$ leaves growth at up to $1/2\varepsilon$ on the zero
    set. Frozen diverges, and a lag of 0.05 in $\tau$ is 6× worse.
  - So Eq. (44) never goes through Neko's scalar solver, whose velocity is an external field.
- **No analytic seed.**
  - Never pair `psi_init = "exact"` with events.
  - With `"redistance"` the analytic distance is never evaluated.
- **The extent is a pseudo-time budget, not a mask.**
  - An event runs `ceiling(band·H/Δτ)` steps. Information travels at speed $\le1$, so $\psi$ is
    rebuilt out to $\text{band}\cdot H$, and the seed stays beyond it.
  - The compression reads $\mathbf n$ out to $9.21\varepsilon$, so at band 2.5 the built field covers
    it only if $N\ge3.68\xi$. That is a quality bound on the build (band mean 1.259 → 1.027 at $N=3$
    with a wider band), not a stability fix.
- **There is no convergence test**, only a fixed step count.
  - Judge a build by the residual on the compression band ($\phi(1-\phi)>10^{-4}$), not the build
    band.
  - The band *minimum* diagnoses nothing: it is $7.5\times10^{-14}$ for the exact $\psi$ in a clean
    run.
- **The integrator is not a suspect.** BDF2, BDF3 and SSP-RK3 agree within 2.2% on §4.4 (measured at
  sign width $H/N$). The one exception: at $N_{svv}=N/4$, $H=1/5$, $N=7$, only RK3's split damping
  held the apex.
- **Walls.** From a flat seed, the reseed front can capture a natural-BC wall at $N\ge7$ (1D replica,
  `examples/redistance_circles/archive/README_process_2026-09.md` §10). The coupled cases are
  periodic; a walled case needs that check.

## 5. How the code solves Eq. (44)

### 5.1 Two clocks

Neko's scalar framework never solves Eq. (44). Its advecting velocity is the fluid's, and here
$\mathbf w$ depends on $\psi$. The user files run their own pseudo-time loop.

| | physical time $t$ | pseudo-time $\tau$ |
|---|---|---|
| equation | Eq. (43): $\psi$ carried by the flow | Eq. (44): $\psi$ relaxed, zero set held |
| solved by | Neko's `scalar_pnpn` for the advection; the SVV is ours, a user source term | `redistance` (coupled cases) or `redistance_standalone` (`redistance_circles`) |
| scheme | BDF$k$ with Neko's extrapolation, $k$ = `time_order` (3 everywhere); dealiased advection | coupled: SSP-RK3, then one implicit SVV step. Circles: BDF2/EXT2 with the Eq. (31) SVV, then the guard |
| step | `case.time.timestep` | `redistance.dtau` if set, else `cfl`$\times h_{\text{GLL,min}}$ |
| length | to `end_time` | a fixed number of steps |

**Order of calls in a coupled case** (Neko's `simulation.f90` calls `initialize` once before the
loop and `compute` after each step):

```
user%initialize    -> redistance(..., "psi_init build")     ! when psi_init = "redistance"
do physical steps:
   user%preprocess -> prescribe the velocity at t_n (Zalesak, Rider–Kothe)
   scalars%step    -> Neko: phi and psi transport (svv_psi enters as a source term)
   user%compute    -> if an event is due: redistance(...), then restart the scalar history
   output
end do
```

In `redistance_circles` the whole solve runs in `initialize`, and the scalar `psi` is only a
container that Neko writes out.

### 5.2 Settings

`case.cdi` in the coupled files:

| key | default | meaning |
|---|---|---|
| `psi_init` | `"exact"` | or `"redistance"`: build $\psi_0$ by Eq. (44) (CDI_METHOD §3) |
| `grad_floor` | $10^{-6}$ | floor on $\lvert\nabla\psi\rvert$ in `unit_normal` |
| `redistance.enabled` | false | periodic events |
| `redistance.dt_tls` | 0.5 | event interval; the first event is one interval in |
| `redistance.seed` | `"phi"` | `"phi"`: Eq. (47) reseed, then relax; `"psi"`: relax in place |
| `redistance.r_f` | 0.1 | Eq. (47) |
| `redistance.band` | 2.5 | pseudo-time extent, in units of $H$ |
| `redistance.cfl` | 0.1 | $\Delta\tau=$ `cfl`$\times h_{\text{GLL,min}}$, in $(0,1]$ |
| `redistance.dtau` | unset | an explicit $\Delta\tau$, overriding `cfl` (e.g. Saini's $H/(N{+}1)$) |
| `redistance.svv.{c0, nsvv_ratio}` | 0, 4 | the Eq. (44) SVV, always implicit. **$c_0=0$ is off**, so a case that enables events must set it (the Zalesak variants set 2.0) |

The sign width of Eq. (46) is `case.cdi.epsilon`, the phase-field width. `redistance_circles` reads
`case.cdi.epsilon` as the sign width (0.25), plus `redistance.dtau` (required),
`redistance.tau_end` (6) and `ic` (`"skewed"` or `"exact"`).

### 5.3 The right-hand side: `rd_sgn`, `rd_rhs`

**`rd_sgn`** is $\operatorname{sgn}(\psi)=\tanh(\psi/2\varepsilon_s)$. It runs pointwise on the
host, with no device kernel. $|\operatorname{sgn}\psi|$ is also the characteristic speed
$|\mathbf w|$, and therefore the weight $\mathbf D_\mu$ of the printed Eq. (31).

**`rd_rhs`** (the coupled path) uses the identity
$\mathbf w\cdot\nabla\psi=\operatorname{sgn}(\psi)|\nabla\psi|$ and evaluates
$\operatorname{sgn}(\psi)(1-|\nabla\psi|)$.
- **The order.** It averages the element gradients at shared nodes (gather–scatter, then
  `coef%mult`) and *then* takes the length.
- **Equivalence.** With one averaged gradient in both places, it equals the non-dealiased GLL
  Galerkin $\mathbf B^{-1}\mathbf C(\mathbf w)\psi$ to $10^{-14}$.
- **What it is.** The equation is Hamilton–Jacobi. Its characteristics move with $\mathbf w$, away
  from the zero set, where $\operatorname{sgn}(0)=0$ holds $\psi$ fixed.
- **Dealiasing breaks the identity.** $\mathbf n$ is then dotted with a gradient evaluated elsewhere,
  so zero-set nodes can move. `redistance_circles` does not use `rd_rhs`.

### 5.4 The pseudo-time loops

**`redistance_circles`: Saini's BDF2/EXT2** (Eqs. 34–35), with BDF1 on step 1.
- The dealiased $\mathbf C(\mathbf w)\psi$ is Neko's `adv_dealias_t` on $\lfloor3(N+1)/2\rfloor$
  points.
- The implicit Eq. (31) SVV comes next, then the guard.
- `rd_niter = tau_end/dtau`.
- The step is written out in
  [`examples/redistance_circles/README.md`](examples/redistance_circles/README.md) §3, with a table
  mapping each printed equation to its routine. This is the loop the coupled files are to adopt
  (`NEXT_SESSION.md`).

**Coupled cases: SSP-RK3, then one backward-Euler SVV step** (a Lie split).
- Each of the three stages re-evaluates `rd_rhs` from the latest $\psi$, so $\mathbf w$ is never
  stale.
- Explicit Euler would not do: the linearised operator is advection, whose collocation spectrum is
  imaginary. SSP-RK3's stability region covers an interval of the imaginary axis.
- If `seed = "phi"`, the loop starts from the reseed $\psi\leftarrow r_f(\phi-\tfrac12)$.

**Saini's $\Delta\tau=H/(N{+}1)$.** It is a pseudo-CFL of 0.9–2.75 at $|\mathbf w|=1$. At width
0.25 the pseudo-velocity is below about 0.1 near the interface, and on Rider–Kothe his
configuration ran 16 events cleanly at this $\Delta\tau$. His code does not subcycle: it caps each
pseudo-step at a Nek CFL of 1 (`lvlSet.f`, `if cfl > 1: dt = dt/cfl`), so 150 steps of $H/(N{+}1)$
are an upper bound on his pseudo-time. His logs print the CFL only at step 1 (0.20), so whether the
cap binds later is not recorded.

### 5.5 The SVV steps

The shared operator is $\mathbf S=\tilde D^T\nu G\tilde D$ with $\tilde D=FD$ and
$F=V\operatorname{diag}((k/N)^{N_{svv}})V^{-1}$ (Saini's Eqs. 24, 27, 28).
- **What it sees.** Only the highest polynomial modes in each element. It is zero on linear
  functions.
- **$\nu$.** `svv_local` builds it in that order, with $\nu=c_0|\mathbf c|H/N$ constant inside the
  bilinear form. $|\mathbf c|$ is the module's `u_max`.
- **The two implicit steps.** Both are solved by CG to $10^{-12}$:
  - **`svv_step_imp`** (shared; the coupled cases):
    $(\mathbf B+\Delta\tau\,\mathbf S)\psi^{new}=\mathbf B\psi^{old}$. It conserves mass and acts on
    the zero set, so in Eq. (44) it can move the interface.
  - **`svv_step_eq31`** (`redistance_circles`): the printed Eq. (31), with
    $\mathbf D_\mu=\operatorname{diag}|\operatorname{sgn}\psi^n|$ left of the assembled operator,
    $(\mathbf B+\tfrac{\Delta\tau}{b_0}\mathbf D_\mu\mathbf S)\psi^{n+1}=\mathbf B\hat\psi$.
    - It is zero on the interface and does not conserve mass, which a distance does not need.
    - The operator is not symmetric. It is solved for $z$ with $\psi^{n+1}=\hat\psi+\mathbf D_\mu z$,
      which gives a symmetric system.

### 5.6 `unit_normal` and the reseed

**`unit_normal`** builds the compression term's $\mathbf n$ from the same averaged gradient as
`rd_rhs`, divided by $\max(|\nabla\psi|,$ `grad_floor`$)$ (CDI_METHOD §3).

**The reseed** (`seed = "phi"`) is $\psi\leftarrow r_f(\phi-\tfrac12)$ before the loop.
`redistance_circles` has no phase field and starts from Eq. (83)'s $\psi_0$.

### 5.7 A weakness: kinks on element faces

At a node shared by elements whose slopes are $-1$ and $+1$ (a V: a cone apex, or the crease between
two circles), the averaged gradient is 0. `rd_rhs` then reads "too shallow" and deepens the V at
unit rate; only the SVV holds it back.
- **The exact solution is not a discrete steady state there.** On Eq. (83)'s $\psi_e$ ($H=1/10$,
  $N=3$) the averaged $|\nabla\psi|$ reads 0.40 on the valley $y=0$, 0.87 on the ridge $x=0$, and
  0.15 at the apexes.
- **Length per element, then average, is no better.** It overshoots at vertex apexes (1.32), and by
  the triangle inequality it reads "too steep" wherever the element gradients merely disagree.
- **No fix yet.** No form has been found that is right at a face V, at a vertex apex, and unbiased
  elsewhere. `rd_rhs` is one of the 17 shared routines and is unchanged.

### 5.8 What was checked

**A numpy replica** of `redistance_circles` reproduced the Fortran node by node, to $10^{-13}$, from
the same $\psi_0$. Against it:
- the shortcut equals the explicit $\operatorname{sgn}-\mathbf w\cdot\nabla\psi$ (1e-15) and the GLL
  Galerkin form (1e-14);
- a finite-difference Jacobian confirms that a $\psi$-dependent $\mathbf w$ adds no hidden term;
- the SVV operator has Eq. (24)'s transfer function, is symmetric and positive semi-definite, and
  gives $\mathbf S\mathbf 1=0$;
- SSP-RK3 is third order.

**An independent audit (2026-09-28) and a rebuild from the paper with Neko's operators** found no
implementation error (`examples/redistance_circles/archive/README_process_2026-09.md` §6). It
recorded two facts:
- $\mathbf D_\mu(\psi^n)$ makes the step first order in $\Delta\tau$, though extrapolating it moves
  $E_r(6)$ by at most 0.31%;
- Neko's `time_order` 2 extrapolates with a modified EXT3, not EXT2.
