> **Carried over verbatim from `neko-multiphase/references/saini_2026_test_cases.md`
> as of the session that created this repo (see `../CDI_METHOD.md`).** This repo
> does not maintain its own copy of the investigation — it was written from,
> and for, the sibling `neko-multiphase` repo, so every relative path below
> (`../examples/...`, `../CDI_IN_SEM.md`) and every "we already have this case" /
> "our cases run at..." statement refers to **that** repo's structure and
> findings, not this one's. This repo's own case set and status live in
> `../CDI_METHOD.md` and each case's own `README.md`. If the notes below and
> `neko-multiphase`'s current copy ever disagree, that repo is the source of
> truth — update or re-copy from there rather than editing this one in place.

# Saini et al. (2026) — test cases, setup, and purpose

N. Saini, A. Tomboulides, D.R. Shaver et al., *"A high order continuous
Galerkin spectrally stabilized level-set approach for incompressible two-phase
flows"*, **J. Comput. Phys. 561 (2026) 114961**, doi:10.1016/j.jcp.2026.114961.
Implemented in Nek5000. PDF: `saini_tomboulides_2026.pdf` (untracked, see
`.gitignore`).

These are reading notes on the paper's **numerical experiments**, written to
answer one question: *which of their tests would expose the grid-scale
dispersion problem we measured in `interface_relaxation_multiphase/`, better
than interface relaxation does?* Notes on the method itself are kept to the
minimum needed to read the test setups.

---

## 1. Why this paper is relevant to us

Same discretisation family (high-order **continuous** Galerkin SEM, GLL nodal
basis, tensor-product elements), same interface model family (conservative
level set with a compression/diffusion balance), and they hit the same class of
problem we did — then fix it with a spectral high-pass operator.

### 1.0 Symbols: theirs are transposed against ours

**Saini call the CLS phase field `psi` and the signed distance `phi`. This
project uses the opposite: `phi` is the phase field, `psi` is the signed
distance.** Their equation numbers therefore read transposed against our code,
and this is the single easiest way to misread the paper. Everything quoted below
is rewritten in *our* symbols unless it is a direct quotation, which is marked as
such. Jain's (2022) algebraic surrogate for the distance is called `jain` here
and never `psi`, since it is a pointwise transform of `phi` rather than a
transported field.

### 1.1 Their equations vs our CDI

Their conservative level-set (CLS) field uses the same profile convention as
ours, so `epsilon` is directly comparable:

```
their symbols:  psi = 0.5*(tanh(phi/(2*eps)) + 1)      (their Eq. 36)
our symbols:    phi = 0.5*(tanh(psi/(2*eps)) + 1)
                s   = 0.5*(1 + tanh(psi/(2*eps)))      (interface_relaxation.f90)
```

But the equations are **split** where ours are **fused**:

| | Saini (CLS) | Ours (CDI) |
|---|---|---|
| Transport | `dphi/dt + v.grad(phi) = S_vv(phi)` (their Eq. 38, our symbols) — pure advection, no compression | one equation: `dphi/dt + u.grad(phi) = div(eps*gamma*grad(phi)) + gamma*div(-phi*(1-phi)*n)` |
| Re-sharpening | separate pseudo-time equation `dphi/dtau + div(phi(1-phi)n) = div(eps*(grad(phi).n)n)` (Eq. 39), run every `dt_cls` for a few iterations | none — the compression term acts continuously, every step, rate set by `gamma` |
| Interface normal `n` | from a **separate, smooth signed-distance field** `psi`, transported alongside and periodically re-distanced (Eqs. 42–47, spelled out in §1.5) | directly `grad(phi)/|grad(phi)|` from the sharp field |
| Stabilisation | SVV in *every* equation above | none |

The compression/diffusion balance in their Eq. (39) is the same balance our
source term implements, so the physics of the sharpening is common. The
differences are all *numerical hygiene*, and two of them speak directly to what
we measured.

### 1.2 They independently document our `n_x` failure

Section 3.1, on why they do **not** compute `n` from the CLS field (quoted
verbatim, so `psi` here is *their* CLS field, i.e. our `phi`):

> "the computation of interface normal using the above equation is
> ill-conditioned due to the sharp profile of `psi`. Minor numerical
> oscillations in `psi` may cause spurious direction switching behavior of `n`
> on numerical integration between subsequent pseudo time steps. Desjardins et
> al. showed that this may lead to the spontaneous formation of spurious
> interfaces in the domain."

That is exactly the mechanism `../examples/interface_relaxation_multiphase/visualize.ipynb`
measured: a tiny grid-scale error in `phi` amplified into a corrupted `n_x`
over ~30 steps. Their answer is two-pronged — compute `n` from a smooth
auxiliary field *and* apply SVV. We currently do neither.

### 1.3 They run interfaces 2–6x sharper than we do

Their smearing width is `eps = xi*H`, `H` = element edge length, with
`xi = {0.5, 1, 1.5}/N`. Since our mean GLL spacing is `dx_GLL = H/N`, this is
exactly

```
eps/dx_GLL = 0.5, 1.0, or 1.5
```

Our resolution criterion (from `neko-multiphase-latex`) is `eps/dx_GLL >= 2.8`,
and our cases run at 2.8 (`zalesak`) to 4.0 (`advecting_drop`,
`interface_relaxation`). **Saini runs 2–6x sharper than our sharpest case** and
still reports 7% over/undershoots at discontinuities — which is the regime SVV
buys them. Two readings, both worth keeping in mind:

- their tests are a harsher probe of the same failure mode than anything we
  currently run, and
- a like-for-like comparison must match `eps/dx_GLL`, not `eps`.

### 1.4 What SVV is, in one paragraph

An artificial diffusion term `S_vv(u)` added to the right-hand side of every
transport equation. It is built by applying a **high-pass filter to the
discrete derivative matrices** in reference-element space, rather than by
setting a local diffusion coefficient (which is what artificial-viscosity
methods do). Concretely, `D_tilde = B*Q*B^-1*D`, where `B` is the 1D Legendre
Vandermonde matrix, `D` the GLL differentiation matrix, and `Q` a diagonal
power kernel

```
Q_kk = (k/N)^N_svv,     k = 0..N          (their Eq. 24)
```

so mode 0 is untouched and mode `N` passes at full strength — diffusion applied
only to the top of the spectrum. Two user parameters: `N_svv` (filter
strength; larger = sharper cutoff = *weaker* overall diffusion; they use
`N/2` by default, `N/4` or `N/6` for the hardest problems) and `c0` (viscosity
scale, `mu = c0*|c|*H/N`, so the cell Peclet number depends only on `N`). The
operator is provably symmetric positive definite, mesh-independent (computed
once on the reference element), and costs no more than the usual stiffness
matrix. **Note the machinery is the same modal transform we already built** in
`../examples/interface_relaxation_multiphase/spectral_dispersion_analysis.ipynb`
— `B` there is `coef.v`, `B^-1` is
`coef.vinv`.

**The detail that matters, and that this note originally left out: the kernel
appears on *both* derivative matrices.** Their Eq. (28) is

```
S_evv = [D_tilde]^T G [D_tilde]          (their Eq. 28)
```

i.e. the same construction as the standard Laplacian stiffness matrix, with
`D_tilde` in place of `D`. Under GLL quadrature `B^T M B` is diagonal, so this
equals a one-sided form (kernel applied once, inside the stiffness matrix) with
`Q` replaced by `Q^2` — and `((k/N)^p)^2 = (k/N)^(2p)`. **Their `N_svv` is a
one-sided form's `N_svv/2`**, so the two conventions must not be mixed when
comparing `c0` against their published values. See the Appendix of
`../CDI_IN_SEM.md`, "The `c0` discrepancy with Saini, resolved".

Three further specifics worth having written down:

- **`H` is local, not global.** `H = 2*J_e^(1/d)` from the element Jacobian, `d`
  the spatial dimension (Eq. 29) — the element edge on a uniform mesh.
- **`mu` is a diagonal matrix, not a scalar.** It is evaluated per GLL point from
  the **local** velocity (Eq. 31) and left-multiplies the assembled element
  operator, `D_mu S_evv` (Eq. 33). Note that placing it on the left is not
  conservative; putting the coefficient inside the bilinear form instead is, and
  is what our implementation does.
- **SVV is treated implicitly** — folded into the Helmholtz operator,
  `H = (b0/dt) M + D_mu S_vv` (Eq. 35), while advection and the non-linear terms
  stay explicit. This is not optional bookkeeping: the operator is stiff, and it
  is why they pay no timestep penalty for it.

### 1.5 The signed-distance (TLS) field, in full

Their §3.1 and Algorithm 1. This is the part of their method that supplies the
normal, and the reason our `n` and theirs are different objects. Two equations
and a re-seeding rule.

**Transport (their Eq. 43).** The distance field `psi` is advected alongside the
CLS field, and that is all it does — no compression, no source:

```
dpsi/dt + v.grad(psi) = S_vv(psi)
```

**Re-distancing (Eq. 44), in pseudo-time `tau`.** Transport does not preserve
`|grad psi| = 1`, so this restores it:

```
dpsi/dtau + w.grad(psi) = sgn(psi) + S_vv(psi)
w   = sgn(psi)*n,   n = grad(psi)/|grad(psi)|      (Eqs. 45, 42)
sgn(psi) = tanh(psi/(2*eps))                       (Eq. 46, smoothed)
```

Note `w.grad(psi) = sgn(psi)*n.grad(psi) = sgn(psi)*|grad psi|`, so Eq. (44) is
algebraically

```
dpsi/dtau = sgn(psi)*(1 - |grad psi|) + S_vv(psi)
```

which needs only `grad`, its magnitude and the SVV operator — no second
convective operator. That is the form worth implementing; it is Sussman's
classic re-distancing with their SVV bolted on.

**Re-seeding (Eq. 47).** Each re-distancing **discards** the transported `psi` and
restarts the pseudo-time iteration from the CLS field:

```
psi_0 = r_f*(phi - 0.5),   r_f = 0.1 in every test in the paper
```

`r_f` is called a regularisation factor and its stated job is to reduce the
gradients in the starting solution. Worth being explicit about what this means
for us: **`phi`'s noise does re-enter `psi` here.** But it enters through a seed
whose gradient is `0.1x` `psi`'s, is then reshaped by ~`O(100)` pseudo-time
iterations of an equation whose fixed point is `|grad d| = 1`, and is regularised
by that equation's SVV throughout. That is a different proposition from reading
`grad(psi)` pointwise every step, which is what `../CDI_IN_SEM.md` §6.2 parked.

**Banded, not global.** The iteration is stopped once the characteristic wave has
travelled `2.5H` from the interface. They are emphatic that an accurate *global*
distance field is unnecessary and not attempted — the requirement is only that
the gradient be "smooth and consistent" in a narrow band, enough to give a smooth
`n` for Eq. (39) and a smooth curvature for surface tension. Their Fig. 16b shows
the recovered `psi` matching the exact distance function only inside that band.

**Cadence.** Both re-initialisations are periodic in simulation time, and TLS
runs far more sparingly than CLS: `dt_tls ~ 10*dt_cls`. Actual values used range
`dt_cls = 0.001-0.05`, `dt_tls = 0.01-0.5`. In Algorithm 1 both happen *after*
the Navier-Stokes solve; TLS re-initialisation comes first, then CLS.

**SVV, per equation.** From §1.6 below: transport of *both* fields at
`N_svv = N/2, c0 = 0.1`; TLS re-distancing at `N_svv = N/4, c0 = 2` (`N/6` at
their hardest), i.e. by far the strongest setting anywhere in the paper. Their
§4.4 states flatly that "precluding spurious oscillations is not feasible without
the stabilizing SVV diffusion operator" for this equation. The reason is
geometric: the fixed point of Eq. (44) has kinks wherever the distance function
does — at medial axes and at any gradient discontinuity of the interface.

**The pseudo-timestep, and a trap of the same kind as the `N_svv` convention.**
They automate the iteration count (their §3.4):

```
dtau_cls = 0.1*H/(N+1),   N_cls = eps/dtau_cls
dtau_tls =     H/(N+1),   N_tls = 2.5*H/dtau_tls
```

**Do not adopt `dtau_tls = H/(N+1)` as written.** At our 1D slab's `N = 10`,
`H = 0.1` it is `9.1e-3` against `h_gll_min = 3.3e-3` — a pseudo-CFL of **2.75**,
unusable for any explicit treatment of a unit-speed characteristic. Their own
§4.4 standalone study of this equation used `dtau = 5e-4` at `H = 1/20, N = 8`,
which is **11x smaller** than their formula gives for those parameters, so the
paper is not self-consistent here. Set `dtau` from `h_gll_min` instead and
recover `N_tls = 2.5*H/dtau` from it, which is what preserves their band
distance — the quantity that actually matters.

Note also that Algorithm 1 line 2 *initialises* `psi` by solving Eq. (44) from
Eq. (47) rather than from an analytic distance, even when one is available.

### 1.6 Their SVV parameters, per test

Both parameters are set **per equation**, not per simulation — transport,
CLS re-initialisation, TLS re-distancing and momentum each get their own.

| test | equation | `N_svv` | `c0` |
|---|---|---|---|
| §4.1 sine wave | transport | `N/2`, `N/4` | `0.01, 0.1, 1` (swept) |
| §4.2 composite advection | transport | `N/2`, `N/4` | `0.1` |
| **§4.3 Zalesak** | **transport** | **`N/2`** | **`1.0`** |
| §4.4 intersecting circles | TLS re-distancing | `N/6` | `2.0` |
| §4.5 Rider-Kothe | transport + CLS re-init | `N/2` | `0.1` |
| §4.5 Rider-Kothe | TLS re-distancing | `N/4` | `2.0` |
| §5.1 stationary bubble | transport | `N/2` | `0.1` |
| §5.1 stationary bubble | CLS re-init | `N/4` | `1.0` |
| §5.1 stationary bubble | TLS re-init | `N/6` | `1.0` |
| §5.3 dam break | momentum | `N/4` | `1.0` |

The §5.1 row is stated in the paper to be the **default for every subsequent
dynamic test**; §5.2 and §5.3 then override only the TLS re-distancing equation
(`N/4, c0 = 2`; `N/6` at `N = 9` for the dam break). Momentum SVV is switched off
entirely for §5.1 so it cannot mask the spurious currents being measured.

The pattern: `N/2, c0 = 0.1-1` is the default everywhere, and only the
re-distancing equation — the one with kinks in its solution — needs more. **The
row that is ours is §4.3**, and it is their weakest setting.

---

## 2. The error norms they use

Defined for the Zalesak case (Eqs. 79–81) and reused everywhere. `psi_e` is the
exact solution, `psi_0` the initial condition (equal for cases that return to
their starting shape).

| Norm | Definition | What it measures |
|---|---|---|
| `E_r` | `int|psi - psi_e| / int psi_e` | relative L1 — total pointwise deviation, **dominated by over/undershoots near the interface** |
| `E_v` | `(int psi - int psi_0) / int psi_e` | volume (mass) conservation error |
| `E_s` | `int|H(psi-0.5) - H(psi_e-0.5)| / (2L int H(psi_e-0.5))`, `L` = exact perimeter | shape error (Herrmann) — how much enclosed area crossed to the wrong side of the interface |
| `l_avg` | `int I(psi) / int delta(psi)`, `I = 1` for `0.05 <= psi <= 0.95` | mean interface thickness, tracked in time; should stay at `l_0` |
| `||psi||_inf` | max over domain | boundedness — how far `psi` escapes `[0,1]` |

The key finding they repeat throughout: **`E_r` and `E_v` decouple.** `E_v`
converges spectrally (down to `1e-11` at `N=7`), while `E_r` converges slowly
because of Gibbs over/undershoots at the sharp `tanh`. The 0.5 isocontour — the
interface *location* — is transported accurately even when the profile around
it is ringing. This is a directly useful framing for us: `E_r` is the norm that
sees our dispersion problem; `E_v` will look fine regardless.

`l_avg` and `||psi||_inf` are the two cheap *time-series* diagnostics, and are
the closest thing in the paper to what we have been plotting.

---

## 3. Kinematic tests (Section 4) — prescribed velocity, no Navier-Stokes

These are the interesting ones for us: no flow solver, no surface tension, no
density ratio. Cheap, and they isolate the interface-capturing discretisation.

### 4.1 Transport of a sine wave — *"does the fix break spectral convergence?"*

| | |
|---|---|
| Domain | `[-1, 1]`, 1D, periodic |
| IC | `u(x) = sin(pi*x)` |
| Velocity | `c = 1` (constant) |
| Mesh | 4 elements, `N = 4..9` |
| Timestep | `dt = 1e-5`, run to `t = 1` (one flow-through) |
| Metric | L2 error vs exact, plotted against `N` |

**Purpose:** a control. A stabilisation operator must not destroy the
convergence rate on a *smooth* solution. Convergence slopes reported: no
diffusion ≈ −2.9, SVV `N_svv=N/2` ≈ −2, `N_svv=N/4` ≈ −1.5, AVM ≈ −2.3.

**The detail that matters to us:** even with a pure sine wave and no
stabilisation, *"the case with no diffusion shows flattening of the error decay
curve at the highest polynomial order, likely due to dispersion errors."*

Their attribution to dispersion is worth flagging rather than adopting. Our own
measurement of the CG SEM advection operator
(`../examples/saini_benchmarks/operator_analysis.ipynb`) finds the physical
branch accurate to `~1e-9` at `N >= 7`, so a smooth sine should transport
essentially exactly and dispersion is an unlikely explanation for a flattening
error curve — time discretisation or the CLS machinery around it are more
plausible in their setup. Either way it is *their* measurement in *their* code,
and it is **not** the phenomenon we are chasing: ours is the compression term
amplifying grid-scale content, not the transport scheme losing accuracy. See
`../CDI_IN_SEM.md`.

**For us:** trivially cheap, and the cleanest possible baseline — it needs no
level set at all, just scalar advection. Would be the first thing to reproduce.

### 4.2 Linear advection of a discontinuous profile — *"how much does it smear?"*

| | |
|---|---|
| Domain | `[0, 1]`, 1D, periodic |
| IC | composite profile (Lu et al.): square wave, triangle, smooth bump — discontinuities in both the solution and its gradient |
| Velocity | constant |
| Mesh | 10 elements, `N = 10` and `N = 20` |
| Timestep | `dt = 1e-4`, snapshots at `t = 1` (1 rotation) and `t = 100` (100 rotations) |
| Metric | qualitative profile comparison + over/undershoot magnitude |

**Purpose:** measure the *diffusivity* of the stabilisation. The exact solution
is known at all times and the discontinuities are **not self-steepening**, so
any smearing is purely numerical — unlike a shock, nothing regenerates the
gradient. The `t = 100` snapshot is a long-time-integration stress test.

**Findings:** ~7% over/undershoots persist at the square-wave edges even with
SVV; SVV is *less* diffusive than AVM (sharper edges, less attenuation of the
left wave), attributed to SVV acting on derivative matrices rather than
prescribing a local diffusivity. At `N=20`, `N_svv = N/2` and `N/4` give
visually identical results — the diffusion saturates, which they read as
evidence of the "spectrally vanishing" property.

**For us:** a very cheap, very sensitive dispersion probe, and it is the
closest 1D analogue of what our elements are doing at the interface. Long
integration (100 flow-throughs) is exactly the amplification regime we found.

### 4.3 Zalesak's disk — *"can it transport a shape with corners?"* **(we already have this case)**

| | |
|---|---|
| Domain | `[0, 1]^2` |
| IC | slotted disk, `r = 0.15` centred `(0.5, 0.75)`; slot half-width `0.025`, slot top at `y = 0.85`. Signed distance via `phi = -max(phi_0, min(phi_1, phi_2))` (Eq. 77), then `psi` from Eq. (36) |
| Velocity | solid-body rotation `u = pi(0.5 - y)`, `v = pi(x - 0.5)` — period 2 |
| Mesh | `H = {1/50, 1/100, 1/200}`, `N = {3, 5, 7, 9}` |
| Interface | `xi = {0.5, 1, 1.5}/N`, i.e. **`eps/dx_GLL = 0.5, 1.0, 1.5`** |
| Timestep | `dt = {1, 0.5, 0.25}e-4` (matched CFL across meshes) |
| Output | `t = 2` (1 rotation) and `t = 20` (**10 rotations**) |
| Metrics | `E_r`, `E_v`, `E_s`; plus the `psi` profile along `y = 0.75` |

**Purpose:** the standard benchmark for interface transport of a shape with
**sharp corners**, which no polynomial basis can represent. Errors concentrate
at the notch corners; the rest of the isocontour is indistinguishable from the
initial one at `N >= 5`. Note the CLS re-initialisation equation is **not**
active here (no re-distancing), so this test isolates transport.

**Key results:**
- `|E_v|` decays spectrally with `N` — `2.95e-11` at `H=1/50, N=7`, vs
  `4.04e-7` for `H=1/100, N=3` **at identical total GLL count** (Table 1).
  4 orders of magnitude for the same DOF.
- `E_r` decays slowly, *"since the formal order of accuracy decreases in the
  vicinity of the sharp tanh profile, which is associated to the well-known
  Gibbs phenomenon in high-order methods for non-smooth solutions."*
- Lower `xi` (sharper interface) gives visible over/undershoots but does **not**
  move the 0.5 isocontour.
- Appendix B: with `eps` held *constant* under h-refinement instead of scaling
  as `xi*H`, `E_r` converges much faster while `|E_v|` is nearly unchanged —
  i.e. `E_r`'s slow convergence is an interface-sharpness effect, not a
  transport-accuracy effect.

**For us — the single highest-value comparison.** We already run this case
(`examples/zalesak_disk_multiphase/`, `[0,1]^2`, 40x40, `N=7`, `eps=0.01`,
one rotation). Differences to close: their mesh is finer (`1/50` vs our
`1/40`), their interface much sharper (`eps/dx_GLL` 0.5–1.5 vs our 2.8), and
they run **ten** rotations, not one. Their Table 1 and Figs. 7–8 are published
numbers we can plot our `E_r/E_v/E_s` against directly.

### 4.4 Re-distancing around intersecting circles — *"the hardest scalar problem"*

| | |
|---|---|
| Domain | `[-2, 2]^2` |
| IC | `phi_0 = ((x-1)^2 + (y-1)^2 + 0.1) * phi_e` — the exact distance field to two intersecting circles (`r = 1`, centres `x = +/-a`, `a = 0.7`), deliberately **skewed** by a multiplicative factor (Eq. 83) |
| Equation | TLS re-distancing (Eq. 44) **only** — no transport, no flow |
| Mesh | `H = {1/5, 1/10, 1/20}`, `N = 3..8` |
| Timestep | `dtau = 5e-4` (CFL 0.28 at finest), run to `tau = 6` (quasi-steady) |
| SVV | strongest of all cases: `N_svv = N/6`, `c0 = 2` |
| Metric | `E_r(phi) = int|phi - phi_e| / int phi_e` |

**Purpose:** the paper's rigour test. Three difficulties stacked: highly
non-uniform gradients in the IC, **kinks (discontinuous gradients) along the
intersection line** of the two circles, and a non-linear convective operator
*with* a non-linear source term. They state flatly: *"Precluding spurious
oscillations is not feasible without the stabilizing SVV diffusion operator."*

**Findings:** error concentrates in the overlap region between the circles and
decreases with refinement, but p-refinement helps only modestly (Table 2) —
attributed to the gradient discontinuities and the non-linearity.

**For us:** this used to be filed as not applicable, on the grounds that CDI
computes `n` straight from the sharp field. `examples/saini_benchmarks/distance_normal/`
now implements Eq. (44), so this section is its **validation case and its warning
label**: it is the paper's own rigour test for the equation we are adding, and
the one place they say the SVV is not optional. Its geometric idea — **a kink,
i.e. a gradient discontinuity, sitting inside an element** — is also what our
high-order operator handles worst, and a distance field has kinks by
construction wherever the medial axis falls.

### 4.5 Rider-Kothe vortex — *"can it survive extreme stretching?"*

| | |
|---|---|
| Domain | `[0, 1]^2` |
| IC | disk `r = 0.15` centred `(0.5, 0.75)` (same as Zalesak, without the slot) |
| Velocity | `u = sin^2(pi x) sin(2 pi y) cos(pi t/T)`, `v = -sin(2 pi x) sin^2(pi y) cos(pi t/T)`, `T = 8` |
| Mesh | `H = {1/32, 1/64, 1/128}`, `N = {4, 5, 6}` |
| Interface | `xi = 1/N` (`eps/dx_GLL = 1`) |
| Timestep | `dt = {8, 4, 2}e-4` (CFL ~0.4); BDF3/EXT3 at the finest mesh |
| Re-init | `dt_tls = 0.5`, `dt_cls = 0.05` |
| Metrics | `E_r`/`E_L1`, `E_v`, `E_s` at `t = 8`; plus `l_avg(t)` and `||psi||_inf(t)` |

**Purpose:** the disk is stretched into a thin spiral filament by `t = 4`, then
the cosine reverses the flow and it must return to its initial shape at
`t = 8`. This is a **non-uniform** strain field (unlike Zalesak's rigid-body
rotation), so the interface genuinely thins and the re-initialisation machinery
is load-bearing. It exercises the *complete* algorithm — CLS transport, TLS
transport, both re-initialisations.

**Findings:** `l_avg` stays within a few percent of `l_0` throughout, including
at maximum stretching; `||psi||_inf` overshoots decrease with h-refinement and
do not grow in time (energy stability from SVV). `|E_v|` converges more slowly
with `N` than Zalesak — attributed to the TLS re-distancing, inactive in
Zalesak.

**For us — the best "reversibility" test.** Same structural idea as
`advecting_drop_multiphase` (return to the exact starting shape, so the error
is unambiguous) but with real deformation instead of pure translation. It also
gives the two time-series diagnostics (`l_avg`, `||psi||_inf`) that would show
CDI thickness drift and boundedness loss directly, and it needs no Navier-Stokes
— only a prescribed velocity field, which our `.f90` files already do for
Zalesak.

---

## 4. Dynamic tests (Section 5) — coupled to Navier-Stokes

These validate the *two-phase flow* implementation (density ratio, surface
tension, pressure splitting), not the interface discretisation in isolation.
Listed for completeness; all are considerably more expensive, and all confound
the dispersion question with flow-solver behaviour.

### 5.1 Stationary bubble — spurious currents

`[-2, 2]^2`, unit-diameter bubble at the centre, symmetric BCs, **density and
viscosity ratios both 1** so only surface tension acts. `We = 1`, `Re` varied
to sweep `La = Re^2/We`. `H = {0.2, 0.1, 0.05}`, `N = {3, 5, 7, 9}`,
`dt = 5e-4`. Metric: `v_rms(t)` — the exact solution is zero velocity, so all
motion is error. SVV is turned **off** in the momentum equation here, to avoid
masking the currents.

*Purpose:* isolate curvature/surface-tension error. Time-averaged `v_rms` drops
from `1.42e-3` (`N=3`) to `3.1e-5` (`N=9`). They note visible jumps in `v_rms`
at each TLS re-initialisation, and that no special curvature treatment (height
function, least-squares, fast marching) was needed.

*For us:* the natural successor to `couette_surface_tension_multiphase/`, and
the standard test if surface tension ever becomes the focus. Not a dispersion
test.

### 5.2 Rayleigh-Taylor instability

`[L, 4L]`, `rho2/rho1 = 3`, `At = 0.5`, `Re = 3000`, no surface tension.
Interface perturbed as `y = 2L + 0.1L cos(2 pi x/L)`. `H = {1/25, 1/50, 1/100}L`,
`N = {3, 4, 5}`, `dt = 1.25e-4`. Metrics: interface position at `x = 0` and
`x = L/2` vs published data (Guermond, Chiu), `|E_v|(t)`, `l_avg(t)`,
`||psi||_inf(t)`.

*Purpose:* topological change under a real flow. `l_avg` stays within 4%;
overshoots stay under 2% until the interface breaks up (`t > 1.5`), then reach
~10%.

*Note:* the initial condition is **a cosine-perturbed interface** — the same
geometry as our `interface_relaxation` case, but with the flow switched on.

### 5.3 Dam break

`[0, 5L] x [0, 1.25L]`, water column `L x L`. Real air/water:
`rho_l/rho_g = 833.33`, `mu_l/mu_g = 55.55`, `Re ~ 42791`, `We ~ 534`, `Fr = 1`.
`H = 1/32` (160x40 in 2D, 160x40x40 in 3D), `N = {3, 5, 7, 9}` in 2D, `N = 5`
in 3D. `xi = 1.5/N`. Metrics: surge-front position and dam height vs Rodriguez
et al.

*Purpose:* validation against experiment at a realistic density ratio, plus the
pressure-extrapolation study (P1 vs P2 vs unsplit). They needed to explicitly
cut off the top two modes of the extrapolated pressure to keep P2 stable —
another instance of high-order modal filtering being load-bearing.

### 5.4 Rising bubbles

**5.4.1 (pressure extrapolation):** air/water, `Re = 217`, `We = 121`, `Fr = 1`,
`rho_l/rho_g = 826.45`, `mu_l/mu_g = 4081.63`. Triply periodic
`[0,2D]^2 x [0,4D]`, 32x32x64 elements, `N = 5`, `dt = 1e-4`. Validated against
Dodd et al. Table 4 is the performance result: P2 extrapolation converges the
pressure Poisson solve in 8.35 iterations vs 190.86 unsplit — 6x faster
wall-clock per step.

**5.4.2 (complex deformation):** `z in [0,20]`, `x,y in [-10,10]`, unit-radius
bubble at `(0,0,8)`, `rho_l/rho_g = 1000`, `mu_l/mu_g = 100`, ~1.46M elements at
`N = 5`. Two conditions from Tripathi et al.: `Ga = 2.316, Eo = 29` (terminal
spherical cap, validated against Bhaga & Weber experiments) and
`Ga = 70.7, Eo = 200` (skirt breakup).

*For us:* out of scope for now — 3D, large, and testing the flow solver rather
than the interface discretisation.

---

## 5. What is worth adopting, ranked

Ordered by (evidence about our dispersion problem) / (effort). The first three
need **no Navier-Stokes and no new physics** — only a prescribed velocity field
and post-processing, both of which our existing cases already do.

1. **§4.1 sine wave** — 1D, 4 elements, seconds to run. Establishes the
   baseline claim ("unstabilised high-order SEM loses convergence to dispersion
   error") in the simplest possible setting, with a published slope to check
   against. Needs only scalar advection, no level set.
2. **§4.2 discontinuous advection** — 1D, 10 elements, and 100 flow-throughs is
   still cheap. The most sensitive dispersion probe in the paper, and directly
   comparable to their Figs. 3–4.
3. **§4.3 Zalesak, extended** — we already have the case. Add: the `E_r/E_v/E_s`
   norms, `xi` matched to theirs (`eps/dx_GLL = 0.5..1.5`), an `N` sweep at
   fixed DOF (their Table 1), and **ten rotations** instead of one. Gives a
   direct numeric comparison against published data.
4. **§4.5 Rider-Kothe** — new case, but structurally close to
   `advecting_drop_multiphase` (prescribed velocity, returns to its initial
   shape). Adds real interface stretching and the `l_avg`/`||psi||_inf`
   time-series diagnostics.
5. **§4.4 intersecting circles** — only if we pursue an auxiliary distance
   field. Its kink geometry is transplantable on its own, though.
6. **§5.x dynamic tests** — deferred. They validate two-phase flow, not the
   interface discretisation, and would confound the dispersion question.

Two cross-cutting things to take regardless of which case we run:

- **The error norms of §2.** `E_r` (sees ringing) and `E_v` (does not) is the
  decomposition that makes the argument legible, and their decoupling is itself
  the finding.
- **`eps/dx_GLL` as the controlled variable**, not `eps`. Their `xi = c/N`
  convention *is* `eps/dx_GLL = c`, and it is the only way to compare across
  `N` and `H`.
