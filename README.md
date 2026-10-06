# neko-multiphase-psi-transport

Neko, plus three demo cases showcasing Neko's CDI (Conservative Diffuse
Interface) method for two-phase flow once its compression term is given a
usable interface normal: computed from a **separately transported
signed-distance field** (`psi`) rather than straight from the sharp phase
field's own gradient.

**All three validated results use transport alone.** `psi` is seeded once from the
exact signed distance and then simply advected — no reinitialization, no
redistancing. The method supports periodic redistancing, and `zalesak_disk`
carries an ablation that switches it on, but nothing in the validated
configurations needs it. See `CDI_METHOD.md` §3–§4, which keeps
*redistancing* (pseudo-time relaxation of an existing field) and
*reinitialization* (discarding it and rebuilding from `phi`) distinct. See `UPSTREAM_README.md` for Neko's own
README (citations, publications, general build docs).

**Start with [`CDI_METHOD.md`](CDI_METHOD.md)** for the equations, the
naming convention, and why this repo is named the way it is (not after any
one paper — see that file's opening note). This repo is a clean showcase of
a working result; the forensic investigation that found it — why the naive
`phi`-gradient normal fails on a spectral element method, and the evidence
behind every design decision here — lives in the sibling repo
`../neko-multiphase/` (`CDI_IN_SEM.md` is its standing conclusion).

[`REDISTANCING.md`](REDISTANCING.md) is the companion on redistancing: why the
validated results transport `psi` without it, why a `psi` *built* by Saini's
Eq. (44) is the path that carries over to geometries with no analytic distance,
and, in its §9, how the Fortran solves Eq. (44) routine by routine.

## Demo cases

Three of the four are the showcase. `redistance_circles` is listed with them
because it shares the machinery, but it is a **reproduced published benchmark**,
not a showcase result. Its [`README.md`](examples/redistance_circles/README.md)
describes it: configuration, the paper against the code, results and limitations.

| Directory | What it shows | Status |
|---|---|---|
| `examples/advecting_slab_1d/` | 1D advection. Swapping the normal from $\nabla\phi$ to $\nabla\psi$ takes $E_r$ from $\approx1.2$ to 0.00006 (at $\xi=0.5$ the $\phi$-normal diverges). $\xi$/$\gamma$ sweeps locate the exactly-bounded envelope, and a p-sweep — where there is no geometry to under-resolve — shows ψ converging 27× while φ saturates at $E_r\approx1.2$. | Validated; swept in $\xi$, $\gamma$ and $N$ |
| `examples/zalesak_disk/` | The main 2D result so far: slotted disk, ten rotations, $E_r = 0.022$ with **zero** boundedness violations. Single-variable $\phi\to\psi$ is **1.07 → 0.022**. SVV on ψ is **what makes $p$-refinement work**: with it $E_r$ falls 25× from $N=3$ to 7, without it $E_r$ *rises* — the on/off gap grows 9× → 47× → **256×**. Under refinement ψ improves and stays exactly bounded while φ degrades. A $5\times4$ $(\xi,\gamma)$ map locates the stability boundary between $\xi=1$ and $\xi=0.5$ and shows $\gamma$ is **not** neutral on an isometry. Carries the repo's one sanctioned ablation, two variants that add redistancing; their recorded numbers predate the current configuration and wait for a re-run. | Validated; SVV, resolution and $(\xi,\gamma)$ studies complete |
| `examples/redistance_circles/` | **A benchmark reproduction, not a showcase result.** Saini et al. §4.4, their own rigour test for the Eq. (44) redistancing PDE, solved standalone, with the authors' configuration from their public case: sign-function $\varepsilon=0.25$ in Eq. (46), dealiased $\mathbf C(\mathbf w)$, their sign guard, BDF2/EXT2 and the Eq. (31) SVV. 12 of the 18 Fig. 12 cells are within 0.1% of their values, four more within 2.2%. Two cells ($H{=}1/5$, $N{=}4$ and 8) fail on our mesh through zero-set round-off; started from their zero-set $\psi_0$ they match exactly. | Reproduced (2026-09-29) |
| `examples/rider_kothe/` | Vortex-in-a-box — the only case with real strain, and so the only one that can settle the open questions. $\|\nabla\psi\|$ drifts to about 8× at maximum stretch and **returns to 1.008**, while $\psi$'s normal stays within 0.6–1.3° of the exact interface. Saini's periodic reseed, run with his settings, completes but points 4–8° off, so redistancing has not paid here yet (open, `NEXT_SESSION.md`). Strain **inverts** the $\xi$ recommendation: it trades shape accuracy against boundedness. $h$-refinement converges at about second order ($E_r$ 0.0452 → 0.0104 over $H=1/64$ → 1/128); at $H=1/64$ it is level with Saini's own code (0.041). | Validated, re-run 2026-10-02 after the diffusion fix and 2026-10-05/06 after the velocity fix; $\xi$/$\gamma$ cross and $h$-series complete |

The $\phi$-normal runs predate the `grad_floor` fix (`CDI_METHOD.md` §4.1) and were not re-run; read them for direction only.

Each case is self-contained: one `.f90` user file, its `.case` configs, a
`run.sh`, and a `README.md` with its exact parameters and results.
Each adds a `visualize.ipynb` and an `evidence/` folder holding the animations
and figures that notebook produces — all from this repo's own runs, not carried
over. The notebooks ship executed, so the numbers are visible without rerunning
them, except `rider_kothe`'s, whose outputs predate the 2026-10-02 fixes
(regenerating them is in `NEXT_SESSION.md`).

The showcase cases carry more than one `.case` on purpose:
`advecting_slab_1d` ships a $\xi$ and $\gamma$ sweep (the operating envelope is
part of what the case demonstrates — see its README), `zalesak_disk` ships
its single-variable $\phi$-normal baseline plus the sanctioned redistancing
ablation, and `rider_kothe` its $\xi$/$\gamma$ cross and $h$-series. Files are named for the parameter they vary. See `CLAUDE.md` for the
conventions, and `CDI_METHOD.md` §7 for what each status claim is (and isn't)
backed by.

**$\xi$, used throughout:** $\xi = \varepsilon N/H$ — the interface width
$\varepsilon$ in units of $H/N$ (element edge over polynomial order), so
$\varepsilon = \xi H/N$. It is **$N$ times Saini et al.'s $\xi$**, and coincides
with the numerator of their $\xi = c/N$. `CDI_METHOD.md` §1 is the canonical
definition and carries the full comparison; never compare the bare symbol across
the two.

## Recommended settings

**Regime-dependent, and deliberately not a single number.**

- **No strain** (rigid rotation, uniform advection): $\xi \approx 1.5$ and the
  **low end of $\gamma$**. A full $(\xi,\gamma)$ map on Zalesak (one rotation,
  $N=7$; see that case's README) shows every $\xi \ge 1.5$ cell bounded to
  round-off across $\gamma \in [0.25, 2]$ — the good region is broad — while
  $\xi=1$ is bounded only to $\sim10^{-5}$ and $\xi=0.5$ fails outright.
  Within that region $\gamma$ is **not** neutral, contrary to what an isometry
  argument predicts: $E_r$ rises monotonically with $\gamma$, about $2\times$
  from 0.25 to 2, because $\gamma$ scales the whole CDI relaxation and a rigid
  rotation's exact solution contains no relaxation at all. That row was re-run
  at **fixed $\Delta t$** to rule out the $\Delta t \propto 1/\gamma$ coupling:
  $\Delta t$ moves $E_r$ by at most 5% over a $6.8\times$ change, $\gamma$ by
  $2.06\times$. Best exactly-bounded
  cell: $\xi=1.5,\ \gamma=0.25$, at $2.1\times$ better $E_r$ than
  $\xi=2.8,\ \gamma=1$. (Over *ten* rotations at $N=5$, $\xi=1.5$ and
  $\xi=2.8$ tie at $E_r \approx 0.022$, both exactly bounded — so the
  $\xi$ preference is settled and the $\gamma$ one is so far a one-rotation
  result.)
- **With strain**: $\xi$ **trades shape accuracy against boundedness**.
  - $\xi \approx 1.0$ recovers the shape best ($E_r$ 0.0452 on Rider–Kothe at $H=1/64$) at a
    violation of $2.4$–$3.0\times10^{-3}$, which does not shrink under $h$-refinement.
  - $\xi \gtrsim 1.5$ is exactly bounded, at about 40% worse $E_r$ (0.063).
  - The reason: a filament $d$ thick holds a core of $1-e^{-d/2\varepsilon}$ at CDI equilibrium,
    so a smaller $\varepsilon$ keeps more of it above 0.5. That is the opposite of the
    strain-free conclusion.
- **SVV on $\psi$ is not optional, and the code enforces it.** A $\psi$-normal run with
  `svv_psi` off stops at startup. It is the largest single effect measured here: on Zalesak,
  without it $E_r$ is 9–256× worse and, decisively, *$p$-refinement reverses*, so a higher-order
  run becomes a worse one. In 1D it buys seven orders of boundedness. Use $c_0=0.1$, $N/2$,
  Saini's TLS-transport value. There is no SVV on $\phi$; SVV is a $\psi$-only knob. See
  `CDI_METHOD.md` §5.
- **Refine if you can.** At fixed $\xi$, $h$-refinement shrinks $\varepsilon$ and
  the mesh together and improves shape, boundedness and filament representation
  simultaneously. It is the only lever with no downside.

Boundedness is the first criterion, not one of several: an unbounded phase field
is unphysical and, coupled to a flow solver, propagates into densities and
viscosities where the damage is neither local nor recoverable. **Speed ranks
last.** See `CDI_METHOD.md`.

## Building Neko

CPU build (in-tree):

```bash
source setup-env.sh
./regen.sh             # only needed once, before the first configure
./configure --prefix=$NEKO_PSI_TRANSPORT_PREFIX
make -j install
```

GPU build (CUDA, out-of-tree — this workstation has an RTX 3090, compute
capability 8.6/`sm_86`; adjust `CUDA_ARCH` for a different GPU, e.g. via
`nvidia-smi --query-gpu=compute_cap --format=csv,noheader`):

```bash
source setup-env-cuda.sh
mkdir -p build-cuda && cd build-cuda
../configure --prefix=$NEKO_PSI_TRANSPORT_CUDA_PREFIX --with-cuda=/usr CUDA_ARCH=-arch=sm_86
make -j install
```

`configure` refuses to run in-tree while an out-of-tree build (or vice versa)
already has generated Makefiles there — if you need to reconfigure the
in-tree CPU build after building CUDA (or vice versa), `make distclean` in
the repo root first. The two builds install to separate prefixes and don't
otherwise conflict (see `setup-env.sh` / `setup-env-cuda.sh` for the exact prefixes).
The user files are backend-agnostic and the 2D cases run unmodified on either
build — whichever `neko` is first on `PATH` (i.e. whichever `setup-env*.sh` you
sourced) is the one that runs. CPU and GPU were checked against each other on
`zalesak_disk` and agree to 1.1e-9 in `phi` after 4000 steps, with GPU ~8.8x
faster.

**`advecting_slab_1d` is CPU-only**, and that is a Neko limit rather than
anything about this repo: its CUDA kernels (`opr_dudxyz`, `opr_cfl`,
`opr_conv1`) cap at `lx = 10`, i.e. `polynomial_order` <= 9, and that case is
settled at `N = 10`.

## Running a case

```bash
source setup-env.sh        # CPU, or setup-env-cuda.sh for GPU — once per shell
cd examples/<case>
genmeshbox <args...>  # once per case; exact command commented in run.sh
./run.sh
```

## Visualizing a result

```bash
jupyter notebook examples/<case>/visualize.ipynb
```

Uses the shared `venv-neko` kernel (pysemtools/matplotlib/mpi4py) at
`/lscratch/sieburgh/local/venv-neko` — same as `../neko-multiphase/`. Open
via an actual Jupyter server, not by opening the `.ipynb` directly in a
browser.

## If you're picking this up fresh

Read `CDI_METHOD.md` first. [`NEXT_SESSION.md`](NEXT_SESSION.md) is the Rider–Kothe
redistancing plan, with Saini §4.5's reference values and the items parked for the
other cases.

**For the Saini §4.4 reproduction, read
[`examples/redistance_circles/README.md`](examples/redistance_circles/README.md)
first.** It is the single record of that work.
