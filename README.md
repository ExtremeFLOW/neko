# neko-multiphase-psi-transport

Neko, plus four cases for its Conservative Diffuse Interface (CDI) method with the compression
term's interface normal taken from a **separately transported signed-distance field** $\psi$, not
from the phase field $\phi$'s own gradient.
- **What comes from where.** The $\psi$ machinery (transport with SVV; re-initialization by
  Eqs. 44–47) follows Saini & Tomboulides (JCP 2026). The $\phi$ equation is Neko's fused CDI.
- **What the cases do.** They replicate Saini §4.2–§4.5 with that combination.
- **The sibling repo.** Why the $\phi$-gradient normal fails on a spectral element method is the
  investigation in `../neko-multiphase/`. Neko's own README is
  [`UPSTREAM_README.md`](UPSTREAM_README.md).

## Results

All validated results transport $\psi$ from the exact distance, with no re-initialization.

| Saini | case | configuration | $E_r$ | worst violation | Saini, our normalisation | status |
|---|---|---|---|---|---|---|
| §4.2 | [`advecting_slab_1d`](examples/advecting_slab_1d/README.md) | $\psi$-normal, $\xi=1$, $\gamma=1$, $N=10$, $t=20$ | 0.00006 | $1.2\times10^{-10}$ | no $E_r$ in §4.2 | validated |
| | | the same at $\xi=1.5$ | 0.00001 | 0 | | |
| | | $\phi$-normal, $\xi=1$ | $\approx1.2$ | $3.1\times10^{-2}$ | | direction only |
| §4.3 | [`zalesak_disk`](examples/zalesak_disk/README.md) | $\xi=1$, $H=1/50$, $N=3/5/7$, ten rotations | 0.1012 / 0.0208 / 0.0040 | $\le9.5\times10^{-6}$ | not compared yet | validated |
| | | shipped: $\xi=2.8$, $N=5$, ten rotations | 0.0221 | 0 | | validated |
| | | the same, $\phi$-normal | 1.07 | $6.3\times10^{-3}$ | | direction only |
| §4.4 | [`redistance_circles`](examples/redistance_circles/README.md) | Eq. (44) alone, his configuration, Table 2's four cells, $\tau=6$ | $1.3$–$6.6\times10^{-3}$ | — | ours/his 0.978–1.011 | reproduced |
| §4.5 | [`rider_kothe`](examples/rider_kothe/README.md) | transport only, $\xi=1$, $H=1/64$, $N=5$, $t=8$ | 0.0452 | $2.6\times10^{-3}$ | 0.0410 | validated |
| | | the same at $H=1/128$ | 0.0104 | $2.4\times10^{-3}$ | 0.0258 (his $N=3$) | validated |
| | | re-initialization in his configuration, our $\Delta\tau$ (scratch file) | 0.0402 | $7.1\times10^{-4}$ | 0.0410 | not in the code yet |

**The norms.** These are Saini's Eqs. (79)–(81), computed by `examples/norms.py`; `redistance_circles`
uses his Eq. (84), $\int|\psi-\psi_e|/\int\psi_e$.
- $E_r=\int|\phi-\phi_e|\,/\int\phi_e$ is the relative $L^1$ error.
- $E_v=(\int\phi-\int\phi_0)/\int\phi_e$ is the volume error.
- $E_s$ is the shape error: the area that crossed the 0.5 contour, over twice the exact perimeter
  times the enclosed area.

**Worst violation** is $\max(\phi-1,-\phi,0)$ over the output frames: every 0.05 on the slab and the
Zalesak primary, and every 0.08 on Rider–Kothe (0.04 at $H=1/96$ and 1/128). Zalesak's $\xi=1$ row
is from the solver's step log (every 200 steps). The primary's 0 holds in the frames; its step log
has a $2.3\times10^{-8}$ transient at step 1.

**Saini's Rider–Kothe values** are his own code's $t=8$ dumps, rescored with `norms.E_r`: his code
divides by the area *outside* the disk. His 0.0410 includes his phase-field re-sharpening (Eq. 39);
without it his code gives 0.208.

**The $\phi$-normal rows** mostly predate the `grad_floor` fix. Read them for direction only.

Runs in other configurations are in each case README's appendix.

## Recommended settings

Boundedness comes first: an unbounded phase field, coupled to a flow solver, corrupts densities and
viscosities. Speed comes last.
- **Without strain** (rigid rotation, uniform advection): $\xi\approx1.5$ and a low $\gamma$
  (0.25–0.5).
  - $\xi\ge1.5$ is bounded to round-off at every $\gamma$.
  - $\gamma$ is not neutral: on Zalesak $E_r$ rises 1.5–2× from $\gamma=0.25$ to 2
    (`examples/zalesak_disk/README.md` §3.5).
- **With strain**, $\xi$ trades shape against boundedness (`examples/rider_kothe/README.md` §3.2).
  - $\xi\approx1$ gives the best shape at a violation of $2$–$3\times10^{-3}$.
  - $\xi\ge1.5$ is bounded to round-off, at about 40% worse $E_r$.
- **SVV on $\psi$ is always on**: $c_0=0.1$, $N_{svv}=N/2$, Saini's value for his $\psi$ transport.
  Without it $p$-refinement reverses (`CDI_METHOD.md` §5). There is no SVV on $\phi$.
- **Refine if you can.** At fixed $\xi$, $h$-refinement improves shape and filament retention
  together.

$\xi=\varepsilon N/H$ is the interface width in units of $H/N$. It is **$N$ times Saini's $\xi$**
(`CDI_METHOD.md` §1).

## The docs

| file | holds |
|---|---|
| [`CDI_METHOD.md`](CDI_METHOD.md) | the method: notation, the $\phi$ equation, $\psi$ and the normal, SVV, and **§6, Saini's configuration against ours** |
| [`REDISTANCING.md`](REDISTANCING.md) | maintaining $\psi$: the terms, when it is needed, what is known, and how the code solves Eq. (44) |
| `examples/<case>/README.md` | each case: problem, configuration, results, limitations, running, evidence, other runs |
| [`NEXT_SESSION.md`](NEXT_SESSION.md) | what is open: Rider–Kothe in Saini's configuration first |
| [`CLAUDE.md`](CLAUDE.md) | the rules this repo is maintained by |

Each case is one `.f90` user file, its `.case` files, a `run.sh`, a `visualize.ipynb` (shipped
executed) and an `evidence/` folder from this repo's own runs.

## Building Neko

CPU build (in-tree):

```bash
source setup-env.sh
./regen.sh             # only needed once, before the first configure
./configure --prefix=$NEKO_PSI_TRANSPORT_PREFIX
make -j install
```

GPU build (CUDA, out-of-tree; `sm_86` is this workstation's RTX 3090, so check yours with
`nvidia-smi --query-gpu=compute_cap --format=csv,noheader`):

```bash
source setup-env-cuda.sh
mkdir -p build-cuda && cd build-cuda
../configure --prefix=$NEKO_PSI_TRANSPORT_CUDA_PREFIX --with-cuda=/usr CUDA_ARCH=-arch=sm_86
make -j install
```

- **Switching between the two.** `configure` refuses to run in-tree while an out-of-tree build has
  generated Makefiles there, and vice versa. Run `make distclean` in the repo root before switching.
- **Which `neko` runs.** The two builds install to separate prefixes. Whichever `setup-env*.sh` you
  sourced decides which `neko` runs.
- **`advecting_slab_1d` is CPU-only.** Neko's CUDA kernels cap at `polynomial_order` $\le9$.

## Running a case

```bash
source setup-env.sh        # or setup-env-cuda.sh, once per shell
cd examples/<case>
genmeshbox <args...>       # once per case; the command is in run.sh and the case README
./run.sh                   # or ./run.sh <case> for one .case
```

## Visualizing a result

```bash
jupyter notebook examples/<case>/visualize.ipynb
```

It uses the `venv-neko` kernel at `/lscratch/sieburgh/local/venv-neko` (pysemtools, matplotlib,
mpi4py). Open it through a Jupyter server, not as a file in a browser.
