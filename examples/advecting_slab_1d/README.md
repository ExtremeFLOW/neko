# advecting_slab_1d

1D linear advection of a $\tanh$ slab through a periodic domain for twenty flow-throughs, under
Neko's CDI equation. It isolates one change, where the compression normal comes from, with
everything else held fixed. It is also the four-minute regression testbed for the machinery shared
with the 2D cases.

**Status: validated.**
- **`svv_psi` is on in every $\psi$-normal case.** The 17 sweep cases use $c_0=0.1$, $N/2$, and the
  tables are from them. `psi_svvon` and the five `psi_rd_*` testbed cases use $c_0=1$.
- **The $\phi$-normal rows** mostly predate the `grad_floor` fix (`../../CDI_METHOD.md` §3); read
  them for direction only.
- **Uniform advection is an isometry**, so this case cannot test whether $\psi$ needs maintaining
  (`../../REDISTANCING.md` §2).

$\xi=\varepsilon N/H$ (`../../CDI_METHOD.md` §1). Paths under `logs/` are gitignored and local.

## 1. The problem

The 1D analogue of Saini §4.2, which advects a composite profile of discontinuities with SVV. It
uses the same 10-element periodic mesh on $[0,1]$, $N=10$ and $\Delta t=10^{-4}$, but here a tanh
slab carried by CDI. He reports no $E_r$ there.
- **Normal.** In 1D the normal is just $\operatorname{sign}(\partial_x\psi)$, so only $\psi$'s zero
  set matters.
- **No drift.** Uniform advection conserves $|\nabla\psi|$ exactly, so $\psi$ has nothing to drift
  from.

## 2. Configuration

- **Mesh and time.** `genmeshbox 0 1 0 0.1 0 0.1 10 1 1 .true. .true. .true.`, $N=10$, so
  $H/N=0.01$ and $\varepsilon=0.01\,\xi$. $\Delta t=10^{-4}$ (the $p$-sweep scales it with $N$),
  `end_time` 20, and a uniform $u=1$ (`case.fluid.freeze = true`).
- **28 `.case` files, named for what they vary:**
  - **$\xi$ sweep:** `psi_xi05` … `psi_xi28` (7).
  - **$\gamma$ sweep:** `gamma*_xi10`, `gamma*_xi15` (7, including the $\gamma=0$ null).
  - **$p$-sweep:** `psi_p6/p8/p12` and `phi_p6/p8/p12` (6).
  - **$\phi$-normal:** `phi_xi05`, `phi_xi10` (2).
  - **`psi_svvon`** (1).
  - **The redistancing testbed:** `psi_rd_*` (5).
- **CPU only.** Neko's CUDA kernels cap at `lx = 10`, i.e. `polynomial_order` $\le9$.

## 3. Results

**How the numbers are taken.**
- **At $t=20$, against the exact slab at the frame's own time.** Neko writes its last frame one step
  past `end_time`, which is worth up to about 9× in $E_r$ (at $\xi=1$; about 2× at 0.5 and 2.8).
- **"Worst violation"** is $\max(\phi-1,-\phi,0)$ over every frame, at 0.05 spacing.
- **Sources.** $\psi$ rows: `logs/svv_off_2026-10-02/slab_table.txt` (made by
  `scripts/slab_table.py`). $\phi$ rows: the shipped outputs, measured the same way on 2026-10-05.

### 3.1 $\psi$-normal against $\phi$-normal

| `.case` | $\xi$ | $\gamma$ | normal | $E_r$ | worst violation |
|---|---|---|---|---|---|
| `gamma0_xi10` | 1.0 | 0 | — | 0.01486 | $1.7\times10^{-2}$ |
| `gamma025_xi10` | 1.0 | 0.25 | $\psi$ | 0.00007 | $6.0\times10^{-6}$ |
| `psi_xi10` | 1.0 | 1 | $\psi$ | **0.00006** | $1.2\times10^{-10}$ |
| `psi_xi15` | 1.5 | 1 | $\psi$ | 0.00001 | **0** |
| `psi_xi28` | 2.8 | 1 | $\psi$ | 0.00042 | **0** |
| `psi_xi05` | 0.5 | 1 | $\psi$ | **0.00044** | $1.3\times10^{-2}$ |
| `phi_xi10` | 1.0 | 1 | $\phi$ | **≈1.2** | $3.1\times10^{-2}$ |
| `phi_xi05` | 0.5 | 1 | $\phi$ | **diverged** at $t=6.80$ | — |

- **Read the $\phi$-normal $E_r$ to one significant figure.** It is an instability, not a converged
  answer: a solver-tolerance change moves it between 1.11 and 1.32.
- **It also ends inside $[0,1]$** after violating on the way, which is why the table reports the
  worst violation over the run.
- **$\gamma=0$ is the null.** Without compression and its balancing diffusion the slab only advects,
  and the method beats it by 250× in $E_r$ and $10^8$ in boundedness.

### 3.2 Where boundedness turns on: the $\xi$ and $\gamma$ sweeps

$\xi$ at $\gamma=1$:

| $\xi$ | 0.5 | 0.6 | 0.8 | 1.0 | **1.5** | 2.0 | 2.8 |
|---|---|---|---|---|---|---|---|
| $E_r$ | 0.00044 | 0.00046 | 0.00023 | 0.00006 | **0.00001** | 0.00005 | 0.00042 |
| worst violation | 1.3e-2 | 6.5e-3 | 3.3e-4 | 1.2e-10 | **0** | 0 | 0 |

$\gamma$ at $\xi=1.0$ and at $\xi=1.5$:

| $\gamma$ | 0 | 0.25 | 0.5 | 1.0 | 1.5 |
|---|---|---|---|---|---|
| $E_r$, $\xi=1$ | 0.01486 | 0.00007 | 0.00006 | 0.00006 | 0.00006 |
| worst violation, $\xi=1$ | 1.7e-2 | 6.0e-6 | 1.8e-5 | 1.2e-10 | **0** |
| $E_r$, $\xi=1.5$ | | **0.00001** | 0.00001 | 0.00001 | 0.00001 |
| worst violation, $\xi=1.5$ | | **0** | 0 | 0 | 0 |
| $\Delta t$ headroom, $0.05/C_{\text{comp}}$, $\xi=1.5$ | | **6.6×** | 3.3× | 1.65× | 1.1× |

- **Exact boundedness switches on between $\xi=1.0$ and 1.5**, and $\xi=1.5$ is also the $E_r$
  optimum.
- **$\gamma$ buys boundedness, not accuracy.** $E_r$ is flat for $\gamma\ge0.25$.
  - At $\xi=1$, exact boundedness needs $\gamma=1.5$, close to the compression ceiling
    ($C_{\text{comp}}=0.0455$).
  - At $\xi=1.5$ it holds at every $\gamma\ge0.25$, with up to 6.6× the time step.
- **What it picks out**, on a flow without strain: $\xi\approx1.5$, $\gamma\approx0.25$–0.5
  (`../../README.md`, "Recommended settings").
- **The shipped primaries stay at $\xi=1$ and 0.5.** They are the points the $\phi$-vs-$\psi$ result
  is quoted at, and 0.5 is Saini's sharpest width.

### 3.3 $p$-refinement isolates the operator

At fixed $\xi=1$, only $N$ changes (with $\varepsilon=H/N$ and $\Delta t$ scaled to match), and 1D
has no geometry to under-resolve. The $\psi$ rows come from `logs/svv_off_2026-10-02/slab_table.txt`;
the $\phi$ $E_r$ were recomputed with each case's own $\varepsilon$
(`logs/audit_2026-10-09/slab_eval_all.txt`; the earlier `slab_rows_phi.txt` used $\varepsilon=0.01$
throughout).

| $N$ | $\psi$ $E_r$ | $\psi$ worst violation | $\phi$ $E_r$ | $\phi$ worst violation |
|---|---|---|---|---|
| 6 | 0.00027 | **0** | 0.722 | $8.7\times10^{-3}$ |
| 8 | 0.00017 | **0** | 1.194 | $2.7\times10^{-2}$ |
| 10 | 0.00006 | $1.2\times10^{-10}$ | 1.214 | $3.1\times10^{-2}$ |
| 12 | **0.00001** | $4.0\times10^{-6}$ | 1.200 | $3.4\times10^{-2}$ |

**$\psi$ converges 29× from $N=6$ to 12, while $\phi$ saturates.** With no geometry in play, the
$\phi$-normal's failure is a property of the operator: a higher order gives the grid-scale growth
more modes. `../zalesak_disk` shows the same in 2D.

### 3.4 $\psi$ without maintenance

`psi_xi10`, every frame, band $\phi(1-\phi)>10^{-4}$:

| quantity | value |
|---|---|
| band $\lvert\nabla\psi\rvert$, range | 0.65–1.48 |
| band $\lvert\nabla\psi\rvert$, mean | stays in $[0.994, 1.009]$ |
| $\max\lvert\psi-\psi_{\text{exact}}\rvert$, in the band / outside | 0.00067 / 0.0056 |
| $\max\lvert\phi-\phi_{\text{exact}}\rvert$ | 0.00084 |

$\psi$'s error is 8× larger outside the band, where the flux $-\phi(1-\phi)\mathbf n$ cannot carry
it to $\phi$.

### 3.5 The redistancing testbed

These cases run the coupled files' **committed events path**, which is not Saini's configuration
(`../../CDI_METHOD.md` §6):
- SSP-RK3, the uniform-$|\mathbf c|$ SVV and sign width $\varepsilon$;
- `band` 2.5, with `svv_psi` at $c_0=1$.

They are the cheap check for changes to that machinery. $\xi=1$, $\gamma=1$, every frame, with the
event's history restart (re-run 2026-10-05):

| $\psi$ at $t=0$ | build $\Delta\tau$ | events | $E_r(t{=}20)$ | worst violation |
|---|---|---|---|---|
| exact | — | off | **0.00006** | $1.2\times10^{-10}$ |
| built by Eq. (44) (`psi_rd_off_cfl`) | $0.1\,h_{\text{GLL,min}}$ | off | **0.00006** | $1.2\times10^{-10}$ |
| built (`psi_rd_on_cfl`) | $0.1\,h_{\text{GLL,min}}$ | 40, `seed="phi"` | **0.00006** | $1.4\times10^{-10}$ |
| built (`psi_rd_off`) | $H/(N{+}1)$ | off | 0.06458 | $2.6\times10^{-3}$ |
| built (`psi_rd_on`) | $H/(N{+}1)$ | 40 | **1.283** | $3.4\times10^{-3}$ ($E_v$ 0.40) |

- **A built $\psi$ costs nothing when the build is resolved.** It matches the analytic distance to
  five figures in $E_r$ and in boundedness, and it never evaluates an analytic distance.
- **$H/(N{+}1)$ fails here because of the sign width.**
  - With the $\varepsilon$-width sign function it is a pseudo-CFL of 2.755. The build is already
    wrong ($|\nabla\psi|$ 0.067–3.165 against 0.999–1.001), and 40 events drive it past $10^{30}$.
  - At $0.1\,h_{\text{GLL,min}}$ each event is nearly a no-op: $\lVert dn\rVert$ is $9\times10^{-8}$ at the
    first event and $7.5\times10^{-5}$ at the 40th.
  - At Saini's width 0.25 his $H/(N{+}1)$ is safe (`../rider_kothe/README.md` §3.5).
- **Neither half of Saini's operation costs anything on an isometry** (40 events, everything else
  fixed):

  | | $E_r$ | $E_v$ | worst violation |
  |---|---|---|---|
  | no events | 0.00006 | $1.7\times10^{-12}$ | $1.2\times10^{-10}$ |
  | `seed="psi"`, relax in place (`psi_rd_on_psi`) | 0.00006 | $-8.1\times10^{-10}$ | $1.2\times10^{-10}$ |
  | `seed="phi"`, Eq. (47) reseed, then relax | 0.00006 | $-3.1\times10^{-7}$ | $1.4\times10^{-10}$ |

  Under strain, and with corners, they do (`../../REDISTANCING.md` §3).

## 4. Limitations and open items

- **There is no strain**, and the normal carries only a sign. What this case picks out has to be
  checked on Rider–Kothe, where $\xi$ trades shape against boundedness.
- **SVV at $c_0=0.1$ raises $E_r$ at large $\xi$** (0.00005 at $\xi=2$, 0.00042 at 2.8; appendix).
  The mechanism is not measured.
- **The evidence predates `svv_psi`** (§6) and is to be regenerated (`../../NEXT_SESSION.md`).

## 5. Running

```bash
genmeshbox 0 1 0 0.1 0 0.1 10 1 1 .true. .true. .true.   # box.nmsh
./run.sh            # seven cases: phi/psi at xi 1 and 0.5, gamma 0 and 0.25, psi_xi28
./run.sh psi_xi15   # any other case by name
```

On the CPU each case gets one rank and they run together, about 13 min each.

## 6. Evidence

Made by `visualize.ipynb`, which ships executed. **These predate 2026-10-02:** the $\psi$-normal
panels are SVV-off runs (archived in `logs/svv_off_2026-10-02/`), and the notebook samples every 10th
frame, which aliases with the elements.
- **`slab_1d_methods.mp4`:** one panel per method at $\xi=1$: $\gamma=0$, $\phi$-normal, and
  $\psi$-normal without and with SVV. The exact slab is dashed.
- **`slab_1d_cdi_off_vs_phi.mp4`:** the first two of those panels.
- **`slab_1d_grad_psi.png`:** band $|\nabla\psi|$ against $t$ for both SVV states. The run is correct
  either way.
- **`slab_1d_redistancing.png`:** the testbed of §3.5. It shows $\psi$ at $t=20$, $E_r(t)$, and the
  band $|\nabla\psi|$, where the coarse $\Delta\tau$ climbs to $10^{30}$.
- **`slab_1d_psi_field.mp4`:** $\phi$, $\psi$ and $\psi$'s error ×500, SVV off and on. The error is
  grid-scale and spikes at the two kinks, outside the band.

## Appendix: other configurations

| run (source) | how it differs from Saini | numbers | note |
|---|---|---|---|
| SVV off, $\xi=1$ (`logs/svv_off_2026-10-02/`) | no `svv_psi` | $E_r$ 0.00006; violation $1.5\times10^{-3}$ (3872 nodes); band $\lvert\nabla\psi\rvert$ 0.001–3.6 | archived run of 2026-09 |
| SVV off, $\xi=2$ and 2.8 | no `svv_psi` | $E_r$ 0.00001 and 0.00003 (with $c_0=0.1$: 0.00005 and 0.00042) | |
| SVV off, $\xi=1.5$, $\gamma=0.25$ | no `svv_psi` | violation $1.1\times10^{-6}$ (0 with SVV) | |
| `psi_svvon` | `svv_psi` $c_0=1$, his phase-field value | $E_r$ 0.00006, violation $1.2\times10^{-10}$ (363 nodes), band $\lvert\nabla\psi\rvert$ 0.975–1.02 | |
