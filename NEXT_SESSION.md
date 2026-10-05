# What is open

The focus is redistancing on Rider–Kothe: whether a $\psi$ maintained by Eq. (44) can give the
compression term a better normal than the transported one. Settled results live in the case
READMEs and `CDI_METHOD.md`, standing rules in `CLAUDE.md`. When an item is finished, its result
goes there and the item leaves this file. The records with every prediction and verdict so far are
`examples/rider_kothe/logs/d5/PREREGISTERED.txt` and
`examples/zalesak_disk/logs/armC_2026-10-02/PREREGISTERED.txt` (gitignored, local).

## Rider–Kothe redistancing

1. **Fix the velocity and $u_{\max}$ lag.**
   - `rider_kothe.f90` (`compute`: the `field_cmult2` of `u0`/`v0`, and `u_max`) and
     `zalesak_disk.f90` (the velocity, set once in `compute`) prescribe the velocity in `compute`.
     Neko calls `compute` after the scalar step (`src/simulation.f90:179`; `preprocess` is called
     at `:150`).
   - So step 1 of both runs with $u=0$, and Rider–Kothe advects with $u(t_{n-1})$. The compression
     and diffusion use $u_{\max}(t_{n-1})$, consistently, so the equilibrium width is unaffected.
   - Fix: move it to the `preprocess` hook. Ask before the Fortran change (`CLAUDE.md`).
   - This changes results. Measure it on `rider_kothe_xi10` first, then decide whether the
     Rider–Kothe tables need re-running.
2. **Run the committed events path on Rider–Kothe.** It has not run since the 2026-10-02 fixes.
   - Configuration: `psi_init = "redistance"`, SSP-RK3, sign-function $\varepsilon$ equal to the
     phase field's, band $2.5H$, `cfl` 0.1, $\Delta t_{tls}=0.5$.
   - One run at $H=1/64$, $N=5$, $\xi=1$, $\gamma=1$. Pre-register it.
   - Compare against `rk_Af` (transport only: $E_r$ 0.0446, normal 0.6–0.8° off;
     `examples/rider_kothe/README.md`).
   - Judge by the normal angle in the compression band (`examples/rider_kothe/logs/d5/rd_quality.py`, gitignored) and by
     $E_r(8)$.
3. **D5: one knob at a time**, from the two measured causes in `CDI_METHOD.md` §4.1d.
   1. The sign-function width: 0.25 against $H/N$-scaled.
   2. Dealiased $\mathbf C(\mathbf w)$ plus the sign guard.
   3. D2: BDF2/EXT2 with the printed Eq. (31) SVV (`svv_step_eq31`) in place of SSP-RK3 with
      `svv_step_imp`. BDF2 can fail at a steep apex near a medial axis where RK3's Lie split held
      (`examples/redistance_circles/archive/README_process_2026-09.md` §7.5).
4. **The thin tail.** Where the filament is thinner than $2\varepsilon$ (5–6% of its length at
   maximum stretch), $\phi$ has no 0.5 contour, so no reseed setting can rebuild $\psi$ there. Each
   candidate is a change to the method, not a knob:
   - a compression flux masked to the interface band;
   - the monotone transform $\psi \leftarrow L\tanh(\psi/L)$ (`REDISTANCING.md` §8).
5. **A convergence measure for the $\tau$ solve:** a residual on the compression band, not the
   build band (`CDI_METHOD.md` §4.3).
6. **Housekeeping.**
   - Regenerate `examples/rider_kothe/evidence/`: all of it predates the fixes. Use kthviz style
     and put a GIF beside every MP4. Working script: `examples/rider_kothe/logs/anim/anim_rk.py`
     (gitignored, local).
   - `rider_h192.case` ships but has never run: run it or remove it.
   - Run $\xi=0.75$ under strain (~30 min) to see whether the cross turns over.

## Saini §4.5, for reference

From the paper (rendered p. 21, checked 2026-10-05). Their field names are transposed: their
$\phi$ is our $\psi$ (`CDI_METHOD.md` §1).

| | Saini §4.5 |
|---|---|
| velocity, time | Eqs. (85)–(86), $T=8$ |
| domain, disk | $\Omega=[0,1]^2$, centre $(0.5,0.75)$, $r=0.15$ |
| $\xi$ | their $1/N$, i.e. our $\xi=1$ |
| mesh, order | $H\in\{1/32,1/64,1/128\}$, $N\in\{4,5,6\}$ |
| $\Delta t$ | $\{8,4,2\}\times10^{-4}$ (CFL $\approx0.4$ at $N=6$); BDF3/EXT3 at $H=1/128$ |
| SVV, CLS and TLS advection, CLS re-initialization Eq. (39) | $N_{svv}=N/2$, $c_0=0.1$ |
| SVV, TLS re-distancing Eq. (44) | $N_{svv}=N/4$, $c_0=2.0$ |
| $\Delta t_{tls}$ | 0.5, "including at the initial step" |
| $\Delta t_{cls}$ | 0.05; no counterpart here (CDI is fused) |
| re-distancing extent | $2.5H$ from the interface (his public code uses $25H$) |
| pseudo-CFL | $\approx0.24$ at $N=6$ |

What his figures measure:
- Figs. 17–18: $E_r$, $|E_v|$, $E_s$ at $t=8$ under $h$- and $p$-refinement.
- Fig. 19: mean interface thickness $l_{avg}/l_0$ against $t$ (Eq. 87).
- Fig. 20: $L_\infty$ of the CLS (boundedness) against $t$.
- Fig. 21: the 0.5 isocontours at $t=8$ against the exact one.
- Fig. 22: $|E_v|$ against $t$.
- Table 3: pairs at equal GLL count, $N=3$ against $N=7$.

His $E_r$ divides by the area outside the disk (`CLAUDE.md`); recompute it with `norms.E_r`
before comparing. He remarks (p. 21) that TLS re-distancing "slows down the convergence rate of the
coupled algorithm", so a shallower slope with redistancing is expected.

## Parked (not Rider–Kothe)

- **Zalesak $\phi$-contour fragmentation under reseeding:** find where the first spurious piece
  appears, on the existing arm C outputs at $t=0.5$–1; the slot corners are the suspect
  (`CDI_METHOD.md` §4.1c).
- **Zalesak `svv_psi` $c_0$:** 1.0 is Saini's CLS value, 0.1 his TLS value; moving it means
  re-running Zalesak's tables.
- **SVV cost at large $\xi$ in 1D:** $E_r$ 0.00001 → 0.00005 at $\xi=2$ and 0.00003 → 0.00042 at
  $\xi=2.8$ (`examples/advecting_slab_1d/README.md`); the mechanism is not measured.
- **Regenerate the slab and Zalesak evidence.** The slab's $\psi$ panels are archived SVV-off runs,
  and its notebook samples every 10th frame, which aliases with the elements: sample every frame.
  Zalesak needs the style only.
- **Re-run the two Zalesak redistancing variants**, now `psi_init = "redistance"` with the timer.
  Their README numbers are from the analytic-$\psi$ and grad-trigger configuration and are not
  citable. Note that $\xi=2.8$ at $N=5$ violates $N\ge3.68\xi$ (`CDI_METHOD.md` §4.1b).
- **Audit every comparison against Saini's $E_r$** for the outside-area denominator: the Zalesak
  §4.3 numbers; `redistance_circles` is likely unaffected.
- **The compression guard's 0.05 is not a measured boundary** (`CDI_METHOD.md` §2).
- **Queued run:** ten rotations on $(\xi,\gamma)=(1.5,0.25)$ and $(2,0.25)$ against $(2.8,1)$
  (~4 h).
- **Open questions:** is a small boundedness violation acceptable for better shape ($\xi=1$,
  $\gamma=0.5$: $E_r$ 0.00240 at a $5.9\times10^{-6}$ violation, `examples/zalesak_disk/README.md`)?
  Should `zalesak.case` stay at $\xi=2.8$?
- **Why the build does not converge at $N=3$** (24% residual). Untested cause: `svv_rd` at
  $N_{svv}=0.75$ gives mode 1 a kernel weight of 0.44.
- **A possible second email to the authors:** the guard's round-off sensitivity at low $N$; his
  reversed with/without-guard pair; $25H$ against the printed $2.5H$; the outside-area
  normalization of `ls_relerr`.
