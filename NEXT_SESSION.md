# What is open

Settled results live in the case READMEs and the method docs, and standing rules in `CLAUDE.md`. When
an item here is finished, its result goes there and the item leaves this file.

**Records** (gitignored, local). These hold every prediction and verdict so far:
- `examples/rider_kothe/logs/s45_2026-10-12/PREREGISTERED.txt` (Saini's configuration, the knobs, his grid);
- `examples/rider_kothe/logs/d5/PREREGISTERED.txt`;
- `examples/rider_kothe/logs/vel_fix_2026-10-05/PREREGISTERED.txt`;
- `examples/zalesak_disk/logs/armC_2026-10-02/PREREGISTERED.txt`.

The measurement scripts are in `examples/rider_kothe/logs/`:
- `d5/rd_quality.py` and `d5/gamma_sweep.py`;
- `d5_2026-10-06/{full_measures,frag_first,phi_normal_exact}.py`;
- `d6_2026-10-07/measure.sh`;
- `s45_2026-10-12/`:
  - `measure.sh` (every measure on every frame) and `summarize.py` (one table per batch);
  - `es_t.py` ($E_s$ at every frame against the exact interface) and `probe_quality.py`;
  - `figs_s45.py` (the three metrics, $|E_v|(t)$, contours with his code);
- `d5/saini/mk_saini_grid.sh` (his code over his grid);
- `figs/figs_rk.py` and `anim/anim_rk.py` (which now also reads his dumps, `saini=<dir>`).

Run them with `/lscratch/sieburgh/local/venv-neko/bin/python`.

## Rider–Kothe: re-initialization as Saini does it

The aim is to replicate Saini §4.5 with our CDI $\phi$ and his $\psi$ machinery, unchanged.
`examples/rider_kothe/README.md` §2 has his settings, and §3.4–§3.5 what is already known.

1. **Put Saini's configuration in the coupled files.**
   - **The checklist** is the "coupled files" column of `CDI_METHOD.md` §6. Every row that differs
     from his becomes his:
     - BDF2/EXT2 in pseudo-time with the printed Eq. (31) SVV (`svv_step_eq31`);
     - dealiased $\mathbf C(\mathbf w)\psi$ and his sign guard against $\phi-\tfrac12$, always on;
     - sign width 0.25, extent $25H$, $\Delta\tau=H/(N{+}1)$, written to full precision. A $\Delta\tau$
       a hair below $H/(N{+}1)$ rounds the step count up by one;
     - local $|\mathbf u|$, left-multiplied, for `svv_psi`;
     - BDF2 in physical time. The paper prints BDF3/EXT3 at $H=1/128$, but his shipped
       $H=1/128$ case uses BDF2 too;
     - the $t=0$ build (`psi_init = "redistance"`) and events every 0.5.
   - **It all runs on the scratch file** `examples/rider_kothe/logs/s45_2026-10-12/rider_kothe_s45.f90`
     (results in the case README §3.4–§3.5).
     - It is the D6 file plus a `redistance.max_events` key (`rider_kothe_s45.diff`). Its diffs against
       the committed file: `logs/d5_2026-10-06/rider_kothe_0c.diff`, then
       `logs/d6_2026-10-07/rider_kothe_d6.diff`.
     - The keys: `case.cdi.redistance.{sgn_eps, band, scheme, svv_form, dealias, guard, dtau}` and
       `case.cdi.svv_psi_local`.
     - With every key at its default it reproduces the committed events run byte for byte in frames
       f00000–f00007.
     - EXT3 in pseudo-time (`s45_2026-10-12/ext3/`) gives the same $\psi$ to 0.05° through $t=6$, and
       $t=8$ errors at the edge of the noise band. Either extrapolation may be ported.
   - **Decide with the user:**
     - $\gamma$ for the Rider–Kothe cases. In Saini's configuration $\gamma=0.5$ beats 1 by about 30%
       in $E_r$ and $E_s$ at $H=1/64$ and 1/128.
     - `svv_psi` explicit or implicit. His is implicit.
   - **Caveats.**
     - **The settings land as the code, not as options** (`CLAUDE.md`). They change the
       shared-routine set, so the user decides, and the change is checked on `advecting_slab_1d`
       first.
     - **$\Delta\tau=H/(N{+}1)$ goes only with width 0.25.** At width $\varepsilon$ it fails
       (`examples/advecting_slab_1d/README.md` §3.5).
     - **Saini scales the sign width and the seed by $L$, the smallest domain extent.** Our meshes
       are 0.1 thick in $z$, so $L$ must be the in-plane extent (1 here).
     - **The Rider–Kothe `.case` files have no `redistance` block**, and `svv_rd`'s $c_0$ defaults to
       0 (off). Adding the block is a committed-`.case` change, so ask first.
     - **Local $|\mathbf c|$ must left-multiply the assembled operator.** Never make $\nu$ vary
       inside `svv_local`.
     - **BDF2 at $N_{svv}=N/4$ once lost an apex that RK3 held** (§4.4, $H=1/5$, $N=7$; archive §7.5).
       §4.5 uses $N/4$ for Eq. (44), so watch for it.
2. **The shape gap to Saini.** His $E_s$ is 2–2.8× lower than ours at $H=1/128$.
   - **Start here (the user's observation).** At $H=1/128$ our $\phi$ does not capture the thin tail;
     his does.
     - At $t=4$ our 0.5 contour stops $5.7H$ ($\gamma=0.5$) and $7.8H$ ($\gamma=1$) short of the exact
       tip; his stops $0.8H$ short.
     - Sources: `examples/rider_kothe/logs/s45_2026-10-12/tables/tail_tip_t4_h128.txt` and
       `examples/rider_kothe/evidence/rider_saini45_tail_h128.png`.
     - The missing tail may be what drives our late $E_s$.
     - The user has ideas for a fix; think them through together before running anything.
     - The $H=1/128$ outputs are kept in `logs/s45_2026-10-12/output_{s2_h128,s8_g05_h128,gd128n*_g05}`.
   - Part of the gap is there before the first event, so the phase field itself contributes.
   - Part is the event coupling: at $t=2$–4 the rebuilt $\psi$ is worse than its seed, and CDI reads
     it at once. $\gamma$ already moves this.
   - **Open:** what else on the $\phi$ side reduces it, without SVV on $\phi$.
3. **$\Delta t_{tls}=0.25$ diverges** at $t=5.69$ after 22 events, inside his range of 0.01–0.5
   (case README appendix). Find the mechanism of this repeated-event instability before any interval
   below 0.5 is used.
4. **$N=6$ in Saini's configuration loses boundedness** (case README §3.7).
   - At $H=1/64$ and 1/128 it is no better than $N=5$; at $H=1/32$ it diverges.
   - No single change fixes it at $H=1/32$: EXT3, $\gamma=1$, half the $\Delta t$, and the Eq. (44) SVV
     at $N/6$ (which only delays it). The explicit terms are far inside their stability limits.
   - In every run the violation passes $10^{-3}$ at $t=1.28$–1.36. Two late frames put the worst values
     within about 1.5$H$ of the interface; the onset frames have not been located.
   - Not yet looked at:
     - how CDI's compression term behaves with a sixth-order $\psi$ normal;
     - the $\psi$ build at $N=6$;
     - whether the first event (at $t=0.5$, $1.0$) is where it starts.
   - Records: `s45_2026-10-12/PREREGISTERED.txt` ("N = 6 diagnosis").
   - Done:
     - his whole grid (§3.6);
     - ours at $\gamma=0.5$ (§3.7, figure `s45_2026-10-12/figs/s45_grid.png`).
   - Still to do:
     - Fig. 16's line plot ($H=1/128$, $N=4$–6);
     - his Table 3 pairs, $(1/64,N{=}3)$ against $(1/32,N{=}7)$ and $(1/128,N{=}3)$ against $(1/64,N{=}7)$.
5. **A convergence measure for the $\tau$ solve:** a residual on the compression band, not the build
   band.
6. **Housekeeping.**
   - `rider_h192.case` and `zalesak_h150.case` ship but have never run: run or remove them.
   - `rider_kothe_xi10` ends at $t=8.00008$, one step past the reversal: decide whether to stop
     exactly at 8.

## Parked

- **Zalesak.**
  - `svv_psi` $c_0$ is 1.0, Saini's phase-field value; his $\psi$ value is 0.1. Moving it means
    re-running the tables.
  - Re-run the two redistancing variants once the coupled files carry Saini's path. The old arm C
    used the committed path.
  - Should `zalesak.case` move from $\xi=2.8$ to 1.5?
- **Audit every comparison against Saini's $E_r$** for his outside-area denominator. That means
  Zalesak's §4.3 numbers; `redistance_circles` is likely unaffected.
- **Regenerate the slab and Zalesak evidence.** The slab's $\psi$ panels are SVV-off runs, and its
  notebook samples every 10th frame, which aliases with the elements. Zalesak needs the style only.
- **The compression guard's 0.05** is not a measured stability boundary (`CDI_METHOD.md` §2).
- **The CDI side.**
  - Queued: ten rotations on $(\xi,\gamma)=(1.5,0.25)$ and $(2,0.25)$ against $(2.8,1)$, about 4 h.
  - Run $\xi=0.75$ under strain, to see whether the Rider–Kothe cross turns over.
  - The SVV cost at large $\xi$ in 1D: the mechanism is not measured.
- **The $\psi$ build at $N=3$ stops short of a distance** (24% residual, Zalesak). Saini's Table 3
  compares $N=3$ with $N=7$. An untested cause: `svv_rd` at $N_{svv}=0.75$ gives mode 1 a weight of
  0.44.
- **A possible second email to the authors:**
  - the guard's round-off sensitivity at low $N$;
  - his reversed with/without-guard pair;
  - $25H$ against the printed $2.5H$;
  - the outside-area normalisation of `ls_relerr`.
