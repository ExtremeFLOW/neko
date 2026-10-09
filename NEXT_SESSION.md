# What is open

Settled results live in the case READMEs and the method docs, and standing rules in `CLAUDE.md`. When
an item here is finished, its result goes there and the item leaves this file.

**Records** (gitignored, local). These hold every prediction and verdict so far:
- `examples/rider_kothe/logs/d5/PREREGISTERED.txt`;
- `examples/rider_kothe/logs/vel_fix_2026-10-05/PREREGISTERED.txt`;
- `examples/zalesak_disk/logs/armC_2026-10-02/PREREGISTERED.txt`.

The measurement scripts are in `examples/rider_kothe/logs/`:
- `d5/rd_quality.py` and `d5/gamma_sweep.py`;
- `d5_2026-10-06/{full_measures,frag_first,phi_normal_exact}.py`;
- `d6_2026-10-07/measure.sh`;
- `figs/figs_rk.py` and `anim/anim_rk.py`.

Run them with `/lscratch/sieburgh/local/venv-neko/bin/python`.

## Rider–Kothe: re-initialization as Saini does it

The aim is to replicate Saini §4.5 with our CDI $\phi$ and his $\psi$ machinery, unchanged.
`examples/rider_kothe/README.md` §2 has his settings, and §3.4–§3.5 what is already known.

1. **Put Saini's configuration in the coupled files.**
   - **The checklist** is the "coupled files" column of `CDI_METHOD.md` §6. Every row that differs
     from his becomes his:
     - BDF2/EXT2 in pseudo-time with the printed Eq. (31) SVV (`svv_step_eq31`);
     - dealiased $\mathbf C(\mathbf w)\psi$ and his sign guard against $\phi-\tfrac12$;
     - sign width 0.25, extent $25H$, $\Delta\tau=H/(N{+}1)$;
     - local $|\mathbf u|$, left-multiplied, for `svv_psi`;
     - BDF2 in physical time. The paper prints BDF3/EXT3 at $H=1/128$, but his shipped
       $H=1/128$ case uses BDF2 too;
     - the $t=0$ build (`psi_init = "redistance"`) and events every 0.5.
   - **Most of it exists already** on the scratch file `examples/rider_kothe/logs/d6_2026-10-07/rider_kothe_d6.f90`
     Its diffs: `logs/d5_2026-10-06/rider_kothe_0c.diff` against the committed file, then
     `logs/d6_2026-10-07/rider_kothe_d6.diff` on top.
     - The keys: `case.cdi.redistance.{sgn_eps, band, scheme, svv_form, dealias, guard}` and
       `case.cdi.svv_psi_local`. `dealias` and `guard` need `scheme = "bdf2"`.
     - With every key at its default it reproduces the committed events run byte for byte in frames
       f00000–f00007 (to $t=0.56$; against `logs/events_2026-10-05/output_rk_events`). With `bdf2`
       and `eq31` set it reproduces `logs/d5_2026-10-06/output_bdf2_e2post` over all 14 frames.
     - Run B (`logs/d6_2026-10-07/runB_s025_25H_dg_local.case`: dealiasing, guard and local
       $|\mathbf c|$) is prepared but has never run.
   - **Caveats.**
     - **The settings land as the code, not as options** (`CLAUDE.md`). They change the
       shared-routine set, so the user decides, and the change is checked on `advecting_slab_1d`
       first.
     - **$\Delta\tau=H/(N{+}1)$ goes only with width 0.25.** At width $\varepsilon$ it fails
       (`examples/advecting_slab_1d/README.md` §3.5).
     - **Saini scales the sign width and the seed by $L$, the smallest domain extent.** Our meshes
       are 0.1 thick in $z$, so $L$ must be the in-plane extent (1 here).
     - **Neko's `time_order: 2` is BDF2 with a modified EXT3.** His code also extrapolates with
       Nek's EXT3 coefficients.
     - **The Rider–Kothe `.case` files have no `redistance` block**, and `svv_rd`'s $c_0$ defaults to
       0 (off). Adding the block is a committed-`.case` change, so ask first.
     - **Local $|\mathbf c|$ must left-multiply the assembled operator.** Never make $\nu$ vary
       inside `svv_local`.
     - **BDF2 at $N_{svv}=N/4$ once lost an apex that RK3 held** (§4.4, $H=1/5$, $N=7$; archive §7.5).
       §4.5 uses $N/4$ for Eq. (44), so watch for it.
2. **If the rebuilt $\psi$ is not good enough**, check in this order:
   - **(a) Pseudo-time.** The extent must be $25H$ at width 0.25. The step count over that extent
     barely matters (150 or 2129 steps give 0.0407 and 0.0402).
   - **(b) Re-initialization frequency.** It has never been varied. Try $\Delta t_{tls}=1.0$, a single
     event at 0.5 and none after, and 0.25.
     - A single event needs a new key (e.g. `redistance.max_events`); nothing limits the count today.
     - Measure every frame with `rd_quality.py`, and compare with transport (0.0452) and the
       Saini-configuration run (0.0402).
   - **(c) Dealiasing against the guard.** They have only ever run as a pair, which costs about 18%.
3. **Replicate his Figs. 17–18:** $E_r$ (Fig. 17), and the $L^1$ norm, $|E_v|$ and $E_s$ (Fig. 18),
   at $t=8$ over $H\in\{1/32,1/64,1/128\}\times N\in\{4,5,6\}$. Expect a shallower slope with re-initialization
   (his p. 21).
4. **A convergence measure for the $\tau$ solve:** a residual on the compression band, not the build
   band.
5. **Housekeeping.**
   - `rider_h192.case` and `zalesak_h150.case` ship but have never run: run or remove them.
   - `rider_kothe_xi10` ends at $t=8.00008$, one step past the reversal: decide whether to stop
     exactly at 8.
   - When item 1 lands, mark `examples/rider_kothe/logs/d6_2026-10-07/HANDOFF.md` as done.

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
  - Is a small violation acceptable for a better shape? On Zalesak, $\xi=1$, $\gamma=0.5$ gives
    $E_r$ 0.00240 at $5.9\times10^{-6}$.
- **The $\psi$ build at $N=3$ stops short of a distance** (24% residual, Zalesak). Saini's Table 3
  compares $N=3$ with $N=7$. An untested cause: `svv_rd` at $N_{svv}=0.75$ gives mode 1 a weight of
  0.44.
- **A possible second email to the authors:**
  - the guard's round-off sensitivity at low $N$;
  - his reversed with/without-guard pair;
  - $25H$ against the printed $2.5H$;
  - the outside-area normalisation of `ls_relerr`.
