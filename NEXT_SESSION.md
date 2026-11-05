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
  - `gamma/`: the $\gamma$ table, the doc check and the observed orders (item 2);
  - `measure.sh` (every measure on every frame) and `summarize.py` (one table per batch);
  - `es_t.py` ($E_s$ at every frame against the exact interface) and `probe_quality.py`;
  - `figs_s45.py` (the three metrics, $|E_v|(t)$, contours with his code);
  - the rate study (case README §3.8): `drift_split.py`, `persist.py` with `persist_stats.py`,
    `mass_disp.py`, `mass_fine.py` with `mass_tests.py`, `hb_predict.py`, and `rate_tables.py`, which
    prints the window rates, drift and retention;
- `d5/saini/mk_saini_grid.sh` (his code over his grid) and `d5/saini/mk_s45_rate.sh` (his repeats);
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
2. **The shape gap to Saini.** His $E_s$ falls 7.1–8.4× per doubling over his repeats, ours 3.5–4×.
   At $H=1/128$ his is 2–2.8× lower.
   - **What is known** (case README §3.4 and §3.8).
     - Per interval, $\phi$ moves away from $\psi$'s zero set by about the same amount in both codes, at the same
       rate.
     - Ours keeps about 70% of the earlier error at every $H$ (less at lower $\gamma$). His keeps 45%, then
       23%.
     - The error patterns are deterministic, and ours locks in more with refinement.
     - At $H=1/128$ our $t\le2$ error pattern is still there at the mirrored times $t=6$–8 (correlation
       +0.70); his much less (+0.27). At $H=1/64$ both keep it (+0.39, +0.32).
     - The tail tip is a small part of $E_s$.
     - In both codes the non-material mass flux keeps its sign through a phase of the flow and reverses
       with it; the mass moved by the late peak ($t=6.48$) is the same at $H=1/128$.
   - **Ruled out** (§3.8): $\psi$'s quality, the curvature–strain drift, noise, a strain-tied deviation,
     his error being profile shape at $t\ge7.4$, and his flux flipping sign from one interval to the next.
     At $t=8$ both codes' contour error is mostly mass moved non-materially.
   - **Tried and stopped: $\phi$'s diffusion form** (case README appendix and §3.8; records in
     `logs/s45_2026-10-12/PREREGISTERED.txt`, "thin tail", steps 1–4, and "Band-projected diffusion: E_s
     rate"; code and tables in
     `logs/s45_2026-10-12/cdi_forms/`).
     - The projected Eq. (39) diffusion matches his $E_s$ to $t\approx1$ and lowers retention from 0.67 to
       0.56 (means over the events at $t=1$–6). But $\phi$ then leaves $[0,1]$ (or diverges), and $t=8$ is
       worse.
     - Projected only inside the band: the same. Its rate is ours, not his: window rates 3.8–4.2× from
       $H=1/64$ to $1/128$, and $E_s(8)$ 4.65× only because its $H=1/64$ value is 51% worse than isotropic.
     - SVV on $\phi$ at his strength: no effect.
     - Untried: an isotropic blend $(1-\alpha)\mathbf n\mathbf n+\alpha\mathbf I$.
   - **Open: why his irreversible remainder falls with $H$ and ours does not.** The remainder is the part
     of the moved mass not undone by $t=8$, $\sum|M(8)|/\sum|M(6.48)|$ (`mass_fine.py`, `mass_tests.py`).
     - His: 0.42 at $H=1/64$, 0.23 at $1/128$. Ours at $\gamma=0.5$: 0.34, then 0.32; at $\gamma=1$: 0.46,
       then 0.53. At $H=1/64$ ours is the smaller at $\gamma\le0.5$.
     - It grows with $\gamma$ at $H=1/64$: 0.27 at 0.25, 0.54 at 2.
     - Look at the remainder, not the forward motion. `massfine_*_s32.npz` holds only $t\le1.44$ and
       $t\ge6.48$, so it cannot say when the irreversible part forms. Running `mass_fine.py` at
       $t=2$–6 needs its sub-point check (c2 in `PREREGISTERED.txt`) to pass there first; the carried
       profile is most compressed mid-run.
   - **$\gamma$: where it stands.** The runs of the $\gamma\propto\varepsilon$ idea are rows in the case
     README appendix: $\gamma=0.25$ isotropic and band-projected, his dose, $N=4$ and 6, $\gamma\propto\sqrt\varepsilon$.
     - **Which $\Gamma$ is the user's call.** With $\beta=(\Gamma/\varepsilon)/S_{\max}$ and $S_{\max}=\pi u_{\max}$
       on Rider–Kothe:
       - $\Gamma=\gamma u_{\max}$ gives $\beta=\gamma N/(\pi\xi H)$;
       - $\Gamma=\gamma(\varepsilon/H)u_{\max}$ gives $\beta=\gamma/(\pi H)$, with no $N$ or $\xi$;
       - $\Gamma=\beta\varepsilon S_{\max}$ holds $\beta$ fixed.
       At fixed $N$ and $H$ the first two are the same $(\xi,\gamma)$ plane in two coordinates. Only a second
       $H$ separates them, by which contours stay put. Rider–Kothe cannot separate a $u_{\max}$ form from a
       strain form (its global $S_{\max}$ is exactly $\pi u_{\max}$), nor the first two over $N=4$–6 at $H=1/64$.
     - **$\phi$'s solver tolerance: it is the $z$-variation, not the in-plane error.** At absolute $10^{-6}$
       $\phi$'s CG/Jacobi stops after 2 iterations. The Jacobi diagonal differs between $z$-face and interior
       nodes, so each step leaves a $z$-dependent error of about $10^{-5}$ (test `zt_tol{6,12}` in
       `logs/s45_2026-10-12/`), which the tail amplifies (0.42 in `s12_g0177_h128`).
       - Rerun at $10^{-12}$ to $t=8$ (`s12t_g0177_h128`): $z$-spread $4\times10^{-5}$; on one $z$-layer $E_r$
         is unchanged and $E_s(8)$ moves 0.5% (windows at most 2.6%); $\lvert E_v\rvert$ falls from $10^{-6}$
         to $10^{-11}$; 56% more wall time.
       - So all-layer numbers carry the $z$-spread: +15% in `s12`'s $E_r$ and $E_s(8)$, within 5% in the
         other $\gamma$ runs. `es_t.py` reads $z$-layer 0; `full_measures.py` and `measure_gamma.py` read
         all layers. Changing the committed cases' tolerance is the user's call.
     - **The measurement** is in `logs/s45_2026-10-12/gamma/`:
       - `runs.py` registers every run (families, repeats); its `.case` gives $\gamma$, $\varepsilon$, $\Delta t$.
       - `measure_gamma.py` writes one table, `gamma_table.{txt,csv}`: $E_r$, $E_s$, $|E_v|$, the worst
         violation, $l/l_0$, the $\phi$ and $\psi$ pieces, the remainder, and $E_r$ split into position and
         shape.
       - `doc_claims.py` checks the 118 $\gamma$ numbers in the docs against the dumps; all reproduce.
       - `fit_orders.py` gives weighted observed orders, with noise from the repeats, plus local orders
         (`orders.txt`). The noise is measured at one cell only ($H=1/64$, $\gamma=1$: 3.7% in $E_r$, 2.7% in
         $E_s$), so the orders there are not citable yet.
     - **Planned, not launched or preregistered.** `gen_gamma_cases.py` generates every case; about 15.5 h of
       GPU in all.
       - A third level at $H=1/32$ for each family.
       - Repeats ($\Delta t\times0.99875$, $r_f=0.1001$) at all three $H$.
       - $\xi=0.75$ and 1.5 at $H=1/64$, under the first two $\Gamma$.
       - Mirjalili, Ivey & Mani's shear drop ($T=4$, $\xi=1$, $\gamma=1$). Its scratch binary
         `mirjalili/neko` reads `case.cdi.rk_T`; its bit-identity check against `output_s1_h64` has not run.
       - Saini's $1/32$ repeats (`d5/saini/mk_s45_rate32.sh`).
       - Preregister in `PREREGISTERED.txt` after the solver gate.
   - **Outputs kept.**
     - The $H=1/128$ outputs: `logs/s45_2026-10-12/output_{s2_h128,s8_g05_h128,gd128n*_g05}` and
       `cdi_forms/output_cdi_*_h128`, plus `cdi_forms/output_cdi_proj_band_h64`.
     - His dumps every 0.04 at both $H$: `logs/d5/saini/cv{64,128}n5_d04`.
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
   - Undocumented scheme facts in `rider_kothe.f90`.
     - The conductivity is set at $t_{n+1}$ (line 887) while the compression uses $u_{\max}$ at $t_n$.
     - Every event restarts $\phi$'s BDF/EXT too, through the shared `ext_bdf` (line 359).
     - `rider_kothe_s45.f90`'s header says $2.5(N+1)$ pseudo-steps; it runs $25(N+1)$.
     Document them or align them.

## CDI's knobs: a $(\xi,\gamma)$ map at one $N$

The user's idea: keep $\Gamma=\gamma u_{\max}$ and map boundedness and $E_s$ over $(\xi,\gamma)$ at one $N$.
It is for a new session, and not on Rider–Kothe.
- **Order of work.**
  1. Pass the $\phi$ solver-tolerance gate (item 2).
  2. Find the exact steady profile of a flat interface under uniform normal strain, squeezing and
     stretching. It is a 1D ODE in $\beta$: solve it with `solve_bvp`, and check it with an independent
     1D solver.
  3. Find the discrete maximum-principle region of the 1D operator at the chosen $N$, from its sign
     structure.
  4. Map a log-spaced $(\xi,\gamma)$ grid on a strained test with a known answer, with repeats near the
     decision boundary.
  5. Repeat the decision cells at a second $H$, to see which contours collapse: constant $\gamma$ or
     constant $\beta$.
- **Preregister** the predictions and the decision rule first.
- **External checks:** Mirjalili, Ivey & Mani's boundedness derivation and order (JCP 401, 2020); mass to
  round-off; $z$-uniformity.
- **Not usable as a test of $\gamma$:** the slab and Zalesak are isometries, so they carry no strain.

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
  - The SVV cost at large $\xi$ in 1D: the mechanism is not measured.
- **The $\psi$ build at $N=3$ stops short of a distance** (24% residual, Zalesak). Saini's Table 3
  compares $N=3$ with $N=7$. An untested cause: `svv_rd` at $N_{svv}=0.75$ gives mode 1 a weight of
  0.44.
- **A possible second email to the authors:**
  - the guard's round-off sensitivity at low $N$;
  - his reversed with/without-guard pair;
  - $25H$ against the printed $2.5H$;
  - the outside-area normalisation of `ls_relerr`.
