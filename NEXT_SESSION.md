# What is open

The plan for `advecting_slab_1d`, `zalesak_disk` and `rider_kothe`, plus the doc and
code fixes left from the 2026-10-05 review (last section). The Rider–Kothe details are in
`NEXT_SESSION_RIDER_KOTHE.md`. Settled results live in
the case READMEs and `CDI_METHOD.md`, standing rules in `CLAUDE.md`. When an item
is finished, its result goes there and the item leaves this file.

## Open — D5: redistancing in the coupled cases (updated 2026-10-05)

Settled results:
- `examples/rider_kothe/README.md` and `examples/advecting_slab_1d/README.md`;
- `CDI_METHOD.md` §2 (how it is stepped), §4–§4.3 and §5.

Records, with every prediction and verdict:
- `examples/rider_kothe/logs/d5/PREREGISTERED.txt`;
- `examples/zalesak_disk/logs/armC_2026-10-02/PREREGISTERED.txt`.

**Open, in order:**
1. **Why $\phi$'s 0.5 contour fragments under periodic reseeding on Zalesak.**
   - At $N=5$ it is in 2 pieces at $t=1$ and 8 at $t=20$. $N=3$ and $N=7$ diverge (§4.1c).
   - Find where the first spurious piece appears, on the existing arm C outputs, frames
     $t=0.5$–1.
   - Candidates: the slot's corners; the $r_f(\phi-\tfrac12)$ seed there; the $2.5H$ build's plateau.
   - Measure before running anything.
2. **Can a reseeded $\psi$ ever beat the transported one on Rider–Kothe?**
   - Saini's configuration completes but points 4–8° off (transport: 0.6–0.8°). Two measured causes:
     - His sign-function width 0.25 keeps the solve from reaching a distance in the compression
       band.
     - The thin tail ($d<2\varepsilon$, 5–6% of the length at maximum stretch) has no $\phi$ contour
       to rebuild from. No redistancing setting fixes that.
   - Cheapest next test: the *committed* events path (SSP-RK3, sign-function width equal to the
     phase field's, $2.5H$, `cfl` 0.1), not run on Rider–Kothe since both fixes. One run;
     pre-register.
   - Until then, redistancing stays off in every validated configuration.
3. **Zalesak's `svv_psi` $c_0$.** 1.0 is Saini's §4.3 value for his CLS (our $\phi$); his TLS value
   is 0.1. Decide whether to move Zalesak; that means re-running its tables.
4. **SVV on $\psi$ costs shape at large $\xi$ in 1D:** 7× at $\xi=2$, 16× at $\xi=2.8$ (unrounded; slab
   README).
   - The mechanism is not measured. The hypothesis: the compression band ($9.21\varepsilon$)
     reaches $\psi$'s medial-axis kinks, and SVV smooths them.
   - Zalesak's primary is $\xi=2.8$ at $c_0=1$, where SVV was a large gain, so check there before
     generalising.
5. **Regenerate `evidence/` from the notebooks, in kthviz style, with a GIF beside every MP4.**
   - **`rider_kothe`:** all of its evidence predates the fixes (frozen diffusion, no SVV).
   - **`advecting_slab_1d`:** the $\psi$ panels are the archived SVV-off runs. Its notebook
     samples violations every 10th frame, which aliases with the element spacing, so make it every
     frame.
   - **`zalesak_disk`:** style only.
   - Working animation script: `examples/rider_kothe/logs/anim/anim_rk.py` (gitignored, local).
6. **Re-run the shipped Zalesak redistancing variants** (`zalesak_redistance_phi`/`_psi`,
   $\xi=0.5$). Both crashed with the history bug, and the README and §7 still quote that.
7. **Audit every comparison against Saini's $E_r$** for the outside-area denominator: Zalesak
   (§4.3 numbers quoted anywhere) and `redistance_circles` (likely unaffected: a different case and
   script).
8. **The compression guard's 0.05 is not a measured boundary** (`CDI_METHOD.md` §2). Measuring where
   the explicit compression actually goes unstable is what could buy a larger $\Delta t$.

Files (gitignored):
- `examples/rider_kothe/logs/`: `d5/`, `anim/`, `frozen_lambda_2026-10-01/`,
  `lambda_fix_nosvv_2026-10-02/`;
- `examples/zalesak_disk/logs/armC_2026-10-02/`;
- `examples/advecting_slab_1d/logs/svv_off_2026-10-02/`.

**D2: whether the coupled cases move to the printed Eq. (31) SVV** (`svv_step_eq31`).
At the earlier sign-function $\varepsilon=H/N$, $N_{svv}=N/4$, $H{=}1/5$, $N{=}7$, BDF2 failed at the
§4.4 apex at every $\Delta\tau$ tried, and RK3's Lie split held it only at $\Delta\tau\ge2.5\times10^{-4}$ (`examples/redistance_circles/archive/README_process_2026-09.md` §7.5), so
moving a coupled case to BDF2 could expose it where its seed is steep near a medial axis.
On Rider–Kothe (D5 rung c) it changed nothing visible, but that run had the history bug. The history
restart is now in the code, so D2 can be decided on runs without it.

**The authors.** The ε question is answered: the fixed sign-function $\varepsilon$ is intentional
(his reply, 2026-10-01). A possible second email, not sent, would cover:
- the guard's round-off sensitivity at low $N$;
- his reversed with/without-guard pair;
- the $25H$ budget against the printed $2.5H$;
- that his `ls_relerr` normalizes by the area outside the disk in circVortex (is Eq. (79) meant that way?).

## Open — the $\phi$-normal re-measurement

The `grad_floor` bug (`CDI_METHOD.md` §4.1) gave a random unit normal wherever
the normal's source field is flat, and the $\phi$ field is flat over most of the
domain. The $\phi$-normal baselines ran before the fix (`zalesak_phi_normal` on
2026-08-29; the slab's `phi_xi10` log was rewritten on 2026-09-08, probably still
pre-fix), so the headline "49×" in `examples/zalesak_disk/README.md` is the
pre-fix number: quote its direction, not its factor, until re-run. The seven
cases: `examples/advecting_slab_1d/phi_{xi05,xi10,p6,p8,p12}.case`,
`examples/zalesak_disk/zalesak_phi_normal.case`, `zalesak_p7_phi.case`. After the
re-run, update those READMEs and `evidence/` (`slab_1d_methods.mp4`,
`slab_1d_psi_field.mp4`, `slab_1d_grad_psi.png`, `zalesak_methods.mp4`), `CDI_METHOD.md` §7 and the root README row.

Until the re-run, the READMEs need a "pre-`grad_floor`-fix" label:
- `examples/zalesak_disk/README.md:26,121-130,190-203`: the "49×", the ~55 000 nodes and the $\phi$-normal
  resolution rows;
- the $\phi$-normal rows of `examples/advecting_slab_1d/README.md`.

Saini's Figs. 9–10 cells (the $\phi$ profile along $y=0.75$, $\xi\in\{0.5,1,1.5\}$,
$N\in\{5,7\}$): never run here; no `f910_*` output or log exists.

## Queued

1. **Ten rotations on the best few $(\xi,\gamma)$ cells** — $(\xi{=}1.5,\gamma{=}0.25)$
   and $(\xi{=}2,\gamma{=}0.25)$ against $(\xi{=}2.8,\gamma{=}1)$ (~4 h). Promotes the
   "2.1× better" one-rotation claim to a real recommendation, or kills it.
2. **Session B** — the full SVV-off $(\xi,\gamma)$ map at $N=7$ (~10 h). Its
   $\gamma=0.25$ row must use A's $\Delta t = 8.8028\times10^{-5}$ so the
   $E_r(\text{off})/E_r(\text{on})$ ratio map compares like with like.

## Open questions

- **Is a small boundedness violation acceptable for better shape error?** An
  application question. $\xi{=}1,\gamma{=}0.5$ gives the grid's best $E_r$
  (0.00240, 3.2× better than $\xi{=}2.8,\gamma{=}1$) at a $5.9\times10^{-6}$
  violation.
- **$\xi<1$ under strain.** After the diffusion fix the Rider–Kothe cross still improves down to
  $\xi=1$ ($E_r$ 0.0446), at a $2.5\times10^{-3}$ violation. $\xi=0.75$ would say whether it turns
  over (~30 min).
- **Would a band-masked compression term rescue Algorithm 1?** The plateau
  (§4.1) only hurts because $\mathbf{n}$ is read where $\phi(1-\phi)$ is
  negligible anyway. Masking the compression flux to the interface band is a
  change to the CDI discretisation, not a redistancing knob.
- **Should the primary `zalesak.case` still be $\xi=2.8$?** It is measurably too
  thick. Queue item 2 would give the evidence to change it.
- **No convergence criterion on the pseudo-time solve** (`CDI_METHOD.md` §4.3): a fixed
  step count, the residual never checked. On Rider–Kothe with Saini's settings it measurably does not
  converge in the compression band (item 2).
- **`rider_h192.case` ships but has never been run.** The $h$-series stops at $H=1/128$. Run it
  or remove it.
- **Why the build does not converge at $N=3$.** There the Eq. (44) relaxation does not reach
  $\lvert\nabla\psi\rvert=1$ (24% build residual against 1.2–1.3% at $N=5,7$), and in
  arm C each event moved the band mean away from 1 (1.184 → 1.265, 1.250 → 1.328,
  1.207 → 1.327). Untested cause: `svv_rd` at $N_{svv}=N/4=0.75$ gives mode 1 a kernel
  weight of 0.44.

## Doc and code fixes (review of 2026-10-05)

Small content changes, found while committing the branch. Each needs only an edit.

**Numbers that disagree with their source**
- `examples/zalesak_disk/README.md`:
  - `:124,356`: the $\phi$-normal worst violation is 5.5e-3; the notebook gives 6.25e-3.
  - `:455`: the `seed="phi"` worst violation is 9.1e-2; the notebook gives 1.06e-1.
  - `:71`: "401 frames" for `zalesak.case`; `output_zalesak/` has 101.
  - `:463`: $E_r$ 1.667 "at $t=3$"; the notebook's last frame is $t=3.20$, at 1.688.
- `examples/rider_kothe/README.md`:
  - `:~165`: "band mean 1.60–1.67 after each event" holds for $t=2$–5 only. At $t=6$ and 8 it is
    1.409 and 1.484 (`logs/d5/ev/rd_quality_eLf_Af.txt`).
  - `:28,160`: "ran 2026-10-01 with the diffusion fixed" needs a clause saying that the scratch file
    had the fix. Otherwise it reads as contradicting "every number before 2026-10-02".
  - The SVV-off $E_r$ 0.0498 ($\xi=1$) has no traceable source.
- No source found, though not shown wrong: circles README `:~137` (9 of 18 zero-set signs) and
  `:~145` (h5n8 0.6% without the guard); slab README `:~137` (4-rank 0.00023).

**`examples/redistance_circles/README.md` against the paper**
- `:46`: "BDF1 on step 1, $\mathbf D_\mu$ from $\psi^n$" is labelled "paper". Neither is printed.
- `:91`: §4.4 starts on p. 16, not p. 18.
- `:49`: $\varepsilon=\xi H/N$ is our $\xi$. The paper prints $\varepsilon=\xi H$ with $\xi=\{1,1.5\}/N$, so
  add the `CDI_METHOD.md` §1 pointer.
- `:26`: "one element thick in $z$, walls" is our setup, but it sits in the row for the paper's $\Omega$.
- `eps_1d/README.md:86` divides by $\int|\psi_e|$ without saying why: $\int\psi_e=0$ there.

**References to files that are deleted or local**
- `CDI_METHOD.md:1084` and circles `README.md:209-210` cite the deleted
  `evidence/archive_2026-09-28_epsHN/`.
- `CDI_METHOD.md:1082` and circles `README.md:206-207` cite `evidence/talk/`, which is gitignored.
- `CLAUDE.md:186` says the circles evidence comes from `logs/mkfigs.py`. It comes from
  `logs/saini_case/figs44.py evidence` and `logs/saini_case/anim_new/`.
- `eps_1d/README.md:5,26` and slab `README.md:298-299` cite local files without a "(gitignored,
  local)" mark.

**Code and scripts**
- `rider_kothe.f90:385` and `zalesak_disk.f90:389` set the prescribed velocity in `compute`, which runs
  after the scalar step. Rider–Kothe advects with $u(t_{n-1})$, and step 1 of both runs with $u=0$.
  This is $O(\Delta t)$ and negligible at the shipped $\Delta t$. Fix: the `preprocess` hook. It is a
  Fortran change, so ask first.
- Stale header comments: `advecting_slab_1d.f90:11-12` ("`svv_psi` defaults to off") and
  `redistance_circles.f90:8-11` ("Rider–Kothe diverges at $t=5.76$").
- The docstring of `examples/redistance_circles/scripts/gen_cases.py:5-7` says `HN:1`. The committed
  cells carry 0.25.
- `eps_1d/eps_mechanism_1d.py:17` hard-codes a workstation path for `sem1d`.
- The "off" column of `advecting_slab_1d/scripts/slab_table.py` needs the local
  `logs/svv_off_2026-10-02/` outputs.
- Small fixes:
  - `examples/rider_kothe/README.md:9` is missing a blank line.
  - The slab `.gitignore` lists `logs/` twice.
  - `contrib/lint_format/lint.sh` was never run, because `flint` is not installed.
