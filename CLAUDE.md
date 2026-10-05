# CLAUDE.md

This repo is a **showcase**, not an investigation: Neko plus demo cases for CDI
with the compression normal taken from a separately transported signed-distance
field $\psi$. The forensic work — why CDI needs a `psi`-normal at all, what was
tried and failed, every measured number behind the decisions in `CDI_METHOD.md`
— lives in the sibling repo `../neko-multiphase/`. Link to it; don't re-derive
it here. Read `README.md` and `CDI_METHOD.md` before touching a case's `.f90`.

## Where things are

| file | holds | read it when |
|---|---|---|
| `README.md` | the cases, their status table, recommended settings, build/run | always, first |
| `examples/redistance_circles/README.md` | the §4.4 case on its own: the problem, its configuration and where each setting comes from, how the code works (paper-against-code table), results against Saini, limitations, running, evidence | **before any §4.4 work, first** |
| `examples/redistance_circles/archive/README_process_2026-09.md` | the archived process record: the audit, the ε = H/N work (drift, capture, apex, N/4 arms), eliminated hypotheses, the study of the authors' code, dated history | when a question is about how or why, or needs an old number |
| `CDI_METHOD.md` | the method: equations, naming, every design decision. §1 $\xi$ and the notation clash with Saini, §4 redistancing, §4.4 the method lessons of the circles case, §5 SVV, §7 what is actually validated | before touching any `.f90` or making a method claim |
| `REDISTANCING.md` | when to redistance, and **§9: what the Fortran computes for Eq. (44), routine by routine** | before touching `rd_*`, `unit_normal`, the SVV steps |
| `NEXT_SESSION.md` | the Rider–Kothe redistancing plan (velocity lag, the events path, D5/D2, the thin tail), Saini §4.5's reference values, and the parked items for the other cases | working on the coupled cases |
| `examples/<case>/README.md` | exact parameters and results of that case | before quoting any number |
| `references/` | the JCP PDF (gitignored) and a frozen snapshot of the test-case notes | see "Writing docs" |
| `../neko-multiphase/references/` | the older ANL report (PDF) and its reading notes, the source of truth for the notes | only for implementation detail |

Plans hold plans; settled results go into the case READMEs and `CDI_METHOD.md`;
standing rules go here. **Saini §4.4 facts go into `examples/redistance_circles/README.md`
first**; `CDI_METHOD.md` §4.4, `REDISTANCING.md` §9 and the root README summarise it
and point there. The README describes the case as it is; process and history go to its
archive, not into the README. Don't grow parallel copies. When an item in a plan is finished, move it out rather
than letting the plan grow.

## Method invariants — settled, don't relitigate

- **SVV is a `psi`-only knob. There is no `svv_phi`: the code has none, and the
  three coupled files stop at startup if `case.cdi.svv_phi` is present.** The architecture is fixed: **`phi` carries the CDI
  equation** (transport + compression + the physical diffusion
  `eps*gamma*u_max`); **`psi` carries transport plus SVV**. **Every run that
  transports `psi` for the normal (`normal = "psi"`) has `svv_psi` on, `c0 > 0`**
  (user, 2026-10-02). SVV-off is not a configuration of this method, so don't
  run it, not even as a baseline row. Saini does the same: SVV on every scalar
  equation, his TLS transport Eq. (43) included. The three coupled files stop at
  startup if `normal = "psi"` and `svv_psi` is off. The value is $N_{svv}=N/2$ with
  `c0` 0.1, his TLS-transport value (§4.5, his §5 default), used by Rider–Kothe and
  the slab's sweep cases (its SVV and redistancing variants keep 1). Zalesak still ships `c0` 1, which comes from his §4.3 setting for the
  *CLS* transport (our $\phi$), not the TLS; whether to change that is open
  (`NEXT_SESSION.md`). The only SVV knob that ever varies is
  `case.cdi.svv_psi.c0`. Don't propose an `svv_phi` arm —
  not as an ablation, and not when comparing against a paper that puts SVV on
  the phase field (Saini et al.'s naming is transposed, so their SVV sits where
  `svv_phi` would; we do not adopt that placement). Nor may SVV replace or merge
  with `phi`'s diffusion: that term is equilibrium-defining, not numerical
  seasoning. See `CDI_METHOD.md` §5.
- **Where the SVV viscosity sits is not a detail for Eq. (44).** `svv_local`
  applies `nu` *inside* the bilinear form, with $|\mathbf c|=1$. Saini's
  Eqs. (7), (31) and (33) instead left-multiply the *assembled* operator by a
  pointwise $\mathbf D_\mu$. For Eq. (44) that is $|\operatorname{sgn}\psi|$, which
  vanishes on the zero set. The two forms agree only for constant $\nu$. The
  printed form is what stops the SVV dragging the interface in
  `redistance_circles` (`CDI_METHOD.md` §4.4).
  - Don't make `nu` vary inside `svv_local` and call it Eq. (31).
  - `redistance_circles` uses the printed form and nothing else, in its own
    `svv_step_eq31`, built on the 17 shared routines without editing them. It is
    how that equation is solved, not a `.case` option; don't reintroduce a
    toggle for the old $|\mathbf c|=1$ form. Moving the coupled cases to it is a
    separate decision (`NEXT_SESSION.md`, D2).
  - The sign convention of an IC (negative vs positive inside) cannot matter to
    Eq. (44): the equation and the scheme are exactly odd in $\psi$, checked bit
    for bit.
- **In Eq. (44), $\mathbf w=\operatorname{sgn}(\psi)\nabla\psi/|\nabla\psi|$ is a function of
  $\psi$, not a velocity field.**
  - Re-form it from each time level's own $\psi$, and apply it only to that $\psi$. In
    BDF/EXT the lags hold the products $\operatorname{sgn}(\psi^m)-\mathbf C(\mathbf w^m)\psi^m$, never $\mathbf w$.
  - Its sign term cancels the source's $\operatorname{sgn}'(\psi)\delta$ on the zero set. A frozen,
    lagged, extrapolated or subcycled $\mathbf w$ leaves growth at up to $1/2\varepsilon$: frozen
    diverges, and a lag of 0.05 in $\tau$ is 6× worse.
  - So Eq. (44) never goes through Neko's scalar solver, whose velocity is an external
    field. In `redistance_circles` (D4, 2026-09-29) $\mathbf C(\mathbf w)$ is **dealiased** on
    $\lfloor 3(N+1)/2\rfloor$ Gauss points (`adv_dealias_t`), as in Saini's code (`convect_new`),
    together with his sign guard: a node whose $\psi^n$ disagrees in sign with $\psi_0$ keeps
    $\psi^n$. Dealiasing moves zero-set nodes and breaks a flat interface's exact steady state.
    The coupled cases keep the non-dealiased GLL form (`conv1`, then B, gather-scatter, Binv,
    equal to $\operatorname{sgn}(\psi)|\nabla\psi|$ to 1e-14) until D5 decides.
  - `redistance_circles` integrates it with Saini's own BDF2/EXT2 (Eqs. 34–35), its
    only scheme since 2026-09-26. At $\varepsilon=H/N$ (2026-09-26) BDF2 and BDF3 agree with RK3
    within 1.9% in the Table 2 cells to $\tau=24$
    (`examples/redistance_circles/archive/README_process_2026-09.md` §7.4), and BDF2 within
    2.2% over the whole Fig. 12 grid (§7.2). So the
    integrator is not a suspect. The exception is the weaker
    $N_{svv}=N/4$ at $H{=}1/5$, $N{=}7$: only RK3's Lie-split damping held the apex
    there (§7.5). The coupled cases still use SSP-RK3 with `svv_step_imp`.
- **Never seed `psi` from an analytic distance in any run that redistances.**
  `case.cdi.psi_init` has two values: `"exact"` (the analytic periodic distance,
  the default, and what the validated no-redistancing results use) and
  `"redistance"` (Saini's Algorithm 1 line 2 — build `psi` by the Eq. (44) solve
  from the `r_f(phi-0.5)` seed). Redistancing on top of an analytic `psi` is the
  incoherent pairing that confounded the first ablation (`CDI_METHOD.md` §4);
  more importantly, **the analytic distance is a crutch that does not exist in
  the applications this method is for**, so anything built on it does not carry
  over. With `psi_init = "redistance"` the analytic distance function must not
  be *evaluated at all* for `psi` — in all three coupled user files its one call
  site sits in the `else` branch of that check, inside `initial_conditions`.
  Keep it that way; don't "helpfully" seed the analytic field first and
  overwrite it.
- **`redistance.band` coverage is a quality knob, not a stability one.** The
  compression term reads the normal out to $9.21\varepsilon = 9.21\xi H/N$ while
  redistancing rebuilds `psi` only within `band`$\times H$, so `band >= 9.21*xi/N`
  ($N >= 3.68\xi$ at the default 2.5) is what makes the built field cover the
  interface — it measurably improves the build (band mean 1.259 -> 1.027 at
  `N=3`). **It does not fix arm C.** With the history bug (2026-09), raising the band
  made `N=3` diverge *sooner*, and `N=5` diverged while satisfying the bound; the
  band-4.0 run has not been repeated since. (An earlier version of this
  rule claimed otherwise.) Re-run with the history restart (2026-10-02): `N=3`
  still diverges ($t=11.2$), `N=5` completes 17× worse than arm B, and `N=7`
  diverges. The cause is $\phi$'s contour fragmenting, with each reseed sustaining
  the fragments (`CDI_METHOD.md` §4.1c).
- **The `unit_normal` gradient floor must stay above the round-off gradient of a
  flat field (~1e-19); 1e-6 is load-bearing.** A banded `psi` is flat over most
  of the domain, and a floor below round-off makes `unit_normal` hand the
  compression term a *random unit vector* there, which diverges the run. It was
  once 1e-30 and its own comment claimed to guard exactly this case
  (`CDI_METHOD.md` §4.1). The three coupled files read it as
  `case.cdi.grad_floor` (default 1e-6); `redistance_circles` fixes it as a
  parameter. Don't lower it, and don't assume `|grad psi|` is order 1 just
  because it is inside the interface band.
- **"Redistancing" and "reinitialization" are not the same operation** — and
  **Saini's usage differs from this repo's**. For them, Eq. (44) is the
  "re-distancing *equation*" and the *operation* is "re-initialization", which
  is **always** the Eq. (47) reseed followed by the Eq. (44) solve; Algorithm 1
  has no in-place mode. This repo splits them: `seed="phi"` is Saini's
  operation, `seed="psi"` (relax in place, no reseed) is our own. With the
  event's history restart, neither costs anything on the 1D slab. The old "the
  reseed is the entire cost" was the history bug. On Zalesak the reseed sustains
  $\phi$'s contour fragments (`CDI_METHOD.md` §4, §4.1c). Use the two words
  precisely, and say whose sense you mean when it could matter.
- **Attribution:** don't describe this method as "Saini's approach". Neko's CDI
  equation is structurally different from Saini et al.'s split CLS scheme
  (`CDI_METHOD.md` §2), and even the transported-signed-distance-normal idea the
  two share is, by their own account, adopted from earlier work (Salami et al.,
  classical TLS redistancing). Cite Saini et al. (2026) only for the specific
  equation forms/parameters this repo's implementation was validated against.

## Repo conventions — keep the showcase minimal

- No diagnostic CSV output, no raw/ungathered diagnostic fields, no comparison
  variants beyond `zalesak_disk`'s one sanctioned ablation, no session-log
  markdown files. If you're tempted to add one, it belongs in `neko-multiphase`
  (or a fresh `neko-multiphase-*` sandbox). The one deliberate exception per
  case is its `visualize.ipynb` — a short pysemtools notebook (load two
  snapshots, one plot). Don't grow it into a multi-part diagnostic notebook.
- **A validated fix becomes *the* code, not an option.** Validate it against the
  reference (temporary bit-identity scaffolding is fine), then remove the
  switch rather than defaulting it, and clean up after it: docs, stale files,
  the notebook, dead toggles. Old results survive only as history in the docs.
- Generated sweep or diagnostic `.case` files go in the session scratchpad (or
  a `scripts/gen_cases.py` output dir), never next to the shipped cases.
- `examples/zalesak_disk/` is the one exception to "one settled `.case` per
  case": a small, labelled ablation set whose **primary is redistancing off**
  (the validated configuration — don't invert this, an earlier draft of these
  docs did), plus two variants that *add* redistancing (`seed="psi"` in place,
  `seed="phi"` reinit+redistance). The two variants build $\psi$ by Eq. (44)
  (`psi_init = "redistance"`) and fire on the timer; their recorded numbers predate
  that configuration and are not citable until re-run. Understanding that distinction
  (`CDI_METHOD.md` §4) is part of what this repo is for. Don't let the pattern
  spread, and don't grow the set beyond those three without a specific reason.
  (`advecting_slab_1d` ships a $\xi$/$\gamma$/$N$ sweep because its operating
  envelope is what it demonstrates; see its README.)
- `examples/redistance_circles/` is a benchmark reproduction, not a showcase
  case, so the one-`.case` rule does not apply: it ships Saini Table 2's four
  cells, and `scripts/gen_cases.py` generates the rest of the Fig. 12 grid on demand
  (ε rule `0.25`). Its committed cases carry $\varepsilon = 0.25$ in Eq. (46), the value
  in Saini's own code (`signls`, `nandu90/nekLS_Examples@jcp`); the paper's $\xi H/N$ is
  the phase field's width, not the sign function's. Don't change it back. Its status is
  **reproduced** (2026-09-29): with his configuration the code gives his `plot.py` values,
  and every remaining difference is zero-set round-off under his guard (runs seeded with
  his zero-set $\psi_0$ match exactly). Its README must say plainly that low-N cells can
  fail on our mesh for that reason. That README describes the case on its own; its §3 table
  maps every printed equation to a routine. Keep it in step with the code. How we got here is
  archived in `archive/README_process_2026-09.md`.
- `examples/rider_kothe/` is complete with redistancing off. Its numbers are from the
  2026-10-02 re-runs with the diffusion fix; every earlier Rider–Kothe number had the
  frozen diffusion, so don't quote one. With both fixes, Saini's full periodic-reseed
  configuration (scratch user file) completes ($E_r$ 0.0407 vs 0.0446 transport only).
  But its $\psi$ normal is 4–8° off where transport's is 0.6–0.8°. So it is not a
  validated alternative; say so plainly.
- **Saini's Rider–Kothe $E_r$ is not on our scale:** his code divides by the area
  *outside* the disk (his CLS is 1 outside), ~13× ours. Recompute on his dumps with
  `norms.E_r` before comparing (`examples/rider_kothe/README.md`).
- Each case's `evidence/` holds the figures and animations produced by this
  repo's own runs: by its `visualize.ipynb`, or for `redistance_circles` by its
  gitignored `logs/saini_case/figs44.py evidence` and `logs/saini_case/anim_new/`. They are meant to be committed, since they are
  the demonstration this repo exists to show. Don't add a placeholder
  `evidence/` to a case that has produced no figures (`CDI_METHOD.md` §8).
- `references/saini_2026_test_cases.md` is a carried-over snapshot of reading
  notes from `neko-multiphase`, not this repo's living document. See its
  disclaimer before editing. If it drifts, `neko-multiphase` is the source of
  truth.

## The four user files

- **The shared machinery exists four times**: `svv_t` and its ~375 lines (the
  power iteration in `svv_init`, the implicit CG, `svv_report`) plus the
  redistancing primitives (`rd_sgn`, `rd_rhs`, `unit_normal`) are in
  `zalesak_disk.f90`, `rider_kothe.f90`, `advecting_slab_1d.f90` and
  `redistance_circles.f90`. The event-driven wrapper (`redistance`,
  `band_grad_stats`) is in the three *coupled* cases only; `redistance_circles`
  solves Eq. (44) standalone with its own `redistance_standalone`/`grad_stats`.
  The duplication is deliberate: it keeps each case to one self-contained user
  file. But **a fix to a shared routine must be made in all four**. Each file
  says so in a comment. **Seventeen** routines are byte-identical; check with
  `python3 examples/tools/check_shared_routines.py <the four .f90 files>`.
  A project hook runs this after every edit to one of the four and names any
  shared routine that has drifted (see "Agents and hooks").
- An **explicit** SVV instance is initialised lazily by its own source-term
  hook, so a scalar with no `"source_terms": [{"type": "user"}]` entry in the
  `.case` silently never gets its SVV term, while the startup header still
  reports it as on (the header reports what was *requested*). All three coupled
  files raise a hard error after step 1 if an explicit `svv_psi` was requested
  but its source hook never fired. Keep that guard in any new case built from
  them.
- Inline comments in the `.f90` files should be rare and explain only
  genuinely non-obvious physics/numerics (e.g. a sign convention, why SSP-RK3
  is required for the redistancing pseudo-time step), not restate the code.
- Fortran: `use neko` brings a lot of names into scope. Avoid module-level
  names that collide with it (`file`, `field`, `vector`, `math`, `space`,
  `csv_file`, …).
- Neko gotchas worth not rediscovering:
  - `field_vdot3` is **broken**: it declares its result `intent(out)`, which
    deallocates the field's storage on entry. Nothing in Neko calls it, so it is
    unexercised upstream. Use `field_col3` + two `field_addcol3`.
  - The CUDA kernels cap at `lx = 10` (`opr_dudxyz`, `opr_cfl`, `opr_conv1`),
    i.e. `polynomial_order` <= 9, so `advecting_slab_1d` ($N=10$) is CPU-only.
  - `makeneko` fails with "Text file busy" if a running `neko` still holds the
    binary: move it aside first.
  - `LOG_SIZE` is 79 characters. A longer `write(mess, ...)` dies at runtime
    with "End of record", possibly thousands of steps in.
  - `neko_log%message` indents by three spaces, so grep log lines as `'^ *tag'`.
  - **A time-varying scalar conductivity must also fill `<name>_lambda_tot`.** The
    solve reads `lambda_tot`, which Neko copies from `lambda` only at
    initialisation unless a turbulence model is set. Filling only `s_lambda` froze
    Rider–Kothe's CDI diffusion at its $t=0$ value until 2026-10-02 and dissolved
    the filament (`CDI_METHOD.md` §4.1d).
  - **A user hook that replaces a scalar field must restart its time history.**
    `compute()` runs after the scalar step's `slag%update()`, so the BDF lags still
    hold the old field, and BDF3 settles at old + 11/6 (new − old). Set
    `neko_user_access%case%fluid%ext_bdf%nadv = 0` and `%ndiff = 0` after the
    replacement (next step BDF1/EXT1, as Saini's `ireset_ls`). The redistancing events
    of the coupled cases lacked this until 2026-10-02; all three now have it
    (`CDI_METHOD.md` §4.1d). A small change
    (an implicit SVV sub-step) is an O(Δt) split and needs no restart.

## Running things

- The `configure`/`make install` build is run manually by the user, not
  automated here. `source setup-env.sh` (CPU) or `setup-env-cuda.sh` (GPU)
  selects which `neko` runs. Each case's `run.sh` runs `makeneko`, then its
  cases.
- **One GPU.** Check it is free before launching
  (`nvidia-smi --query-compute-apps=pid`; don't count `neko` processes, since a
  CPU run of `advecting_slab_1d` is one). Launch long chains detached:
  `setsid nohup ./run.sh ... > chain.log 2>&1 < /dev/null &`.
- To stop a run, kill by PID from
  `ps -eo comm,args --no-headers | awk '$1=="neko"{print}'`. A `pkill -f`
  pattern matches the calling shell and kills the wrong thing.
- **`advecting_slab_1d` is the 4-minute regression testbed.** Check any change
  to the shared machinery there before spending hours on Zalesak or
  Rider–Kothe.
- **Run logs are huge**: `run_*.log` reaches 620 MB (7.6M lines). Never `cat` or
  `Read` one. Grep anchored tags with `-m`/`tail`, or hand a multi-log question
  to the `log-probe` agent.
- Python is `python3` (there is no `python`). Notebooks and figure scripts use
  the `venv-neko` kernel at `/lscratch/sieburgh/local/venv-neko`.
  `pdftotext`/`pdftoppm` are installed.

## Research practice — each of these was paid for once

- **Vary one thing at a time.** Three conclusions were drawn and retracted on
  2026-09-08 because a run changed two variables at once.
- **A one-shot test does not reveal repeated-event instability.** That is how a
  pseudo-timestep that is fine for one build and unstable under 40 events
  nearly became the default (`CDI_METHOD.md` §4.2).
- **The band minimum of $|\nabla\psi|$ diagnoses nothing.** It is
  $7.5\times10^{-14}$ for the exact analytic $\psi$ in a run that is clean.
  Watch the band mean.
- **Don't judge by one number at one time.** Run long enough to see whether a
  minimum is transient; for `redistance_circles` judge by the curve shapes and
  the Fig. 13 map split, not only $E_r(6)$.
- **Compare two implementations from the same $\psi_0$** before calling a
  late-time difference a bug. Round-off on zero-set nodes alone moves results.
- **Neko is not a suspect** when a result doesn't reproduce a Nek5000 paper: the
  two are near-identical. Look at our user file's scheme choices.
- **Treat any "Done"/status claim in these docs as provisional.** An earlier
  pass overclaimed by trusting prior summaries over the case files and notebook
  *outputs* (`CDI_METHOD.md` verification note and §7). When a status claim
  matters, check the primary source (or run `claim-auditor`) before repeating
  it.
- **Ask the user before:** changing Neko itself (`src/`); changing any of the
  17 shared routines; changing a committed `.case` file's parameters or the
  Fortran's log line formats (the notebooks and `logs/*.py` parse them);
  reading Nek5000's source; contacting the paper's authors.

## Writing docs

- Write equations in markdown as LaTeX (`$...$` inline, `$$...$$` block), not
  ASCII pseudo-equations. Exception: `references/saini_2026_test_cases.md` is a
  frozen verbatim snapshot; don't reformat it in place.
- **$\xi$ means $\varepsilon N/H$ here, which is $N$ times Saini's $\xi$**
  (`CDI_METHOD.md` §1). Their field names are also transposed: their $\psi$ is
  our $\phi$. Never compare the bare symbols across the two.
- Quote $E_r$ for `redistance_circles` with Eq. (84)'s printed denominator
  $\int\psi_e$ only. The $\int|\psi_e|$ reading is dropped.
- Use absolute dates. Say where a number comes from (file, log, notebook
  output).

## Agents and hooks

Local to this checkout (upstream Neko's `.gitignore` ignores `.claude/`).

- A **hook** runs `examples/tools/check_shared_routines.py` after every edit to one of the four user files
  and names any of the 17 shared routines that has drifted. Silent when all is
  well; mid-propagation it is the list of files still to update.
- **Delegate reading, not thinking.** PDF lookups go to `paper-lookup`,
  multi-log summaries to `log-probe`, status-claim checks to `claim-auditor`.
  Keep planning, experiment design and Fortran edits in the main thread: they
  need the plan's context and the constraints above. When an equation's exact form
  matters, ask `paper-lookup` to read the **rendered** PDF pages, not the text
  cache, which garbles math.
