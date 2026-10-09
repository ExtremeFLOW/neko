# CLAUDE.md

This repo is a **showcase**: Neko plus four cases for CDI with the compression normal taken from a
separately transported signed distance $\psi$.
- **What comes from where.** The $\psi$ machinery is meant to be Saini & Tomboulides' (JCP 2026) as
  published. The $\phi$ equation is Neko's CDI, which is the one intended difference.
- **Attribution.** It is not "Saini's approach". The transported-distance normal predates him
  (Al-Salami et al.), and Eq. (44) is classical TLS re-distancing. Cite him for the equation forms
  and parameters we validate against.
- **The investigation lives elsewhere.** Why CDI needs a $\psi$-normal at all, and what failed, is in
  `../neko-multiphase/`. Link to it; don't re-derive it.

Read `README.md` and `CDI_METHOD.md` before touching any `.f90`.

## Where things are

| file | holds | read it when |
|---|---|---|
| `README.md` | the cases, the results table, recommended settings, build/run | always, first |
| `CDI_METHOD.md` | the method; §5 SVV; **§6 Saini's configuration against ours** | before a method claim or a `.f90` change |
| `REDISTANCING.md` | maintaining $\psi$: terms, what is known, rules; §5 the code for Eq. (44), routine by routine | before touching `rd_*`, `unit_normal`, the SVV steps or the events path |
| `NEXT_SESSION.md` | what is open, including the Rider–Kothe porting spec | before working on the coupled cases |
| `examples/<case>/README.md` | exact parameters, results, and an appendix of runs in other configurations | before quoting any number |
| `examples/redistance_circles/README.md` | the §4.4 record, with its paper-against-code table | before any §4.4 work |
| `examples/redistance_circles/archive/README_process_2026-09.md` | the frozen §4.4 process record | for an old number, or a how or why |
| `examples/<case>/logs/` | gitignored, local: run records (`PREREGISTERED.txt`), scratch user files, measurement and figure scripts | when a number's source matters |
| `examples/tools/check_shared_routines.py`, `examples/norms.py` | the shared-routine checker; the error norms | |
| `references/` | the JCP PDF and Saini's email of 2026-09 (both gitignored) | through `paper-lookup` |
| `../neko-multiphase/references/` | the older ANL report and the reading notes on Saini's test cases (don't copy them here) | for implementation detail |

## Writing and placing docs

- **Results live in the case READMEs only.** §4.4 facts go in the circles README first.
  `README.md` summarises, and the method docs point. Don't grow parallel copies.
- **Plans live in `NEXT_SESSION.md`.** A finished item leaves it.
- **Delete superseded text**; don't archive or date it. Git is the record, and the circles archive is
  frozen.
- **Non-Saini runs go in the appendix.** A run in a configuration Saini doesn't use becomes one row
  of its case README's appendix table, with no discussion.
- **No dated narrative.** Dates appear only as the provenance of a number or a status.
  - Say where every number comes from: file, log, or notebook output.
  - Name configurations by content ("Saini's settings, our $\Delta\tau$"), not by run label.
- **Write equations in LaTeX** (`$...$`, `$$...$$`), not ASCII.
- **$\xi$ means $\varepsilon N/H$ here**, which is $N$ times Saini's $\xi$. His field names are
  transposed: his $\psi$ is our $\phi$. Never compare the bare symbols.
- **$E_r$.**
  - `redistance_circles` uses Eq. (84) with $\int\psi_e$ only.
  - Saini's Rider–Kothe $E_r$ divides by the area *outside* the disk, about 13× ours. Recompute on
    his dumps with `norms.E_r` before comparing.
- **Don't quote Rider–Kothe numbers from before 2026-10-02** (frozen diffusion). Mark those from
  before 2026-10-05 as "before the velocity fix".

## Method rules — settled, don't relitigate

- **SVV is a $\psi$-only knob.** $\phi$ carries the CDI equation: transport, compression, and the
  physical diffusion $\varepsilon\gamma u_{\max}$. $\psi$ carries transport plus SVV.
  - There is no `svv_phi`; the coupled files stop at startup if `case.cdi.svv_phi` is present. Don't
    propose one, not as an ablation, and not when comparing with Saini's SVV on his phase field.
  - SVV never replaces or merges with $\phi$'s diffusion.
- **`svv_psi` is always on for `normal = "psi"`.** $c_0>0$ is enforced at startup.
  - The value is $N/2$, $c_0=0.1$, Saini's $\psi$-transport value. Zalesak still ships 1.0, his
    phase-field value; that is open.
  - **Never run SVV-off**; old SVV-off results are appendix rows.
  - $c_0$ is the per-case knob. $|\mathbf c|$ is the flow's $u_{\max}$ in the committed files and the
    local $|\mathbf u|$ in Saini's.
- **Eq. (44)'s SVV is the printed Eq. (31) form.** A pointwise
  $\mathbf D_\mu=|\operatorname{sgn}\psi|$ left-multiplies the assembled operator (`svv_step_eq31`).
  - Don't make $\nu$ vary inside `svv_local` and call it Eq. (31).
  - No toggle back to the uniform form. The coupled files still use `svv_step_imp`; moving them is
    `NEXT_SESSION.md` item 1.
- **$\mathbf w=\operatorname{sgn}(\psi)\nabla\psi/|\nabla\psi|$ is a function of $\psi$, not a
  velocity.**
  - Re-form it from each level's own $\psi$. BDF/EXT lags hold the products, never $\mathbf w$.
  - So Eq. (44) never goes through Neko's scalar solver.
  - An IC's sign convention cannot matter: Eq. (44) and the scheme are odd in $\psi$.
- **The integrator is not a suspect.** BDF2, BDF3 and RK3 agree within 2.2% on §4.4. The exception:
  at $N_{svv}=N/4$ only RK3's split damping held an apex.
- **No analytic seed when redistancing.** `psi_init` is `"exact"` (analytic; the validated runs) or
  `"redistance"` (Algorithm 1, line 2).
  - With `"redistance"` the analytic distance is never evaluated. Its one call site stays in the
    `else` branch of `initial_conditions`; don't seed it first and overwrite it.
  - Never pair `"exact"` with events.
- **`grad_floor` is $10^{-6}$ and load-bearing.** It must stay above a flat field's round-off
  gradient, about $10^{-19}$. The coupled files read `case.cdi.grad_floor`; `redistance_circles` fixes
  it as a parameter.
  - Don't lower it.
  - Don't assume $|\nabla\psi|\sim1$ inside the band.
- **Band coverage $N\ge3.68\xi$** (at band 2.5) is a build-quality bound, not a stability fix.
- **Redistancing and reinitialization are different words.**
  - Saini's "re-distancing" names Eq. (44); his operation, "re-initialization", is the Eq. (47)
    reseed plus Eq. (44), and he has no in-place mode.
  - Ours: `seed = "phi"` is his operation, `seed = "psi"` our in-place relaxation.
  - Say whose sense you mean.
- **The $\psi$ target is Saini's configuration** (`CDI_METHOD.md` §6). A deviation is a work item, not
  an alternative.

## Repo conventions — keep the showcase minimal

- **No diagnostics in the repo.** No diagnostic CSV output, no raw or ungathered diagnostic fields,
  no session-log markdown. Those belong in `neko-multiphase` or a scratch sandbox.
- **One short `visualize.ipynb` per case.** Don't grow it into a diagnostic notebook.
- **A validated fix becomes the code, not an option.**
  1. Validate it; temporary bit-identity scaffolding is fine.
  2. Remove the switch.
  3. Clean up the docs, stale files, the notebook and dead toggles.
- **Generated `.case` files** go in the session scratchpad or the case's gitignored `logs/`, never
  beside the shipped cases.
- **Shipped `.case` sets.** Don't grow them without a specific reason.
  - **Zalesak, 8:** the primary with redistancing **off** (don't invert this), its $\phi$-normal
    baseline, four resolution cases, and two variants that *add* redistancing. The variants' old
    numbers are not citable.
  - **Rider–Kothe, 8:** the $\xi$/$\gamma$ cross and the $h$-series.
  - **Slab, 28:** its $\xi$/$\gamma$/$N$ envelope is what it demonstrates, plus the redistancing
    testbed.
  - **Circles, 4:** Table 2. `scripts/gen_cases.py` generates the rest of the Fig. 12 grid.
- **`redistance_circles` is a reproduction.**
  - Its cases use $\varepsilon=0.25$ in Eq. (46), as Saini's code does (`signls`); don't change it
    back.
  - Its README must say that low-$N$ cells can fail on our mesh, through zero-set round-off under the
    guard.
  - Keep its §3 paper-to-code table in step with the code.
- **`rider_kothe` re-initialization** configurations exist only as scratch user files under `logs/`.
  None is in `rider_kothe.f90`.
- **`evidence/`** holds figures and animations from this repo's own runs, made by the notebook or,
  for circles and Rider–Kothe, by local scripts under `logs/`. Commit them.
  - Use the kthviz style (Figtree, the KTH palette, `cmap('phase')` for $\phi$), and ship every
    animation as a GIF too.
  - Don't add a placeholder `evidence/`.

## The four user files

- **The shared machinery exists four times**: `zalesak_disk.f90`, `rider_kothe.f90`,
  `advecting_slab_1d.f90`, `redistance_circles.f90`. It keeps each case to one self-contained file.
  - **17 routines are byte-identical, in-body comments included**, so a fix to one is made in all
    four.
  - A hook runs `examples/tools/check_shared_routines.py` after every edit to one of these files and
    names any routine that has drifted. Mid-propagation it lists the files still to update.
  - The `!>` doc blocks above those routines are not checked. Keep them identical by hand, except
    where `redistance_circles` uses a routine differently.
  - `redistance`, `band_grad_stats` and `initialize` are identical across the three coupled files,
    also by hand.
- **An explicit SVV instance is initialised lazily** by its own source-term hook. A scalar without
  `"source_terms": [{"type": "user"}]` never gets its SVV term, while the header still reports it as
  on. The coupled files stop after step 1 when that happens; keep the guard in any new case.
- **The knobs are read from `case.cdi`**, not from `case.scalar(s)`.
- **Inline comments are rare.** Only non-obvious physics or numerics; no restating the code, no
  history.
- **`use neko` brings many names into scope.** Avoid module-level names such as `file`, `field`,
  `vector`, `math`, `space` and `csv_file`.

**Neko gotchas worth not rediscovering:**
- **`field_vdot3` is broken.** Its result is `intent(out)`, which deallocates the field's storage. Use
  `field_col3` and two `field_addcol3`.
- **The CUDA kernels cap at `lx = 10`** (`polynomial_order` $\le9$), so `advecting_slab_1d`
  ($N=10$) is CPU-only.
- **`makeneko` fails with "Text file busy"** if a running `neko` holds the binary. Move it aside
  first.
- **`LOG_SIZE` is 79.** A longer `write(mess, ...)` dies at runtime with "End of record".
- **`neko_log%message` indents by three spaces**, so grep log lines as `'^ *tag'`.
- **A time-varying scalar conductivity must also fill `<name>_lambda_tot`.** Neko copies it from
  `lambda` only at initialisation unless a turbulence model is set.
- **A prescribed velocity goes in `preprocess`, at `time%tlag(1)`** ($t_n$). `compute` runs after the
  scalar step, and `time%t` in `preprocess` is already $t_{n+1}$.
- **A hook that replaces a scalar field must restart its history:**
  `neko_user_access%case%fluid%ext_bdf%nadv = 0` and `%ndiff = 0`, so the next step is BDF1/EXT1, as
  Saini's `ireset_ls`.
  - Without it BDF3 settles at old + 11/6 (new − old).
  - A small change, such as an implicit SVV sub-step, needs no restart.

## Running things

- **The `configure`/`make install` build is run by the user.** `source setup-env.sh` (CPU) or
  `setup-env-cuda.sh` (GPU) picks which `neko` runs. Each case's `run.sh` runs `makeneko`, then its
  cases.
- **There is one GPU.** Check it is free first (`nvidia-smi --query-compute-apps=pid`). Don't count
  `neko` processes: a CPU slab run is one.
- **Launch long chains detached:** `setsid nohup ./run.sh ... > chain.log 2>&1 < /dev/null &`.
- **Stop a run by PID**, from `ps -eo comm,args --no-headers | awk '$1=="neko"{print}'`. A `pkill -f`
  pattern also matches the calling shell.
- **`advecting_slab_1d` is the 4-minute regression testbed.** Check shared-machinery changes there
  before Zalesak or Rider–Kothe.
- **Run logs reach 620 MB.** Never `cat` or `Read` one. Grep anchored tags with `-m` or `tail`, or use
  `log-probe`.
- **Python is `python3`.** Notebooks and figure scripts use the `venv-neko` kernel at
  `/lscratch/sieburgh/local/venv-neko`. `pdftotext` and `pdftoppm` are installed.

## Research practice — each of these was paid for once

- **Vary one thing at a time.**
- **A one-shot test does not reveal repeated-event instability.**
- **The band minimum of $|\nabla\psi|$ diagnoses nothing.** Watch the band mean.
- **Don't judge by one number at one time.** For `redistance_circles`, judge by the curve shapes and
  the Fig. 13 split, not only $E_r(6)$.
- **Compare two implementations from the same $\psi_0$** before calling a late-time difference a bug.
  Zero-set round-off alone moves results.
- **Neko is not a suspect** when a Nek5000 result doesn't reproduce. Look at our user file's scheme
  choices.
- **Treat any "done" or status claim in the docs as provisional.** Check the primary source, or run
  `claim-auditor`, before repeating one.
- **Ask the user before:**
  - changing Neko itself (`src/`);
  - changing any of the 17 shared routines;
  - changing a committed `.case` file's parameters, or the Fortran's log line formats (the notebooks
    and `logs/*.py` parse them);
  - reading Nek5000's source (Saini's fork, `nandu90/Nek5000@nekLS` and `nandu90/nekLS_Examples@jcp`,
    is allowed);
  - contacting the paper's authors.

## Agents and hooks

These are local to this checkout; upstream Neko's `.gitignore` ignores `.claude/`.
- **The hook** runs the shared-routine check after every edit to one of the four user files.
- **Delegate reading, not thinking.**
  - PDF lookups go to `paper-lookup`. Ask it to read the rendered pages when an equation's form
    matters, not the text cache, which garbles math.
  - Multi-log summaries go to `log-probe`, and status-claim checks to `claim-auditor`.
  - Planning, experiment design and Fortran edits stay in the main thread.
