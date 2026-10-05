# redistance_circles — Saini et al. (2026) §4.4 in Neko

Saini & Tomboulides (JCP 2026), §4.4 "Re-distancing around intersecting circles": the
re-distancing equation, Eq. (44), relaxed on its own from a skewed initial condition toward a
known signed distance function, measured against their Table 2 and Figs. 12–13.

**Status: reproduced.** With the authors' configuration, taken from their public case
(`nandu90/nekLS_Examples`, branch `jcp`, `intersectingCircles/`), this case gives their
published errors: 12 of the 18 Fig. 12 cells within 0.1% (eight of them to 4–5 figures), four
more within 2.2%. Two cells fail on our mesh through round-off; §5 explains why.

How we got here, with every intermediate result, is archived in
[`archive/README_process_2026-09.md`](archive/README_process_2026-09.md).

Naming: Saini's distance field is called $\phi$; this repo calls it $\psi$ (and the phase field
$\phi$). Their $E_r(\phi)$ is our $E_r(\psi)$ (`../../CDI_METHOD.md` §1).

## 1. The problem

Eq. (44), standalone: no flow, no phase field. The whole solve runs in the user `initialize` hook.

$$\frac{\partial\psi}{\partial\tau} + \mathbf{w}\cdot\nabla\psi = \operatorname{sgn}(\psi) + S_{vv}(\psi), \qquad \mathbf{w} = \operatorname{sgn}(\psi)\,\frac{\nabla\psi}{|\nabla\psi|}, \qquad \operatorname{sgn}(\psi) = \tanh\!\Big(\frac{\psi}{2\varepsilon}\Big)$$

| | |
|---|---|
| domain | $\Omega=[-2,2]^2$, one element thick in $z$, walls on all four sides |
| exact solution, Eq. (83) | signed distance to two circles, $r=1$, centres $(\pm0.7,0)$, negative inside (`circle_distance`) |
| initial condition, Eq. (83) | $\psi_0=((x-1)^2+(y-1)^2+0.1)\,\psi_e$ (`skewed_ic`) |
| meshes | $H\in\{1/5,1/10,1/20\}$ (`box20/40/80.nmsh`), $N=3$…8 |
| pseudo time | $\Delta\tau=5\times10^{-4}$ to $\tau=6$, 12000 steps |
| metric, Eq. (84) | $E_r=\int\lvert\psi-\psi_e\rvert\,d\Omega\,/\int\psi_e\,d\Omega$ (`error_report`) |
| targets | Table 2 (four cells), Fig. 12 ($E_r(6)$ over the grid), Fig. 13 ($\lvert e\rvert$ at $N{=}8$) |

## 2. Configuration

Everything below is either printed in the paper or taken from the authors' code and email
(Nadish Saini, 2026-09). Three settings are not in the paper; they are marked.

| item | value | source |
|---|---|---|
| sign function width, Eq. (46) | **$\varepsilon=0.25$, fixed** (not in the paper) | their `lvlSet.f`, `signls`: `tanh(phi/(2*0.25))` |
| $\mathbf C(\mathbf w)\psi$ | **dealiased**, $\lfloor3(N+1)/2\rfloor$ Gauss points per direction (not in the paper) | their `SIZE` (`lxd = lx1*3/2`), `convect_new`; email: "active for all my cases" |
| sign guard | **a node whose $\psi^n$ disagrees in sign with $\psi_0$ keeps $\psi^n$** (not in the paper) | their `constrainTLSR`; email |
| $\mathbf n$ | $\nabla\psi$ averaged at shared nodes, then normalised | Eq. (42); email |
| SVV, Eqs. (24), (29)–(33) | $\hat Q=(k/N)^{N_{svv}}$, $N_{svv}=N/6$, $c_0=2$; $\mathbf D_\mu=\lvert\operatorname{sgn}\psi\,\mathbf n\rvert\,c_0H/N$ left-multiplying the assembled operator | paper; their `usrdat3`, `svv.f`; email |
| time scheme, Eqs. (34)–(35) | BDF2/EXT2, BDF1 on step 1, $\mathbf D_\mu$ from $\psi^n$ | paper (their code extrapolates with Nek's EXT3, `settime_cls`: `NAB = 3`; equal to EXT2 within $5\times10^{-5}$ here) |
| boundary condition | zero Neumann (natural) | email |

The phase-field width of the method ($\varepsilon=\xi H/N$, Eq. (36)) is a different quantity:
it sets the CDI interface, not the sign function of Eq. (46).

## 3. How it works

Neko's scalar solver is not used: its advecting velocity is an external field, and here
$\mathbf w$ is a function of $\psi$ itself (`../../REDISTANCING.md` §9.0). The user file runs its
own pseudo-time loop, `redistance_standalone`. One step $n\to n+1$ is Eqs. (34)–(35):

$$F^m=\operatorname{sgn}(\psi^m)-\mathbf B^{-1}\mathbf C(\mathbf w^m)\psi^m,\qquad \hat\psi=\tfrac{1}{b_0}\Big(\textstyle\sum_j\beta_j\psi^{n+1-j}+\Delta\tau\sum_j a_jF^{n+1-j}\Big),\qquad \big(\mathbf B+\tfrac{\Delta\tau}{b_0}\mathbf D_\mu(\psi^n)\mathbf S_{vv}\big)\psi^{n+1}=\mathbf B\hat\psi$$

with $b_0=3/2$, $\beta=(2,-\tfrac12)$, $a=(2,-1)$, then the guard.

- **$\mathbf w$ is re-formed from each level's own $\psi$** and applied only to it; the EXT history
  holds the products $F^m$, never $\mathbf w$. A frozen, lagged or extrapolated $\mathbf w$ leaves
  growth at up to $1/2\varepsilon$ on the zero set.
- **$\mathbf n$:** the gradient is averaged at shared nodes (`unit_normal`, floor $10^{-6}$);
  $\mathbf w=\operatorname{sgn}(\psi)\mathbf n$ (`rd_sgn`).
- **$\mathbf C(\mathbf w)\psi$:** Neko's `adv_dealias_t` on $\lfloor3(N+1)/2\rfloor$ Gauss points,
  then gather-scatter and $\mathbf B^{-1}$. Dealiasing evaluates $\mathbf w$ between nodes, so
  nodes on the zero set move and can change sign.
- **The implicit SVV step** (`svv_step_eq31`): $\mathbf D_\mu=\lvert\operatorname{sgn}\psi^n\rvert$
  vanishes on the zero set and makes the operator non-symmetric. The substitution
  $\psi^{n+1}=\hat\psi+\mathbf D_\mu z$ gives a symmetric system, solved by CG to $10^{-12}$
  relative.
- **The guard:** after the step, a node whose $\psi^n$ has the opposite sign to $\psi_0$ keeps
  $\psi^n$. Their code reads the sign off the phase field; in this problem that is the sign of
  $\psi_0$.

The paper's equations against the code:

| paper | what it prints | code |
|---|---|---|
| Eq. (7), (33), p. 4, 7 | $\mathbf M\dot u=-\mathbf Cu+\mathbf Mf-\mathbf D_\mu\mathbf S_{vv}u$, $\mathbf D_\mu$ per element then assembled | `svv_step_eq31`: $\mathbf D_\mu$ left-multiplies the assembled $\mathbf S_{vv}$ |
| Eq. (18), p. 5 | $(G_{lm})=\rho_i\rho_j\rho_k\{\tilde G_{lm}J_e\}$ | Neko `G11`, `G22` in `svv_local` |
| Eq. (22)–(28), p. 6–7 | Legendre basis; $\hat Q_{ii}=(k/N)^{N_{svv}}$; $\tilde D=\hat B\hat Q\hat B^{-1}\hat D$; $\mathbf S^e_{vv}=\tilde D^TG\tilde D$ | `svv_init`, `svv_local` |
| Eq. (29), (31), p. 7 | $\mu=c_0\lvert\mathbf c\rvert\mathbb H/N$, $\mathbb H=2J_e^{1/d}$ | $\nu=c_0H/N$ times $\mathbf D_\mu=\lvert\operatorname{sgn}\psi\rvert$ |
| Eq. (34)–(35), p. 8 | BDF2/EXT2 | `redistance_standalone` |
| Eq. (42), p. 8 | $\mathbf n=\nabla\phi/\lvert\nabla\phi\rvert$ | `unit_normal` |
| Eq. (44)–(46), p. 9 | as §1 | `rd_sgn`, `redistance_standalone` |
| Eq. (83), p. 18 | exact solution and skewed IC | `circle_distance`, `skewed_ic` |
| Eq. (84), p. 19 | the global relative error | `ensure_psie`, `error_report` |
| §4.4, p. 18–19 | $\Omega$, $r$, $a$, $H$, $N$, $\Delta\tau$, $\tau$, $N_{svv}$, $c_0$; "the entire domain" | the four committed `.case` files; no band |

Routine by routine: `../../REDISTANCING.md` §9. The seventeen routines shared with the coupled
cases are unchanged.

## 4. Results

"Theirs" is `plot.py` L83–85 of their public case, the values behind Table 2 and Fig. 12.

**Table 2**, $E_r(\tau=6)$:

| cell | theirs | ours | ours / theirs |
|---|---|---|---|
| $H{=}1/10$, $N{=}3$ | 6.7613e-3 | 6.6156e-3 | 0.978 |
| $H{=}1/5$, $N{=}7$ | 3.8735e-3 | 3.8735e-3 | 1.000 |
| $H{=}1/10$, $N{=}7$ | 1.3484e-3 | 1.3486e-3 | 1.000 |
| $H{=}1/20$, $N{=}3$ | 1.8575e-3 | 1.8784e-3 | 1.011 |

**Fig. 12**, ours / theirs over the grid ($E_r(6)$ values in `run_circles_h*_n*.log`):

| | $N{=}3$ | $4$ | $5$ | $6$ | $7$ | $8$ |
|---|---|---|---|---|---|---|
| $H{=}1/5$ | 1.000 | fails | 1.000 | 1.000 | 1.000 | fails |
| $H{=}1/10$ | 0.978 | 0.995 | 1.000 | 0.999 | 1.000 | 1.000 |
| $H{=}1/20$ | 1.011 | 1.007 | 1.000 | 1.000 | 1.000 | 1.000 |

**Fig. 13**, $\lvert e\rvert$ at $\tau=6$, $N{=}8$, split into the "diamond" (the corner branch of
Eq. (83), vertices $(\pm0.7,0)$, $(0,\pm0.714)$) and the rest. Theirs is their code's own field,
run here:

| $N{=}8$ | $E_r$ | diamond | outside | outside L/R, T/B |
|---|---|---|---|---|
| $H{=}1/5$, theirs / ours (seeded, §5) | 1.6800e-3 / 1.6800e-3 | 3.960e-3 / 3.959e-3 | 1.784e-3 / 1.785e-3 | 1.00, 0.95 / 1.00, 0.95 |
| $H{=}1/10$, theirs / ours | 1.3150e-3 / 1.3149e-3 | 4.091e-3 / 4.091e-3 | 4.058e-4 / 4.051e-4 | 0.85, 1.15 / 0.82, 1.18 |
| $H{=}1/20$, theirs / ours | 5.1163e-4 / 5.1163e-4 | 1.584e-3 / 1.584e-3 | 1.658e-4 / 1.658e-4 | 0.55, 1.79 / 0.55, 1.79 |

**Beyond $\tau=6$.** The paper calls $\tau=6$ quasi-steady. It is not, in their code or ours:
$E_r$ has a minimum between $\tau\approx5$ and 6 and grows after it. From $\tau=6$ to 24 their
code's $E_r$ rises 1.9× at $H{=}1/10$ $N{=}3$, 1.7× at $H{=}1/5$ $N{=}7$, 2.7× at $H{=}1/10$
$N{=}7$ and 1.5× at $H{=}1/20$ $N{=}3$. Ours at $\tau=24$ equals theirs at $H{=}1/5$ $N{=}7$, is within
3% at the two $N{=}3$ cells, and blows up at $H{=}1/10$ $N{=}7$ (§5; `evidence/er_vs_tau.png`).

## 5. Limitations

- **Round-off decides some low-$N$ cells.** Dealiasing lets zero-set nodes move, and the guard
  then freezes nodes by the sign of $\psi_0$ there, which is $\pm10^{-15}$. Our mesh generator
  (genmeshbox) and theirs (genbox) differ at that level, and so do 9 of the 18 zero-set signs at
  $H{=}1/10$, $N{=}3$.
  - On our mesh $H{=}1/5$ $N{=}4$ departs from $\tau\approx4.25$ and $N{=}8$ from $\tau\approx5.25$;
    $H{=}1/10$ $N{=}7$ fails at $\tau\approx8.2$ (after the paper's $\tau=6$); $H{=}1/10$ $N{=}3$ is
    2.2% low.
  - Started from their $\psi_0$ at the zero-set nodes (the `seed` runs in `logs/saini_case/`), our
    code gives their values exactly in every one of these; at $H{=}1/10$, $N{=}7$ it follows their
    run to $\tau=24$. Their code, run here, holds in all of them.
  - Without the guard (dealiasing only) $H{=}1/5$ $N{=}4$ and $N{=}8$ hold, 9% and 0.6% from their
    values.
- **The field drifts.** Every signed distance function is a steady state of Eq. (44); nothing
  restores an interface's position, and the SVV beside the zero set moves it slowly. Dealiasing
  makes even a straight interface drift (without it, a tilted line is an exact discrete steady
  state).
- **The fixed $\varepsilon=0.25$ is intentional** (the author's reply to our 2026-09-29 question):
  the Eq. (47) seed already carries the $H$ scaling, and an $H$-free sign function is meant to give
  a mesh-independent re-distancing distance. He recommends it in general. Their shared `lvlSet.f`
  divides $\psi$ by the smallest domain extent before the sign function, so their other cases use
  $0.25L$; this case overrides that to $L=1$ (`gfac = 1.0`). Tested in 1D (`eps_1d/README.md`):
  0.25 slows every near-interface rate 3.75–40× rather than speeding anything up, and it delays the
  drift past $\tau=6$. In the coupled solve it relies on the code's $25H$ budget, not the paper's
  $2.5H$.
- **1D re-initialization** (the Eq. (47) seed, then Eq. (44)): at $N{=}8$, $\psi$ at the walls
  collapses to 0.06 after the front arrives (`evidence/anim4a_1d_reinit.mp4`). The coupled cases
  are periodic.

## 6. Running it

```bash
source ../../setup-env-cuda.sh
genmeshbox -2 2 -2 2 -0.1 0 20 20 1 .false. .false. .true.  && mv box.nmsh box20.nmsh
genmeshbox -2 2 -2 2 -0.1 0 40 40 1 .false. .false. .true.  && mv box.nmsh box40.nmsh
genmeshbox -2 2 -2 2 -0.1 0 80 80 1 .false. .false. .true.  && mv box.nmsh box80.nmsh
./run.sh                      # Table 2's four cells, GPU-guarded
./run.sh circles_h10_n3       # or one cell
```

- The committed `.case` files are Table 2's four cells.
- `scripts/gen_cases.py <outdir> 0.25 <Hden>:<N> …` emits the rest of the Fig. 12 grid.
- `"ic": "exact"` under `case.cdi` starts from $\psi_e$ instead of the skewed IC.
- Wall time to $\tau=6$ on the RTX 3090: $H{=}1/5$ 0.4–2 min, $H{=}1/10$ 1–8.5 min, $H{=}1/20$
  4–42 min.
- `redistance_circles.f90` is the user file; seventeen of its routines are shared byte-identically
  with the coupled cases (`../../CLAUDE.md`).
- `logs/saini_case/` (gitignored) holds the study of their code: its build and runs (`nek/`), the
  seeded runs, `results.txt`, `RUNS.txt`, the replica with their configuration
  (`replica/saini_scheme.py`) and the figure script (`figs44.py evidence`).

## 7. Evidence

Figures, from the shipped runs; Fig. 13 at $H{=}1/5$ uses the seeded run:
- `fig11_ic_and_exact.png`: their Fig. 11, the IC and the exact solution.
- `fig12_error_decay.png`: their Fig. 12, $E_r(6)$ against $N$ and $H$; theirs dashed.
- `fig13_error_maps.png`: their Fig. 13, $\lvert e\rvert$ at $N{=}8$ on their colour scale.
- `fig13_vs_saini.png`: their fields above ours.
- `er_vs_tau.png`: $E_r(\tau)$ to 24 in the Table 2 cells, against their code.

Animations (scripts `logs/saini_case/anim_new/`, gitignored, local):
- `anim1_1d.mp4`: 1D, the skewed IC relaxes onto $\lvert x\rvert-1$ ($N=3$, 8).
- `anim4a_1d_reinit.mp4`: 1D re-initialization; the wall collapse at $N{=}8$.
- `anim2_flat_vs_curved.mp4`: from $\psi_e$, a tilted line, one circle and the two circles drift,
  the line least.
- `anim3_fig13_over_tau.mp4`: Fig. 13 over $\tau$ at $N{=}8$ on the three meshes, both ICs, with
  $H{=}1/5$ on our mesh and seeded side by side.
- `anim3b_apex.mp4`: the steep-side cone apex at $H{=}1/20$, $N{=}3$ holds from both ICs.
- `anim4b_blowup_h10_n7.mp4`: $H{=}1/10$, $N{=}7$ after $\tau=6$: our mesh blows up, the seeded runs
  follow their code.
- `anim4c_guard_n3.mp4`, `anim4c_guard_n8_apex.mp4`: the guard against none.
- `talk/`: the same stories for a presentation, one message per clip, with stills of the key
  frames in `talk/snapshots/`.

`archive_2026-09-28_epsHN/` holds the earlier figures and animations, made with $\varepsilon=H/N$
before the authors' configuration was known.

Their raster figures are not reproduced here (CC BY-NC-ND).
