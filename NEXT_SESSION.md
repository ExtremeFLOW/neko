# What is open

The focus is redistancing on Rider–Kothe: whether a $\psi$ maintained by Eq. (44) can give the
compression term a better normal than the transported one. Settled results live in the case
READMEs and `CDI_METHOD.md`, standing rules in `CLAUDE.md`. When an item is finished, its result
goes there and the item leaves this file. The records with every prediction and verdict so far are
`examples/rider_kothe/logs/d5/PREREGISTERED.txt` and
`examples/zalesak_disk/logs/armC_2026-10-02/PREREGISTERED.txt` (gitignored, local).

## Rider–Kothe redistancing

D5 and D6 are done (`CDI_METHOD.md` §4.1d items 5–6; settings, tables and findings in
`examples/rider_kothe/README.md`, "Redistancing"). They ran on the scratch user files
`examples/rider_kothe/logs/d5_2026-10-06/rider_kothe_0c.f90` and
`examples/rider_kothe/logs/d6_2026-10-07/rider_kothe_d6.f90`. With their keys at the defaults both
reproduce the committed file's output byte for byte to $t=0.56$ (the README says what was compared).

1. **Which events path should the coupled files carry?** A decision for the user; it touches the
   shared routines.
   - **BDF2/EXT2 with the printed Eq. (31) SVV** (`svv_step_eq31`) in place of SSP-RK3 with
     `svv_step_imp`. It removes the event-made zero sets ($E_r(8)$ 0.835 → 0.0960). Adopting it
     makes `svv_step_eq31` and the BDF2 loop part of the three coupled files, so the shared-routine
     set changes. Check it on `advecting_slab_1d` first, then on Zalesak arm C (§4.1c).
   - **Saini's sign width 0.25 with $25H$** on top of that: 0.0342 against transport's 0.0452
     ($H=1/64$), 0.0085 against 0.0104 ($H=1/128$). But $\psi$ is then not a distance in the band;
     it is a smoothed copy of $\phi$'s interface, refreshed every 0.5. That is a method decision,
     not a knob.
   - **Saini's dealiasing with his sign guard** on top of that reproduces his settings (0.0402,
     against `rk_eLf`'s 0.0407 on the older build) and costs 18% against 0.0342, with $\psi$'s normal 2.4–4.2× further
     off at $t=2$–5. Being faithful to his code is worse here. Which of the two carries the cost is
     not measured; they ran as a pair.
2. **Why the rebuilt $\psi$'s normal is worse than transport's** against the exact interface, in
   every events configuration (with 0.25 and $25H$: 1.0–2.6° at $t=2$–5, 29° at $t=8$; transport
   0.6–0.8° and 1.3°). Against $\phi$'s own contour it is closer than transport's, and $\phi$'s own
   normal is off the exact interface by about as much (1.4–2.3° at $t=3$–5, 27° at $t=8$;
   transport-only $\phi$ 0.9–1.0°, 2.0°). So the question is why $\phi$ drifts from the exact
   interface once $\psi$ follows it. Where, measured by region: the tail tip, where $\phi$'s contour
   retreats behind the exact one (visible in `evidence/rider_redistancing_events.mp4`)? The
   end-time rise, which every 0.25/$25H$ run and Saini's own code share?
3. **The thin tail.** Where the filament is thinner than $2\varepsilon$ (5–6% of its length at
   maximum stretch), $\phi$ has no 0.5 contour, so no reseed setting can rebuild $\psi$ there. Each
   candidate is a change to the method, not a knob:
   - a compression flux masked to the interface band;
   - the monotone transform $\psi \leftarrow L\tanh(\psi/L)$ (`REDISTANCING.md` §8).
4. **A convergence measure for the $\tau$ solve:** a residual on the compression band, not the
   build band (`CDI_METHOD.md` §4.3).
5. **Housekeeping.**
   - `rider_h192.case` ships but has never run: run it or remove it.
   - Run $\xi=0.75$ under strain (~30 min) to see whether the cross turns over.
   - The $H=1/64$ runs end at $t=8.00008$ (100001 steps; Neko's summed time is a round-off below 8
     after 100000), so $E_r(8)$ is one step past the reversal. Decide whether to stop exactly at 8.

## Saini §4.5, for reference

From the paper (rendered p. 21, checked 2026-10-05). Their field names are transposed: their
$\phi$ is our $\psi$ (`CDI_METHOD.md` §1).

| | Saini §4.5 | his `circVortex` at $H=1/64$, $N=5$ (`logs/d5/saini/cv64n5`) |
|---|---|---|
| velocity, time | Eqs. (85)–(86), $T=8$ | same |
| domain, disk | $\Omega=[0,1]^2$, centre $(0.5,0.75)$, $r=0.15$ | same |
| $\xi$ | their $1/N$, i.e. our $\xi=1$ | same |
| mesh, order | $H\in\{1/32,1/64,1/128\}$, $N\in\{4,5,6\}$ | |
| $\Delta t$ | $\{8,4,2\}\times10^{-4}$ (CFL $\approx0.4$ at $N=6$); BDF3/EXT3 at $H=1/128$, else BDF2/EXT2 | $4\times10^{-4}$, BDF2 |
| SVV, CLS and TLS advection, CLS re-initialization Eq. (39) | $N_{svv}=N/2$, $c_0=0.1$; $|\mathbf c|$ not stated per equation | same; $|\mathbf c|$ the local $|\mathbf u|$, and $|\mathbf n|=1$ for Eq. (39) |
| SVV, TLS re-distancing Eq. (44) | $N_{svv}=N/4$, $c_0=2.0$ | same; $|\mathbf c|=|\operatorname{sgn}\psi|$ at each node |
| $\Delta t_{tls}$ | 0.5, "including at the initial step" | same |
| $\Delta\tau_{tls}$, steps | $H/(N{+}1)$, $N_{tls}=2.5H/\Delta\tau_{tls}$ (§3.4) | $H/(N{+}1)$, 150 steps |
| re-distancing extent | $2.5H$ from the interface | $25H$ |
| $\Delta t_{cls}$ | 0.05; no counterpart here (CDI is fused) | 0.05, also at $t=0$ |
| $\Delta\tau_{cls}$, steps | $0.1H/(N{+}1)$, $N_{cls}=\varepsilon/\Delta\tau_{cls}$ | same, 12 steps |
| pseudo-CFL | $\approx0.24$ at $N=6$, one number for both equations | Nek's CFL of the Eq. (39) step at $|\mathbf c|=1$: 0.2006 at $N=5$, which scales to 0.24 at $N=6$; the Eq. (44) step at $|\mathbf c|=1$ would be 10× that |

His §5.1 defaults differ from §4.5: CLS re-initialization $N/4$, $c_0=1.0$; TLS re-distancing $N/6$,
$c_0=1.0$.

What his figures measure:
- Fig. 16: $\phi$ (his CLS) and $\psi$ (his TLS) along $y=0.75$ at $t=8$, $H=1/128$, $N=4$–6.
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

- **Zalesak arm C with the printed Eq. (31) SVV.** Its first spurious pieces are $\psi$ loops
  that the $t=0$ build leaves beside the slot's top corners (`CDI_METHOD.md` §4.1c). Whether the
  printed form removes them there, as it does on Rider–Kothe, is not measured.
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
