# Rider–Kothe: the redistancing study (not run)

**Superseded (2026-10-05).** Every number in this file predates the two 2026-10-02 fixes:
the frozen CDI diffusion and the redistancing events' broken time history
(`CDI_METHOD.md` §4.1d).
- The current state is in `examples/rider_kothe/README.md`.
- What is still open is item 2 of `NEXT_SESSION.md`.
- This file is history only.

**Status (2026-09-30).** Everything with redistancing off is done and in
`examples/rider_kothe/README.md`: the baseline, the $(\xi,\gamma)$ cross, the
$h$-series and the SVV runs. The redistancing study below has **not been run**. The
one earlier periodic-reseed run, `rk_h64_n4` (2026-09-08, built $\psi$ + Eq. (47)
reseeds), diverged at $t=5.76$ of 8, the band mean $|\nabla\psi|$ after each event
collapsing 1.277 → 0.058 over eight events. That is not Zalesak arm C's signature,
where the band mean drifts *up* (`CDI_METHOD.md` §4.1c).

**Its redistancing settings wait for D5** (`NEXT_SESSION.md`). Saini's own settings were
run here on 2026-10-01, at $H=1/64$, $N=5$, $\xi=1$ (`examples/rider_kothe/README.md`). Every
events arm diverged. Two causes are now known: a stale-BDF-history bug in our events, and $\phi$
losing the filament at $t\approx3.8$ (`CDI_METHOD.md` §4.1d; open items in `NEXT_SESSION.md` D5).

## Why this case

`advecting_slab_1d` and `zalesak_disk` are isometries: $|\nabla\psi|$ is provably
conserved, so redistancing has nothing to correct there. `rider_kothe` is the only
case with real strain (normal strain rates to $\pm2.5$). With redistancing off, the
band $|\nabla\psi|$ swings three decades and returns:

| $t$ | band min | mean | max |
|---|---|---|---|
| 0.08 | 0.699 | 1.020 | 1.279 |
| **4.0** (max stretch) | **0.0043** | **6.10** | **21.8** |
| **8.0** (returned) | **0.9995** | **1.000** | **1.001** |

## What Saini actually do here — their §4.5

Their Rider–Kothe setup, read from the paper rather than from summaries:

| | Saini §4.5 | ours |
|---|---|---|
| velocity, domain, disk | Eqs. (85)–(86), $T=8$, centre $(0.5,0.75)$, $r=0.15$ | **identical** |
| $\xi$ | $1/N$ — i.e. **our $\xi = 1.0$** | `rider_kothe_xi10` |
| mesh, $N$ | $H\in\{1/32,1/64,1/128\}$, $N\in\{4,5,6\}$ | $H=1/64$, $N=5$ |
| $\Delta t$ | $\{8,4,2\}\times10^{-4}$, CFL $\approx0.4$ at $N=6$ | ours is set by the CDI compression guard instead |
| SVV, **advection** (both fields) | $N_{svv}=N/2$, **$c_0 = 0.1$** | `svv_psi` off / 0.1 / 1.0 run with redistancing off (case README) |
| SVV, TLS re-distancing | $N_{svv}=N/4$, $c_0 = 2.0$ | matches `svv_rd` default |
| $\Delta t_{tls}$ | **0.5, "including at the initial step"** | `psi_init="redistance"` + periodic events |
| $\Delta t_{cls}$ | 0.05 | **no counterpart** — CDI is fused, their Eq. (39) was never reintroduced |
| band | $2.5H$ | same |

**Two things to carry from this.** First, **their $c_0$ here is 0.1, not the 1.0
we used on Zalesak** — 1.0 is the §4.3 value and does not transfer. Second,
unlike §4.3 they *do* redistance here, and say why: *"Re-initialization of the
CLS field, using Eq. (39), is essential to maintain a sharp tanh profile"*
under stretching.

Remember the notation transposition (`CDI_METHOD.md` §1): their $\phi$ is our
$\psi$, their $\psi$ is our $\phi$.

## The runs

A convergence study of one method, varying only $H$ and $N$ as Saini's §4.5 does.
$\Delta t$ is set by the CDI compression guard
$\gamma u_{\max}\Delta t/h_{\text{GLL,min}}\le0.05$ (10% margin), not Saini's
$\{8,4,2\}\times10^{-4}$. Generated `.case` files go in the scratchpad; they don't ship.

| cell | $H$ | $N$ | $\Delta t$ | steps | est. | why this order |
|---|---|---|---|---|---|---|
| `rk_h64_n3` | 1/64 | 3 | 1.94e-04 | 41 168 | **3 min** | repo convention; + Saini Table 3 |
| `rk_h64_n4` | 1/64 | 4 | 1.21e-04 | 65 893 | 11 min | Saini §4.5 range |
| `rk_h64_n5` | 1/64 | **5** | 8.26e-05 | 96 856 | 28 min | settled base, shared with $h$-sweep |
| `rk_h64_n6` | 1/64 | 6 | 5.97e-05 | 134 035 | 60 min | Saini §4.5 range |
| `rk_h64_n7` | 1/64 | 7 | 4.51e-05 | 177 419 | 119 min | repo convention; + Saini Table 3 |
| | | | | | **3.7 h** | |
| `rk_h64_n5_nord` | 1/64 | 5 | 8.26e-05 | 96 856 | 28 min | the redistancing control |
| `rk_h96_n5` | 1/96 | 5 | 5.51e-05 | 145 283 | 93 min | $h$-refinement |
| `rk_h128_n5` | 1/128 | 5 | 4.13e-05 | 193 714 | 220 min | $h$-refinement, **optional** |
| `rk_h128_n3` | 1/128 | 3 | 9.70e-05 | 82 474 | 28 min | Saini Table 3, **optional** |

- Run the $p$-sweep cheapest first. $\{4,5,6\}$ is their Fig. 17, $\{3,5,7\}$ the
  Zalesak convention, and $\{3,7\}$ at equal GLL point counts (`rk_h64_n7`,
  `rk_h128_n3`: 262 144 each) is their Table 3.
- `rk_h64_n5_nord` is `rk_h64_n5` with `redistance.enabled = false`: the only matched
  pair that says what redistancing buys under strain.
- **Recheck the band bound** $9.21\,\xi H/N \le \text{band}\cdot H$ for every cell
  (`CDI_METHOD.md` §4.1b). At band 2.5 that is $N\ge3.68\xi$, so $N=3$ at $\xi=1$ fails it.
- If SVV on $\psi$ looks over-diffused (a low $\phi$ peak at $t=4$, $E_r$ that stops
  improving with $N$), re-run `rk_h64_n5` at Saini's $c_0=0.1$.

## State the predictions before running

- **The case README's claim that redistancing at $t=4$ would be "actively wrong"
  — that it "would reset $|\nabla\psi|$ at exactly the moment the field is
  supposed to be strained, destroying the information that carries it home" — is
  probably overstated.** The normal depends only on $\psi$'s *level sets*
  (`REDISTANCING.md` §1), and `seed="phi"` re-registers the zero contour to
  $\phi = 0.5$ at every event. What redistancing discards is the far-field level
  *spacing*, which the compression term never reads. The `rk_h64_n5` vs
  `rk_h64_n5_nord` pair tests this directly, and the README must be corrected
  either way.
- If redistancing helps anywhere in this repo, it is here, at $t \approx 4$,
  where the band minimum reaches 0.0043 and the mean 6.10.
- **Errors should fall with both $h$ and $p$.** Saini see that, but also note
  that TLS redistancing *slows the convergence rate* relative to their Zalesak
  case — their §4.5 says the $|E_v|$ decrease with $N$ "has a slower rate as
  compared to the Zalesak problem ... likely due to the errors introduced by the
  TLS re-distancing equation, which is inactive for the Zalesak case." So a
  shallower slope here than in `zalesak_disk` is **expected**, not a defect.

## What to measure — Saini's §4.5 figures, recreated

$E_r$ is quoted **only at $t=8$**: the flow reverses, so the exact solution there
is the initial condition. At intermediate times the shape is a spiral with no
closed form. Everything else below needs no reference and is valid throughout.

Their figures, and what each maps to here:

| Saini | quantity | ours |
|---|---|---|
| Fig. 17 | $E_r$, $E_v$, $E_s$ vs **$p$** at $t=8$ | `rk_h64_n{4,5,6}` |
| Fig. 18 | same vs **$h$** (they add $E_{L_1}$ for literature comparison) | `rk_h{64,96,128}_n5` |
| Table 3 | error at **identical GLL points**, $N=\{3,7\}$ | optional extra pair |
| Fig. 19 | $l_{avg}/l_0$, mean interface thickness vs $t$ (their Eq. 87, bounds $\phi\in[0.05,0.95]$) | already a column in this case's README |
| Fig. 20 | $L_\infty$ of $\phi$ vs $t$ — boundedness | from `bnd` log lines |
| Fig. 21 | the 0.5 isocontour at $t=8$, all cases, against exact | contour plot |
| Fig. 22 | $\lvert E_v\rvert$ vs $t$ | per-frame |

Plus the two diagnostics specific to this repo's question:

1. **band $|\nabla\psi|$ (min/mean/max) vs $t$** — already logged every
   `band_report_every` steps. The go/no-go number and the thing redistancing acts
   on. Compare the redistancing pair through $t=4$ and back.
2. **$\phi$ peak at $t=4$** — the filament-resolution measure. At $\xi=1.0$,
   $H=1/64$ it is 0.344 in the old no-SVV baseline; anything that raises it is
   resolving the filament better. This is the number that showed the $\xi$
   trade-off, and it should improve with $h$-refinement (0.344 → 0.433 → 0.506).

Wall time is worth recording as `rk_h64_n5` minus `rk_h64_n5_nord` — both carry
identical diagnostics, so the difference is what the 16 events actually cost.
