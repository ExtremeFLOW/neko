"""Two animations for README.md section 6, each with one message:
eps = H/N (left, blue) against Saini's fixed 0.25 (right, orange).

  anim_eps_A_standalone.mp4  same answer, but 0.25 gets there later
                             (psi_e = |x|-1, skewed IC, H=1/10, N=5, tau 0..8)
  anim_eps_B_seeded.mp4      from the Eq. (47) seed, 0.25 moves ~10x slower, so the
                             code needs 25H, not the paper's 2.5H (H=1/80, N=5)

python anim_eps_1d.py [A|B]   (both by default; about 30 s)
"""
import os, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from eps_mechanism_1d import SEM1D, solve, seed, HERE

INK, INK2, MUTED, GRID, AXIS, SURF = ("#0b0b0b", "#52514e", "#898781", "#e1e0d9",
                                      "#c3c2b7", "#fcfcfb")
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURF,
                     "figure.facecolor": SURF, "text.color": INK})


def style(ax):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def write(fig, update, frames, name, fps=20):
    w = FFMpegWriter(fps=fps, bitrate=2000)
    with w.saving(fig, os.path.join(HERE, name), dpi=110):
        for k in frames:
            update(k); w.grab_frame()
    plt.close(fig)


def header(fig, title, sub):
    fig.text(0.06, 0.955, title, fontsize=14, va="top")
    fig.text(0.06, 0.895, sub, fontsize=10.5, color=INK2, va="top")


# ---------------------------------------------------------------- A
def anim_a():
    Hden, N = 10, 5
    sem = SEM1D(4*Hden, N); x = sem.xg; H = sem.H
    psie = np.abs(x) - 1; ic = ((x - 1)**2 + 0.1)*psie
    den = float(sem.Bg @ np.abs(psie))
    runs = []
    for lab, eps, c in (("H/N", H/N, BLUE), ("0.25", 0.25, ORANGE)):
        _, snaps, _ = solve(sem, ic, eps, 5e-4, 8.0, snap_every=100)
        er = np.array([float(sem.Bg @ np.abs(p - psie))/den for _, p in snaps])
        runs.append((lab, c, snaps, er))
    taus = np.array([t for t, _ in runs[0][2]])
    m = (x >= 0) & (x <= 2); xs = x[m]

    fig = plt.figure(figsize=(11, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 0.8], hspace=0.38, wspace=0.12,
                          left=0.07, right=0.98, top=0.80, bottom=0.09)
    header(fig, r"Same answer, but $\varepsilon=0.25$ gets there later",
           r"Eq. (44) from a skewed start, $H=1/10$, $N=5$. Only $\varepsilon$ in "
           r"$\mathrm{sgn}(\psi)=\tanh(\psi/2\varepsilon)$ differs. A wider sign layer means a "
           "lower speed\n" r"$|w|=|\mathrm{sgn}\,\psi|$ near the interface ($x=1$).")
    lines, readouts = [], []
    for j, (lab, c, snaps, er) in enumerate(runs):
        a = fig.add_subplot(gs[0, j]); style(a)
        a.plot(xs, psie[m], ":", color=MUTED, lw=1.8, label="exact distance")
        ln, = a.plot(xs, snaps[0][1][m], color=c, lw=2.4, label=r"$\psi$")
        a.set_xlim(0, 2); a.set_ylim(-1.1, 1.1); a.set_xlabel("$x$")
        a.set_title(r"$\varepsilon=H/N=0.02$" if lab == "H/N" else r"$\varepsilon=0.25$",
                    color=c, fontsize=12)
        readouts.append(a.text(0.03, 0.95, "", transform=a.transAxes, va="top", color=c))
        if j == 0:
            a.set_ylabel(r"$\psi$"); a.legend(frameon=False, loc="lower right")
        lines.append(ln)
    e = fig.add_subplot(gs[1, :]); style(e)
    trails = []
    for lab, c, snaps, er in runs:
        e.semilogy(taus, er, color=c, lw=1, alpha=0.25)
        trails.append(e.semilogy([], [], color=c, lw=2.4, label=fr"$\varepsilon={lab}$")[0])
    e.axvline(6, color=MUTED, ls="--", lw=1)
    e.text(6.08, 0.3, r"$\tau=6$ (paper)", color=MUTED, fontsize=9)
    e.set_xlim(0, 8); e.set_ylim(1e-5, 3)
    e.set_xlabel(r"pseudo-time $\tau$"); e.set_ylabel("error $E_r$")
    e.legend(frameon=False, loc="upper right")
    clock = fig.text(0.98, 0.955, "", fontsize=14, ha="right", va="top")

    def update(k):
        for (lab, c, snaps, er), ln, rd, tr in zip(runs, lines, readouts, trails):
            ln.set_ydata(snaps[k][1][m])
            rd.set_text(f"$E_r$ = {er[k]:.1e}")
            tr.set_data(taus[:k+1], er[:k+1])
        clock.set_text(fr"$\tau={taus[k]:.2f}$")
    n = len(taus)
    write(fig, update, list(range(n)) + [n - 1]*40, "anim_eps_A_standalone.mp4")


# ---------------------------------------------------------------- B
def anim_b():
    Hden, N = 80, 5
    sem = SEM1D(4*Hden, N); x = sem.xg; H = sem.H; T = 25*H
    phi0 = seed(sem)
    dt = 0.05*sem.hmin; nst = int(np.ceil(T/dt)); dt = T/nst
    every = max(1, nst//200)
    i1 = np.argmin(np.abs(x - 1))
    runs = []
    for lab, eps, c in (("H/N", H/N, BLUE), ("0.25", 0.25, ORANGE)):
        _, snaps, _ = solve(sem, phi0, eps, dt, T, snap_every=every)
        runs.append((lab, c, snaps))
    th = np.array([t for t, _ in runs[0][2]])/H          # tau in units of H
    s = (x - 1)/H; m = (s > -6) & (s < 30); ss = s[m]

    fig = plt.figure(figsize=(11, 5.0))
    gs = fig.add_gridspec(1, 2, wspace=0.12, left=0.07, right=0.98, top=0.74, bottom=0.13)
    header(fig, r"From the seed, $\varepsilon=0.25$ moves about 10$\times$ slower",
           r"One re-distancing event from Eq. (47), $\phi_0=0.1(\psi_{CLS}-0.5)$, so "
           r"$|\phi_0|\leq0.05$ and, at $\varepsilon=0.25$, $|\mathrm{sgn}\,\phi_0|\leq0.1$."
           "\n" r"The paper budgets $2.5H$ of pseudo-time assuming speed 1; Saini's code "
           r"runs $25H$. $H=1/80$, $N=5$; both axes in units of $H$.")
    lines, readouts = [], []
    for j, (lab, c, snaps) in enumerate(runs):
        a = fig.add_subplot(gs[0, j]); style(a)
        a.plot(ss, ss, ":", color=MUTED, lw=1.8, label="exact distance")
        a.plot(ss, phi0[m]/H, "--", color=MUTED, lw=1.1, label=r"seed $\phi_0$")
        for k, t in ((2.5, "2.5H\npaper"), (25, "25H\ncode")):
            a.axvline(k, color=AXIS, lw=1)
            a.text(k + 0.4, -5.3, t, color=MUTED, fontsize=8.5)
        ln, = a.plot(ss, snaps[0][1][m]/H, color=c, lw=2.4, label=r"$\phi$")
        a.set_xlim(-6, 30); a.set_ylim(-6, 30)
        a.set_xlabel("distance from interface / $H$")
        a.set_title(r"$\varepsilon=H/N=0.0025$" if lab == "H/N" else r"$\varepsilon=0.25$",
                    color=c, fontsize=12)
        readouts.append(a.text(0.03, 0.95, "", transform=a.transAxes, va="top", color=c))
        if j == 0:
            a.set_ylabel(r"$\phi\,/\,H$"); a.legend(frameon=False, loc="lower right")
        lines.append(ln)
    clock = fig.text(0.98, 0.955, "", fontsize=14, ha="right", va="top")

    def update(k):
        for (lab, c, snaps), ln, rd in zip(runs, lines, readouts):
            p = snaps[k][1]
            ln.set_ydata(p[m]/H)
            rd.set_text(fr"slope at interface $|\phi_x|$ = {abs(sem.grad_avg(p))[i1]:.2f}"
                        "\n(target 1)")
        tag = "   paper's budget" if abs(th[k] - 2.5) < 0.07 else \
              "   code's budget" if k == len(th) - 1 else ""
        clock.set_text(fr"$\tau={th[k]:4.1f}\,H$" + tag)
    k25 = int(np.argmin(np.abs(th - 2.5))); n = len(th)
    frames = list(range(k25 + 1)) + [k25]*40 + list(range(k25 + 1, n)) + [n - 1]*50
    write(fig, update, frames, "anim_eps_B_seeded.mp4")


if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B"]
    if "A" in which: anim_a()
    if "B" in which: anim_b()
