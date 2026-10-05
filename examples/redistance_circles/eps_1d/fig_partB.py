"""Part B figures from results.txt (T = 2.5H), results_b25.txt and b25.pkl
(T = 25H, the budget in Saini's core lvlSet.f)."""
import os, re, pickle, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
INK2, MUTED, GRID, AXIS, SURF = "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURF,
                     "figure.facecolor": SURF})


def style(ax):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def parse(fn):
    pat = re.compile(r"N=(\d) H=1/(\d+)\s+eps=(\S+)\s.*?reach=([\d.]+) H")
    out = {}
    for line in open(os.path.join(HERE, fn)):
        m = pat.search(line)
        if m:
            out[(int(m[1]), int(m[2]), m[3])] = float(m[4])
    return out

# (d) reach against H
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6))
for ax, (fn, bud) in zip(axs, (("results.txt", 2.5), ("results_b25.txt", 25.0))):
    style(ax); r = parse(fn)
    for N, c in zip((3, 5, 8), RAMP[1:]):
        for lab, mk, ls in (("H/N", "o", "-"), ("0.25", "s", "--")):
            Hd = [h for h in (10, 20, 40, 80) if (N, h, lab) in r]
            ax.plot([1/h for h in Hd], [r[(N, h, lab)] for h in Hd], ls, marker=mk,
                    color=c, lw=1.6, label=fr"$N={N}$, $\varepsilon={lab}$")
    ax.axhline(bud, color=MUTED, ls=":", lw=1.2)
    ax.text(0.1, bud, f"budget {bud:g}H", color=MUTED, va="bottom", ha="right", fontsize=8.5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("$H$")
    ax.set_title(fr"pseudo-time budget $T={bud:g}H$", color=INK2)
    if bud == 25.0:
        ax.axvspan(1/22, 0.11, color="#f0efec", zorder=0)
        ax.text(0.068, 1.6, "front reaches\nthe wall ($25H>1$)", color=MUTED,
                fontsize=8, ha="center", va="bottom")
axs[0].set_ylabel("reach of a characteristic, in units of $H$")
h, l = axs[1].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=8.5, ncol=6, loc="lower center")
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(os.path.join(HERE, "partB_reach.png"), dpi=150)
plt.close(fig)

# (e) profiles after 25H, N=5
d = pickle.load(open(os.path.join(HERE, "b25.pkl"), "rb"))
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
for ax, lab in zip(axs, ("H/N", "0.25")):
    style(ax)
    for Hd, c in zip((10, 20, 40, 80), RAMP):
        phi0, phi, x, H, eps = d[(5, Hd, lab)]
        m = (x > 0.7) & (x < 1.5)
        ax.plot(x[m] - 1, phi[m], color=c, lw=1.8, label=fr"$H=1/{Hd}$")
    xx = np.linspace(-0.3, 0.5, 3)
    ax.plot(xx, xx, ":", color=MUTED, lw=1.2, label="distance")
    ax.plot(x[m] - 1, phi0[m], "--", color=MUTED, lw=1, label=r"$\phi_0$ ($H=1/80$)")
    ax.set_title(fr"$\varepsilon={lab}$, $N=5$, after $T=25H$", color=INK2)
    ax.set_xlabel("distance from interface $x-1$")
axs[0].set_ylabel(r"$\phi(T)$"); axs[0].set_ylim(-0.35, 0.5)
axs[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "partB_profiles_25H.png"), dpi=150)
plt.close(fig)
