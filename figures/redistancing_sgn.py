"""REDISTANCING.md §9.2: what sgn(psi) = tanh(psi/2eps) is, and what it does in
Eq. (44). Panel (b) is a fine-grid first-order upwind (Godunov) solve of the 1D
equation, used only to show the equation's behaviour -- it is not our scheme."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.polynomial import legendre as L

INK, INK2, MUTED, GRID, AXIS, SURF, BAND = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#f0efec"
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
H, N = 0.1, 3
eps = H/N
sgn = lambda p: np.tanh(p/(2*eps))

plt.rcParams.update({"font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.facecolor": SURF, "figure.facecolor": SURF, "text.color": INK})
fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.6))
for ax in (a, b):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

# (a) the function
p = np.linspace(-0.2, 0.2, 801)
a.axvspan(-2*eps, 2*eps, color=BAND, zorder=0)
a.step([-0.2, 0, 0.2], [-1, 1, 1], where="post", color=MUTED, lw=1.5, ls="--", label="sharp sign")
a.plot(p, sgn(p), color=BLUE, lw=2, label=r"sgn$(\psi)=\tanh(\psi/2\varepsilon)$")
c = np.zeros(N+1); c[N] = 1
z = np.concatenate(([-1.0], np.sort(L.legroots(L.legder(c))), [1.0]))
faces = -0.16 + H*np.arange(4)
nodes = np.unique(np.round(np.concatenate([f + (z + 1)*H/2 for f in faces]), 12))
nodes = nodes[(nodes >= -0.2) & (nodes <= 0.2)]
a.plot(nodes, sgn(nodes), "o", ms=8, color=ORANGE, mec=SURF, mew=2, zorder=5,
       label="GLL nodes ($H=1/10$, $N=3$)")
inside = np.sum(np.abs(nodes) < 2*eps)
a.annotate("interface, $\\psi=0$:\nsgn = 0, so $\\psi_\\tau=0$ there\n(the interface cannot move)",
           xy=(0, 0), xytext=(0.035, -0.62), color=INK2, fontsize=9,
           arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
a.text(-0.195, -1.15, "inside, $\\psi<0$: sgn \u2192 \u22121", color=INK2, fontsize=9)
a.text(0.195, 0.60, "outside, $\\psi>0$: sgn \u2192 +1", color=INK2, fontsize=9, ha="right")
a.text(0, 1.1, "sign layer $|\\psi|<2\\varepsilon$: %d nodes" % inside, color=INK2, fontsize=9, ha="center")
a.set_xlim(-0.2, 0.2); a.set_ylim(-1.2, 1.25)
a.set_xlabel(r"$\psi$  (for a distance function: the distance from the interface)")
a.set_ylabel(r"sgn$(\psi)$")
a.set_title(r"(a) the smoothed sign, $\varepsilon=H/N=1/30$: slope $1/2\varepsilon=15$ at $\psi=0$",
            loc="left", fontsize=10, color=INK)
a.legend(loc="upper left", frameon=False, fontsize=8.5, bbox_to_anchor=(0.0, 0.93))

# (b) what it does: psi_t = sgn(psi)(1 - |psi_x|) from psi_0 = 3x, fine-grid Godunov solve
h = 2.5e-4; x = np.arange(-0.4, 0.4 + h/2, h); u = 3*x; dt = 0.5*h
taus = (0.0, 0.04, 0.08, 0.12); out = {0.0: u.copy()}
for it in range(1, int(round(taus[-1]/dt)) + 1):
    am = np.empty_like(u); bp = np.empty_like(u)
    am[1:] = (u[1:] - u[:-1])/h; am[0] = am[1]; bp[:-1] = am[1:]; bp[-1] = bp[-2]
    s = sgn(u)
    Gp = np.sqrt(np.maximum(np.maximum(am, 0)**2, np.minimum(bp, 0)**2))
    Gm = np.sqrt(np.maximum(np.minimum(am, 0)**2, np.maximum(bp, 0)**2))
    u = u - dt*s*(np.where(s > 0, Gp, Gm) - 1)
    u[0] = 2*u[1] - u[2]; u[-1] = 2*u[-2] - u[-3]
    for t in taus:
        if abs(it*dt - t) < dt/2: out[t] = u.copy()
m = np.abs(x) <= 0.2
b.plot(x[m], x[m], color=MUTED, lw=1.5, ls="--", label=r"target: distance $\psi=x$ (slope 1)")
for t, col in zip(taus, RAMP):
    b.plot(x[m], out[t][m], color=col, lw=2, label=r"$\tau=%.2f$" % t)
    b.text(0.203, out[t][m][-1], r"$\tau=%.2f$" % t, color=INK2, fontsize=8.5, va="center")
b.plot([0], [0], "o", ms=8, color=ORANGE, mec=SURF, mew=2, zorder=5)
b.annotate("", xy=(0.15, 0.20), xytext=(0.15, 0.43), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.2))
b.text(0.143, 0.47, "sgn = +1, slope 3 > 1:\npulled down", color=INK2, fontsize=9, ha="right")
b.annotate("", xy=(-0.15, -0.20), xytext=(-0.15, -0.43), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.2))
b.text(-0.143, -0.53, "sgn = \u22121: pulled up", color=INK2, fontsize=9)
b.text(0.012, -0.075, "sgn = 0: stays at 0", color=INK2, fontsize=9)
b.set_xlim(-0.2, 0.235); b.set_ylim(-0.62, 0.62)
b.set_xlabel(r"$x$"); b.set_ylabel(r"$\psi$")
b.set_title(r"(b) what it does in $\psi_\tau=$ sgn$(\psi)\,(1-|\psi_x|)$, from a too-steep $\psi_0=3x$",
            loc="left", fontsize=10, color=INK)
b.legend(loc="upper left", frameon=False, fontsize=8.5)
fig.tight_layout()
fig.savefig(__file__.replace(".py", ".png"), dpi=160)
print("nodes in the sign layer:", inside, "| psi(0.2) at tau:", {t: round(float(out[t][m][-1]), 4) for t in taus})
