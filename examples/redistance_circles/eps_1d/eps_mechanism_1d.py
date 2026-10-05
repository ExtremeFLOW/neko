"""What the sign-function width eps in Saini's Eq. (46), sgn(psi) = tanh(psi/2eps),
does to Eq. (44) in 1D: eps = H/N (the phase-field width) against the fixed 0.25 of
their code (lvlSet.f, signls). See README.md in this directory.

Part A: the standalone Sec. 4.4 analogue, psi_e = |x| - 1 on [-2, 2], skewed IC.
Part B: the coupled-solve seed, Eq. (47), phi0 = r_f (psi_CLS - 0.5), run for the
        automated budget T = 2.5H of Sec. 3.4 ("b25": the 25H of their core lvlSet.f;
        "fig16": their Rider-Kothe mesh, H = 1/128, against Figs. 15/16b).

Scheme: SSP-RK3 on sgn(psi)(1 - |grad psi|) (averaged gradient, no dealiasing),
then one implicit SVV step (B + dtau D_mu S) psi = B psi*, with D_mu = |sgn(psi^n)|
left-multiplying the assembled operator (Eq. 31 as printed, as in
redistance_circles). svv="const" sets D_mu = 1, which is sem1d.run().
The 1D SEM pieces come from the forensic replica (sem1d.py), not a copy.
"""
import sys, os, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
HERE = os.path.dirname(os.path.abspath(__file__))
# sem1d.py lives in the sibling repo neko-multiphase, next to this checkout
# (repo root is three levels up); NEKO_MULTIPHASE overrides its location.
_NM = os.environ.get("NEKO_MULTIPHASE",
                     os.path.join(HERE, "..", "..", "..", "..", "neko-multiphase"))
sys.path.insert(0, os.path.join(_NM, "examples", "saini_benchmarks", "redistance_eq44"))
from sem1d import SEM1D, svv_matrix, run as sem1d_run


def solve(sem, psi0, eps, dtau, tau_end, svv="dmu", c0=2.0, ratio=6.0,
          snap_every=None, psie=None):
    """Integrate Eq. (44) from psi0. Returns (psi, snapshots, trace)."""
    S = sp.csr_matrix(svv_matrix(sem, c0, ratio)) if svv != "off" else None
    B = sp.diags(sem.Bg)
    lu_const = spla.splu((B + dtau*S).tocsc()) if svv == "const" else None
    den = float(sem.Bg @ np.abs(psie)) if psie is not None else None

    def rhs(u):
        return np.tanh(u/(2*eps))*(1.0 - np.abs(sem.grad_avg(u)))

    psi = psi0.copy(); nstep = int(round(tau_end/dtau))
    snaps, trace = [(0.0, psi.copy())], []
    for it in range(1, nstep + 1):
        u0 = psi.copy()
        psi = psi + dtau*rhs(psi)
        psi = 0.75*u0 + 0.25*(psi + dtau*rhs(psi))
        psi = (1/3)*u0 + (2/3)*(psi + dtau*rhs(psi))
        if svv == "const":
            psi = lu_const.solve(sem.Bg*psi)
        elif svv == "dmu":
            Dmu = sp.diags(np.abs(np.tanh(u0/(2*eps))))
            psi = spla.spsolve((B + dtau*(Dmu @ S)).tocsc(), sem.Bg*psi)
        if not np.all(np.isfinite(psi)) or np.max(np.abs(psi)) > 1e6:
            trace.append((it*dtau, np.inf)); return psi, snaps, trace
        if snap_every and it % snap_every == 0:
            snaps.append((it*dtau, psi.copy()))
            if den:
                trace.append((it*dtau, float(sem.Bg @ np.abs(psi - psie))/den))
    return psi, snaps, trace


# ---------------------------------------------------------------- Part A
def part_a(out):
    out.write("# Part A: standalone, psi_e=|x|-1, skewed IC, dtau=5e-4\n")
    # verification against the replica (constant nu, eps=H/N)
    tr, _ = sem1d_run(40, 3, eps_xi=1.0, tau_end=6.0)
    sem = SEM1D(40, 3); x = sem.xg; psie = np.abs(x) - 1
    _, _, t2 = solve(sem, ((x-1)**2+0.1)*psie, sem.H/3, 5e-4, 6.0, svv="const",
                     snap_every=12000, psie=psie)
    out.write(f"verify: sem1d.run E_r(6)={tr[-1][1]:.6e}  this={t2[-1][1]:.6e}\n")

    res = {}
    for Hden in (10, 20):
        for N in (3, 5, 8):
            sem = SEM1D(4*Hden, N); x = sem.xg; psie = np.abs(x) - 1
            ic = ((x-1)**2 + 0.1)*psie
            for lab, eps in (("H/N", sem.H/N), ("0.25", 0.25)):
                _, snaps, tr = solve(sem, ic, eps, 5e-4, 24.0, snap_every=100,
                                     psie=psie)
                res[(Hden, N, lab)] = (sem, snaps, tr)
                d = dict(tr)
                out.write(f"H=1/{Hden} N={N} eps={lab:4s} ({eps:.4f})  "
                          + "  ".join(f"Er({t:g})={d[t]:.3e}" for t in
                                      (1.0, 2.0, 3.0, 6.0, 12.0, 24.0)) + "\n")
    # arrival of the correct slope at distance d outside x=1: first tau with
    # |psi - psi_e| < 1e-3 at every node in (1, 1+d]
    out.write("# arrival tau(d) (|psi-psi_e|<1e-3 on (1,1+d]) vs "
              "2eps ln[sinh(d/2eps)/sinh(d0/2eps)], d0=H/N\n")
    for (Hden, N, lab), (sem, snaps, tr) in res.items():
        if N != 5: continue
        x = sem.xg; psie = np.abs(x) - 1
        eps = sem.H/N if lab == "H/N" else 0.25
        for d in (0.1, 0.25, 0.5, 0.75):
            m = (x > 1) & (x <= 1 + d + 1e-12)
            ta = next((t for t, p in snaps if np.max(np.abs(p[m]-psie[m])) < 1e-3),
                      np.nan)
            d0 = sem.H/N
            tf = 2*eps*np.log(np.sinh(d/(2*eps))/np.sinh(d0/(2*eps)))
            out.write(f"H=1/{Hden} N={N} eps={lab:4s} d={d:.2f}  "
                      f"tau_arrive={ta:.2f}  formula={tf:.2f}\n")
    return res


# ---------------------------------------------------------------- Part B
def seed(sem, rf=0.1, xi=1.0):
    """Eq. (36) CLS (1 inside |x|<1) and Eq. (47); sign flipped so phi<0 inside
    like psi_e (Eq. 44 is odd in phi)."""
    epsH = xi*sem.H/sem.N
    cls = 0.5*(1 - np.tanh((np.abs(sem.xg) - 1)/(2*epsH)))
    return -rf*(cls - 0.5)


def part_b(out, budget=2.5, Hdens=(10, 20, 40, 80)):
    out.write(f"\n# Part B: Eq. (47) seed, rf=0.1, eps_H=H/N, budget T={budget}H\n")
    out.write("# s = (x-1)/H. reach = distance (in H) a characteristic "
              "dx/dtau=|sgn phi| starting at s=1/N travels in T\n")
    res = {}
    for N in (3, 5, 8):
        for Hden in Hdens:
            sem = SEM1D(4*Hden, N); x = sem.xg; H = sem.H
            phi0 = seed(sem); T = budget*H
            for lab, eps in (("H/N", H/N), ("0.25", 0.25)):
                hmin = sem.hmin
                dt = 0.05*hmin
                nst = int(np.ceil(T/dt)); dt = T/nst
                phi, snaps, _ = solve(sem, phi0, eps, dt, T,
                                      snap_every=max(1, nst//200))
                # characteristic tracing on the evolving field, outer side x>1
                xc = 1 + H/N; tprev = 0.0
                for t, p in snaps[1:]:
                    w = np.abs(np.tanh(np.interp(xc, x, p)/(2*eps)))
                    xc += w*(t - tprev); tprev = t
                reach = (xc - 1)/H - 1/N
                # Saini's own step, dtau = H/(N+1), N_tls = 2.5(N+1)
                dS = H/(N + 1)
                phiS, _, trS = solve(sem, phi0, eps, dS, T)
                okS = "blows up" if trS and not np.isfinite(trS[-1][1]) else \
                      f"max|phi-phi_fine|={np.max(np.abs(phiS-phi)):.2e}"
                m1 = np.argmin(np.abs(x - (1 + H)))
                gs = np.abs(sem.grad_avg(phi))
                gat = "  ".join(f"|grad|(s={k:g})={gs[np.argmin(np.abs(x-(1+k*H)))]:.3f}"
                                for k in (2.5, 5, 10))
                g = np.abs(sem.grad_avg(phi)); g0 = np.abs(sem.grad_avg(phi0))
                iface = np.argmin(np.abs(x - 1))
                res[(N, Hden, lab)] = (sem, phi0, phi, eps)
                out.write(
                    f"N={N} H=1/{Hden:<2d} eps={lab:4s} sgn(phi0) max="
                    f"{np.max(np.abs(np.tanh(phi0/(2*eps)))):.3f}  "
                    f"|w|(s=1)={abs(np.tanh(phi[m1]/(2*eps))):.3f}  "
                    f"reach={reach:.3f} H  |grad|@iface {g0[iface]:.3f}->"
                    f"{g[iface]:.3f}  |grad|(s=1) {g0[m1]:.3f}->{g[m1]:.3f}  "
                    f"max|dphi|/H={np.max(np.abs(phi-phi0))/H:.3f}  "
                    f"Saini dtau: {okS}  {gat}\n")
                out.flush()
    return res


# ---------------------------------------------------------------- figures
def figures(resA, resB):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    INK2, MUTED, GRID, AXIS, SURF = "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
    BLUE, ORANGE = "#2a78d6", "#eb6834"
    RAMP = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": AXIS,
                         "axes.labelcolor": INK2, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.facecolor": SURF,
                         "figure.facecolor": SURF})

    def style(ax):
        ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
        for s in ("top", "right"): ax.spines[s].set_visible(False)

    # (a) |w| = |sgn| against distance from the interface
    fig, ax = plt.subplots(figsize=(6.5, 4.2)); style(ax)
    d = np.linspace(0, 1.0, 801)
    for (Hd, N), c in zip(((5, 3), (10, 3), (20, 8)), RAMP[1:]):
        e = 1/(Hd*N)
        ax.plot(d, np.tanh(d/(2*e)), color=c, lw=2,
                label=fr"$\varepsilon=H/N$, $H=1/{Hd}$, $N={N}$ ({e:.3f})")
    ax.plot(d, np.tanh(d/0.5), color=ORANGE, lw=2.5, label=r"$\varepsilon=0.25$")
    ax.set_xlabel(r"distance from interface $|\psi|$")
    ax.set_ylabel(r"characteristic speed $|\mathbf{w}|=|\mathrm{sgn}\,\psi|$")
    ax.legend(frameon=False, fontsize=8.5); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "eps_speed.png"), dpi=150); plt.close(fig)

    # (b) Part A: E_r(tau)
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
    for ax, Hd in zip(axs, (10, 20)):
        style(ax)
        for N, c in zip((3, 5, 8), RAMP[1:]):
            for lab, ls in (("H/N", "-"), ("0.25", "--")):
                tr = resA[(Hd, N, lab)][2]
                t, e = np.array(tr).T
                ax.semilogy(t, e, ls, color=c, lw=1.8,
                            label=fr"$N={N}$, $\varepsilon={lab}$")
        ax.set_title(fr"$H=1/{Hd}$", color=INK2); ax.set_xlabel(r"$\tau$")
    axs[0].set_ylabel(r"$E_r=\int|\psi-\psi_e|/\int|\psi_e|$")
    axs[1].legend(frameon=False, fontsize=8, ncol=2); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "partA_er_tau.png"), dpi=150); plt.close(fig)

    # (c) Part B: phi(T) against s = (x-1)/H, N=5, both eps
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
    for ax, lab in zip(axs, ("H/N", "0.25")):
        style(ax)
        for Hd, c in zip((10, 20, 40, 80), RAMP):
            sem, phi0, phi, eps = resB[(5, Hd, lab)]
            s = (sem.xg - 1)/sem.H; m = (s > -4) & (s < 6)
            ax.plot(s[m], phi[m]/sem.H, color=c, lw=1.8, label=fr"$H=1/{Hd}$")
        ax.plot(s[m], s[m], ":", color=MUTED, lw=1.2, label="distance")
        ax.plot(s[m], phi0[m]/sem.H, "--", color=MUTED, lw=1, label=r"$\phi_0/H$ ($H=1/80$)")
        ax.set_title(fr"$\varepsilon={lab}$, $N=5$, after $T=2.5H$", color=INK2)
        ax.set_xlabel(r"$s=(x-1)/H$")
    axs[0].set_ylabel(r"$\phi(T)/H$"); axs[0].set_ylim(-4, 6)
    axs[1].legend(frameon=False, fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(HERE, "partB_profiles.png"), dpi=150); plt.close(fig)


if __name__ == "__main__" and sys.argv[1:] == ["b25"]:
    # Saini's core lvlSet.f: dt_tls = dxmin/lx1, nsteps_tls = floor(25*dxave/dt_tls)
    with open(os.path.join(HERE, "results_b25.txt"), "w") as out:
        rB = part_b(out, budget=25.0)
    import pickle
    pickle.dump({k: (v[1], v[2], v[0].xg, v[0].H, v[3]) for k, v in rB.items()},
                open(os.path.join(HERE, "b25.pkl"), "wb"))
elif __name__ == "__main__" and sys.argv[1:] == ["fig16"]:
    with open(os.path.join(HERE, "results_fig16.txt"), "w") as out:
        out.write("# H=1/128, eps=0.25, T=25H, dtau=H/(N+1); Fig. 16b: plateau 0.074, "
                  "slope ~1\n")
        for N in (4, 5, 6):
            sem = SEM1D(4*128, N); x = sem.xg; H = sem.H; T = 25*H
            phi0 = seed(sem); i0 = np.argmin(np.abs(x - 1))
            y0 = abs(sem.grad_avg(phi0))[i0]
            ylog = 1/(1 - (1 - 1/y0)*np.exp(-T/(2*0.25)))
            for svv, c0 in (("dmu", 2.0), ("off", 0.0), ("const", 2.0),
                            ("dmu", 20.0), ("const", 20.0)):
                n = int(round(T/(H/(N + 1))))
                phi, _, _ = solve(sem, phi0, 0.25, T/n, T, svv=svv, c0=c0)
                g = np.abs(sem.grad_avg(phi))
                out.write(f"N={N} svv={svv:5s} c0={c0:4.1f} steps={n} plateau="
                          f"{phi[np.argmin(np.abs(x-1.3))]:.4f} (0.05e^(50H)="
                          f"{0.05*np.exp(50*H):.4f})  |grad|iface {y0:.2f}->"
                          f"{g[i0]:.3f} (logistic {ylog:.3f})\n")
elif __name__ == "__main__":
    with open(os.path.join(HERE, "results.txt"), "w") as out:
        rA = part_a(out); out.flush()
        rB = part_b(out)
    figures(rA, rB)
