#!/usr/bin/env python3
"""Regenerate every figure used in the manuscript and README.

    python scripts/make_figures.py          # writes figures/*.pdf and figures/*.png
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from qtime import clock_design as cd  # noqa: E402
from qtime import clocks as qc  # noqa: E402
from qtime import estimates as est  # noqa: E402
from qtime import interferometer as ifm  # noqa: E402

OUT = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3, "savefig.bbox": "tight"})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=200)
    plt.close(fig)
    print("wrote", name)


# ----------------------------------------------------------------------------
# Fig. 1  Exact bipartite conditioning: sharp vs unsharp readout vs diffusive clock
# ----------------------------------------------------------------------------
def fig_conditioning():
    N, r, sigma = 256, 4, 0.15
    psi = qc.pw_constraint_state(N, [0, r], [1, 1])
    th = qc.phase_grid(N)
    p_sharp, rho_sharp = qc.condition_on_clock(psi, qc.sharp_phase_povm(N))
    p_uns, rho_uns = qc.condition_on_clock(psi, qc.smeared_phase_povm(N, sigma))

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    ax[0].plot(th, p_sharp * N, "k-", lw=1.5, label="sharp readout")
    ax[0].plot(th, p_uns * N, "C1--", lw=1.5, label="unsharp readout")
    ax[0].set_ylim(0, 1.6)
    ax[0].set_xlabel(r"clock reading $\theta$")
    ax[0].set_ylabel(r"$N\,p(\theta_k)$  (raw clock marginal)")
    ax[0].set_title("(a) exact constraint state: uniform marginal")
    ax[0].legend(loc="lower right", fontsize=8)

    ax[1].plot(th, np.real(rho_sharp[:, 0, 1]), "k-", lw=1.2, label=r"sharp: Re $\rho_{01}(\theta)$")
    ax[1].plot(th, np.real(rho_uns[:, 0, 1]), "C1-", lw=1.8, label=r"unsharp ($\sigma_\theta=0.15$)")
    amp = 0.5 * qc.coherence_suppression_static(r, sigma)
    ax[1].axhline(amp, color="C2", ls=":", lw=1.5, label=r"$\frac{1}{2} e^{-r^2\sigma_\theta^2/2}$")
    ax[1].axhline(-amp, color="C2", ls=":", lw=1.5)
    ax[1].set_xlabel(r"clock reading $\theta=\omega t$")
    ax[1].set_title("(b) unsharp readout: constant, not growing")
    ax[1].set_ylim(-0.62, 0.75)
    ax[1].legend(loc="upper right", fontsize=7, ncol=3)

    # diffusive clock model vs unsharp readout
    E = np.array([0.0, float(r)])
    rho0 = 0.5 * np.ones((2, 2))
    t = np.linspace(0, 2 * np.pi, 120)
    T_C = 0.01
    diff = qc.diffusive_clock_channel(rho0, E, t, T_C, n_quad=401)
    lind = qc.lindblad_dephasing_solution(rho0, E, t, T_C)
    ax[2].plot(th, np.abs(rho_uns[:, 0, 1]), "C1-", lw=1.8, label="unsharp readout (exact bipartite)")
    ax[2].plot(t, np.abs(diff[:, 0, 1]), "C0o", ms=3, label=r"diffusive clock, $\sigma_t^2=T_C t$ (quadrature)")
    ax[2].plot(t, np.abs(lind[:, 0, 1]), "C3-", lw=1, label="GKSL closed form")
    ax[2].set_xlabel(r"$\omega t$")
    ax[2].set_ylabel(r"$|\rho_{01}|$")
    ax[2].set_ylim(0, 0.55)
    ax[2].set_title("(c) Lindblad decay needs the diffusion assumption")
    ax[2].legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    save(fig, "fig1_conditioning")


# ----------------------------------------------------------------------------
# Fig. 2  Winding-number marginalisation on a cyclic clock
# ----------------------------------------------------------------------------
def fig_winding():
    nu = np.linspace(-0.25, 3.25, 3501)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for M, c in [(1, "0.6"), (5, "C0"), (40, "C3")]:
        vals = np.array([abs(qc.winding_marginal_coherence(x, M)) for x in nu])
        ax.plot(nu, vals, color=c, lw=1.4, label=f"M = {M} windings")
    ax.set_xlabel(r"$\nu=(E_j-E_k)/\hbar\omega$")
    ax.set_ylabel(r"surviving coherence $|\overline{\rho_{jk}}|/|\rho_{jk}|$")
    ax.set_title("Circle clock read without its winding number")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, "fig2_winding")


# ----------------------------------------------------------------------------
# Fig. 3  Relativistic clock interferometer (corrected plot: revivals shown)
# ----------------------------------------------------------------------------
def fig_interferometer():
    p1, p2, omega0, dn = 0.2, 2.5, 8.0, 6.0
    delta = 1 / ifm.lorentz_gamma(p1) - 1 / ifm.lorentz_gamma(p2)
    pops = ifm.gaussian_populations(50, dn)
    T = ifm.revival_period(omega0, delta)
    dt = ifm.suggested_time_step(dn, omega0, delta)
    t = np.arange(0, 4.2 * T, dt)
    V = ifm.visibility(t, pops, omega0, delta)
    env = ifm.gaussian_envelope(t, dn, omega0, delta)
    fig, ax = plt.subplots(figsize=(7.5, 3.3))
    ax.plot(t, V, "C0-", lw=1.2, label=r"$V(t)=|\sum_n p_n e^{-in\omega_0\Delta t}|$")
    ax.plot(t, env, "C3--", lw=1.2, label="Gaussian envelope (valid for $t\\ll T_{rev}$)")
    for k in range(1, 5):
        ax.axvline(k * T, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel("coordinate time $t$ (units with $m=c=1$)")
    ax.set_ylabel("visibility")
    ax.set_title(rf"$p_1={p1},\ p_2={p2},\ \Delta n={dn:g}$: revivals every $T_{{rev}}=2\pi/(\omega_0|\Delta|)={T:.3f}$")
    ax.set_ylim(-0.03, 1.08)
    ax.legend(fontsize=8, loc="upper center")
    fig.tight_layout()
    save(fig, "fig3_interferometer")


# ----------------------------------------------------------------------------
# Fig. 4  Circuit-QED estimates
# ----------------------------------------------------------------------------
def fig_cqed():
    nbar = np.logspace(1, 4, 100)
    rows = est.circuit_qed_table(nbar_values=nbar)
    g_chi = np.array([r["gamma_chi"] for r in rows])
    g_q = np.array([r["gamma_qubit"] for r in rows])
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.loglog(nbar, g_q, "C3-", lw=2, label=r"$\Delta E=\hbar\omega_s$ (qubit splitting, 5 GHz)")
    ax.loglog(nbar, g_chi, "C0-", lw=2, label=r"$\Delta E=\hbar\chi$ (manuscript's choice, 1 MHz)")
    ax.axhspan(1e3, 1e4, color="0.8", alpha=0.6, label=r"typical transmon $1/T_2$ (0.1–1 ms)")
    ax.set_xlabel(r"mean cavity photon number $\bar n$")
    ax.set_ylabel(r"model rate $\Gamma=\Delta E^2/(2\hbar\Delta E_C)$  [s$^{-1}$]")
    ax.axvspan(250, 1.2e4, color="C3", alpha=0.08, label=r"$\bar n>n_{\rm crit}\approx250$ (non-dispersive)")
    ax.set_title(r"Circuit-QED model rate, $\Delta E_C=\hbar\omega_c\sqrt{\bar n}$, $\omega_c/2\pi=6$ GHz")
    ax.legend(fontsize=7.5, loc="center right")
    fig.tight_layout()
    save(fig, "fig4_cqed")


# ----------------------------------------------------------------------------
# Fig. 5  Finite-clock scaling: the N^-4 law for three-level spectra
# ----------------------------------------------------------------------------
def fig_clock_scaling():
    phi = (1 + 5**0.5) / 2
    specs = {
        "golden  (0,1,phi)": cd.Spectrum([0, 1, phi]),
        "silver  (0,1,1+sqrt2)": cd.Spectrum([0, 1, 1 + 2**0.5]),
        "sqrt2   (0,1,sqrt2)": cd.Spectrum([0, 1, 2**0.5]),
        "pi-1    (0,1,pi-1)": cd.Spectrum([0, 1, np.pi - 1]),
        "e-1     (0,1,e-1)": cd.Spectrum([0, 1, np.e - 1]),
    }
    nmax, pmax = 200_000, 2000
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for name, sp in specs.items():
        N, V = cd.scaling_curve(sp, nmax)
        ax[0].loglog(N, V, lw=1.2, label=name)
    guide = np.array([30.0, nmax])
    v_ref = cd.three_level_optimum(specs["golden  (0,1,phi)"], 10)[0]
    ax[0].loglog(guide, v_ref * 10**4 * guide**-4, "k--", lw=1, label=r"$\propto N^{-4}$")
    ax[0].set(xlabel=r"clock levels $N$", ylabel=r"$V^\star_N$",
              title="optimal mismatch variance (staircase)")
    ax[0].legend(fontsize=7, ncol=2)

    r_phi = phi
    v_e = 2 * (r_phi**2 - r_phi + 1) / 9
    asymptote = r_phi**4 / (135 * v_e)
    for name, sp in specs.items():
        r = sp.energies[2] / sp.energies[1]
        conv = [(pp, qq) for pp, qq in cd.convergents(r, 40) if 8 <= pp <= pmax]
        pp = np.array([c[0] for c in conv], float)
        vv = np.array([cd.mismatch_variance(sp, (0, qq, pp_)) for pp_, qq in conv])
        ax[1].loglog(pp, vv * pp**4, "o-", ms=3.5, lw=1, label=name)
    ax[1].axhline(asymptote, color="k", ls="--", lw=1,
                  label=r"golden constant $\varphi^4/(135V_E)$")
    ax[1].set(xlabel=r"largest label $p$ (convergent labels, $p\leq2000$)",
              ylabel=r"$p^4\,V^\star_p$", title=r"flat = $N^{-4}$ law")
    ax[1].legend(fontsize=7, ncol=2)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "fig5_clock_scaling")


# ----------------------------------------------------------------------------
# Fig. 6  The decay exponent depends on the spectrum, not the dimension
# ----------------------------------------------------------------------------
def fig_clock_dimension():
    phi = (1 + 5**0.5) / 2
    cases = {
        "d=3 (0,1,phi)": cd.Spectrum([0, 1, phi]),
        "d=4 phi powers": cd.Spectrum([0, 1, phi, phi**2]),
        "d=5 phi powers": cd.Spectrum([0, 1, phi, phi**2, phi**3]),
        "d=4 (1,sqrt2,sqrt3)": cd.Spectrum([0, 1, 2**0.5, 3**0.5]),
        "d=3 (0,1,sqrt2)": cd.Spectrum([0, 1, 2**0.5]),
    }
    M = 200_000
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    curves = {}
    for name, sp in cases.items():
        N, V = cd.proportional_curve_stable(sp, M)
        curves[name] = (N, V)
        good = V > 0
        ax.loglog(N[good], V[good], lw=1.3, label=name)
    g = np.array([30.0, M])
    ax.loglog(g, 0.114 * 30.0**4 * g**-4, "k--", lw=1, label=r"$\propto N^{-4}$")
    N0, V0 = curves["d=4 (1,sqrt2,sqrt3)"]
    v0 = V0[np.searchsorted(N0, 30.0)]
    ax.loglog(g, v0 * (30.0 / g) ** 3.1, "k:", lw=1, label=r"$\propto N^{-3.1}$")
    ax.set(xlabel=r"clock levels $N$", ylabel=r"$V^\star_N$ (proportional family)",
           title="the exponent depends on the spectrum, not the dimension")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    save(fig, "fig6_clock_dimension")


# ----------------------------------------------------------------------------
# Fig. 2  The clock face has no date: covering, relabeling, and a seam
# ----------------------------------------------------------------------------
def fig_clock_cover():
    """Daylight saving time as a literal illustration of winding marginalization.

    The reading map t -> t mod 24 h is a covering of R onto the circle. A global
    shift is a rotation of the face (physically inert); a shift applied only
    after a seam is not -- it cuts the cover and deletes a label segment.
    Everything asserted in the caption is computed here.
    """
    t = np.arange(0.0, 96.0, 1.0 / 60.0)
    reading = np.mod(t, 24.0)

    # (a) one reading, many instants
    hits = np.flatnonzero(np.isclose(reading, 2.0, atol=1e-9))
    n_hits = hits.size
    assert n_hits == 4, n_hits                      # 2 AM occurs once per day

    # (b) three maps
    rotation = np.mod(t + 5.0, 24.0)                # global relabeling
    seam_t = 50.0                                   # clocks go forward at 50 h
    dst = np.where(t < seam_t, np.mod(t, 24.0), np.mod(t + 1.0, 24.0))
    step = np.diff(dst)
    wraps = np.abs(step + 24.0) < 1e-6              # benign: the covering closes
    cuts = step > 0.5                               # the seam itself (a forward jump)
    assert cuts.sum() == 1, cuts.sum()
    assert np.all(np.abs(step[wraps] + 24.0) < 1e-6)
    assert step[cuts][0] > 0.9                      # it skips a whole hour
    seam_day = (t >= 48.0) & (t < 72.0)
    labels_seam_day = {int(np.floor(x)) % 24 for x in dst[seam_day]}
    labels_normal = {int(np.floor(x)) % 24 for x in reading[seam_day]}
    assert len(labels_normal) == 24
    assert len(labels_seam_day) == 23 and 2 not in labels_seam_day

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 3.9))

    ax[0].plot(t, reading, "C0-", lw=1.4)
    ax[0].axhline(2.0, color="0.4", ls="--", lw=0.9)
    ax[0].plot(t[hits], np.full(n_hits, 2.0), "o", color="C3", ms=5, zorder=5)
    for k, i in enumerate(hits):
        ax[0].annotate(rf"$m={k}$", (t[i], 2.0), textcoords="offset points",
                       xytext=(2, 7), fontsize=8, color="C3")
    ax[0].set(xlabel="physical time  [h]", ylabel=r"clock reading  $\theta/2\pi$  [h]",
              title="(a) one reading, four instants: the face has no date",
              ylim=(-0.5, 25), xlim=(0, 96))
    ax[0].set_yticks([0, 6, 12, 18, 24])

    # The clean way to compare: the relabeling function itself, g(t) = label(t) - t.
    g_rot = (rotation - reading + 12.0) % 24.0 - 12.0      # global shift: constant
    g_dst = (dst - reading + 12.0) % 24.0 - 12.0           # DST: a step
    ax[1].plot(t, np.zeros_like(t), color="0.72", lw=1.2, label="no relabeling:  $g=0$")
    ax[1].step(t, g_rot, "k--", lw=1.3, where="post", label="global shift:  $g=5$ h")
    ax[1].step(t, g_dst, "C3-", lw=1.6, where="post", label="DST:  a step in $t$")
    ax[1].axvline(seam_t, color="C3", lw=0.8, ls=":")
    ax[1].annotate("seam", (seam_t, 4.3), textcoords="offset points", xytext=(4, 0),
                   fontsize=8, color="C3")
    ax[1].annotate("a constant offset is a rotation of the\nface: every label is still "
                   "displayed\nsomewhere (inert, no $m$ lost)",
                   (2.0, 5.35), fontsize=7.5, color="k")
    ax[1].annotate("a step cuts the cover: the seam hour\nis skipped, so 23 labels are "
                   "displayed\nthat day, not 24",
                   (seam_t + 0.8, 0.55), textcoords="offset points", xytext=(6, 22),
                   fontsize=7.5, color="C3",
                   arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
    ax[1].set(xlabel="physical time  [h]", ylabel=r"relabeling  $g(t)=\theta/2\pi-t$  [h]",
              title="(b) a rotation of the face is inert; a cut is not",
              ylim=(-1.2, 7.0), xlim=(0, 96))
    ax[1].set_yticks([0, 2, 4, 6])
    ax[1].legend(fontsize=7.5, loc="lower right")
    ax[1].legend(fontsize=7.5, loc="lower right", ncol=1)
    for a in ax:
        a.grid(alpha=0.25)
    fig.tight_layout()
    save(fig, "fig_cover")


# ----------------------------------------------------------------------------
# Fig. 3  The spiral staircase: why the line and the circle are one object
# ----------------------------------------------------------------------------
def fig_staircase():
    """Three views of one clock.

    A uniform staircase is a helix. Along it the height is the unrolled time
    (the universal cover); from above it is a circle (the clock face, traversed
    again and again, with no way to say which lap you are on); from the side it
    is a sinusoid (the oscillation). A staircase that narrows as it rises is a
    conical helix: from above it is a logarithmic spiral, and the turn rate
    diverges at the axis, so the walker completes infinitely many turns before
    arriving and simply fades out of the projection.
    """
    th = np.linspace(0.0, 6.0 * np.pi, 900)
    h = th / (2.0 * np.pi)                      # height = winding number
    fig = plt.figure(figsize=(13.0, 3.9))

    # ---- (a) the uniform staircase ----
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    ax.plot(np.cos(th), np.sin(th), h, color="C0", lw=1.4)
    walk = np.array([1.0, 2.0]) * np.pi          # the walker
    ax.plot(np.cos(walk), np.sin(walk), walk / (2 * np.pi), "o", color="C3", ms=8)
    ax.text(np.cos(walk[1]), np.sin(walk[1]), walk[1] / (2 * np.pi) + 0.25,
            "walker", color="C3", fontsize=8)
    ax.set(xlabel=r"$x$", ylabel=r"$y$", zlabel="height = unrolled time")
    ax.set_title("(a) the staircase itself", fontsize=10)
    ax.view_init(elev=22, azim=-60)
    ax.tick_params(labelsize=7)

    # ---- (b) the two projections ----
    ax1 = fig.add_subplot(1, 3, 2)
    ax1.plot(th, np.cos(th), color="C0", lw=1.3)
    ax1.set(xlabel=r"clock angle $\theta$", ylabel=r"height $=\theta/2\pi$")
    ax1.plot(th, h, color="0.55", ls="--", lw=1.1)
    ax1.text(0.3, 1.35, "height: keeps counting", fontsize=7.5, color="0.35")
    ax1.text(0.3, -1.45, "side view: a sinusoid", fontsize=7.5, color="C0")
    ax1.set_title("(b) the same object, two projections", fontsize=10)
    ax1.set_ylim(-1.9, 2.1)
    ax1.grid(alpha=0.25)

    # ---- (c) top view: the ambiguity ----
    ax2 = fig.add_subplot(1, 3, 3)
    tt = np.linspace(0, 2 * np.pi, 400)
    ax2.plot(np.cos(tt), np.sin(tt), color="C0", lw=1.3)
    for m in range(4):
        ax2.plot([np.cos(np.pi / 2)], [np.sin(np.pi / 2)], "o", color="C3", ms=7,
                 alpha=1.0 - 0.2 * m)
        ax2.annotate(f"lap {m}", (0.0, 1.0), textcoords="offset points",
                     xytext=(6, 4 - 11 * m), fontsize=7.5, color="C3")
    ax2.set(xlim=(-1.3, 1.5), ylim=(-1.3, 1.3), aspect="equal",
            xlabel=r"$x$", ylabel=r"$y$")
    ax2.set_title("(c) top view: same point, four instants", fontsize=10)
    ax2.annotate("the face cannot\nsay which lap", (0.05, 0.35), fontsize=7.5, color="0.3")
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    save(fig, "fig_staircase")


def fig_staircase_conical():
    """The vanishing point: a narrowing staircase never lets the walker arrive."""
    th = np.linspace(0.0, 26.0 * np.pi, 6000)
    DECAY_PER_TURN = 1.6                        # r shrinks by e^-1.6 per lap
    b = DECAY_PER_TURN / (2.0 * np.pi)          # log-spiral rate
    r = np.exp(-b * th)
    fig = plt.figure(figsize=(9.6, 3.9))

    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.plot(r * np.cos(th), r * np.sin(th), th / (2 * np.pi), color="C3", lw=1.1)
    z = np.linspace(0.0, 13.0, 60)
    rr = np.exp(-b * 2 * np.pi * z)
    for k in range(0, 60, 7):                   # a walker that fades as it climbs
        ax.plot([rr[k] * np.cos(2 * np.pi * z[k])], [rr[k] * np.sin(2 * np.pi * z[k])],
                [z[k]], "o", color="C0", ms=6, alpha=1.0 - z[k] / 15.0)
    ax.set(xlabel=r"$x$", ylabel=r"$y$", zlabel="height")
    ax.set_title("(a) a narrowing staircase", fontsize=10)
    ax.view_init(elev=16, azim=-70)
    ax.tick_params(labelsize=7)

    # (b) the top view is a log spiral, and the turns diverge
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(r * np.cos(th), r * np.sin(th), color="C3", lw=1.2)
    ax2.set(xlim=(-1.05, 1.05), ylim=(-1.05, 1.05), aspect="equal",
            xlabel=r"$x$", ylabel=r"$y$")
    ax2.set_title("(b) top view: a logarithmic spiral", fontsize=10)
    ax2.annotate("self-similar: every lap is the previous\none scaled by "
                 r"$e^{-2\pi b}$", (0.0, 0.62), fontsize=7.5, color="0.3")
    ax2.annotate("the axis: angle undefined,\ntick rate divergent",
                 (0.08, -0.72), fontsize=7.5, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
    ax2.grid(alpha=0.25)

    # Assert the caption's claims. The turns completed by the time the radius has
    # shrunk to r is -ln(r/r0)/(2 pi b); check that against a numerical integration
    # of d theta = -dr/(b r) along the sampled curve rather than against itself.
    turns_from_integral = np.concatenate(
        [[0.0], np.cumsum(-np.diff(np.log(r)) / b)]
    ) / (2.0 * np.pi)
    turns_analytic = -(np.log(r) - np.log(r[0])) / (2.0 * np.pi * b)
    assert np.allclose(turns_from_integral, turns_analytic, rtol=1e-10, atol=1e-12)
    assert np.allclose(turns_analytic, th / (2 * np.pi), rtol=1e-12)
    assert np.exp(-2 * np.pi * b) == pytest.approx(np.exp(-DECAY_PER_TURN), rel=1e-12)
    # the turns to reach the axis diverge: twice the decades of radius, twice the turns
    r_half = np.exp(-b * th[len(th) // 2])
    assert -np.log(r_half**2) / (b) == pytest.approx(2 * (-np.log(r_half) / b), rel=1e-12)
    fig.tight_layout()
    save(fig, "fig_staircase_conical")


if __name__ == "__main__":
    fig_staircase()
    fig_staircase_conical()
    fig_clock_cover()
    fig_conditioning()
    fig_winding()
    fig_interferometer()
    fig_cqed()
    fig_clock_scaling()
    fig_clock_dimension()
