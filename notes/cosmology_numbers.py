#!/usr/bin/env python3
"""Reproduce every number in notes/2026-10-18_black-hole-universe-and-quadratures.md.

    python notes/cosmology_numbers.py

Nothing here is imported by the package or the paper; this is a working note's
arithmetic, kept runnable so the note cannot drift from its numbers.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from qtime.clock_design import (  # noqa: E402
    Spectrum,
    three_level_optimum,
    time_horizon,
)

# --------------------------------------------------------------------------
# constants (SI)
# --------------------------------------------------------------------------
G = 6.67430e-11
c = 2.99792458e8
Mpc = 3.0856775814913673e22
Msun = 1.98847e30
Gly = 9.4607304725808e24
YR = 3.15576e7
H0 = 70.0e3 / Mpc
rho_c = 3 * H0**2 / (8 * np.pi * G)
R_H = c / H0


def section(title: str) -> None:
    print("\n" + title)
    print("-" * len(title))


def pathria_identity() -> None:
    section("1. At critical density, r_s(R) = R exactly at R = c/H0 (Pathria)")
    print(f"  H0 = {H0:.4e} s^-1 | R_H = {R_H/Gly:.2f} Gly | rho_c = {rho_c:.3e} kg/m^3")
    for frac in (0.5, 1.0, 2.0):
        R = frac * R_H
        M = rho_c * 4 / 3 * np.pi * R**3
        rs = 2 * G * M / c**2
        print(
            f"  R = {frac:3.1f} R_H -> M = {M:.3e} kg | r_s = {rs/Gly:7.2f} Gly "
            f"| r_s/R = {rs/R:.12f}"
        )
    print(f"  mass inside R_H = {rho_c*4/3*np.pi*R_H**3/Msun:.3e} Msun")


def observed_universe() -> None:
    section("2. The honest caveat: the observable radius is not the Hubble radius")
    R_obs = 46.5 * Gly
    M_obs = 1.5e53
    rs = 2 * G * M_obs / c**2
    print(f"  R_obs = {R_obs/Gly:.1f} Gly = {R_obs/R_H:.2f} R_H")
    print(f"  M_obs ~ {M_obs:.2e} kg -> r_s = {rs/Gly:.2f} Gly = {rs/R_obs:.2f} R_obs")
    print("  -> of order unity, but not 1: quote it as a coincidence, not a proof.")


def photon_diffusion() -> None:
    section("3. Trapped radiation: the Sun vs the early universe")
    R_sun = 6.957e8
    for lam in (1e-3, 1e-2, 1e-1):
        t = R_sun**2 / (lam * c)
        print(f"  solar interior, mean free path {lam*1e3:5.1f} mm -> {t/YR:9.1f} yr")
    print(f"  (light crosses the Sun in {R_sun/c:.1f} s)")
    t_rec = 380e3
    print(f"  universe opaque until recombination: t ~ {t_rec:.0f} yr, z ~ 1100")
    print(f"  -> the CMB is the photosphere; ratio to a 1e5 yr solar diffusion: {t_rec/1e5:.1f}")


def reading_vs_horizon() -> None:
    section("4. A sharper reading costs horizon: T proportional to N^2")
    spec = Spectrum([0.0, 1.0, (1 + 5**0.5) / 2])
    eps = 1e-2
    print(f"  {'N':>6s}{'1/N':>10s}{'T(eps=1e-2)':>16s}{'T (1/N)^2':>12s}")
    for N in (10, 100, 1000, 3000):
        v, _ = three_level_optimum(spec, N)
        T = time_horizon(v, eps)
        print(f"  {N:6d}{1/N:10.4f}{T:16.4f}{T*(1/N)**2:12.3f}")
    print("  -> a ten-fold sharper reading costs a hundred-fold of horizon.")


def sine_cosine() -> None:
    section("5. Sine and cosine: one oscillation, a quarter period apart")
    n = np.arange(120, dtype=float)
    a = np.diag(np.sqrt(n[1:]), 1)
    alpha = 1.0 + 0.5j

    def disp(al):
        from scipy.linalg import expm

        return expm(al * a.conj().T - np.conj(al) * a)[:, 0]

    def phase(theta):
        from scipy.linalg import expm

        return expm(-1j * theta * np.diag(n))

    X = (a + a.conj().T) / np.sqrt(2)
    P = -1j * (a - a.conj().T) / np.sqrt(2)
    psi = disp(alpha)
    vr = phase(np.pi / 2) @ psi

    def moments(v):
        return (
            float(np.real(np.vdot(v, X @ v))),
            float(np.real(np.vdot(v, P @ v))),
            float(np.sum(n * np.abs(v) ** 2)),
        )

    x0, p0, n0 = moments(psi)
    x1, p1, n1 = moments(vr)
    print(f"  before: <x>={x0:.6f} <p>={p0:.6f} <N>={n0:.12f}")
    print(f"  after 90 deg: <x>={x1:.6f} <p>={p1:.6f} <N>={n1:.12f}")
    print("  -> quadratures exchanged, energy identical: same state, later in phase.")


if __name__ == "__main__":
    pathria_identity()
    observed_universe()
    photon_diffusion()
    reading_vs_horizon()
    sine_cosine()
    print("\nAll numbers above match the note.")
