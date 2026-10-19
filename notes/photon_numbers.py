#!/usr/bin/env python3
"""The numbers a photon visualisation should display, for a given wavelength.

    python notes/photon_numbers.py

Used to audit the readouts in the PHOTON-E screenshots archived under
archive/sketches/photon_app/. Pure bookkeeping: E = hc/lambda, nu = c/lambda,
hbar*omega = E, momentum p = E/c, and the helicity range for a massless spin-1
quantum.
"""
from __future__ import annotations

import numpy as np

H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
HBAR = H_PLANCK / (2 * np.pi)


def photon_readouts(wavelength_nm: float) -> dict[str, float]:
    lam = float(wavelength_nm) * 1e-9
    if lam <= 0.0:
        raise ValueError("wavelength must be positive")
    energy = H_PLANCK * C_LIGHT / lam
    return {
        "wavelength_m": lam,
        "energy_J": energy,
        "frequency_Hz": C_LIGHT / lam,
        "angular_frequency_rad_s": 2 * np.pi * C_LIGHT / lam,
        "momentum_kg_m_s": energy / C_LIGHT,
        "energy_eV": energy / 1.602176634e-19,
    }


def main() -> None:
    print(f"{'lambda (nm)':>12s}{'E = hc/lambda (J)':>20s}{'nu = c/lambda (Hz)':>21s}"
          f"{'omega (rad/s)':>18s}{'E (eV)':>11s}")
    for lam in (280, 370, 507, 720):
        r = photon_readouts(lam)
        print(f"{lam:12d}{r['energy_J']:20.4e}{r['frequency_Hz']:21.4e}"
              f"{r['angular_frequency_rad_s']:18.4e}{r['energy_eV']:11.3f}")

    print("\naudit of the screenshot values:")
    shown = {507: 3.92e-19, 280: 7.09e-19, 720: 2.76e-19, 370: 5.36e-19}
    for lam, value in shown.items():
        exact = photon_readouts(lam)["energy_J"]
        print(f"  {lam:4d} nm: shown {value:.2e} J, hc/lambda = {exact:.4e} J, "
              f"ratio {value/exact:.4f}")

    print("\nhelicity of a single photon: the projection of spin on momentum is +-1 only")
    print("  (massless spin-1 has no zero-helicity state). A fractional readout in")
    print("  [-1, 1] is the expectation <h> = Stokes V/I of a superposition or ensemble,")
    print("  not the result of any single measurement.")


if __name__ == "__main__":
    main()
