"""
Order-of-magnitude estimates used in the review and in the revised manuscript.
SI units throughout (CODATA values from scipy.constants).
"""
from __future__ import annotations

import numpy as np
from scipy import constants as C

hbar, c, G, kB, g_earth = C.hbar, C.c, C.G, C.k, C.g
M_sun = 1.98847e30  # kg


def clock_dephasing_rate(delta_E_sys: float, delta_E_clock: float) -> float:
    """Manuscript's model rate Gamma = (Delta E_sys)^2 / (2 hbar Delta E_C)   [1/s]."""
    return delta_E_sys**2 / (2.0 * hbar * delta_E_clock)


def static_suppression(delta_E_sys: float, delta_E_clock: float) -> float:
    """Exact static readout smearing with sigma_t = hbar/(2 Delta E_C): exp(-(dE)^2/(8 dE_C^2))."""
    return float(np.exp(-(delta_E_sys**2) / (8.0 * delta_E_clock**2)))


def circuit_qed_table(nbar_values=(10, 1e2, 1e3, 1e4), chi_2pi=1.0e6, omega_c_2pi=6.0e9, omega_s_2pi=5.0e9):
    """
    Two possible identifications of the 'system energy gap' in Sec. VI:
      (a) the paper's choice dE = hbar*chi (dispersive shift),
      (b) the physical qubit splitting dE = hbar*omega_s.
    Clock energy spread: coherent state, dE_C = hbar*omega_c*sqrt(nbar).
    """
    chi, wc, ws = 2 * np.pi * chi_2pi, 2 * np.pi * omega_c_2pi, 2 * np.pi * omega_s_2pi
    rows = []
    for nbar in nbar_values:
        dEC = hbar * wc * np.sqrt(nbar)
        g_chi = clock_dephasing_rate(hbar * chi, dEC)
        g_s = clock_dephasing_rate(hbar * ws, dEC)
        rows.append(
            dict(
                nbar=nbar,
                gamma_chi=g_chi,
                T_chi=1 / g_chi,
                gamma_qubit=g_s,
                T_qubit=1 / g_s,
                static_chi=static_suppression(hbar * chi, dEC),
                static_qubit=static_suppression(hbar * ws, dEC),
            )
        )
    return rows


def critical_photon_number(detuning_2pi: float = 1.0e9, chi_2pi: float = 1.0e6) -> float:
    """
    Dispersive-regime limit n_crit = Delta^2/(4 g^2) with chi ~ g^2/Delta (two-level estimate),
    i.e. n_crit ~ Delta/(4 chi). Above ~n_crit the dispersive Hamiltonian hbar*chi*a^dag a|1><1|
    used in Sec. VI is not valid. Transmon multilevel corrections change the prefactor only.
    """
    return detuning_2pi / (4.0 * chi_2pi)


def chi_for_dispersive_sweep(nbar_max: float, detuning_2pi: float = 1.0e9, margin: float = 10.0) -> float:
    """Largest chi/2pi [Hz] keeping nbar_max <= n_crit/margin."""
    return detuning_2pi / (4.0 * margin * nbar_max)


def redshift_difference_lab(height_m: float = 1.0) -> float:
    """Weak-field Delta N_12 ~ g*dx/c^2 for a vertical separation on Earth."""
    return g_earth * height_m / c**2


def optimal_bandwidth(delta_E_sys: float, delta_N: float) -> float:
    """Delta E_C* = |dE| / (2^{1/4} sqrt|dN|)  (heuristic sum-of-rates model)."""
    return abs(delta_E_sys) / (2 ** 0.25 * np.sqrt(abs(delta_N)))


def minimal_total_rate(delta_E_sys: float, delta_N: float) -> float:
    """Gamma_min = 2 sqrt(a b) = 2^{1/4} |dE| sqrt|dN| / hbar   (corrected)."""
    return 2 ** 0.25 * abs(delta_E_sys) * np.sqrt(abs(delta_N)) / hbar


def hawking_temperature(mass_kg: float = M_sun) -> float:
    return hbar * c**3 / (8 * np.pi * G * mass_kg * kB)


def schwarzschild_lapse(r: np.ndarray, r_s: float) -> np.ndarray:
    return np.sqrt(1.0 - r_s / r)


def proper_radial_distance(r: np.ndarray, r_s: float) -> np.ndarray:
    """rho(r) = int_{r_s}^{r} dr'/sqrt(1 - r_s/r') (closed form)."""
    r = np.asarray(r, dtype=float)
    return np.sqrt(r * (r - r_s)) + r_s * np.log((np.sqrt(r) + np.sqrt(r - r_s)) / np.sqrt(r_s))


if __name__ == "__main__":
    print("Circuit-QED (Sec. VI) -- model rate Gamma = dE^2/(2 hbar dE_C), dE_C = hbar wc sqrt(nbar)")
    for row in circuit_qed_table():
        print(
            f"  nbar={row['nbar']:>7.0f}:  dE=hbar*chi -> Gamma={row['gamma_chi']:.3g}/s (T={row['T_chi']*1e3:.3g} ms);"
            f"  dE=hbar*omega_s -> Gamma={row['gamma_qubit']:.3g}/s (T={row['T_qubit']*1e9:.3g} ns);"
            f"  static factor (chi) = 1-{1-row['static_chi']:.2e}"
        )
    nc = critical_photon_number()
    chi_ok = chi_for_dispersive_sweep(1e4)
    print(f"Dispersive limit for Delta/2pi=1 GHz, chi/2pi=1 MHz: n_crit ~ {nc:.0f}")
    print(f"chi/2pi allowing nbar<=1e4 with 10x margin: {chi_ok/1e3:.1f} kHz;"
          f" model rate (dE=hbar*chi) then {circuit_qed_table((10, 1e4), chi_2pi=chi_ok)[0]['gamma_chi']:.2e}"
          f" to {circuit_qed_table((10, 1e4), chi_2pi=chi_ok)[1]['gamma_chi']:.2e} /s")
    dN = redshift_difference_lab(1.0)
    dE = hbar * 2 * np.pi * 5e9
    print(f"Lab redshift difference for 1 m: dN = {dN:.3e}")
    print(f"Optimal bandwidth for a 5 GHz gap: dE_C* = {optimal_bandwidth(dE, dN):.3e} J = {optimal_bandwidth(dE, dN)/C.e:.3e} eV"
          f"  (= {optimal_bandwidth(dE, dN)/dE:.3e} x gap)")
    print(f"Minimal total rate: {minimal_total_rate(dE, dN):.3e} /s")
    print(f"Hawking temperature, 1 solar mass: {hawking_temperature():.3e} K")
