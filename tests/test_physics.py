import numpy as np
import pytest

from qtime import estimates as est
from qtime import interferometer as ifm
from qtime.modular import reduced_states, tomita_delta


# ---------------------------------------------------------------- modular flow
def test_modular_operator_is_rhoC_inverse_tensor_rhoS():
    rng = np.random.default_rng(1)
    d = 3
    alpha = rng.random(d) + 0.2
    alpha /= np.linalg.norm(alpha)
    # random local unitaries so the test is not basis-trivial
    Uc, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    Us, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    psi = Uc @ np.diag(alpha) @ Us.T
    rho_C, rho_S = reduced_states(psi)
    Delta = tomita_delta(psi)
    correct = np.kron(np.linalg.inv(rho_C), rho_S)
    as_in_original = np.kron(rho_C, np.linalg.inv(rho_S))
    assert np.allclose(Delta, correct, atol=1e-9)
    assert not np.allclose(Delta, as_in_original, atol=1e-3)


def test_gibbs_schmidt_spectrum_modular_flow_is_time_evolution():
    E = np.array([0.0, 0.7, 1.9])
    beta = 1.3
    p = np.exp(-beta * E)
    p /= p.sum()
    rho_S = np.diag(p)
    A = np.array([[0, 1, 2], [1, 0, 1j], [2, -1j, 0]], dtype=complex)
    s = 0.37
    mod = np.diag(p ** (1j * s)) @ A @ np.diag(p ** (-1j * s))
    t = -beta * s  # hbar = 1
    heis = np.diag(np.exp(1j * E * t)) @ A @ np.diag(np.exp(-1j * E * t))
    assert np.allclose(mod, heis)
    assert np.allclose(rho_S, np.diag(p))


# ---------------------------------------------------------------- interferometer
def test_visibility_is_exactly_periodic_and_envelope_valid_at_short_times():
    g1, g2 = ifm.lorentz_gamma(0.2), ifm.lorentz_gamma(2.5)
    delta = 1 / g1 - 1 / g2
    pops = ifm.gaussian_populations(50, 6.0)
    T = ifm.revival_period(8.0, delta)
    t = np.linspace(0, 0.05, 50)
    V = ifm.visibility(t, pops, 8.0, delta)
    assert np.allclose(ifm.visibility(t + T, pops, 8.0, delta), V)
    assert np.allclose(V[:10], ifm.gaussian_envelope(t[:10], 6.0, 8.0, delta), atol=2e-3)


def test_original_window_is_dominated_by_full_revivals_not_decay():
    # Original script: t in [0, 15]. The internal spectrum is equally spaced, so the
    # "visibility loss" is reversible dephasing with exact revivals; the Gaussian
    # envelope plotted alongside it is only valid before the first revival.
    g1, g2 = ifm.lorentz_gamma(0.2), ifm.lorentz_gamma(2.5)
    delta = 1 / g1 - 1 / g2
    pops = ifm.gaussian_populations(50, 6.0)
    T = ifm.revival_period(8.0, delta)
    assert 15.0 / T > 11  # more than eleven full revivals inside the plotted window
    V_rev = ifm.visibility(np.array([3 * T]), pops, 8.0, delta)[0]
    env_rev = ifm.gaussian_envelope(np.array([3 * T]), 6.0, 8.0, delta)[0]
    assert V_rev == pytest.approx(1.0) and env_rev < 1e-100


# ---------------------------------------------------------------- estimates
def test_optimal_bandwidth_minimises_sum_and_corrected_minimum():
    dE, dN = 1.0e-24, 1e-3
    a = dE**2 / (2 * est.hbar)
    b = abs(dN) / (np.sqrt(2) * est.hbar)
    x = np.logspace(-30, -15, 20001)
    tot = a / x + b * x
    x_star = est.optimal_bandwidth(dE, dN)
    assert x[np.argmin(tot)] == pytest.approx(x_star, rel=3e-3)
    assert est.minimal_total_rate(dE, dN) == pytest.approx(2 * np.sqrt(a * b), rel=1e-12)
    # original text had |dN| instead of sqrt|dN|
    original = 2 ** 0.25 * dE * abs(dN) / est.hbar
    assert not np.isclose(original, 2 * np.sqrt(a * b), rtol=1e-2)


def test_hawking_temperature_solar_mass():
    assert est.hawking_temperature() == pytest.approx(6.17e-8, rel=2e-3)


def test_near_horizon_lapse_is_linear_in_proper_distance():
    r_s = 1.0
    r = r_s * (1 + np.logspace(-8, -4, 5))
    rho = est.proper_radial_distance(r, r_s)
    N = est.schwarzschild_lapse(r, r_s)
    assert np.allclose(N, rho / (2 * r_s), rtol=1e-3)  # N ~ kappa*rho/c^2 with kappa/c^2 = 1/(2 r_s)
    # hence Gamma ~ 1/N(r_C) scales as rho^-1, i.e. as (r - r_s)^(-1/2) in the coordinate distance
    slope = np.polyfit(np.log(rho), np.log(1 / N), 1)[0]
    assert slope == pytest.approx(-1.0, abs=1e-3)
    slope_coord = np.polyfit(np.log(r - r_s), np.log(1 / N), 1)[0]
    assert slope_coord == pytest.approx(-0.5, abs=1e-3)


def test_circuit_qed_identifications():
    rows = est.circuit_qed_table()
    # paper's choice (dE = hbar*chi): rates between ~5 and ~170 per second
    assert 4 < rows[-1]["gamma_chi"] < 6 and 150 < rows[0]["gamma_chi"] < 180
    # physical qubit splitting: >1e8 /s for all nbar <= 1e4, i.e. T2 < 10 ns -- excluded by experiment
    assert all(r["gamma_qubit"] > 1e8 for r in rows)
    # the exact static effect with sigma_t = hbar/(2 dE_C) is unobservably small
    assert all(1 - r["static_chi"] < 1e-9 for r in rows)


def test_dispersive_limit_restricts_photon_sweep():
    from qtime.estimates import critical_photon_number, chi_for_dispersive_sweep, circuit_qed_table
    nc = critical_photon_number(1e9, 1e6)
    assert abs(nc - 250) < 1e-9          # the paper's parameters: n_crit ~ 250 << 1e4
    chi_ok = chi_for_dispersive_sweep(1e4)  # 2.5 kHz with a 10x margin
    assert abs(chi_ok - 2.5e3) < 1e-9
    g = circuit_qed_table((10, 1e4), chi_2pi=chi_ok)
    assert g[0]["gamma_chi"] < 2e-3 and g[1]["gamma_chi"] < 1e-4   # unmeasurably small rates
