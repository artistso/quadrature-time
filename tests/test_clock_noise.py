"""Tests for the derived clock-diffusion time T_C (v14).

The point of these tests is that T_C is no longer a free parameter: it follows from
the clock cavity's linewidth kappa and photon number nbar, and the resulting
dephasing lands where experiments live.
"""
import numpy as np
import pytest

from qtime.clock_noise import (
    CAVITY_6GHZ,
    PHOTON_LOSS_1MHZ,
    QUBIT_5GHZ,
    circuit_qed_noise_table,
    coherence_time,
    dephasing_from_T_C,
    dephasing_rate,
    diffusion_to_clock_time,
    schawlow_townes_diffusion,
)


def test_schawlow_townes_scales_inversely_with_photon_number():
    d1 = schawlow_townes_diffusion(PHOTON_LOSS_1MHZ, 1.0)
    for n in (10.0, 100.0, 1e4):
        assert schawlow_townes_diffusion(PHOTON_LOSS_1MHZ, n) == pytest.approx(d1 / n)


def test_schawlow_townes_is_linear_in_the_linewidth():
    d1 = schawlow_townes_diffusion(PHOTON_LOSS_1MHZ, 1e3)
    d2 = schawlow_townes_diffusion(3 * PHOTON_LOSS_1MHZ, 1e3)
    assert d2 == pytest.approx(3 * d1)


def test_diffusion_matches_the_textbook_linewidth():
    """D = kappa/(2 nbar) is the standard Schawlow-Townes result.

    The laser linewidth is Δν_FWHM = κ/(4π n̄) in Hz, i.e. HWHM = κ/(4 n̄) in angular
    frequency; a phase-diffusing field has correlation e^{-D t/2}, whose Lorentzian
    has HWHM = D/2. Hence D = κ/(2 n̄). Checked here rather than asserted.
    """
    kappa, nbar = PHOTON_LOSS_1MHZ, 137.0
    fwhm_hz = kappa / (4 * np.pi * nbar)        # Schawlow-Townes, in hertz
    hwhm_angular = 2 * np.pi * fwhm_hz / 2.0    # radians per second
    d_phi = schawlow_townes_diffusion(kappa, nbar)
    assert d_phi / 2.0 == pytest.approx(hwhm_angular)


def test_clock_time_is_diffusion_over_omega_squared():
    d = 1.23e4
    for omega in (1e9, 6e9, 3.7e10):
        assert diffusion_to_clock_time(d, omega) == pytest.approx(d / omega**2)


def test_dephasing_rate_matches_the_master_equation_formula():
    """Gamma = T_C omega_jk^2 / 2 with T_C from the cavity."""
    kappa, nbar, omega, wjk = PHOTON_LOSS_1MHZ, 1e3, CAVITY_6GHZ, QUBIT_5GHZ
    t_c = diffusion_to_clock_time(schawlow_townes_diffusion(kappa, nbar), omega)
    assert dephasing_rate(kappa, nbar, omega, wjk) == pytest.approx(
        dephasing_from_T_C(t_c, wjk)
    )
    # and the closed form quoted in the paper
    assert dephasing_rate(kappa, nbar, omega, wjk) == pytest.approx(
        kappa * wjk**2 / (4 * nbar * omega**2)
    )


def test_rate_scales_inversely_with_the_clock_field():
    """A weaker clock is a worse clock: Gamma ∝ 1/nbar."""
    g1 = dephasing_rate(PHOTON_LOSS_1MHZ, 10.0, CAVITY_6GHZ, QUBIT_5GHZ)
    g2 = dephasing_rate(PHOTON_LOSS_1MHZ, 1000.0, CAVITY_6GHZ, QUBIT_5GHZ)
    assert g1 / g2 == pytest.approx(100.0)


def test_the_prediction_lands_in_the_measured_transmon_window():
    """At nbar ~ 1e3 the predicted T_2 is 0.1-10 ms, where real transmons live."""
    rows = {r["nbar"]: r for r in circuit_qed_noise_table()}
    t2 = rows[1e3]["T_2"]
    assert 1e-4 < t2 < 1e-2
    # and it is far too slow to be seen with a strong field
    assert rows[1e4]["T_2"] > 1e-3
    # and comfortably fast with a single photon
    assert rows[1.0]["T_2"] < 1e-5


def test_gap_identification_changes_the_rate_by_the_square_of_the_ratio():
    """Using chi = 1 MHz instead of omega_s = 5 GHz suppresses Gamma by (chi/omega_s)^2."""
    chi = 2 * np.pi * 1e6
    g_s = dephasing_rate(PHOTON_LOSS_1MHZ, 1e3, CAVITY_6GHZ, QUBIT_5GHZ)
    g_chi = dephasing_rate(PHOTON_LOSS_1MHZ, 1e3, CAVITY_6GHZ, chi)
    assert g_s / g_chi == pytest.approx((QUBIT_5GHZ / chi) ** 2)


def test_input_validation():
    with pytest.raises(ValueError):
        schawlow_townes_diffusion(0.0, 1.0)
    with pytest.raises(ValueError):
        schawlow_townes_diffusion(1.0, 0.0)
    with pytest.raises(ValueError):
        diffusion_to_clock_time(1.0, 0.0)
    with pytest.raises(ValueError):
        diffusion_to_clock_time(-1.0, 1.0)
    assert coherence_time(0.0) == float("inf")
