"""Tests for the finite-clock design module (notebooks/finite_clock, v1-v8).

Each test corresponds to a numbered claim in notebooks/finite_clock/README.md.
"""
import itertools

import numpy as np
import pytest

from qtime.clock_design import (
    Spectrum,
    continued_fraction,
    convergents,
    diophantine_ratio,
    exhaustive_optimum,
    fidelity_lower_bound,
    fit_quality,
    mismatch_variance,
    optimal_spacing,
    predicted_exponent,
    scaling_curve,
    three_level_closed_form,
    three_level_optimum,
    time_horizon,
    universal_bound,
    universal_horizon,
)

PHI = (1.0 + 5.0**0.5) / 2.0
SILVER = 1.0 + 2.0**0.5

SPECTRA = {
    "golden": [0.0, 1.0, PHI],
    "silver": [0.0, 1.0, SILVER],
    "sqrt2": [0.0, 1.0, 2.0**0.5],
    "pi": [0.0, 1.0, np.pi - 1.0],
    "e": [0.0, 1.0, np.e - 1.0],
    "generic": [0.0, 1.0, 1.234567],
}


def uniform_energy_variance(E):
    E = np.asarray(E, float)
    return float(np.mean((E - E.mean()) ** 2))


# --------------------------------------------------------------------------
# claim 1: the three-level closed form is exact (uniform weights)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("r", [PHI, SILVER, np.pi - 1.0, np.e - 1.0, 1.234567])
def test_three_level_closed_form_matches_weighted_least_squares(r):
    spec = Spectrum([0.0, 1.0, r])
    for q in range(1, 25):
        for p in range(q + 1, 35):
            exact = mismatch_variance(spec, (0, q, p))
            closed = three_level_closed_form(r, q, p)
            # the covariance form loses digits when V* is tiny; the closed form is the
            # numerically stable one, so compare absolutely at the round-off floor
            assert closed == pytest.approx(exact, rel=1e-7, abs=1e-18)


def test_closed_form_needs_uniform_weights():
    """With non-uniform weights the closed form does not apply (documented)."""
    r = PHI
    spec = Spectrum([0.0, 1.0, r], weights=[0.4, 0.35, 0.25])
    exact = mismatch_variance(spec, (0, 3, 5))
    closed = three_level_closed_form(r, 3, 5)
    assert abs(exact - closed) > 1e-8


# --------------------------------------------------------------------------
# claim 2: the optimal spacing is the weighted regression slope
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["golden", "silver", "generic"])
def test_optimal_spacing_minimises_variance(name):
    spec = Spectrum(SPECTRA[name])
    for labels in [(0, 3, 5), (0, 7, 12), (0, 1, 2)]:
        d_star = optimal_spacing(spec, labels)
        v_star = mismatch_variance(spec, labels)
        for delta in np.linspace(0.2 * d_star, 1.8 * d_star, 37):
            gapped = spec.shifted() - np.asarray(labels) * delta
            assert float(spec.weights @ (gapped - spec.weights @ gapped) ** 2) >= v_star - 1e-15


def test_variance_is_offset_invariant():
    spec = Spectrum([0.0, 1.0, PHI])
    shifted = Spectrum(np.asarray([0.0, 1.0, PHI]) + 17.0)
    assert mismatch_variance(spec, (0, 3, 5)) == pytest.approx(mismatch_variance(shifted, (0, 3, 5)))
    # a common shift of the labels is the same as shifting the energies
    assert mismatch_variance(spec, (0, 3, 5)) == pytest.approx(mismatch_variance(spec, (4, 7, 9)))


def test_fit_quality_is_one_minus_normalised_variance():
    spec = Spectrum(SPECTRA["golden"])
    v = mismatch_variance(spec, (0, 8, 13))
    assert fit_quality(spec, (0, 8, 13)) == pytest.approx(1.0 - v / spec.variance)


# --------------------------------------------------------------------------
# claim 3: the universal bound V*_N <= R^2 / (4(N-1)^2)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", list(SPECTRA))
def test_universal_bound_holds_for_rounding(name):
    spec = Spectrum(SPECTRA[name])
    R = spec.energy_range
    for N in [4, 5, 6, 8, 12, 24, 48, 96, 128]:
        delta = R / (N - 1)
        labels = np.rint(spec.shifted() / delta).astype(int)
        assert labels.max() <= N - 1
        v = mismatch_variance(spec, labels)
        assert v <= universal_bound(R, N) * (1.0 + 1e-12)


@pytest.mark.parametrize("name", list(SPECTRA))
def test_optimum_never_exceeds_universal_bound(name):
    spec = Spectrum(SPECTRA[name])
    R = spec.energy_range
    for N in [6, 10, 20, 40, 80]:
        v, _ = three_level_optimum(spec, N)
        assert v <= universal_bound(R, N) * (1.0 + 1e-12)


def test_universal_bound_rejects_a_one_level_clock():
    with pytest.raises(ValueError):
        universal_bound(1.0, 1)


def test_universal_horizon_matches_the_bound():
    R, N, eps = 2.0, 50, 0.01
    assert universal_horizon(R, N, eps) == pytest.approx(
        time_horizon(universal_bound(R, N), eps), rel=1e-12
    )


# --------------------------------------------------------------------------
# claim 4/5: the O(N) optimiser agrees with exhaustive enumeration
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["golden", "silver", "pi", "e", "generic"])
@pytest.mark.parametrize("N", [4, 6, 10, 16, 24, 40, 60])
def test_three_level_optimum_matches_exhaustive(name, N):
    spec = Spectrum(SPECTRA[name])
    v_ex, lab_ex = exhaustive_optimum(spec, N)
    v_fast, lab_fast = three_level_optimum(spec, N)
    assert v_fast == pytest.approx(v_ex, rel=1e-12, abs=1e-18)
    assert len(lab_fast) == 3


@pytest.mark.parametrize("name", ["golden", "silver", "pi"])
def test_scaling_curve_matches_the_optimiser(name):
    spec = Spectrum(SPECTRA[name])
    N, V = scaling_curve(spec, 200)
    for n in [4, 10, 25, 60, 120, 199]:
        v = three_level_optimum(spec, n)[0]
        assert V[np.searchsorted(N, n)] == pytest.approx(v, rel=1e-8, abs=2e-15)


def test_scaling_curve_is_monotone():
    spec = Spectrum(SPECTRA["golden"])
    N, V = scaling_curve(spec, 5000)
    assert np.all(np.diff(V) <= 1e-30)


# --------------------------------------------------------------------------
# claim 6: the N^-4 law and its exact constant
# --------------------------------------------------------------------------
def test_golden_ratio_decays_as_n_minus_four_with_exact_constant():
    """At Fibonacci labels, V*_N * N^4 -> phi^4 / (135 V_E)."""
    spec = Spectrum([0.0, 1.0, PHI])
    r = PHI
    v_e = 2.0 * (r * r - r + 1.0) / 9.0
    predicted = r**4 / (135.0 * v_e)
    assert predicted == pytest.approx(0.1142350327708281, rel=1e-12)
    for p, q in convergents(PHI, 20):
        if q < 3:
            continue
        v = three_level_closed_form(PHI, q, p)
        # the O(1/q) correction to the denominator |B|^2 = q^2 P(1 + O(1/q))
        # is a 2% effect at p = 5 and below the 1e-3 level from p = 21 on
        tol = 3e-2 if p < 21 else 3e-3
        assert v * p**4 == pytest.approx(predicted, rel=tol)


def test_measured_exponent_is_close_to_minus_four():
    spec = Spectrum([0.0, 1.0, PHI])
    # sample at Fibonacci labels, where the optimum actually re-equilibrates
    fibs = [f for _, f in convergents(PHI, 24)]
    fibs = [n for n in fibs if 10**3 <= n <= 10**6]
    assert len(fibs) >= 4
    exponent = predicted_exponent(spec, fibs)
    assert exponent == pytest.approx(-4.0, abs=0.05)


def test_commensurate_spectrum_is_exact():
    """A rational gap ratio is reproduced exactly once N reaches the denominator."""
    spec = Spectrum([0.0, 1.0, 2.5])
    v, labels = three_level_optimum(spec, 10)
    assert v == pytest.approx(0.0, abs=1e-30)
    assert labels == (0, 2, 5)


# --------------------------------------------------------------------------
# claim 7: continued fractions and the Diophantine constant
# --------------------------------------------------------------------------
def test_convergents_of_the_golden_ratio_are_fibonacci():
    assert convergents(PHI, 8) == [(1, 1), (2, 1), (3, 2), (5, 3), (8, 5), (13, 8), (21, 13), (34, 21)]


def test_diophantine_constants_hit_the_classical_values():
    """Hurwitz: 1/sqrt(5) is the largest possible kappa; it is the golden ratio's."""
    # the estimate converges from below like O(1/q), so 1% is the honest tolerance
    assert diophantine_ratio(PHI) == pytest.approx(1.0 / 5.0**0.5, rel=1.5e-2)
    assert diophantine_ratio(SILVER) == pytest.approx(1.0 / (2.0 * 2.0**0.5), rel=1.5e-2)
    assert diophantine_ratio(PHI) > diophantine_ratio(SILVER)


def test_continued_fraction_of_a_rational_terminates():
    assert continued_fraction(2.5, 10) == [2, 2]
    assert continued_fraction(PHI, 5) == [1, 1, 1, 1, 1]


# --------------------------------------------------------------------------
# claim 8: operational wrappers
# --------------------------------------------------------------------------
def test_fidelity_bound_and_horizon_are_consistent():
    var, eps = 1e-6, 0.05
    T = time_horizon(var, eps)
    assert fidelity_lower_bound(T, var) == pytest.approx(1.0 - eps)
    assert fidelity_lower_bound(0.5 * T, var) > 1.0 - eps
    assert time_horizon(0.0, eps) == float("inf")


def test_validation_errors():
    with pytest.raises(ValueError):
        Spectrum([0.0, 1.0, 0.5])            # not increasing
    with pytest.raises(ValueError):
        Spectrum([0.0, 1.0], weights=[1.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        Spectrum([0.0, 1.0], weights=[0.0, 0.0])
    with pytest.raises(ValueError):
        mismatch_variance(Spectrum([0.0, 1.0, 2.0]), (0, 0, 1))   # repeated label
    with pytest.raises(ValueError):
        three_level_optimum(Spectrum([0.0, 1.0, 2.0, 3.0]), 8)    # wrong size
    with pytest.raises(ValueError):
        three_level_closed_form(PHI, 5, 5)                        # p must exceed q


# ==========================================================================
# v9: exact arithmetic and the general (d >= 4) optimiser
# ==========================================================================
from fractions import Fraction

from qtime.clock_design import (
    decay_exponent,
    mismatch_variance_stable,
    exact_variance,
    fitted_design,
    gap_ratio_fraction,
    legalize_labels,
    nearest_labels,
    proportional_curve,
    proportional_curve_stable,
    search_optimum,
    self_consistent,
    three_level_closed_form_exact,
)


def exhaustive_general(spec, N):
    """Exact optimum, scored with the same stable residual form the optimiser reports."""
    d = spec.size
    best, arg = spec.variance, tuple(range(d))
    for tail in itertools.combinations(range(1, N), d - 1):
        lab = (0,) + tail
        v = mismatch_variance_stable(spec, lab)
        if v < best:
            best, arg = v, lab
    return best, arg


# ---- exact arithmetic ----------------------------------------------------
def test_exact_variance_matches_the_float_formula():
    spec = Spectrum([0.0, 1.0, PHI])
    E = [Fraction(0), Fraction(1), gap_ratio_fraction(spec)]
    w = [Fraction(1, 3)] * 3
    for labels in [(0, 3, 5), (0, 8, 13)]:
        exact = exact_variance(E, w, labels)
        assert float(exact) == pytest.approx(mismatch_variance(spec, labels), rel=1e-9)


def test_exact_closed_form_matches_the_float_closed_form():
    spec = Spectrum([0.0, 1.0, PHI])
    r = gap_ratio_fraction(spec)
    for q, p in [(3, 5), (13, 21), (233, 377)]:
        assert float(three_level_closed_form_exact(r, q, p)) == pytest.approx(
            three_level_closed_form(PHI, q, p), rel=1e-12
        )


@pytest.mark.parametrize("labels", [(0, 2, 5), (0, 11, 27)])
def test_commensurate_spectrum_has_exactly_zero_variance(labels):
    """A rational spectrum is reproduced exactly -- an exact statement, not a float one."""
    spec = Spectrum([0.0, 1.0, 2.5])
    r = gap_ratio_fraction(spec)
    assert exact_variance([Fraction(0), Fraction(1), r], [Fraction(1, 3)] * 3, (0, 2, 5)) == 0
    # and V* itself is a rational number of denominator dividing the labels' geometry
    v = exact_variance([Fraction(0), Fraction(1), r], [Fraction(1, 3)] * 3, labels)
    assert isinstance(v, Fraction)


def test_exact_arithmetic_survives_where_the_covariance_form_does_not():
    """The covariance form cancels two numbers of order V_E; exact arithmetic cannot.

    At p = 4181 the residual is ~1e-4 while V_E is of order 1, so
    ``mismatch_variance`` has lost the answer entirely; the exact value still
    pins the asymptotic constant, up to the truncation of r itself.
    """
    spec = Spectrum([0.0, 1.0, PHI])
    r = gap_ratio_fraction(spec)
    v_e = 2 * (PHI**2 - PHI + 1) / 9
    predicted = PHI**4 / (135 * v_e)
    q, p = 2584, 4181
    exact = float(three_level_closed_form_exact(r, q, p)) * p**4
    assert exact == pytest.approx(predicted, rel=1e-6)

    floaty = mismatch_variance(spec, (0, q, p)) * p**4
    assert abs(floaty / predicted - 1.0) > 1e-3          # hopeless at this size

    # and the truncation of r, not the arithmetic, is what limits us next
    q2, p2 = 46368, 75025
    exact2 = float(three_level_closed_form_exact(r, q2, p2)) * p2**4
    assert exact2 == pytest.approx(predicted, rel=1e-4)


# ---- the general optimiser ----------------------------------------------
@pytest.mark.parametrize("N", [10, 20, 40, 80, 120])
def test_search_optimum_matches_the_exact_three_level_optimiser(N):
    spec = Spectrum([0.0, 1.0, PHI])
    op = search_optimum(spec, N)
    assert op.variance == pytest.approx(three_level_optimum(spec, N)[0], rel=1e-8, abs=1e-20)


@pytest.mark.parametrize("N", [10, 24, 60])
def test_search_optimum_handles_weighted_three_level_spectra(N):
    spec = Spectrum([0.0, 1.0, PHI], weights=[0.4, 0.35, 0.25])
    ex, _ = exhaustive_general(spec, N)
    assert search_optimum(spec, N).variance == pytest.approx(ex, rel=1e-8, abs=1e-20)


@pytest.mark.parametrize("N", [8, 12, 16, 20, 24, 30])
def test_search_optimum_matches_exhaustive_for_four_levels(N):
    spec = Spectrum([0.0, 1.0, PHI, PHI**2])
    ex, _ = exhaustive_general(spec, N)
    assert search_optimum(spec, N).variance == pytest.approx(ex, rel=1e-8, abs=1e-20)


def test_search_optimum_matches_exhaustive_for_other_four_level_spectra():
    for E in ([0.0, 1.0, 2.0**0.5, 3.0**0.5], [0.0, 1.0, 2.0, 3.0], [0.0, 0.73, 1.41, 2.26]):
        spec = Spectrum(E)
        for N in (12, 20, 30):
            ex, _ = exhaustive_general(spec, N)
            assert search_optimum(spec, N).variance == pytest.approx(ex, rel=1e-8, abs=1e-20)


def test_search_optimum_respects_the_universal_bound():
    for E in ([0.0, 1.0, PHI, PHI**2], [0.0, 1.0, 2.0**0.5, 3.0**0.5], [0.0, 0.73, 1.41, 2.26, 3.3]):
        spec = Spectrum(E)
        for N in (10, 40, 200, 2000):
            op = search_optimum(spec, N)
            assert op.variance <= universal_bound(spec.energy_range, N) * (1 + 1e-12)
            assert len(op.labels) == spec.size
            assert list(op.labels) == sorted(set(op.labels))
            assert 0 <= op.labels[0] and op.labels[-1] <= N - 1


def test_proportional_family_finds_the_four_level_phi_optimum():
    spec = Spectrum([0.0, 1.0, PHI, PHI**2])
    op = search_optimum(spec, 30)
    assert op.labels == (0, 8, 13, 21)          # labels proportional to (1, phi, phi^2)
    assert op.self_consistent


def test_legalize_and_nearest_labels_are_legal():
    spec = Spectrum([0.0, 1.0, PHI])
    for N in (4, 7, 33):
        for delta in (0.01, 0.3, 5.0):
            for offset in (-1.0, 0.0, 0.7):
                lab = nearest_labels(spec, delta, offset, N)
                assert len(lab) == 3
                assert list(lab) == sorted(set(lab))
                assert 0 <= lab[0] and lab[-1] <= N - 1
    with pytest.raises(ValueError):
        legalize_labels((0, 1, 2), 2)


def test_fitted_design_minimises_the_joint_least_squares():
    spec = Spectrum([0.0, 1.0, PHI])
    labels = (0, 8, 13)
    delta, offset = fitted_design(spec, labels)
    def loss(d, a):
        r = spec.energies - a - np.asarray(labels) * d
        return float(spec.weights @ r**2)
    best = loss(delta, offset)
    for dd in np.linspace(0.5 * delta, 1.5 * delta, 11):
        for aa in np.linspace(-0.3, 0.3, 11):
            assert loss(dd, aa) >= best - 1e-15


def test_self_consistency_flag_is_honest():
    spec = Spectrum([0.0, 1.0, PHI])
    best = search_optimum(spec, 40)
    assert self_consistent(spec, best.labels, 40) == best.self_consistent


# ---- what the exponent really is ---------------------------------------
def test_exponent_is_four_for_phi_power_spectra_in_every_dimension():
    for E, M in ([[0.0, 1.0, PHI], 4000],
                 [[0.0, 1.0, PHI, PHI**2], 4000],
                 [[0.0, 1.0, PHI, PHI**2, PHI**3], 4000]):
        N, V = proportional_curve(Spectrum(E), M)
        assert decay_exponent(N, V) == pytest.approx(-4.0, abs=0.12)


def test_exponent_is_worse_for_independent_irrational_ratios():
    """(1, sqrt2, sqrt3) needs simultaneous approximation in two dimensions."""
    N, V = proportional_curve(Spectrum([0.0, 1.0, 2.0**0.5, 3.0**0.5]), 4000)
    exponent = decay_exponent(N, V)
    assert -4.0 < exponent < -2.5


def test_commensurate_spectrum_reaches_exactly_zero():
    N, V = proportional_curve(Spectrum([0.0, 1.0, 2.5]), 20)
    assert V[-1] == 0.0


# ---- the stable proportional curve ------------------------------------
def test_stable_proportional_curve_agrees_with_exact_arithmetic():
    for E in ([0.0, 1.0, PHI], [0.0, 1.0, PHI, PHI**2]):
        spec = Spectrum(E)
        N_exact, V_exact = proportional_curve(spec, 1500)
        N_fast, V_fast = proportional_curve_stable(spec, 1500)
        idx = np.searchsorted(N_fast, N_exact)
        ok = idx < N_fast.size
        rel = np.abs(V_fast[idx[ok]] - V_exact[ok]) / np.maximum(V_exact[ok], 1e-300)
        assert rel.max() < 1e-6


def test_stable_proportional_curve_is_monotone_and_positive():
    N, V = proportional_curve_stable(Spectrum([0.0, 1.0, PHI]), 20000)
    assert np.all(np.diff(V) <= 0)
    assert np.all(V >= 0)
    assert N[0] >= 3


def test_proportional_stable_supports_six_levels():
    spec = Spectrum([0.0, 1.0, PHI, PHI**2, PHI**3, PHI**4])
    N, V = proportional_curve_stable(spec, 5000)
    assert decay_exponent(N, V) == pytest.approx(-4.0, abs=0.15)


# ==========================================================================
# v10: commensurability in an expanding universe
# ==========================================================================
from qtime.clock_design import (
    closed_mode_frequency_s3,
    commensurability_ratio,
    drift_over_efolds,
    drift_rate,
    first_unprotected_efolds,
    scale_factor,
    winding_envelope,
)


def test_scale_factor_is_normalised_and_increasing():
    assert scale_factor(0.0) == pytest.approx(1.0)
    ts = np.linspace(0.0, 3.0, 20)
    a = np.array([scale_factor(t) for t in ts])
    assert np.all(np.diff(a) > 0)
    # de Sitter by default: a = e^{H0 t}
    assert scale_factor(1.0) == pytest.approx(np.e)


def test_scale_factor_matter_plus_lambda_matches_the_ode():
    H0, om, ol = 1.0, 0.3, 0.7
    a = lambda t: scale_factor(t, H0=H0, omega_m=om, omega_l=ol)
    h = 1e-6
    t = 0.7
    adot = (a(t + h) - a(t - h)) / (2 * h)
    # Friedmann: (adot/a)^2 = H0^2 (om a^-3 + ol)
    assert (adot / a(t)) ** 2 == pytest.approx(H0**2 * (om / a(t) ** 3 + ol), rel=1e-6)
    with pytest.raises(ValueError):
        scale_factor(1.0, omega_m=1.0, omega_l=0.0)


def test_closed_mode_frequency_scales_inversely_with_a():
    base = closed_mode_frequency_s3(3, 1.0)
    for a in (0.5, 1.0, 2.0, 10.0):
        assert closed_mode_frequency_s3(3, a) == pytest.approx(base / a)
    assert closed_mode_frequency_s3(2, 1.0) == pytest.approx(np.sqrt(3.0))
    with pytest.raises(ValueError):
        closed_mode_frequency_s3(1, 1.0)      # the homogeneous mode has zero frequency
    with pytest.raises(ValueError):
        closed_mode_frequency_s3(2, 0.0)


def test_winding_envelope_matches_the_manuscript_closed_form():
    nu = np.linspace(-2, 3, 5001)
    for M in (1, 4, 25):
        got = winding_envelope(nu, M)
        sin_nu = np.sin(np.pi * nu)
        safe = np.abs(sin_nu) >= 1e-12
        ref = np.ones_like(nu)
        ref[safe] = np.abs(np.sin(np.pi * M * nu[safe]) / (M * sin_nu[safe]))
        assert np.allclose(got, ref)
        assert got.max() <= 1.0 + 1e-12
        assert np.allclose(winding_envelope(np.round(nu).astype(float), M), 1.0)


def test_commensurability_ratio_and_drift_are_consistent():
    a, k, dE = 2.0, 3, 5.0
    omega = closed_mode_frequency_s3(k, a)
    assert commensurability_ratio(dE, omega) == pytest.approx(dE / omega)
    # d ln omega / dt = -H  =>  dnu/dt = +nu H
    nu = commensurability_ratio(dE, omega)
    assert drift_rate(dE, k, a, hdot_over_h=0.7) == pytest.approx(0.7 * nu)
    with pytest.raises(ValueError):
        commensurability_ratio(1.0, 0.0)


def test_expansion_breaks_commensurability_in_order_one_efold():
    """omega ~ 1/a drives nu up like a^1, so a protected coherence is lost fast."""
    M = 200
    for nu0 in (1, 2, 3, 10):
        Ns = np.linspace(0.0, 6.0, 60001)
        D = winding_envelope(nu0 * np.exp(Ns), M)
        N_half = Ns[np.argmax(D < 0.5)]
        assert 0.0 < N_half < 0.01           # under a hundredth of an e-fold
    # the fraction travelled is mode-independent: nu -> nu e^N
    for k in (2, 5, 50):
        nu0 = commensurability_ratio(closed_mode_frequency_s3(k, 1.0), closed_mode_frequency_s3(k, 1.0))
        assert drift_over_efolds(nu0, 1.0)[1] == pytest.approx(np.e)
    assert first_unprotected_efolds(1.0, 0.05) == pytest.approx(np.log(1.05))
