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
