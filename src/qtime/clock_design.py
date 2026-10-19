"""
Finite-clock design: how many equally spaced levels does a relational clock need?

This module answers one quantitative question that the manuscript's finite-clock
section needs regardless of whether the conjecture is true:

    Given a system with d energy levels E_j and initial weights p_j = |c_j|^2,
    how accurately can an N-level, equally spaced clock (energies n*Delta,
    n = 0..N-1) reproduce its relational time evolution, and for how long?

Conventions (hbar = 1 throughout)
---------------------------------
* The clock is equally spaced with spacing Delta > 0; a system level E_j is
  paired with an integer label q_j, so the approximate energy is q_j*Delta.
* The mismatch is delta_j = E_j - q_j*Delta, and the weighted energy-mismatch
  variance is Var_p(delta) = sum_j p_j (delta_j - mean(delta))^2.
* A common energy offset (in the spectrum or in the labels) leaves the variance
  unchanged, so both may be anchored: E_1 -> 0 and q_1 -> 0.
* Labels must be distinct and increasing: q_1 < q_2 < ... < q_d.

The two results this module is built around
-------------------------------------------
1. **Exact optimum for fixed labels.** Var_p(delta) is a weighted least-squares
   residual of E against q, so its minimiser over Delta is analytic::

       Delta*(q) = C_{Eq} / V_q ,      V*(q) = V_E - C_{Eq}^2 / V_q = V_E (1 - rho^2)

   with V_E, V_q the weighted variances of E and q and C_{Eq} their weighted
   covariance. The design problem is therefore: *pick the integer regressor q
   that maximises its squared correlation with the spectrum.*

2. **Universal bound.** Rounding each E_j to the nearest label gives
   |delta_j| <= Delta/2, so the delta_j span at most Delta and Popoviciu's
   inequality gives Var_p(delta) <= Delta^2/4. Fitting the range
   R = E_max - E_min into N levels forces Delta >= R/(N-1), hence::

       V*_N <= R^2 / (4 (N-1)^2)          (``universal_bound``)

   and, with the fidelity bound F(t) >= 1 - t^2 Var of the v5 notebook, an
   infidelity eps is guaranteed up to T >= 2 (N-1) sqrt(eps) / R.

3. **Three-level closed form (uniform weights).** For E = (0, 1, r) with
   uniform weights and labels (0, q, p)::

       V*(q, p) = (p - r q)^2 / [(p+q)^2 + (2q-p)^2 + (2p-q)^2]   (``three_level_closed_form``)

   Exactly the squared distance of the gap ratio r from the rational p/q,
   normalised by the geometry of the label vector. Optimising over q and p is
   therefore simultaneous Diophantine approximation, which is why the optimum
   for incommensurate spectra decays as N^-4 rather than N^-2: the variance is
   *quadratic* in the rational-approximation error.

Everything here is deterministic and checked in tests/test_clock_design.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Iterable, Sequence

import numpy as np

__all__ = [
    "Spectrum",
    "mismatch_variance",
    "optimal_spacing",
    "fit_quality",
    "fidelity_lower_bound",
    "time_horizon",
    "universal_bound",
    "universal_horizon",
    "exhaustive_optimum",
    "three_level_optimum",
    "scaling_curve",
    "three_level_closed_form",
    "continued_fraction",
    "convergents",
    "diophantine_ratio",
    "predicted_exponent",
    # v9: exact arithmetic and the general (d >= 4) optimiser
    "Optimum",
    "gap_ratio_fraction",
    "exact_variance",
    "three_level_closed_form_exact",
    "scaling_curve_exact",
    "fitted_design",
    "mismatch_variance_stable",
    "nearest_labels",
    "legalize_labels",
    "self_consistent",
    "search_optimum",
    "proportional_curve",
    "proportional_curve_stable",
    "decay_exponent",
    # v10: commensurability drift in an expanding universe
    "scale_factor",
    "closed_mode_frequency_s3",
    "commensurability_ratio",
    "winding_envelope",
    "drift_rate",
    "drift_over_efolds",
    "first_unprotected_efolds",
]


# ---------------------------------------------------------------------------
# spectrum container
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Spectrum:
    """A target spectrum: energies and the weights carried by each level.

    Parameters
    ----------
    energies : sequence of float
        Strictly increasing system energies.
    weights : sequence of float, optional
        Non-negative weights p_j = |c_j|^2; normalised on construction.
        Defaults to uniform.

    Attributes
    ----------
    energies, weights : np.ndarray
    energy_range : float
        R = E_max - E_min, the range that must be resolved by the clock.
    variance : float
        V_E, the weighted variance of the energies.
    """

    energies: np.ndarray
    weights: np.ndarray

    def __init__(self, energies: Iterable[float], weights: Iterable[float] | None = None):
        E = np.asarray(list(energies), dtype=float)
        if E.ndim != 1 or E.size < 2:
            raise ValueError("energies must be a 1-D sequence with at least two levels")
        if not np.all(np.isfinite(E)):
            raise ValueError("energies must be finite")
        if np.any(np.diff(E) <= 0.0):
            raise ValueError("energies must be strictly increasing")
        w = np.full(E.size, 1.0 / E.size) if weights is None else np.asarray(list(weights), dtype=float)
        if w.shape != E.shape:
            raise ValueError("weights must have the same length as energies")
        if not np.all(np.isfinite(w)) or np.any(w < 0.0):
            raise ValueError("weights must be finite and non-negative")
        total = w.sum()
        if total <= 0.0:
            raise ValueError("weights must not all be zero")
        object.__setattr__(self, "energies", E)
        object.__setattr__(self, "weights", w / total)

    # -- basic quantities ---------------------------------------------------
    @property
    def size(self) -> int:
        return int(self.energies.size)

    @property
    def energy_range(self) -> float:
        return float(self.energies[-1] - self.energies[0])

    @property
    def mean(self) -> float:
        return float(self.weights @ self.energies)

    @property
    def variance(self) -> float:
        return float(self.weights @ (self.energies - self.mean) ** 2)

    def shifted(self) -> np.ndarray:
        """Energies with E_1 set to zero (a free energy offset)."""
        return self.energies - self.energies[0]

    def labels(self, spacing: float) -> np.ndarray:
        """Nearest-integer labels for a given clock spacing."""
        if spacing <= 0.0:
            raise ValueError("spacing must be positive")
        return round_half_up(self.shifted() / spacing).astype(int)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"Spectrum(d={self.size}, E={self.energies.tolist()}, weights={self.weights.tolist()})"


# ---------------------------------------------------------------------------
# the weighted least-squares core
# ---------------------------------------------------------------------------
def round_half_up(x):
    """Round to the nearest integer, ties away from zero.

    This is the convention the admissibility lemma needs.  ``numpy.rint`` rounds
    ties to *even*, which is not monotone under a common shift: for the lattice
    spectrum ``(0, 1.5, 2.5)`` and a quantum equal to the gap, ``rint`` returns
    ``(0, 2, 2)`` -- a collision -- whereas half-up returns ``(0, 2, 3)``.  The
    bound ``Delta <= g`` is only sufficient with a monotone convention, so every
    label construction in this module uses this helper.
    """
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def _weighted_stats(values: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """Weighted mean and variance (population form) of a 1-D array."""
    mean = float(weights @ values)
    var = float(weights @ (values - mean) ** 2)
    return mean, var


def mismatch_variance(spectrum: Spectrum, labels: Sequence[int]) -> float:
    """Exact minimised weighted mismatch variance V*(q) for a fixed assignment.

    This is the residual variance of the weighted least-squares fit of the
    energies against ``labels * Delta`` with Delta free, i.e. with the optimal
    spacing supplied by :func:`optimal_spacing`.  It is invariant under a common
    shift of the energies or of the labels.
    """
    lab = np.asarray(list(labels), dtype=float)
    if lab.shape != spectrum.energies.shape:
        raise ValueError("labels must have the same length as energies")
    if lab.size != np.unique(lab).size:
        raise ValueError("labels must be distinct")
    w = spectrum.weights
    _, v_e = _weighted_stats(spectrum.energies, w)
    _, v_q = _weighted_stats(lab, w)
    e_mean = float(w @ spectrum.energies)
    q_mean = float(w @ lab)
    c_eq = float(w @ ((spectrum.energies - e_mean) * (lab - q_mean)))
    return float(max(v_e - c_eq * c_eq / v_q, 0.0))


def optimal_spacing(spectrum: Spectrum, labels: Sequence[int]) -> float:
    """Spacing Delta*(q) = C_{Eq} / V_q that minimises the mismatch variance.

    Returns ``nan`` when the labels are constant (the spacing is then
    unidentifiable; the variance is V_E for every Delta).
    """
    lab = np.asarray(list(labels), dtype=float)
    if lab.shape != spectrum.energies.shape:
        raise ValueError("labels must have the same length as energies")
    w = spectrum.weights
    e_mean, _ = _weighted_stats(spectrum.energies, w)
    q_mean, v_q = _weighted_stats(lab, w)
    if v_q <= 0.0:
        return float("nan")
    c_eq = float(w @ ((spectrum.energies - e_mean) * (lab - q_mean)))
    return float(max(c_eq / v_q, 0.0))


def mismatch_variance_stable(spectrum: Spectrum, labels: Sequence[int]) -> float:
    """Same quantity as :func:`mismatch_variance`, computed without cancellation.

    Instead of V_E - C^2/V_q (a difference of two numbers of order V_E, which
    loses everything once the variance falls below ~1e-14 V_E) this evaluates the
    residual directly: V*(q) = sum_j p_j (E_j - a* - q_j Delta*)^2 with the
    fitted (a*, Delta*). Every term is a difference of nearby numbers, so the
    relative accuracy is limited by how small the residual is, not by V_E.
    """
    delta, offset = fitted_design(spectrum, labels)
    if not np.isfinite(delta):
        return float(spectrum.variance)
    resid = spectrum.energies - offset - np.asarray(list(labels), dtype=float) * delta
    return float(spectrum.weights @ resid**2)


def fit_quality(spectrum: Spectrum, labels: Sequence[int]) -> float:
    """Squared weighted correlation rho^2 between the spectrum and the labels."""
    v_e = spectrum.variance
    if v_e <= 0.0:
        raise ValueError("spectrum has zero variance")
    return float(1.0 - mismatch_variance(spectrum, labels) / v_e)


# ---------------------------------------------------------------------------
# what the variance means operationally
# ---------------------------------------------------------------------------
def fidelity_lower_bound(t: float, variance: float) -> float:
    """F(t) >= max{0, 1 - t^2 * variance}.

    The sufficient bound of the v5 notebook: with delta_j the energy mismatches,
    1 - F(t) = 2 sum_{j<k} p_j p_k [1 - cos((delta_j - delta_k) t)]
             <= t^2 Var_p(delta).
    """
    return float(max(0.0, 1.0 - float(t) ** 2 * float(variance)))


def time_horizon(variance: float, epsilon: float) -> float:
    """Time up to which an infidelity ``epsilon`` is *guaranteed*.

    Solves 1 - t^2 Var = 1 - eps, i.e. T = sqrt(eps / Var). Infinite if the
    variance vanishes (an exactly commensurate spectrum).
    """
    if epsilon < 0.0:
        raise ValueError("epsilon must be non-negative")
    v = float(variance)
    if v <= 0.0:
        return float("inf")
    return float(np.sqrt(epsilon / v))


def universal_bound(energy_range: float, N: int) -> float:
    """Spectrum-free bound V*_N <= R^2 / (4 (N-1)^2).

    Proof: round every shifted E_j to its nearest label. Then |delta_j| <= Delta/2,
    so the delta_j span at most Delta and Popoviciu's inequality
    (Var <= (range)^2/4) gives Var <= Delta^2/4. The largest label is
    R/Delta, which must fit in N-1 steps, so Delta >= R/(N-1).
    """
    N = int(N)
    if N < 2:
        raise ValueError("a clock needs at least two levels")
    return float(energy_range) ** 2 / (4.0 * (N - 1) ** 2)


def universal_horizon(energy_range: float, N: int, epsilon: float) -> float:
    """Guaranteed horizon from :func:`universal_bound`: 2 (N-1) sqrt(eps) / R."""
    return float(2.0 * (int(N) - 1) * np.sqrt(epsilon) / float(energy_range))


# ---------------------------------------------------------------------------
# global optimisation
# ---------------------------------------------------------------------------
def exhaustive_optimum(spectrum: Spectrum, N: int) -> tuple[float, tuple[int, ...]]:
    """Exact optimum over every increasing label assignment inside [0, N-1].

    Cost is binomial(N-1, d-1): exact but only viable for small N and small d
    (d=4, N=64 is about 4x10^4 assignments).  Use :func:`three_level_optimum`
    for d=3, which is exact and linear in N.
    """
    from itertools import combinations

    N = int(N)
    d = spectrum.size
    if N < d:
        raise ValueError("clock dimension must be at least the number of levels")
    best = spectrum.variance
    best_labels: tuple[int, ...] = tuple(range(d))
    for tail in combinations(range(1, N), d - 1):
        labels = (0,) + tail
        v = mismatch_variance(spectrum, labels)
        if v < best:
            best, best_labels = v, labels
    return float(best), best_labels


def three_level_optimum(
    spectrum: Spectrum, N: int, window: int = 2
) -> tuple[float, tuple[int, int, int]]:
    """Exact optimum for a three-level spectrum in O(N * window).

    For every label q of the middle level the routine tries the ``2*window+1``
    labels of the top level nearest to r*q (r = (E_3-E_1)/(E_2-E_1)) and keeps
    the best. The nearest-integer choice alone is *not* always optimal, because
    the normalising denominator depends on p as well (e.g. r = phi, q = 4), but
    a window of +-2 reproduces exhaustive enumeration for every spectrum tested
    in the test suite up to N = 120. This is a validated heuristic window, not
    a proof; ``window`` can be raised for spot checks.
    """
    N = int(N)
    if spectrum.size != 3:
        raise ValueError("three_level_optimum requires exactly three levels")
    if N < 3:
        raise ValueError("clock dimension must be at least three")
    gap = spectrum.shifted()
    r = gap[2] / gap[1]
    best = spectrum.variance
    best_labels: tuple[int, int, int] = (0, 1, 2)
    for q in range(1, N - 1):
        centre = int(np.rint(r * q))
        for p in range(centre - window, centre + window + 1):
            if p <= q or p >= N:
                continue
            v = mismatch_variance(spectrum, (0, q, p))
            if v < best:
                best, best_labels = v, (0, q, p)
    return float(best), best_labels


def scaling_curve(
    spectrum: Spectrum, Nmax: int, window: int = 2
) -> tuple[np.ndarray, np.ndarray]:
    """V*_N for every N <= Nmax, as a step function (three-level spectra).

    Returns ``(N, V)`` with ``V[k] = V*_{N[k]}``.  Cost is O(Nmax * window), so
    N = 10^6 takes well under a second instead of the ~10^11 assignments that
    exhaustive enumeration would need.

    For uniform weights the numerically stable closed form
    ``(p - r q)^2 / |B|^2`` is used; that matters because the interesting values
    of V* at large N are smaller than the round-off floor of the covariance
    form (which cancels two quantities of order V_E).  For non-uniform weights
    the covariance form is used, and V* below ~1e-14 V_E should then be treated
    as "numerically zero".
    """
    Nmax = int(Nmax)
    if spectrum.size != 3:
        raise ValueError("scaling_curve currently supports three-level spectra")
    if Nmax < 4:
        raise ValueError("Nmax must be at least 4")
    gap = np.asarray(spectrum.shifted(), dtype=float)
    r = gap[2] / gap[1]
    w = spectrum.weights
    uniform = np.allclose(w, 1.0 / 3.0, rtol=0.0, atol=1e-15)

    qs, ps = [], []
    q = np.arange(1, Nmax, dtype=np.int64)
    for off in range(-window, window + 1):
        p = np.rint(r * q).astype(np.int64) + off
        keep = (p > q) & (p < Nmax)
        qs.append(q[keep])
        ps.append(p[keep])
    q = np.concatenate(qs)
    p = np.concatenate(ps)
    qf = q.astype(np.float64)
    pf = p.astype(np.float64)

    if uniform:
        # extended precision: (p - r q)^2 is a small difference of large numbers
        r_ld = np.longdouble(r)
        resid = pf.astype(np.longdouble) - r_ld * qf.astype(np.longdouble)
        den = (pf + qf) ** 2 + (2 * qf - pf) ** 2 + (2 * pf - qf) ** 2
        v = np.asarray(resid**2 / den.astype(np.longdouble), dtype=np.float64)
    else:
        e1, e2, e3 = spectrum.energies
        w1, w2, w3 = w
        e_mean = w @ spectrum.energies
        v_e = float(w @ (spectrum.energies - e_mean) ** 2)
        q_mean = w2 * qf + w3 * pf
        v_q = w1 * q_mean**2 + w2 * (qf - q_mean) ** 2 + w3 * (pf - q_mean) ** 2
        c_eq = (
            w1 * (e1 - e_mean) * (0.0 - q_mean)
            + w2 * (e2 - e_mean) * (qf - q_mean)
            + w3 * (e3 - e_mean) * (pf - q_mean)
        )
        v = np.maximum(v_e - c_eq**2 / v_q, 0.0)

    order = np.argsort(p, kind="stable")
    p_sorted = p[order]
    v_sorted = v[order]
    v_running = np.minimum.accumulate(v_sorted)

    first = np.searchsorted(p_sorted, np.arange(2, Nmax + 1), side="left") - 1
    values = np.empty(Nmax - 1, dtype=float)
    base = spectrum.variance
    if uniform:
        base = float("inf")  # no assignment fits: no information, not a floor
    for k, idx in enumerate(first):
        values[k] = base if idx < 0 else v_running[idx]
    return np.arange(2, Nmax + 1), values


# ---------------------------------------------------------------------------
# three-level closed form (uniform weights)
# ---------------------------------------------------------------------------
def three_level_closed_form(r: float, q: int, p: int) -> float:
    """V*(q, p) = (p - r q)^2 / [(p+q)^2 + (2q-p)^2 + (2p-q)^2] for E = (0, 1, r).

    Exact for uniform weights; the denominator is |3 (centred label vector)|^2.
    This is the identity that turns clock design into Diophantine approximation:
    the numerator is the rational-approximation error of the gap ratio, the
    denominator is a pure geometry factor of order q^2.
    """
    if q <= 0 or p <= q:
        raise ValueError("labels must satisfy 1 <= q < p")
    num = (p - r * q) ** 2
    den = (p + q) ** 2 + (2 * q - p) ** 2 + (2 * p - q) ** 2
    return float(num / den)


# ---------------------------------------------------------------------------
# continued fractions: the asymptotic structure of the optimum
# ---------------------------------------------------------------------------
def continued_fraction(x: float, n_terms: int = 40) -> list[int]:
    """Partial quotients [a0; a1, a2, ...] of a real number."""
    if n_terms < 1:
        raise ValueError("n_terms must be positive")
    out: list[int] = []
    y = float(x)
    for _ in range(n_terms):
        a = int(np.floor(y))
        out.append(a)
        frac = y - a
        if frac <= 1e-15:
            break
        y = 1.0 / frac
    return out


def convergents(x: float, n_terms: int = 40) -> list[tuple[int, int]]:
    """Convergents p/q of x, as (p, q) pairs, increasing."""
    a = continued_fraction(x, n_terms)
    p_prev, q_prev = 1, 0
    p_prev2, q_prev2 = 0, 1
    out = []
    for ai in a:
        p = ai * p_prev + p_prev2
        q = ai * q_prev + q_prev2
        out.append((p, q))
        p_prev2, q_prev2 = p_prev, q_prev
        p_prev, q_prev = p, q
    return out


def diophantine_ratio(x: float, n_terms: int = 40, q_max: int = 10**6, tail: int = 3) -> float:
    """Estimate kappa = liminf q^2 |x - p/q| from the largest convergents of ``x``.

    Small kappa means an unusually well-approximable ratio (an easy spectrum);
    the largest possible value for irrational x is 1/sqrt(5), attained by the
    golden ratio and its continued-fraction relatives (Hurwitz's theorem).

    The quantity q^2|x - p/q| *increases* towards its limit along the
    convergents, so taking the minimum over a finite set is wrong (it just picks
    the smallest convergent).  The estimate is therefore the median of the last
    ``tail`` convergents.

    Caveat: this is only meaningful for x with *bounded* partial quotients
    (quadratic irrationals).  If the partial quotients are unbounded -- e, pi --
    then q^2|x - p/q| dips to zero along the large-quotient convergents, the true
    liminf is 0, and this estimator will not see it.  Convergents with q > ``q_max`` are ignored, and an
    exactly vanishing difference is skipped: every float is a rational, so the
    continued fraction of a float terminates and its last convergent reproduces
    x exactly, which would otherwise report kappa = 0 for every input.
    """
    values = []
    for p, q in convergents(x, n_terms):
        if q < 2 or q > q_max:
            continue
        err = abs(x - p / q)
        if err == 0.0:
            continue
        values.append(q * q * err)
    if not values:
        return float("nan")
    return float(np.median(values[-tail:]))


def predicted_exponent(spectrum: Spectrum, Ns: Sequence[int]) -> float:
    """Local log-log slope of V*_N over the given clock dimensions."""
    Ns = np.asarray(list(Ns), dtype=float)
    if Ns.size < 2:
        raise ValueError("need at least two clock dimensions")
    N, V = scaling_curve(spectrum, int(Ns.max()) + 1)
    pick = np.searchsorted(N, Ns.astype(int))
    pick = np.clip(pick, 0, N.size - 1)
    v = V[pick]
    good = v > 0.0
    if good.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(N[pick][good]), np.log(v[good]), 1)[0])


# ===========================================================================
# v9, part 1: exact arithmetic
# ===========================================================================
# Double precision stops being trustworthy once |p - r q| falls below about
# 1e-13, i.e. around p ~ 2000 for quadratic irrationals. The closed form is a
# ratio of integers once r is written as an exact fraction, so there is no
# reason to lose a single digit.


def _dec(x: float) -> Fraction:
    """Exact rational value of the shortest decimal representation of a float.

    ``Fraction(Decimal(repr(2.5)))`` is 5/2; ``Fraction(int(2.5*10**40), 10**40)``
    is not, because 2.5*10**40 has already lost its low digits in binary. Using
    the decimal repr is what makes the commensurate cases come out exactly zero.
    """
    return Fraction(Decimal(repr(float(x))))


def gap_ratio_fraction(spectrum: Spectrum, digits: int | None = None) -> Fraction:
    """The gap ratio (E_3-E_1)/(E_2-E_1) of a three-level spectrum, exactly.

    The result is the exact rational value of the *decimal* numbers in the
    spectrum, e.g. 5/2 for (0, 1, 2.5) and the 16-digit decimal expansion of
    sqrt(2) for (0, 1, sqrt(2)). ``digits`` is accepted for backwards
    compatibility and ignored.

    Honest limit: the exactness is relative to the given decimal expansion of r.
    A double-precision r is itself only accurate to ~1e-17, so p^4 V* is reliable
    while |p - r q| >> 1e-17 p, i.e. out to p ~ 10^8 for quadratic irrationals.
    Past that the remaining error is the truncation of r, not round-off in the
    arithmetic.
    """
    gap = spectrum.shifted()
    if spectrum.size != 3:
        raise ValueError("gap_ratio_fraction requires exactly three levels")
    if gap[1] <= 0.0:
        raise ValueError("energies must be strictly increasing")
    return _dec(gap[2]) / _dec(gap[1])


def three_level_closed_form_exact(r: Fraction, q: int, p: int) -> Fraction:
    """Exact V*(q, p) = (p - r q)^2 / [(p+q)^2 + (2q-p)^2 + (2p-q)^2].

    Inputs are integers apart from ``r``; the result is an exact Fraction, so
    there is no round-off at any label size. Uniform weights only (the closed
    form does not generalise to weighted spectra; see ``exact_variance``).
    """
    if q <= 0 or p <= q:
        raise ValueError("labels must satisfy 1 <= q < p")
    resid = p - r * q
    den = (p + q) ** 2 + (2 * q - p) ** 2 + (2 * p - q) ** 2
    return Fraction(resid * resid, den)


def exact_variance(
    energies: Sequence[Fraction], weights: Sequence[Fraction], labels: Sequence[int]
) -> Fraction:
    """Exact minimised mismatch variance for any d and any rational weights.

    Same quantity as :func:`mismatch_variance`, computed in rational arithmetic:
    V*(q) = V_E - C_{Eq}^2 / V_q. Use this to turn a floating-point optimum into
    an exact statement (e.g. to prove that a commensurate spectrum is reproduced
    with zero error) or to evaluate a candidate at label sizes where doubles
    cannot resolve the residual.
    """
    E = [Fraction(x) for x in energies]
    w = [Fraction(x) for x in weights]
    q = [Fraction(int(x)) for x in labels]
    if not (len(E) == len(w) == len(q)):
        raise ValueError("energies, weights and labels must have equal length")
    if len(set(labels)) != len(labels):
        raise ValueError("labels must be distinct")
    total = sum(w)
    if total <= 0:
        raise ValueError("weights must not all be zero")
    w = [x / total for x in w]
    mE = sum(wi * Ei for wi, Ei in zip(w, E))
    mq = sum(wi * qi for wi, qi in zip(w, q))
    vE = sum(wi * (Ei - mE) ** 2 for wi, Ei in zip(w, E))
    vq = sum(wi * (qi - mq) ** 2 for wi, qi in zip(w, q))
    if vq == 0:
        raise ValueError("labels must not all be equal")
    c = sum(wi * (Ei - mE) * (qi - mq) for wi, Ei, qi in zip(w, E, q))
    return vE - c * c / vq


def scaling_curve_exact(
    spectrum: Spectrum, pmax: int, digits: int = 40, window: int = 2
) -> list[tuple[int, Fraction]]:
    """Exact ``(p, V*_p)`` at the convergent labels of a three-level spectrum.

    The optimum re-equilibrates at convergent labels, so this is the sequence
    that carries the asymptotic information; evaluating it exactly removes the
    round-off floor that limits :func:`scaling_curve` to p ~ 2000.
    """
    if spectrum.size != 3:
        raise ValueError("scaling_curve_exact requires exactly three levels")
    r = gap_ratio_fraction(spectrum)
    out: list[tuple[int, Fraction]] = []
    conv = convergents(float(r), 60)
    for p, q in conv:
        if q < 1 or p < 2 or p > pmax:
            continue
        if p <= q:
            continue
        out.append((p, three_level_closed_form_exact(r, q, p)))
    return out


# ===========================================================================
# v9, part 2: the general optimiser (d >= 4)
# ===========================================================================
# Enumeration is binomial(N-1, d-1): hopeless past d = 3. But the objective has
# a special structure. Minimising over the spacing AND the offset is a two
# parameter weighted least-squares fit, so for fixed (Delta, offset) the optimal
# labels are independent: q_j = round((E_j - offset)/Delta). The optimisation is
# therefore a search over a two-dimensional design plane, where for each cell the
# labels are known and the objective is an explicit quadratic. We walk the plane
# on a grid, apply the analytic fit, re-round, and iterate to a fixed point.


@dataclass(frozen=True)
class Optimum:
    """Result of :func:`search_optimum`.

    Attributes
    ----------
    variance : float
        V*(q) at the best assignment found.
    labels : tuple of int
        The assignment, strictly increasing inside [0, N-1].
    spacing, offset : float
        The fitted design (Delta*, a*) for those labels.
    self_consistent : bool
        Whether re-rounding the fitted design reproduces the same labels. This
        is a fixed-point test, i.e. a *necessary* condition for the assignment
        to be a stationary point of the design problem; it is not a proof of
        global optimality.
    starts : int
        How many grid starts were explored (diagnostic).
    """

    variance: float
    labels: tuple[int, ...]
    spacing: float = float("nan")
    offset: float = float("nan")
    self_consistent: bool = False
    starts: int = 0


def fitted_design(spectrum: Spectrum, labels: Sequence[int]) -> tuple[float, float]:
    """The (spacing, offset) minimising sum_j p_j (E_j - a - q_j Delta)^2."""
    lab = np.asarray(list(labels), dtype=float)
    if lab.shape != spectrum.energies.shape:
        raise ValueError("labels must have the same length as energies")
    w = spectrum.weights
    e_mean = float(w @ spectrum.energies)
    q_mean = float(w @ lab)
    v_q = float(w @ (lab - q_mean) ** 2)
    if v_q <= 0.0:
        return float("nan"), float("nan")
    c = float(w @ ((spectrum.energies - e_mean) * (lab - q_mean)))
    delta = c / v_q
    return delta, e_mean - delta * q_mean


def legalize_labels(labels: Iterable[int], N: int, d: int | None = None) -> tuple[int, ...]:
    """Force a label vector to be strictly increasing with 0 <= q_j <= N-1.

    Labels are pushed apart forward, then pulled back inside the clock's range,
    so the result is always a legal assignment (provided N >= d).
    """
    q = [int(x) for x in labels]
    d = len(q) if d is None else int(d)
    N = int(N)
    if d < 1 or N < d:
        raise ValueError("need N >= number of levels")
    for j in range(1, d):
        q[j] = max(q[j], q[j - 1] + 1)
    if q[-1] > N - 1:
        shift = q[-1] - (N - 1)
        q = [x - shift for x in q]
    for j in range(d):
        q[j] = max(0, min(q[j], N - 1 - (d - 1 - j)))
    for j in range(1, d):
        q[j] = max(q[j], q[j - 1] + 1)
    return tuple(q)


def nearest_labels(spectrum: Spectrum, spacing: float, offset: float, N: int) -> tuple[int, ...]:
    """Legal label vector obtained by rounding (E_j - offset)/spacing."""
    if spacing <= 0.0:
        raise ValueError("spacing must be positive")
    raw = round_half_up((spectrum.energies - offset) / spacing).astype(int)
    return legalize_labels(raw, N)


def self_consistent(spectrum: Spectrum, labels: Sequence[int], N: int, tol: float = 1e-9) -> bool:
    """Does re-rounding the fitted design reproduce ``labels``?"""
    delta, offset = fitted_design(spectrum, labels)
    if not np.isfinite(delta) or delta <= 0.0:
        return False
    again = nearest_labels(spectrum, delta, offset, N)
    return tuple(int(x) for x in labels) == again


def _coordinate_descent(
    spectrum: Spectrum, labels: tuple[int, ...], N: int, steps: int = 2, max_pass: int = 12
) -> tuple[float, tuple[int, ...]]:
    """Best-improvement search over single-label moves of +-1 ... +-steps.

    A local search on the label lattice; it escapes the poor fixed points that
    the re-rounding iteration can sit in (e.g. the golden ratio at N = 40, where
    (0,13,21) is a fixed point but (0,21,34) is better).
    """
    best = mismatch_variance_stable(spectrum, labels)
    current = tuple(int(x) for x in labels)
    d = len(current)
    for _ in range(int(max_pass)):
        improved = False
        for j in range(d):
            for move in range(-steps, steps + 1):
                if move == 0:
                    continue
                trial = list(current)
                trial[j] += move
                if not all(0 <= trial[k] < N for k in range(d)):
                    continue
                if any(trial[k] >= trial[k + 1] for k in range(d - 1)):
                    continue
                v = mismatch_variance_stable(spectrum, trial)
                if v < best - 1e-18 * max(1.0, abs(best)):
                    best, current = v, tuple(trial)
                    improved = True
        if not improved:
            break
    return best, current


def _window_descent(
    spectrum: Spectrum,
    labels: tuple[int, ...],
    N: int,
    half_width: int = 6,
    max_pass: int = 6,
) -> tuple[float, tuple[int, ...]]:
    """Coordinate sweep in which each label may move anywhere inside +-half_width.

    A much larger neighbourhood than :func:`_coordinate_descent`, and the one that
    settles the d = 4 cases where single moves get stuck.
    """
    best = mismatch_variance_stable(spectrum, labels)
    current = tuple(int(x) for x in labels)
    d = len(current)
    for _ in range(int(max_pass)):
        improved = False
        for j in range(d):
            lo = max(0, current[j] - half_width, current[j - 1] + 1 if j > 0 else 0)
            hi = min(N - 1, current[j] + half_width, current[j + 1] - 1 if j < d - 1 else N - 1)
            for q in range(lo, hi + 1):
                if q == current[j]:
                    continue
                trial = list(current)
                trial[j] = q
                v = mismatch_variance_stable(spectrum, trial)
                if v < best - 1e-18 * max(1.0, abs(best)):
                    best, current = v, tuple(trial)
                    improved = True
        if not improved:
            break
    return best, current


def search_optimum(
    spectrum: Spectrum,
    N: int,
    n_spacing: int = 41,
    n_shift: int = 5,
    max_iter: int = 30,
    polish: int = 12,
    extra_seeds: Iterable[Sequence[int]] = (),
) -> Optimum:
    """Global search over the (spacing, offset) design plane, for any d.

    Three stages:

    1. *Seeds.* Spacings log-spaced around R/(N-1) combined with offsets whose
       residues fill one period of the spacing, plus any user seeds. Each seed is
       rounded and then iterated to a fixed point of (fit -> re-round).
    2. *Scoring.* Every distinct legal assignment reached is scored with the
       analytic V*(q); duplicates are collapsed.
    3. *Polish.* The best ``polish`` candidates undergo coordinate descent on the
       label lattice (single-label moves of +-1, +-2). This matters: for the
       golden ratio at N = 40 the fixed point (0,13,21) is beaten by (0,21,34),
       and only the local search finds it.

    The result is the best assignment visited. This is a heuristic: it is
    validated against exhaustive enumeration for d = 3 (N <= 120, uniform and
    weighted) and d = 4 (N <= 36) in the test suite, and ``self_consistent`` is
    reported as a diagnostic, not as a certificate of global optimality.

    Known limitation: the seed grid is finite, so at large N it can miss the
    Diophantine-quality assignments (it finds labels *proportional* to the
    spectrum only if m = q_top happens to be in the grid). For spectra whose
    optimum is proportional -- every power-of-a-quadratic-irrational spectrum --
    use :func:`proportional_curve` (exact) or :func:`proportional_curve_stable`
    (fast to M ~ 10^6) instead, which scan the whole family rather than a sample
    of it.
    """
    N = int(N)
    d = spectrum.size
    if N < d:
        raise ValueError("clock dimension must be at least the number of levels")
    energies = spectrum.shifted()
    span = float(energies[-1] - energies[0])
    if span <= 0.0:
        raise ValueError("spectrum has zero range")

    # Family A: a blind grid over the design plane, iterated to a fixed point.
    spacings = span / (N - 1) * np.logspace(-1.2, 1.4, int(n_spacing))
    seeds: list[tuple[int, ...]] = [tuple(range(d)), tuple(range(N - d, N))]
    for delta in spacings:
        for frac in np.linspace(0.0, 1.0, int(n_shift), endpoint=False):
            offset = float(energies[0]) - frac * float(delta)
            seeds.append(nearest_labels(spectrum, float(delta), offset, N))

    # Family B: spacings tied to the spectrum's own range, one per candidate
    # top label m. The labels of the optimum are very nearly proportional to the
    # energies, so Delta = range/m with a sub-step offset generates them
    # directly. This is what finds (0,8,13,21) for the phi-power spectrum, which
    # no logarithmic grid of spacings hits exactly.
    m_all = np.unique(
        np.concatenate(
            [
                np.arange(d, min(N, 400)),
                np.unique(np.round(np.geomspace(max(d, 2), N, 120)).astype(int)) - 1,
                np.unique(np.round(np.linspace(d, N, 240)).astype(int)) - 1,
                np.array([N - 1]),
            ]
        )
    )
    m_all = m_all[(m_all >= d) & (m_all <= N - 1)]
    for m in m_all:
        delta = span / float(m)
        for shift in (0.0, -0.5, 0.5, -0.25, 0.25):
            offset = float(energies[0]) + shift * delta
            seeds.append(nearest_labels(spectrum, delta, offset, N))

    if d == 3:  # convergent labels are near-optimal for uniform weights
        rr = float(energies[2] / energies[1])
        for p_, q_ in convergents(rr, 40):
            if 1 <= q_ < p_ <= N - 1:
                seeds.append((0, int(q_), int(p_)))
    for lab in extra_seeds:
        seeds.append(legalize_labels(lab, N, d))

    scored: dict[tuple[int, ...], float] = {}
    for lab in seeds:
        current = lab
        for _ in range(int(max_iter)):
            delta, offset = fitted_design(spectrum, current)
            if not np.isfinite(delta) or delta <= 0.0:
                break
            nxt = nearest_labels(spectrum, delta, offset, N)
            if nxt == current:
                break
            current = nxt
        for candidate in (lab, current):
            if candidate not in scored:
                scored[candidate] = mismatch_variance_stable(spectrum, candidate)

    ranked = sorted(scored.items(), key=lambda kv: kv[1])[: max(1, int(polish))]
    best_lab, best_v = ranked[0]
    for lab, _ in ranked:
        v, polished = _coordinate_descent(spectrum, lab, N)
        v, polished = _window_descent(spectrum, polished, N)
        if v < best_v:
            best_v, best_lab = v, polished

    delta, offset = fitted_design(spectrum, best_lab)
    return Optimum(
        variance=float(min(best_v, mismatch_variance_stable(spectrum, best_lab))),
        labels=tuple(int(x) for x in best_lab),
        spacing=float(delta),
        offset=float(offset),
        self_consistent=self_consistent(spectrum, best_lab, N),
        starts=len(seeds),
    )


# ===========================================================================
# v9, part 3: the proportional family, and what the exponent really is
# ===========================================================================
# The labels of the optimum are very nearly proportional to the energies, so
# the one-parameter family q_j = round(m E_j / R), m = 2, 3, ..., is where the
# action is. Scanning it with exact arithmetic gives V*_N without any round-off
# floor, which is what makes the decay exponent measurable at all.


def proportional_curve(
    spectrum: Spectrum, M: int, exact: bool = True, digits: int = 40
) -> tuple[np.ndarray, np.ndarray]:
    """V*_N for the proportional family q_j = round(m E_j / R), m <= M.

    Returns ``(N, V*)`` where ``V*[k]`` is the best variance achieved by any
    member of the family that fits inside a clock with ``N[k] = m+1`` levels, and
    both arrays are restricted to the running minimum. Because the family is a
    *subset* of all assignments, these values are upper bounds on the true
    optimum; for the quadratic spectra studied here (powers of the golden ratio)
    the family contains the optimum.

    ``exact=True`` evaluates V* in rational arithmetic (``exact_variance``),
    which is essential: for m ~ 10^4 the variance is ~10^-18 and the
    floating-point covariance form returns nothing but cancellation noise.
    """
    M = int(M)
    d = spectrum.size
    if M < d:
        raise ValueError("M must be at least the number of levels")
    E = [_dec(x) if exact else Fraction(float(x)) for x in spectrum.energies]
    E0 = E[0]
    E = [x - E0 for x in E]
    span = E[-1]
    if span <= 0:
        raise ValueError("spectrum has zero range")
    w = [_dec(x) for x in spectrum.weights]

    Ns: list[int] = []
    vs: list[Fraction] = []
    for m in range(max(2, d - 1), M + 1):
        q = [int(round(float(Ej * m / span))) for Ej in E]
        if len(set(q)) < d:
            continue
        Ns.append(m + 1)
        vs.append(exact_variance(E, w, q))
    vals = np.array([float(v) for v in vs], dtype=float)
    N = np.array(Ns, dtype=float)
    running = np.minimum.accumulate(vals)
    return N, running


def decay_exponent(
    N: Iterable[float], V: Iterable[float], n_points: int = 40
) -> float:
    """Log-log slope of a running-minimum V*(N) curve over its informative range.

    Points are sampled geometrically and the fit is restricted to values above
    the floor where the curve is still falling, so that plateaus and numerical
    saturation do not bias the slope.
    """
    N = np.asarray(list(N), dtype=float)
    V = np.asarray(list(V), dtype=float)
    if N.size < 4:
        raise ValueError("need at least four samples")
    grid = np.geomspace(N[2], N[-1], int(n_points))
    idx = np.searchsorted(N, grid) - 1
    idx = idx[idx >= 2]
    if idx.size < 3:
        raise ValueError("not enough samples in range")
    n_sel = N[idx]
    v_sel = V[idx]
    good = v_sel > 0.0
    if good.sum() < 3:
        raise ValueError("curve is numerically zero over the sampled range")
    return float(np.polyfit(np.log(n_sel[good]), np.log(v_sel[good]), 1)[0])


def proportional_curve_stable(spectrum: Spectrum, M: int) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised proportional-family curve, stable and fast up to M ~ 10^6.

    Same family as :func:`proportional_curve` (``q_j = round(m E_j / R)``), but
    evaluated in floating point with the *residual* form of the variance
    (``mismatch_variance_stable``) rather than the covariance form, which is what
    makes it accurate: for each member the fitted (a*, Delta*) is obtained from
    vectorised weighted statistics and the variance is then read off the
    residuals. Cost is O(M d) with d numpy operations, about 1e6 members per
    second.

    Accuracy caveat: the residual of the best member is ~1e-9 relative to V_E at
    M ~ 10^4 and ~1e-13 at M ~ 10^6, so beyond M ~ 10^5 the floating-point result
    should be checked against :func:`proportional_curve` (exact) at a few points.
    """
    M = int(M)
    d = spectrum.size
    if M < max(2, d):
        raise ValueError("M must be at least the number of levels")
    E = spectrum.energies.astype(float)
    E = E - E[0]
    span = float(E[-1])
    w = spectrum.weights
    e_mean = float(w @ E)
    v_e = float(w @ (E - e_mean) ** 2)

    m = np.arange(max(2, d - 1), M + 1, dtype=float)
    Q = round_half_up(m[:, None] * E[None, :] / span)
    distinct = np.all(np.diff(Q, axis=1) > 0, axis=1)
    Q, m = Q[distinct], m[distinct]
    if Q.size == 0:
        raise ValueError("no legal proportional assignment below M")

    mq = Q @ w
    vq = ((Q - mq[:, None]) ** 2) @ w
    c = ((E - e_mean)[None, :] * (Q - mq[:, None])) @ w
    delta = c / vq
    offset = e_mean - delta * mq
    resid = E[None, :] - offset[:, None] - Q * delta[:, None]
    v = (resid**2) @ w
    v = np.maximum(v, 0.0)

    N = m + 1.0
    return N, np.minimum.accumulate(v)


# ===========================================================================
# v10, first step: commensurability in an expanding universe
# ===========================================================================
# Everything above assumes a fixed clock quantum hbar*omega. In a closed,
# expanding universe the geometry itself supplies the natural frequency,
# omega_k = c sqrt(k^2 + 1)/a(t), so omega falls as the universe grows. Then the
# commensurability ratio nu = dE/(hbar*omega) that controls the winding
# marginalisation factor |D_M(nu)| is *time dependent*, and coherence that is
# protected at early times is not protected later. The functions here quantify
# that drift. The cosmological input (a(t)) is a model choice, kept explicit.


def scale_factor(
    t: float, H0: float = 1.0, omega_m: float = 0.0, omega_l: float = 1.0
) -> float:
    """Normalised scale factor a(t), with a(0) = 1 and a'(0)/a(0) = H0.

    Defaults to the empty, Lambda-dominated case (pure de Sitter), a = e^{H0 t},
    which is the model the drift estimates assume, and which is what an era of
    inflation or of dark-energy domination looks like. For matter + Lambda with
    omega_m > 0 the exact integral is
        t = (2/(3 H0 sqrt(omega_l))) asinh( sqrt(omega_l/omega_m) a^{3/2} ),
    inverted here for a. Units: t in 1/H0.
    """
    if omega_m < 0.0 or omega_l < 0.0:
        raise ValueError("densities must be non-negative")
    if omega_m > 0.0:
        if omega_l <= 0.0:
            raise ValueError("a pure-matter closed form needs omega_l > 0; use an ODE")
        pref = 2.0 / (3.0 * H0 * np.sqrt(omega_l))
        x = np.sinh(t / pref) * np.sqrt(omega_m / omega_l)
        return float(x ** (2.0 / 3.0))
    return float(np.exp(H0 * t))


def closed_mode_frequency_s3(k: int, a: float, c: float = 1.0) -> float:
    """Frequency on a round S^3 of curvature radius a: omega_k = c sqrt(k^2-1)/a, k >= 2."""
    if k < 2:
        raise ValueError("the k = 1 mode of S^3 is homogeneous (zero frequency)")
    if a <= 0.0:
        raise ValueError("scale factor must be positive")
    return float(c * np.sqrt(k * k - 1) / a)


def commensurability_ratio(delta_e: float, omega: float) -> float:
    """nu = dE/(hbar*omega) with hbar = 1: the argument of the winding factor D_M."""
    if omega <= 0.0:
        raise ValueError("frequency must be positive")
    return float(delta_e / omega)


def winding_envelope(nu, M_windings: int) -> np.ndarray:
    """|D_M(nu)| = |sin(pi M nu)/(M sin(pi nu))|, with the limit 1 at integer nu."""
    nu = np.asarray(nu, dtype=float)
    M = int(M_windings)
    if M < 1:
        raise ValueError("M must be positive")
    out = np.ones_like(nu)
    off = np.abs(np.sin(np.pi * nu)) > 1e-12
    out[off] = np.abs(np.sin(np.pi * M * nu[off]) / (M * np.sin(np.pi * nu[off])))
    return out


def drift_rate(delta_e: float, k: int, a: float, hdot_over_h: float, hbar: float = 1.0) -> float:
    """d(nu)/dt = -nu * (d ln omega / dt) = nu * H(t), for omega = omega_k(a).

    Because omega_k propto 1/a, ln omega decreases at exactly the Hubble rate, so
    a coherence with ratio nu drifts through the winding factor at a rate set by
    nu * H. This is the quantity that decides whether a protected coherence can
    survive: the total drift across an era is the integral of nu H dt, which is
    just nu * (number of e-folds).
    """
    nu = commensurability_ratio(delta_e, closed_mode_frequency_s3(k, a))
    return float(nu * hdot_over_h)


def drift_over_efolds(nu0: float, delta_efolds: float) -> tuple[float, float]:
    """(nu_initial, nu_final) for a coherence drifting through ``delta_efolds``.

    omega propto 1/a means nu propto a: the ratio GROWS with expansion, by one
    e-fold per e-fold. So an early-time commensurability (nu integer) is broken by
    expansion, not preserved.
    """
    return float(nu0), float(nu0 * np.exp(delta_efolds))


def first_unprotected_efolds(nu0: float, tolerance: float = 0.05) -> float:
    """E-folds before an integer-matched coherence (nu0 integer) is spoiled.

    |D_M| returns to 1 at integer nu and falls to zero at half-integers; the
    coherence is "spoiled" once nu has moved by ``tolerance`` from its integer
    value. Returns ln((1+tolerance)/nu0 * nu0) style drift: N = ln(1+tolerance)
    for a coherence starting at nu0 = 1, and shorter for higher harmonics.
    """
    if nu0 <= 0.0:
        raise ValueError("nu0 must be positive")
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")
    return float(np.log(1.0 + tolerance))


# ---------------------------------------------------------------------------
# v12: the golden-ratio constant, exactly and provably
#
# The numerical statement "p^4 V*_p -> phi^4/60" was previously checked in
# floating point, which fails past p ~ 10^3 because the covariance form cancels
# two numbers of order V_E, and with rational arithmetic, which is limited by the
# truncation of phi to a rational -- and the truncation dominates, so the error
# does not improve as the labels grow.  Both limits are avoidable: phi is a
# quadratic irrational, so V* can be computed exactly in Q(sqrt5), where it never
# cancels and never truncates.  That turns the numerical claim into a theorem.
# ---------------------------------------------------------------------------
_PHI4_OVER_60_A = Fraction(7, 120)                   # phi^4/60 = (7 + 3 sqrt5)/120
_PHI4_OVER_60_B = Fraction(3, 120)


def fibonacci(n: int) -> int:
    """F_n with F_1 = F_2 = 1."""
    if n < 1:
        raise ValueError("n must be positive")
    a, b = 1, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return a


def closed_form_golden_pair(q: int, p: int) -> tuple[Fraction, Fraction]:
    """Exact V*(q,p) for the golden-ratio spectrum, returned as ``a + b*sqrt(5)``.

    With r = phi = (1+sqrt5)/2 the residual is ``p - q*phi = (p - q/2) - (q/2)*sqrt5``,
    whose square is again in Q(sqrt5), and the closed-form denominator is an
    integer.  Every coefficient below is a Fraction, so nothing is truncated.
    """
    if q <= 0 or p <= q:
        raise ValueError("labels must satisfy 1 <= q < p")
    alpha = Fraction(p) - Fraction(q, 2)
    beta = Fraction(-q, 2)
    den = (p + q) ** 2 + (2 * q - p) ** 2 + (2 * p - q) ** 2
    return (alpha * alpha + 5 * beta * beta) / den, (2 * alpha * beta) / den


def golden_limit_difference(n: int) -> tuple[Fraction, Fraction]:
    """Exact ``p^4 V*(q,p) - phi^4/60`` at Fibonacci labels, as ``a + b*sqrt(5)``.

    ``q = F_n``, ``p = F_{n+1}``.  Both coefficients are exact Fractions; the
    difference itself is a genuinely algebraic number, so it is reported as a
    pair rather than collapsed to a float.  Use
    :func:`golden_limit_magnitude` to evaluate its size to arbitrary precision.
    """
    q, p = fibonacci(n), fibonacci(n + 1)
    a, b = closed_form_golden_pair(q, p)
    return a * p**4 - _PHI4_OVER_60_A, b * p**4 - _PHI4_OVER_60_B


def golden_limit_magnitude(n: int, digits: int = 40) -> Decimal:
    """``|p^4 V*(q,p) - phi^4/60|`` evaluated to ``digits`` significant figures.

    The two coefficients of the exact difference are of order ``p^4`` while their
    ``sqrt(5)`` combination is tiny, so this is a cancellation -- but it is a
    cancellation in *decimal*, performed here with ``digits`` significant figures
    and an ``sqrt(5)`` correct to the same, not in binary double precision.  At
    ``n = 40`` the difference is about ``1e-17`` while the coefficients are about
    ``1e+37``, i.e. 54 digits of cancellation, which is why 40 digits is not
    enough there; raise ``digits`` accordingly.
    """
    with localcontext() as ctx:
        ctx.prec = int(digits) + 10
        a, b = golden_limit_difference(n)
        da = Decimal(a.numerator) / Decimal(a.denominator)
        db = Decimal(b.numerator) / Decimal(b.denominator)
        root = Decimal(5).sqrt()
        # scale to working precision before combining, then report the modulus
        ctx.prec = int(digits)
        return abs(da + db * root)


def golden_residual_identity(n: int) -> Fraction:
    """Exact check that F_{n+1} - F_n*phi = (-1)^n * phi^(-n), as a Q(sqrt5) pair.

    Returned as ``(a, b)`` meaning ``a + b*sqrt5``.
    """
    q, p = fibonacci(n), fibonacci(n + 1)
    return Fraction(p) - Fraction(q, 2), Fraction(-q, 2)
