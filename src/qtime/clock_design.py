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

from dataclasses import dataclass
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
        return np.rint(self.shifted() / spacing).astype(int)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"Spectrum(d={self.size}, E={self.energies.tolist()}, weights={self.weights.tolist()})"


# ---------------------------------------------------------------------------
# the weighted least-squares core
# ---------------------------------------------------------------------------
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
