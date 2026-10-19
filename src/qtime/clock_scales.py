"""
Which clocks can be relational clocks, and what does tuning the spacing buy?

``qtime.clock_design`` treats the level spacing ``Delta = hbar*omega`` as a free
design parameter. That is right for an *ideal* clock, and it is what produces the
``N^-4`` law of Proposition 4: choose ``Delta`` at a convergent of the gap ratio
and the mismatch variance falls four powers of ``N``. Real clocks do not offer a
continuum of quanta -- a cavity mode, an atomic line or a Zeeman splitting is
fixed by nature -- so two questions remain:

1. **Admissibility.** Given a quantum ``Delta``, which spectra can it label at
   all? With ``d`` levels spanning ``R`` and smallest gap ``g``:

   * necessary: ``Delta <= R/(d-1)`` (``d`` distinct integers need that span);
   * sufficient: ``Delta <= g`` (then rounding is strictly increasing);
   * cost: ``N >= 1 + R/Delta`` levels to cover the range.

2. **What tuning buys.** A fixed quantum forces the rounding assignment, which is
   exactly the construction behind the universal bound
   ``V*_N <= R^2/(4(N-1)^2)``. So an untunable clock is stuck with the
   universal ``N^-2`` variance, i.e. a horizon growing only like ``N^1``, whereas
   a tunable one achieves ``N^-4`` and a horizon like ``N^2``. The difference is
   not a constant: it grows as ``N^2``.

That second point is the useful one, and it is a *scope condition on Proposition
4* rather than a new effect: the ``N^-4`` law is a statement about clock design,
not about clocks as found.
"""
from __future__ import annotations

import numpy as np

from .clock_design import Spectrum, three_level_optimum, universal_bound

__all__ = [
    "design_rule",
    "labels_for_quantum",
    "admissible",
    "levels_needed",
    "mean_variance_at_quantum",
    "horizon_at_quantum",
    "tuning_gain",
    "quantum_from_wavelength",
    "quantum_from_energy",
    "PHYSICAL_CLOCKS",
]

_H_PLANCK = 6.62607015e-34
_C_LIGHT = 2.99792458e8
_EV = 1.602176634e-19


# ---------------------------------------------------------------------------
# 1. admissibility
# ---------------------------------------------------------------------------
def design_rule(spectrum: Spectrum) -> dict[str, float]:
    """Admissible quanta and the level budget for a spectrum.

    ``range``                  R = E_d - E_1
    ``min_gap``                g = min_j (E_{j+1} - E_j); Delta <= g suffices
    ``quantum_max_necessary``  R/(d-1); Delta <= R/(d-1) is necessary
    ``level_floor``            1 + ceil(R/g); no clock can do it with fewer levels
    """
    E = spectrum.energies
    d = spectrum.size
    rng = float(E[-1] - E[0])
    g = float(np.diff(E).min())
    return {
        "range": rng,
        "min_gap": g,
        "quantum_max_necessary": rng / (d - 1),
        "level_floor": float(1.0 + np.ceil(rng / g)),
    }


def labels_for_quantum(
    spectrum: Spectrum, quantum: float, offset: float | None = None
) -> tuple[int, ...]:
    """Nearest-integer labels for a fixed quantum, with an optional offset."""
    if quantum <= 0.0:
        raise ValueError("quantum must be positive")
    a = float(spectrum.energies[0]) if offset is None else float(offset)
    return tuple(int(x) for x in np.rint((spectrum.energies - a) / quantum).astype(int))


def admissible(spectrum: Spectrum, quantum: float, n_offsets: int = 20001) -> bool:
    """Can some offset place every level in its own bin?

    The sufficient test ``quantum <= min_gap`` is tried first, then the
    necessary test ``quantum <= R/(d-1)``; between those bounds a scan over the
    offset decides the question exactly rather than conservatively.
    """
    if quantum <= 0.0:
        raise ValueError("quantum must be positive")
    rule = design_rule(spectrum)
    if quantum <= rule["min_gap"]:
        return True
    if quantum > rule["quantum_max_necessary"]:
        return False
    d = spectrum.size
    for a in spectrum.energies[0] + np.linspace(-quantum, quantum, int(n_offsets)):
        if len(set(labels_for_quantum(spectrum, quantum, a))) == d:
            return True
    return False


def levels_needed(spectrum: Spectrum, quantum: float) -> int:
    """Fewest levels that cover the range: 1 + ceil(R/Delta)."""
    if quantum <= 0.0:
        raise ValueError("quantum must be positive")
    rng = float(spectrum.energies[-1] - spectrum.energies[0])
    return int(1 + np.ceil(rng / quantum))


# ---------------------------------------------------------------------------
# 2. what a fixed quantum achieves, and what tuning buys
# ---------------------------------------------------------------------------
def _variance_at(spectrum: Spectrum, quantum: float, offset: float) -> float:
    q = labels_for_quantum(spectrum, quantum, offset)
    if len(set(q)) < spectrum.size:
        return float("nan")
    resid = spectrum.energies - np.asarray(q, dtype=float) * quantum
    mean = float(spectrum.weights @ resid)
    return float(spectrum.weights @ (resid - mean) ** 2)


def mean_variance_at_quantum(
    spectrum: Spectrum, quantum: float, n_offsets: int = 4001
) -> float:
    """Mismatch variance for a fixed quantum, averaged over the clock's zero.

    The zero of a clock is not physical, so for a *given* quantum the honest
    figure of merit is the average over offsets (or the worst case). Returns
    ``nan`` when no offset gives distinct labels, i.e. when the quantum is
    inadmissible.
    """
    offs = spectrum.energies[0] + np.linspace(0.0, quantum, int(n_offsets), endpoint=False)
    vals = np.array([_variance_at(spectrum, quantum, a) for a in offs])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan")
    return float(vals.mean())


def horizon_at_quantum(
    spectrum: Spectrum, quantum: float, epsilon: float = 1e-2, n_offsets: int = 4001
) -> float:
    """Guaranteed horizon for a fixed quantum, in the same units as the energies.

    Uses the Lemma of Sec. V.C, ``F >= 1 - t^2 Var/hbar^2``, with the
    offset-averaged variance. Zero for an inadmissible quantum.
    """
    v = mean_variance_at_quantum(spectrum, quantum, n_offsets)
    if not np.isfinite(v):
        return 0.0
    if v <= 0.0:
        return float("inf")
    return float(np.sqrt(float(epsilon) / v))


def tuning_gain(spectrum: Spectrum, N: int, n_offsets: int = 4001) -> float:
    """How much better the *tuned* design is than a fixed, bandwidth-matched one.

    Compares the optimised variance ``V*_N`` (spacing free) with the variance of a
    clock whose quantum is pinned to ``R/(N-1)`` and whose zero is unknown. The
    ratio grows quadratically in ``N``, which is the whole difference between the
    universal ``N^-2`` variance and the designed ``N^-4`` one.
    """
    rng = float(spectrum.energies[-1] - spectrum.energies[0])
    if N < 3:
        raise ValueError("need at least three levels")
    quantum = rng / (N - 1)
    fixed = mean_variance_at_quantum(spectrum, quantum, n_offsets)
    if not np.isfinite(fixed) or fixed <= 0.0:
        return float("nan")
    if spectrum.size == 3:
        tuned = three_level_optimum(spectrum, N)[0]
    else:                       # no exact three-level optimiser: compare to the bound
        tuned = universal_bound(rng, N)
    if tuned <= 0.0:
        return float("inf")
    return float(fixed / tuned)


# ---------------------------------------------------------------------------
# 3. physical quanta
# ---------------------------------------------------------------------------
def quantum_from_energy(energy_joules: float) -> float:
    """Clock quantum in eV from an energy in joules."""
    return float(energy_joules) / _EV


def quantum_from_wavelength(wavelength_nm: float) -> float:
    """Clock quantum in eV for a photon of the given wavelength in nanometres."""
    lam = float(wavelength_nm) * 1e-9
    if lam <= 0.0:
        raise ValueError("wavelength must be positive")
    return float(_H_PLANCK * _C_LIGHT / lam / _EV)


#: Reference quanta for clocks one might actually build, in eV.
PHYSICAL_CLOCKS: dict[str, float] = {
    "optical 507 nm": quantum_from_wavelength(507.0),
    "optical 1550 nm (telecom)": quantum_from_wavelength(1550.0),
    "hydrogen 21 cm line": quantum_from_wavelength(2.1e8),
    "caesium hyperfine 9.19 GHz": 9.192631770e9 * _H_PLANCK / _EV,
    "Zeeman splitting, 1 T": 13.996e9 * _H_PLANCK / _EV,
    "microwave cavity 6 GHz": 6.0e9 * _H_PLANCK / _EV,
    "transmon splitting 5 GHz": 5.0e9 * _H_PLANCK / _EV,
}
