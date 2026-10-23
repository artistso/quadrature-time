"""
Clock interferometers: an internal quantum clock carried in a superposition of
two paths whose proper-time rates differ (special-relativistic: Lorentz
factors; gravitational: lapse/redshift factors).

Visibility of the centre-of-mass interference (Zych et al. 2011; Pikovski et
al. 2015; Smith & Ahmadi 2020):

    V(t) = | sum_n p_n exp(-i n omega0 * Delta * t) |,

with p_n the internal-clock energy populations and Delta the difference of
proper-time rates between the two arms (Delta = 1/gamma_1 - 1/gamma_2, or
Delta = N(r_2) - N(r_1)).  Because the internal spectrum n*omega0 is
equally spaced, V(t) is EXACTLY periodic with revival period
T_rev = 2 pi / (omega0 |Delta|); the Gaussian envelope
exp(-Delta_n^2 omega0^2 Delta^2 t^2 / 2) is valid only for t << T_rev.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "lorentz_gamma",
    "gaussian_populations",
    "visibility",
    "gaussian_envelope",
    "revival_period",
    "suggested_time_step",
]


def lorentz_gamma(p: float, m: float = 1.0, c: float = 1.0) -> float:
    return float(np.sqrt(1.0 + p**2 / (m**2 * c**2)))


def gaussian_populations(n_max: int, delta_n: float, n0: int | None = None) -> np.ndarray:
    n = np.arange(n_max)
    n0 = n_max // 2 if n0 is None else n0
    p = np.exp(-((n - n0) ** 2) / (2.0 * delta_n**2))
    return p / p.sum()


def visibility(t: np.ndarray, populations: np.ndarray, omega0: float, delta_rate: float) -> np.ndarray:
    n = np.arange(len(populations))
    phases = np.exp(-1j * np.outer(t, n) * omega0 * delta_rate)
    return np.abs(phases @ populations)


def gaussian_envelope(t: np.ndarray, delta_n: float, omega0: float, delta_rate: float) -> np.ndarray:
    return np.exp(-0.5 * (delta_n * omega0 * delta_rate * t) ** 2)


def revival_period(omega0: float, delta_rate: float) -> float:
    return float(2.0 * np.pi / (omega0 * abs(delta_rate)))


def suggested_time_step(delta_n: float, omega0: float, delta_rate: float, samples_per_decay: int = 20) -> float:
    """Resolve the decay time 1/(Delta_n omega0 |Delta|) with the given number of samples."""
    return float(1.0 / (delta_n * omega0 * abs(delta_rate)) / samples_per_decay)
