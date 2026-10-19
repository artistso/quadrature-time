"""
Where does the clock-diffusion time T_C come from?

The manuscript's dephasing model

    drho/dt = -(i/hbar)[H, rho] - (T_C/2 hbar^2)[H,[H,rho]],
    Gamma_jk = T_C (E_j - E_k)^2 / (2 hbar^2)                        (eq:lindblad)

contains one undetermined parameter, T_C, and the referee's central objection
(REVIEW.md, M2) is that T_C is *assumed* rather than derived: with T_C free, the
model can be tuned to agree with, or to contradict, any experiment, and therefore
makes no prediction.  This module removes that freedom for the concrete case of a
clock that is an electromagnetic cavity mode.

The argument, in three steps
----------------------------
1. A clock mode with mean photon number nbar and photon-loss rate kappa does not
   hold its phase exactly.  The standard result for a coherent field -- the
   Schawlow-Townes linewidth of a laser far above threshold -- is that the phase
   performs Brownian motion with variance growing at

       D_phi = kappa / (2 nbar)      [rad^2 per unit time].

   The physics is that each spontaneous-emission event into the mode shifts the
   phase by of order 1/sqrt(nbar), and such events occur at the cavity rate.

2. The clock's *reading* is t = phi/omega, so a phase variance D_phi * T becomes a
   time variance D_phi T / omega^2.  Comparing with Var(delta t) = T_C * T gives

       T_C = D_phi / omega^2 = kappa / (2 nbar omega^2).               (derived)

3. Substituting into eq:lindblad:

       Gamma_jk = kappa omega_jk^2 / (4 nbar omega^2).                 (derived)

Both are now fixed by measurable device parameters (kappa, nbar, omega) with no
free knob.  That is the whole point: the framework becomes falsifiable.

What this does and does not claim
---------------------------------
It does NOT claim a new decoherence mechanism.  Schawlow-Townes phase diffusion is
textbook, and a noisy clock dephasing a system is Egusquiza-Garay-Raya (1999).  What
is new here is only that the framework's free parameter is *eliminated*, and that
the resulting rate can be compared with measured coherence times.

It also does not by itself distinguish the framework from standard quantum
mechanics: the nbar^-1 and kappa scalings below are shared with ordinary clock
noise.  The value is that the numbers are forced, so the framework can now be
wrong.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "schawlow_townes_diffusion",
    "diffusion_to_clock_time",
    "dephasing_rate",
    "coherence_time",
    "circuit_qed_noise_table",
    "PHOTON_LOSS_1MHZ",
    "CAVITY_6GHZ",
    "QUBIT_5GHZ",
]

# Reference device numbers, in angular frequencies (rad/s).
PHOTON_LOSS_1MHZ = 2.0 * np.pi * 1.0e6      # kappa/2pi = 1 MHz, Q ~ 6000 at 6 GHz
CAVITY_6GHZ = 2.0 * np.pi * 6.0e9           # clock mode
QUBIT_5GHZ = 2.0 * np.pi * 5.0e9            # system gap (transmon splitting)


# ---------------------------------------------------------------------------
# 1. the phase noise of the clock itself
# ---------------------------------------------------------------------------
def schawlow_townes_diffusion(kappa: float, nbar: float) -> float:
    """Phase-diffusion coefficient D_phi = kappa/(2 nbar), in rad^2 per unit time.

    ``kappa`` is the photon-loss rate (not the field decay rate: for amplitude
    damping ``da/dt = -(kappa/2) a``, so the *energy* decays at ``kappa``) and
    ``nbar`` the mean photon number of the clock field.  This is the Schawlow-
    Townes result for a coherent field: the phase walks, the amplitude does not.
    """
    if kappa <= 0.0:
        raise ValueError("kappa must be positive")
    if nbar <= 0.0:
        raise ValueError("nbar must be positive")
    return float(kappa / (2.0 * nbar))


def diffusion_to_clock_time(d_phi: float, omega: float) -> float:
    """Convert a phase-diffusion rate to the clock-diffusion time T_C = D_phi/omega^2.

    The reading is ``t = phi/omega``, so ``Var(phi) = D_phi T`` becomes
    ``Var(delta t) = (D_phi/omega^2) T``; identifying the bracket with T_C is the
    definition used in eq:lindblad.
    """
    if omega <= 0.0:
        raise ValueError("omega must be positive")
    if d_phi < 0.0:
        raise ValueError("diffusion coefficient cannot be negative")
    return float(d_phi / omega**2)


# ---------------------------------------------------------------------------
# 2. the dephasing this forces on the system
# ---------------------------------------------------------------------------
def dephasing_rate(
    kappa: float, nbar: float, omega: float, omega_jk: float
) -> float:
    """Gamma_jk = kappa omega_jk^2 / (4 nbar omega^2), in inverse time.

    Obtained by substituting T_C = kappa/(2 nbar omega^2) into
    ``Gamma_jk = T_C omega_jk^2 / 2``.
    """
    return float(
        dephasing_from_T_C(
            diffusion_to_clock_time(schawlow_townes_diffusion(kappa, nbar), omega), omega_jk
        )
    )


def dephasing_from_T_C(t_c: float, omega_jk: float) -> float:
    """Gamma_jk = T_C omega_jk^2 / 2, the rate in eq:lindblad."""
    if t_c < 0.0:
        raise ValueError("T_C cannot be negative")
    return float(t_c * omega_jk**2 / 2.0)


def coherence_time(gamma: float) -> float:
    """T_2 = 1/Gamma (the coherence time against this channel alone)."""
    if gamma <= 0.0:
        return float("inf")
    return float(1.0 / gamma)


# ---------------------------------------------------------------------------
# 3. the circuit-QED numbers
# ---------------------------------------------------------------------------
def circuit_qed_noise_table(
    kappa: float = PHOTON_LOSS_1MHZ,
    omega: float = CAVITY_6GHZ,
    omega_jk: float = QUBIT_5GHZ,
    nbars: tuple[float, ...] = (1.0, 10.0, 100.0, 1.0e3, 1.0e4),
) -> list[dict[str, float]]:
    """Predicted dephasing for a transmon timed by a 6 GHz cavity clock.

    Returns one dict per ``nbar`` with the diffusion coefficient, the induced
    T_C, the dephasing rate and the resulting T_2.  Nothing here is fitted: every
    number follows from ``kappa``, ``omega``, ``omega_jk`` and ``nbar``.
    """
    rows = []
    for n in nbars:
        d_phi = schawlow_townes_diffusion(kappa, n)
        t_c = diffusion_to_clock_time(d_phi, omega)
        gamma = dephasing_from_T_C(t_c, omega_jk)
        rows.append(
            {
                "nbar": float(n),
                "D_phi": d_phi,
                "T_C": t_c,
                "Gamma": gamma,
                "T_2": coherence_time(gamma),
            }
        )
    return rows
