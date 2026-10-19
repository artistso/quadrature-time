"""
Finite-dimensional models of relational ("quadrature") time.

Everything here is built from an explicit bipartite state |Psi> in
C^N (clock) (x) C^d (system); conditional system states are obtained by
actual partial traces against clock POVM elements -- no shortcut
convolution of pre-computed "ideal" states.

Conventions (match the manuscript)
----------------------------------
* Clock energy basis |n>, n = 0..N-1, H_C = hbar*omega*N  (hbar = omega = 1).
* Pegg--Barnett phase states  |theta_k> = N^{-1/2} sum_n e^{i n theta_k}|n>,
  theta_k = 2 pi k / N.  They form an orthonormal basis (a PVM).
* Relational constraint  (hbar*omega*N (x) 1 - 1 (x) H_S)|Psi> = 0, so that a
  system level of energy E_k = n_k (integer, in units of hbar*omega) is
  paired with clock level n_k.  Conditioning on the bra (theta| gives
  |Psi(theta)> = sum_k alpha_k e^{-i n_k theta}|E_k>, i.e. forward
  Schroedinger evolution with t = theta/omega.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "phase_grid",
    "pegg_barnett_basis",
    "pw_constraint_state",
    "sharp_phase_povm",
    "smeared_phase_povm",
    "wrapped_gaussian_weights",
    "condition_on_clock",
    "coherence_suppression_static",
    "diffusive_clock_channel",
    "lindblad_dephasing_solution",
    "integrate_lindblad",
    "winding_marginal_coherence",
    "szego_kernel_regularized",
]


# --------------------------------------------------------------------------
# Kinematics
# --------------------------------------------------------------------------
def phase_grid(N: int) -> np.ndarray:
    """theta_k = 2 pi k / N, k = 0..N-1."""
    return 2.0 * np.pi * np.arange(N) / N


def pegg_barnett_basis(N: int) -> np.ndarray:
    """Matrix B with columns |theta_k> in the energy basis: B[n, k] = e^{i n theta_k}/sqrt(N)."""
    n = np.arange(N)[:, None]
    th = phase_grid(N)[None, :]
    return np.exp(1j * n * th) / np.sqrt(N)


def pw_constraint_state(N: int, levels, amplitudes) -> np.ndarray:
    """
    Exact solution of (N (x) 1 - 1 (x) H_S)|Psi> = 0 for a system whose k-th
    energy eigenvalue is the integer levels[k] (units of hbar*omega).

    Returns the coefficient matrix Psi[n, k] of |n>_C (x) |E_k>_S.
    The constraint has a nonzero solution ONLY because levels are integers in
    [0, N): this is the commensurability condition discussed in Sec. II.
    """
    levels = np.asarray(levels, dtype=int)
    amps = np.asarray(amplitudes, dtype=complex)
    if np.any(levels < 0) or np.any(levels >= N):
        raise ValueError("system levels must lie in the clock spectrum {0..N-1}")
    psi = np.zeros((N, len(levels)), dtype=complex)
    psi[levels, np.arange(len(levels))] = amps
    return psi / np.linalg.norm(psi)


# --------------------------------------------------------------------------
# Clock POVMs
# --------------------------------------------------------------------------
def sharp_phase_povm(N: int) -> np.ndarray:
    """Pegg--Barnett PVM: array E[k] = |theta_k><theta_k| (shape N x N x N)."""
    B = pegg_barnett_basis(N)
    return np.einsum("ik,jk->kij", B, B.conj())


def wrapped_gaussian_weights(N: int, sigma: float) -> np.ndarray:
    """Discrete wrapped Gaussian on the phase grid, centred at 0, summing to 1."""
    th = phase_grid(N)
    d = np.angle(np.exp(1j * th))  # map to (-pi, pi]
    if sigma <= 0:
        w = np.zeros(N)
        w[0] = 1.0
        return w
    # wrap a few times for accuracy when sigma is large
    w = np.zeros(N)
    for m in range(-3, 4):
        w += np.exp(-((d + 2 * np.pi * m) ** 2) / (2 * sigma**2))
    return w / w.sum()


def smeared_phase_povm(N: int, sigma: float) -> np.ndarray:
    """
    Unsharp but still shift-covariant clock POVM
        E_sigma(theta_k) = sum_l G_sigma(theta_k - theta_l) |theta_l><theta_l|.
    sum_k E_sigma(theta_k) = 1 exactly.  sigma is the READOUT phase resolution.
    """
    sharp = sharp_phase_povm(N)
    G = wrapped_gaussian_weights(N, sigma)
    idx = (np.arange(N)[:, None] - np.arange(N)[None, :]) % N  # k - l
    W = G[idx]  # W[k, l] = G(theta_k - theta_l)
    return np.einsum("kl,lij->kij", W, sharp)


def condition_on_clock(psi: np.ndarray, povm: np.ndarray):
    """
    For each clock POVM element E_k compute
        p_k     = <Psi| E_k (x) 1 |Psi>
        rho_S,k = Tr_C[(E_k (x) 1)|Psi><Psi|] / p_k
    using Tr_C[(E (x) 1)|Psi><Psi|]_{ss'} = sum_{n n'} E[n', n] Psi[n, s] Psi*[n', s'].
    Returns (p, rho) with p.shape = (K,), rho.shape = (K, d, d).
    """
    unnorm = np.einsum("ns,kmn,mt->kst", psi, povm, psi.conj())
    p = np.real(np.einsum("kss->k", unnorm))
    rho = unnorm / np.where(p > 1e-300, p, 1.0)[:, None, None]
    return p, rho


# --------------------------------------------------------------------------
# Finite-clock corrections
# --------------------------------------------------------------------------
def coherence_suppression_static(gap: float, sigma_theta: float) -> float:
    """Exact Gaussian readout smearing: |rho_jk| -> |rho_jk| exp(-gap^2 sigma^2 / 2)."""
    return float(np.exp(-0.5 * gap**2 * sigma_theta**2))


def diffusive_clock_channel(rho0: np.ndarray, energies, t: np.ndarray, T_C: float, n_quad: int = 801):
    """
    rho_obs(t) = int dtau G_{sigma(t)}(tau) U(t+tau) rho0 U(t+tau)^dag,
    with ASSUMED diffusive variance sigma(t)^2 = T_C * t (hbar = 1).
    Evaluated by direct numerical quadrature over tau (not the closed form).
    """
    E = np.asarray(energies, dtype=float)
    out = np.empty((len(t),) + rho0.shape, dtype=complex)
    for i, ti in enumerate(t):
        s = np.sqrt(max(T_C * ti, 0.0))
        if s == 0.0:
            taus, wts = np.array([0.0]), np.array([1.0])
        else:
            taus = np.linspace(-8 * s, 8 * s, n_quad)
            wts = np.exp(-(taus**2) / (2 * s**2))
            wts /= wts.sum()
        acc = np.zeros_like(rho0, dtype=complex)
        for tau, w in zip(taus, wts):
            ph = np.exp(-1j * E * (ti + tau))
            acc += w * (ph[:, None] * rho0 * ph.conj()[None, :])
        out[i] = acc
    return out


def lindblad_dephasing_solution(rho0: np.ndarray, energies, t: np.ndarray, T_C: float):
    """Closed form of d rho/dt = -i[H,rho] - (T_C/2)[H,[H,rho]] in the energy basis (hbar = 1)."""
    E = np.asarray(energies, dtype=float)
    dE = E[:, None] - E[None, :]
    return rho0[None] * np.exp(-1j * dE[None] * t[:, None, None] - 0.5 * T_C * dE[None] ** 2 * t[:, None, None])


def integrate_lindblad(rho0: np.ndarray, H: np.ndarray, t: np.ndarray, T_C: float):
    """Numerically integrate the GKSL equation with jump operator L = sqrt(T_C) H (hbar = 1)."""
    from scipy.integrate import solve_ivp

    d = H.shape[0]

    def rhs(_t, y):
        r = y.reshape(d, d)
        comm = H @ r - r @ H
        dcomm = H @ comm - comm @ H
        return (-1j * comm - 0.5 * T_C * dcomm).ravel()

    sol = solve_ivp(rhs, (t[0], t[-1]), rho0.astype(complex).ravel(), t_eval=t, rtol=1e-10, atol=1e-12)
    return sol.y.T.reshape(len(t), d, d)


def winding_marginal_coherence(nu: float, M: int) -> complex:
    """
    Factor multiplying rho_jk when a circle clock is read WITHOUT its winding
    number, averaged uniformly over M windings: (1/M) sum_{m<M} e^{-2 pi i m nu},
    nu = (E_j - E_k)/(hbar omega).  -> 1 if nu in Z, -> 0 (as 1/M) otherwise.
    """
    m = np.arange(M)
    return complex(np.mean(np.exp(-2j * np.pi * m * nu)))


def szego_kernel_regularized(x: np.ndarray, r: float) -> np.ndarray:
    """Abel-regularized Szego sum  sum_{n>=0} r^n e^{i n x} = 1/(1 - r e^{i x}),  0 < r < 1."""
    return 1.0 / (1.0 - r * np.exp(1j * x))
