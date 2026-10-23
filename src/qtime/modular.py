"""
Brute-force Tomita--Takesaki modular operator for a pure bipartite state.

For |Psi> in H_C (x) H_S (both C^d) with full Schmidt rank and the algebra
A = 1_C (x) B(H_S), the Tomita operator S: (1 (x) a)|Psi> -> (1 (x) a^dag)|Psi>
is antilinear and Delta = S^* S.  We construct Delta numerically from this
definition (no closed-form input), so the result can be used to check the
closed form  Delta = rho_C^{-1} (x) rho_S  claimed in the revised manuscript.
"""
from __future__ import annotations

import numpy as np

__all__ = ["tomita_delta", "reduced_states"]


def reduced_states(psi_mat: np.ndarray):
    """psi_mat[c, s] = coefficient of |c>_C|s>_S.  Returns (rho_C, rho_S)."""
    rho_C = psi_mat @ psi_mat.conj().T
    rho_S = psi_mat.T @ psi_mat.conj()
    return rho_C, rho_S


def tomita_delta(psi_mat: np.ndarray) -> np.ndarray:
    dC, dS = psi_mat.shape
    psi = psi_mat.reshape(-1)
    V, W = [], []
    for i in range(dS):
        for j in range(dS):
            Eij = np.zeros((dS, dS))
            Eij[i, j] = 1.0
            A = np.kron(np.eye(dC), Eij)
            V.append(A @ psi)
            W.append(A.conj().T @ psi)
    V = np.array(V).T  # columns v_ij = A_ij Psi
    W = np.array(W).T  # columns w_ij = A_ij^dag Psi
    if np.linalg.matrix_rank(V) < V.shape[1]:
        raise ValueError("Psi is not cyclic for 1 (x) B(H_S) (need dC >= dS and full Schmidt rank)")
    # <Sx, Sy> = <y, Delta x>  with x = V a, Sx = W conj(a)  =>  V^dag Delta V = (W^dag W)^T
    M = (W.conj().T @ W).T
    Vinv = np.linalg.pinv(V)
    return Vinv.conj().T @ M @ Vinv
