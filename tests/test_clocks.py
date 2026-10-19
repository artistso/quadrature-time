import numpy as np
import pytest

from qtime import clocks as qc


@pytest.mark.parametrize("N", [8, 33, 64])
def test_pegg_barnett_is_orthonormal_and_complete(N):
    B = qc.pegg_barnett_basis(N)
    assert np.allclose(B.conj().T @ B, np.eye(N))
    assert np.allclose(qc.sharp_phase_povm(N).sum(axis=0), np.eye(N))


@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.3, 1.0])
def test_smeared_povm_is_a_covariant_povm(sigma):
    N = 48
    E = qc.smeared_phase_povm(N, sigma)
    assert np.allclose(E.sum(axis=0), np.eye(N))  # completeness
    for k in (0, 7, 31):  # positivity
        assert np.linalg.eigvalsh(E[k]).min() > -1e-12
    # covariance: U(phi) E(theta_k) U^dag = E(theta_k + phi) with phi = 2 pi / N, U = e^{i phi N}
    U = np.diag(np.exp(1j * 2 * np.pi / N * np.arange(N)))
    for k in (0, 5):
        assert np.allclose(U @ E[k] @ U.conj().T, E[(k + 1) % N])


def test_exact_constraint_state_gives_uniform_marginal_and_exact_unitary_evolution():
    N, levels = 64, [3, 10]  # integer gap r = 7 (commensurate by construction)
    psi = qc.pw_constraint_state(N, levels, [1, 1])
    p, rho = qc.condition_on_clock(psi, qc.sharp_phase_povm(N))
    assert np.allclose(p, 1.0 / N)  # Theorem 1(i): uniform raw marginal
    theta = qc.phase_grid(N)
    for k in range(N):  # Theorem 1(ii): pure, Schroedinger-evolved conditional state
        v = np.exp(-1j * np.array(levels) * theta[k]) / np.sqrt(2)
        assert np.allclose(rho[k], np.outer(v, v.conj()))


def test_constraint_has_no_solution_off_the_clock_spectrum():
    with pytest.raises(ValueError):
        qc.pw_constraint_state(16, [0, 20], [1, 1])


@pytest.mark.parametrize("sigma", [0.05, 0.1, 0.2])
def test_unsharp_readout_gives_static_not_growing_suppression(sigma):
    N, levels = 256, [0, 4]
    psi = qc.pw_constraint_state(N, levels, [1, 1])
    _, rho = qc.condition_on_clock(psi, qc.smeared_phase_povm(N, sigma))
    coh = np.abs(rho[:, 0, 1])
    # constant over the whole circle (no decay with theta) ...
    assert np.ptp(coh) < 1e-12
    # ... and equal to the Gaussian characteristic function
    assert coh[0] == pytest.approx(0.5 * qc.coherence_suppression_static(4, sigma), rel=1e-6)


def test_diffusive_assumption_reproduces_lindblad_exactly():
    E = np.array([0.0, 1.3, 2.1])
    v = np.ones(3) / np.sqrt(3)
    rho0 = np.outer(v, v)
    t = np.linspace(0, 6, 25)
    T_C = 0.2
    quad = qc.diffusive_clock_channel(rho0, E, t, T_C)
    closed = qc.lindblad_dephasing_solution(rho0, E, t, T_C)
    ode = qc.integrate_lindblad(rho0, np.diag(E), t, T_C)
    assert np.allclose(quad, closed, atol=1e-8)
    assert np.allclose(ode, closed, atol=1e-7)


def test_winding_marginalization_kills_incommensurate_coherence():
    assert abs(qc.winding_marginal_coherence(3.0, 1000)) == pytest.approx(1.0)
    nu = (np.sqrt(5) - 1) / 2
    M = 1000
    bound = 1.0 / (M * abs(np.sin(np.pi * nu)))
    assert abs(qc.winding_marginal_coherence(nu, M)) <= bound + 1e-12


def test_szego_kernel_corrected_form():
    # pointwise (x != 0): 1/(1 - e^{ix}) = 1/2 + (i/2) cot(x/2), no extra delta term off x = 0
    x = np.array([0.3, 1.0, 2.5, -1.7])
    lim = qc.szego_kernel_regularized(x, 1 - 1e-9)
    assert np.allclose(lim, 0.5 + 0.5j / np.tan(x / 2), rtol=1e-6)
    # distributional: Re part -> 1/2 + pi*delta(x); test against f(x) = exp(cos x)
    xs = np.linspace(-np.pi, np.pi, 400001)
    f = np.exp(np.cos(xs))
    r = 1 - 1e-3
    val = np.trapezoid(np.real(qc.szego_kernel_regularized(xs, r)) * f, xs)
    expected = np.pi * np.exp(1.0) + 0.5 * np.trapezoid(f, xs)
    assert val == pytest.approx(expected, rel=2e-3)
    # the ORIGINAL formula pi*delta + 1/(1-e^{i(x+i0)}) would give 2*pi*f(0) + ... : off by pi f(0)
    wrong = 2 * np.pi * np.exp(1.0) + 0.5 * np.trapezoid(f, xs)
    assert abs(val - wrong) > 5.0
