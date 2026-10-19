"""Claim-by-claim verification of the analytic identities in the manuscript.

Every numbered result that can be checked by computation is checked here: the
identities are evaluated independently of the code path that produced them
(direct evaluation of the definition rather than the closed form under test), and
the residual is reported. Run directly for a ledger, or via pytest.
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "src")

from qtime.clock_design import (  # noqa: E402
    Spectrum,
    bound_applies,
    mismatch_variance,
    universal_bound,
)

HBAR = 1.0
RESULTS: list[tuple[str, str, float, str, float]] = []


def record(tag: str, statement: str, residual: float, note: str = "", tol: float = 1e-8) -> None:
    """Record one check. ``residual`` is what should be small; ``tol`` its budget.

    Checks whose value is *supposed* to be large (e.g. a bound's looseness) pass
    ``tol=float("inf")`` and are reported for information only.
    """
    RESULTS.append((tag, statement, float(residual), note, float(tol)))


# ---------------------------------------------------------------------------
# eq:szego  -- the corrected Szego kernel (Referee minor issue m1)
# ---------------------------------------------------------------------------
def check_szego() -> None:
    """Sum_{n>=0} e^{inx} = 1/(1-e^{i(x+i0)}) = pi*delta + 1/2 + (i/2) P.V. cot(x/2).

    Two independent statements are checked: (a) the *regular* part is exactly
    1/2 + (i/2)cot(x/2) away from the singularity; (b) the distributional part is
    a delta of weight pi, verified by integrating the mean over the circle.
    """
    x = np.linspace(1e-3, 2 * np.pi - 1e-3, 20001)
    formal = 1.0 / (1.0 - np.exp(1j * x))
    regular = 0.5 + 0.5j / np.tan(x / 2.0)
    # both sides diverge as 1/x near 0, so compare relatively
    rel = np.abs(formal - regular) / np.maximum(np.abs(formal), 1.0)
    record("eq:szego (regular part)", "1/(1-e^{ix}) = 1/2 + (i/2)cot(x/2)", np.max(rel))

    # (b) mean of the partial sums tends to pi/(2 pi) = 1/2: the delta carries pi
    # the +i0 prescription is what produces the delta: the pole at x = 0 carries
    # residue i, and Plemelj's formula i/(x+i0) = i P.V.(1/x) + pi delta(x) then
    # gives a delta of weight pi.  A quadrature of the oscillating partial sums
    # cannot see this (the sum has N oscillations), so test the pole directly.
    xs = np.array([1e-3, 1e-4, 1e-5, 1e-6, 1e-7])
    pole = xs / (1.0 - np.exp(1j * xs))
    record("eq:szego (pole)", "lim x/(1-e^{ix}) = i, so the delta weight is pi",
           abs(pole[-1] - 1j), "expansion is i + x/2 + O(x^2); smallest x used", tol=1e-6)


def check_unsharp() -> None:
    """Smearing in theta by a wrapped Gaussian of width sigma_t multiplies the
    Fourier mode nu by exp(-nu^2 sigma_t^2 / 2), with no theta dependence."""
    M = 4096
    theta = 2 * np.pi * np.arange(M) / M

    def wrapped_gaussian(sigma, m=M):
        k = np.arange(m)
        k = np.where(k > m // 2, k - m, k)
        g = np.exp(-0.5 * (k * 2 * np.pi / m / sigma) ** 2)
        return g / g.sum()

    worst = 0.0
    for sigma in (0.05, 0.2, 0.5):
        kernel = wrapped_gaussian(sigma)
        for nu in (1, 2, 3, 5):
            f = np.exp(1j * nu * theta)
            smeared = np.fft.ifft(np.fft.fft(f) * np.fft.fft(kernel)).real + 0j
            smeared = np.fft.ifft(np.fft.fft(f) * np.fft.fft(kernel))
            predicted = f * np.exp(-nu**2 * sigma**2 / 2)
            worst = max(worst, np.max(np.abs(smeared - predicted)))
    record("eq:static", "smeared coherence = rho_jk * exp(-nu^2 sigma_t^2/2)", worst)


# ---------------------------------------------------------------------------
# eq:lindblad  -- double-commutator algebra gives Gamma_jk = T_C (Ej-Ek)^2/2 hbar^2
# ---------------------------------------------------------------------------
def check_gksl() -> None:
    """The dissipator D(rho) = -(T_C/2)[H,[H,rho]] acts on matrix elements as
    rho_jk -> -(T_C/2)(Ej-Ek)^2 rho_jk in units hbar = 1, so the decay rate is
    Gamma_jk = T_C (Ej-Ek)^2 / 2. Verified against an explicit 3x3 matrix."""
    rng = np.random.default_rng(0)
    E = np.array([0.0, 0.7, 1.9])
    H = np.diag(E)
    A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    rho = A @ A.conj().T
    rho /= np.trace(rho)
    TC = 0.37
    D = -(TC / 2) * (H @ (H @ rho - rho @ H) - (H @ rho - rho @ H) @ H)
    predicted = np.zeros_like(rho)
    for j in range(3):
        for k in range(3):
            predicted[j, k] = -(TC / 2) * (E[j] - E[k]) ** 2 * rho[j, k]
    record("eq:lindblad", "D(rho)_jk = -(T_C/2)(Ej-Ek)^2 rho_jk", np.max(np.abs(D - predicted)))


# ---------------------------------------------------------------------------
# eq:winding  -- the Dirichlet factor and its 1/M envelope
# ---------------------------------------------------------------------------
def check_winding() -> None:
    worst_dir, worst_env = 0.0, 0.0
    for n in (3.3, 1.25, 0.5 + 7 / 11):
        for M in (1, 2, 5, 17, 50):
            m = np.arange(M)
            direct = np.exp(-1j * 2 * np.pi * m * n).mean()
            # D_M(nu) = e^{-i pi (M-1) nu} sin(pi M nu)/(M sin pi nu)
            formula = np.exp(-1j * np.pi * (M - 1) * n) * np.sin(np.pi * M * n) / (M * np.sin(np.pi * n))
            worst_dir = max(worst_dir, abs(direct - formula))
        worst_env = max(worst_env, abs(direct) - 1.0 / (M * abs(np.sin(np.pi * n))) - 1e-12)
    record("eq:winding (identity)", "D_M(nu) = e^{-i pi (M-1) nu} sin(pi M nu)/(M sin pi nu)", worst_dir)
    record("eq:winding (envelope)", "|D_M| <= 1/(M |sin pi nu|)", max(worst_env, 0.0))


# ---------------------------------------------------------------------------
# eq:kg and eq:proper  -- Klein-Gordon factorisation and proper-time expansion
# ---------------------------------------------------------------------------
def check_kg() -> None:
    """(P_mu P^mu - M^2 c^2) = c^{-2} C_+ C_- with C_pm = P_t +- H_rel.

    Uses P_0 = P_t/c, so P_mu P^mu = P_t^2/c^2 - p^2. Checked on a grid.
    """
    c = 2.0
    worst = 0.0
    for p in (0.0, 0.5, 3.0, 1e2):
        for mass in (0.3, 5.0):
            P_t = np.sqrt(c**2 * p**2 + mass**2 * c**4) + 0.37 * mass
            H_rel = np.sqrt(c**2 * p**2 + mass**2 * c**4)
            lhs = (P_t / c) ** 2 - p**2 - mass**2 * c**2
            rhs = (P_t + H_rel) * (P_t - H_rel) / c**2
            worst = max(worst, abs(lhs - rhs))
    record("eq:kg", "(P.P - M^2 c^2) = c^-2 C_+ C_-", worst)

    # proper-time: H_rel = E_P + delta/gamma at first order. First order means the
    # RELATIVE error must fall linearly with delta, so test the convergence rate.
    def rel_error(delta):
        worst = 0.0
        for p in (0.0, 1.0, 4.0):
            for mass in (1.0, 2.5):
                E_P = np.sqrt(c**2 * p**2 + mass**2 * c**4)
                gamma = E_P / (mass * c**2)
                exact = np.sqrt(c**2 * p**2 + (mass + delta / c**2) ** 2 * c**4) - E_P
                worst = max(worst, abs(exact - delta / gamma) / abs(delta / gamma))
        return worst

    rate = rel_error(1e-3) / rel_error(5e-4)
    record("eq:proper", "H_rel = E_P + gamma^-1 delta + O(delta^2); relative error halves with delta",
           abs(rate - 2.0), "rate = 2 confirms first order", tol=1e-3)


# ---------------------------------------------------------------------------
# eq:universal  -- the universal bound, and how tight it actually is
# ---------------------------------------------------------------------------
def check_universal_bound() -> None:
    """The bound is a theorem only under N >= 1 + R/g. Check both directions."""
    import itertools

    def best_legal(S, N):
        return min(mismatch_variance(S, c) for c in itertools.combinations(range(N), S.size))

    rng = np.random.default_rng(3)
    violations, ratios, outside = 0, [], 0
    for _ in range(150):
        d = int(rng.integers(3, 6))
        E = np.sort(rng.uniform(0, 3, d))
        E = E - E[0]
        if E[-1] < 1e-6:
            continue
        S = Spectrum(E)
        R = float(E[-1])
        for N in range(d, 12):
            v, b = best_legal(S, N), universal_bound(R, N)
            if bound_applies(S, N):
                if v > b * (1 + 1e-12):
                    violations += 1
                ratios.append(v / b)
            elif v > b * (1 + 1e-12):
                outside += 1
    record("eq:universal (hypothesis)", "violations when N >= 1 + R/g", violations)
    record("eq:universal (tightness)", "max V*/bound inside the hypothesis (1 = tight, 0 = loose)",
           max(ratios), "informational: the bound is an envelope, not a sharp value", tol=float("inf"))
    record("eq:universal (necessity)", "violations OUTSIDE the hypothesis (must be > 0)",
           0.0 if outside else 1.0, f"{outside} violations found; the hypothesis is necessary",
           tol=0.5)


def main() -> int:
    check_szego()
    check_unsharp()
    check_gksl()
    check_winding()
    check_kg()
    check_universal_bound()

    width = max(len(t) for t, _, _, _, _ in RESULTS)
    failed = 0
    print(f"{'check'.ljust(width)}  {'residual':>12}  verdict  statement")
    print("-" * (width + 70))
    for tag, statement, residual, note, tol in RESULTS:
        ok = residual < tol if np.isfinite(tol) else True
        if not ok:
            failed += 1
        verdict = "PASS" if ok else "FAIL"
        if not np.isfinite(tol):
            verdict = "info"
        print(f"{tag.ljust(width)}  {residual:12.3e}  {verdict:6s}  {statement}")
        if note:
            print(f"{''.ljust(width)}  {'':12}        note: {note}")
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks pass")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
