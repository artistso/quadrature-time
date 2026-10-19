# The finite relational clock series (v1 → v7)

A self-contained, deliberately narrow study of one question:

> **How well can a finite, equally spaced clock spectrum represent a system with
> arbitrary (incommensurate) energies, and for how long?**

This is the quantitative core of the manuscript's finite-clock section, kept separate
from the conjecture so that the numbers stand on their own. Every notebook states its
scope explicitly: these are controlled approximation results about *equally spaced clock
spectra*, not evidence that time is circular and not derivations of intrinsic decoherence.

Conventions: ħ = 1. Clock levels E_n = nΔ with 0 ≤ n ≤ N−1. System energies E_j paired
with integer labels q_j give approximate energies Ẽ_j = q_jΔ and mismatches
δ_j = E_j − Ẽ_j. With weights p_j = |c_j|² the conditioned-state overlap is
A(t) = Σ_j p_j e^{−iδ_j t} and the fidelity is F(t) = |A(t)|².

## Chain

| File | Adds | Status |
|---|---|---|
| `Finite_Dimensional_Relational_Clock_Model_v1.ipynb` | Explicit bipartite state, exact phase POVM, classical readout kernel. Separates readout smearing (a parameter) from physical decoherence (not derived). | complete |
| `Finite_Relational_Clock_Accuracy_Time_Horizon_v4.ipynb` | Grid search over spacing and N; instantaneous fidelity, energy-mismatch phase bound T·max_j|δ_j|, time-horizon trade-off. | complete (superseded by v5) |
| `Finite_Relational_Clock_Fidelity_Bound_v5.ipynb` | **The bound.** 1 − F(t) ≤ t² Var_p(δ), hence F(t) ≥ max{0, 1 − t² Var_p(δ)}. Invariant under a common energy offset; sufficient but not tight. | complete |
| v6 | Grid search refined (not in this repository) | *missing* |
| `Finite_Relational_Clock_Exact_Assignment_Optimization_v7.ipynb` | Replaces the grid: for fixed labels q the optimal spacing is analytic, Δ*(q) = C_{Eq}/V_q, with V*(q) = V_E − C²_{Eq}/V_q. Integer assignments enumerated exhaustively, q_0 anchored to 0. | complete |

### The two-line core

For a fixed assignment the objective is a weighted least-squares fit of energies
against integer labels, so

    Δ*(q) = C_{Eq} / V_q,        V*(q) = V_E − C_{Eq}² / V_q = V_E (1 − ρ²(E, q)).

The design problem is therefore: **choose the integer regressor q, with q_j distinct,
increasing, and in [0, N−1], that maximises the squared correlation with the spectrum.**
Global optimality is claimed only for that finite design problem (v7 states this).

Combining with v5 gives the operational statement:

    F(t) ≥ 1 − t² V*_N   ⟹   a target infidelity ε is guaranteed up to
                             T(ε) = sqrt(ε / V*_N).

So the entire clock-design question reduces to the decay of V*_N with N.

## Where it stands: the open question

`V*_N` is monotone non-increasing in N and decreases in **plateaus** (a larger clock is
only useful when a new label size enables a better simultaneous approximation of the
gaps). Brute-force enumeration in v7 is exact but combinatorial, so it stops at N ≈ 64.

Two things are not yet written down, and they are the natural v8.

**1. A universal upper bound with an explicit constant.** This is elementary but it is
missing from the series, and it makes the linear-in-N horizon unconditional. Rounding
each E_j/Δ to its nearest label gives |δ_j| ≤ Δ/2, so the range of the δ_j is at most Δ,
and Popoviciu's inequality gives Var_p(δ) ≤ Δ²/4. Fitting the spectrum's range
R = E_max − E_min into N levels forces Δ ≥ R/(N−1). Since Var_p(δ) is invariant under a
common energy offset (so the shift E_min → 0 is free),

    V*_N ≤ R² / (4(N−1)²)          and therefore       T(ε) = sqrt(ε/V*_N) ≥ 2(N−1)sqrt(ε)/R.

**Verified numerically: 0 violations for N = 4…128 over six spectra** (uniform weights,
`R` per spectrum), and the rounded
labels stay distinct in every case (distinctness is required by v7's problem statement).
So a clock with N levels always holds ε-fidelity for a time that grows at least linearly
in N.

**2. The true asymptotic of V*_N.** The exact optimum beats that bound by **2–5 orders of
magnitude** (725× at N = 16 for golden-ratio gaps; 6.6×10⁴ for a π-like spectrum),
because Δ is free and the labels need only be *proportional* to the energies, not to
round them. This is precisely a simultaneous Diophantine approximation problem: V* is
quadratic in the residual, so good rational approximants drive it far below N⁻², and the
worst case should be governed by badly approximable spectra with golden-ratio-type gap
ratios. Whether the worst-case constant is the golden ratio's — and whether the exact
optimum can be computed without enumeration (branch and bound, or a closest-vector
formulation) — is the interesting theorem.

## How to pick this up

```bash
python -m pytest                      # 24 tests, ~35 s, unrelated to this series
# then open notebooks/finite_clock/ in order: v1, v4, v5, v7
```

The exact-optimum routine is ~10 lines (`weighted_stats` + `solve_spacing_for_labels` +
`exhaustive_optimize` in v7). Any v8 should (i) state the universal bound above as a
proposition with proof, (ii) compute V*_N exactly without enumeration — a branch-and-bound
or a closest-vector/diophantine search — and (iii) plot V*_N·N² and V*_N·N⁴ against N for
golden, silver, π- and e-like spectra to separate the badly-approximable worst case from
the generic case.

## Honest caveats carried forward

- v7's optimum is global only for: fixed finite target spectrum, fixed weights, equally
  spaced clock, distinct ordered labels, free common offset.
- The fidelity bound is sufficient, not necessary; a clock can do better than the bound
  predicts at particular times because phase errors wrap.
- None of this says anything about whether time is fundamentally circular. It says how
  accurately a finite equally spaced spectrum can act as a clock, which the manuscript
  needs regardless of the conjecture's truth.
