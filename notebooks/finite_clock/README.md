# The finite relational clock series (v1 → v8)

A self-contained, deliberately narrow study of one question:

> **How well can a finite, equally spaced clock spectrum represent a system with
> arbitrary (incommensurate) energies, and for how long?**

This is the quantitative core of the manuscript's finite-clock section, kept separate
from the conjecture so that the numbers stand on their own. Every notebook states its
scope explicitly: these are controlled approximation results about *equally spaced clock
spectra*, not evidence that time is circular and not derivations of intrinsic
decoherence.

Conventions: ħ = 1. Clock levels E_n = nΔ with 0 ≤ n ≤ N−1. System energies E_j paired
with integer labels q_j give approximate energies Ẽ_j = q_jΔ and mismatches
δ_j = E_j − Ẽ_j. With weights p_j = |c_j|² the conditioned-state overlap is
A(t) = Σ_j p_j e^{−iδ_j t} and the fidelity is F(t) = |A(t)|².

## Chain

| File | Adds | Status |
|---|---|---|
| `Finite_Dimensional_Relational_Clock_Model_v1.ipynb` | Explicit bipartite state, exact phase POVM, classical readout kernel. Separates readout smearing (a parameter) from physical decoherence (not derived). | complete |
| `Finite_Relational_Clock_Accuracy_Time_Horizon_v4.ipynb` | Grid search over spacing and N; instantaneous fidelity, energy-mismatch phase bound T·max_j\|δ_j\|, time-horizon trade-off. | complete (superseded by v5) |
| `Finite_Relational_Clock_Fidelity_Bound_v5.ipynb` | **The bound.** 1 − F(t) ≤ t² Var_p(δ), hence F(t) ≥ max{0, 1 − t² Var_p(δ)}. Invariant under a common energy offset; sufficient but not tight. | complete |
| v6 | Grid search refined (not in this repository) | *missing* |
| `Finite_Relational_Clock_Exact_Assignment_Optimization_v7.ipynb` | Replaces the grid: for fixed labels q the optimal spacing is analytic, Δ\*(q) = C_{Eq}/V_q, with V\*(q) = V_E − C²_{Eq}/V_q. Integer assignments enumerated exhaustively, q_0 anchored to 0. | complete |
| `Finite_Relational_Clock_Asymptotics_v8.ipynb` | **The scaling law.** Universal bound, three-level closed form, N⁻⁴ decay with exact constant, the Diophantine reformulation, and an O(N) optimiser that reaches N = 10⁶. | **current** |

Code: `src/qtime/clock_design.py`; tests: `tests/test_clock_design.py` (72 checks, in
addition to the 24 manuscript tests).

## The three results

### 1. The design problem is a weighted least-squares fit

For a fixed assignment q the objective is a weighted fit of energies against integer
labels, so the spacing and the residual are analytic:

    Δ*(q) = C_{Eq} / V_q,        V*(q) = V_E − C_{Eq}² / V_q = V_E (1 − ρ²(E, q)).

**Design = choose the integer regressor q, distinct and increasing in [0, N−1], that
maximises its squared weighted correlation with the spectrum.** Global optimality is
claimed only for this finite design problem (v7 says so explicitly).

### 2. A universal, spectrum-free bound (Proposition 1)

Round each shifted E_j to its nearest label: then |δ_j| ≤ Δ/2, the δ_j span at most Δ,
and Popoviciu's inequality gives Var_p(δ) ≤ Δ²/4. Fitting the range R = E_max − E_min
into N levels forces Δ ≥ R/(N−1). Hence

    V*_N ≤ R² / (4(N−1)²)        and        T(ε) ≥ 2(N−1)√ε / R.

Unconditional: no number theory, any spectrum, any weighting. Verified with 0 violations
for N = 4…128 over five spectra, where the rounding construction never exceeds 0.222 of
the bound. It is loose by two to five orders of magnitude — that gap is the story.

### 3. The three-level closed form, and N⁻⁴ (Proposition 2)

For E = (0, 1, r) with uniform weights and labels (0, q, p),

    V*(q, p) = (p − rq)² / [(p+q)² + (2q−p)² + (2p−q)²].

The numerator is the *rational-approximation error* of the gap ratio; the denominator is
a label-geometry factor equal to 27 V_E up to O(q) corrections. Verified against the
density-matrix computation to 10⁻⁹. Two consequences:

* **Clock design is simultaneous Diophantine approximation.** Because V\* is *quadratic*
  in the error, and |r − p/q| = κ/q² along convergents, the optimum decays as

      V*_N  ≃  κ(r)² r⁴ / (27 V_E) · N⁻⁴,

  four powers of N, not two. The guaranteed horizon therefore grows **quadratically**
  in the number of clock levels: T ≈ N² √(27 V_E ε / r⁴) / κ.

* **The golden ratio is the extremal spectrum.** At Fibonacci labels the constant is
  exactly φ⁴/(135 V_E) = 0.1142350…, matched to 3 × 10⁻⁸ by p = 4181. Hurwitz's theorem
  (κ ≤ 1/√5, equality for φ and its continued-fraction relatives) means no three-level
  spectrum can decay more slowly. Measured exponents at convergent labels between
  10³ and 10⁶: φ, 1+√2, √2 give −4.000; e−1 gives −4.06; π−1 gives −3.70 because its
  convergent 355/113 is exceptional.

## Where it stands

`V*_N` is monotone non-increasing in N and decreases in **plateaus**: a larger clock only
helps when it crosses a convergent, which is why the curve is a staircase and why the
right-hand panel of `figures/fig5_clock_scaling.png` is flat exactly at the convergents.

Done: the closed form makes three-level spectra O(N); `scaling_curve` reaches N = 200 000
in 0.1 s and N = 10⁶ in 0.6 s, against ~10¹¹ assignments for enumeration.

## Open questions — the v9 candidates

1. **d ≥ 4, and non-uniform weights.** The closed form is three-level and uniform-weight
   only. With weights the optimum is a *weighted* least-squares problem, and V_E(1 − ρ²)
   no longer collapses to a determinant of two vectors. The general statement is a
   lattice/Diophantine problem in dimension d−2, and nothing in the series addresses it.
2. **Prove the optimiser.** The ±2 window around nint(rq) reproduces exhaustive
   enumeration for every spectrum tested up to N = 120, and the nearest integer alone is
   *not* always optimal (r = φ, q = 4). No proof of the window, and no proof that the
   optimum is attained at a convergent label rather than merely near one.
3. **The full extremal statement.** Hurwitz extremises κ alone; the constant also contains
   r⁴/V_E(r), and that factor is not monotone in κ. "The golden ratio is the hardest
   three-level spectrum" is therefore supported numerically, not proved.
4. **Precision, not physics, at large N.** For p ≳ 2000 the residual (p − rq)² falls below
   the double-precision floor and the apparent spikes are round-off. A v9 that wants
   N = 10⁹ needs exact integer arithmetic (`fractions.Fraction` or sympy) rather than
   floats — cheap to do, and it would make the asymptotic constants measurable to many
   more digits.
5. **Feedback into the manuscript.** The manuscript's finite-clock section could carry
   Prop. 1 and Prop. 2 as a short subsection: the first is a two-line proof, the second
   is a clean, checkable, genuinely new statement, and together they say something the
   paper currently cannot — that a relational clock's accuracy is governed by how well
   its gap ratios approximate rationals, not by how finely it is spaced.

## How to pick this up

```bash
python -m pytest tests/test_clock_design.py -q     # 72 checks, < 1 s
# then, in order: v1, v4, v5, v7, v8
```

## Honest caveats carried forward

- v7's optimum is global only for: fixed finite target spectrum, fixed weights, equally
  spaced clock, distinct ordered labels, free common offset.
- The fidelity bound is sufficient, not necessary; a clock can do better at particular
  times because phase errors wrap.
- None of this says anything about whether time is fundamentally circular. It says how
  accurately a finite equally spaced spectrum can act as a clock, which the manuscript
  needs regardless of the conjecture's truth.
