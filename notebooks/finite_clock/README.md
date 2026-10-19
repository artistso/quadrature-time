# The finite relational clock series (v1 → v11)

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
| `Finite_Relational_Clock_Asymptotics_v8.ipynb` | **The scaling law.** Universal bound, three-level closed form, N⁻⁴ decay with exact constant, the Diophantine reformulation, and an O(N) optimiser that reaches N = 10⁶. | complete |
| `Finite_Relational_Clock_Exact_Arithmetic_v9.ipynb` | **Exact arithmetic, higher dimensions, and a correction.** Commensurability becomes an exact statement; a general optimiser for d ≥ 4 replaces enumeration; and the exponent turns out to depend on the *spectrum*, not the dimension. | **current** |

Code: `src/qtime/clock_design.py`; tests: `tests/test_clock_design.py` (128 checks, in
addition to the 24 manuscript tests). Propositions 1–2 and the N⁻⁴ law are now Sec. V.C of
the manuscript.

## The results

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

### 4. Exact arithmetic, and d ≥ 4 (v9)

**Exact.** The general objective is a covariance cancellation ($V_E-C^2/V_q$: two numbers of
order $V_E$), so it returns noise once the variance drops below ~10⁻¹⁴$V_E$. Rational
arithmetic fixes that: commensurate spectra give variance **exactly zero** (2.5 is exactly 5/2),
and the golden-ratio constant is confirmed to 3 × 10⁻⁸. What still limits precision beyond
p ~ 10⁵ is the truncation of *r* itself, not the arithmetic — worth being clear about.

**d ≥ 4.** Enumeration is $\binom{N-1}{d-1}$: hopeless. But minimising over spacing *and* offset
is a two-parameter weighted least-squares fit, so for fixed $(\Delta,a)$ the labels decouple:
$q_j=\mathrm{round}((E_j-a)/\Delta)$. `search_optimum` walks that design plane, iterates
fit-then-round to a fixed point, then polishes by coordinate and window descent. It matches
exhaustive enumeration on 300 randomised instances (d = 3 to N = 120, d = 4 to N = 36) and runs
d = 4, N = 10⁵ in 0.1 s. It is a heuristic, and it has a documented large-N limit: its seed grid
can miss convergent-denominator assignments, which is why the proportional family
(`proportional_curve`, `proportional_curve_stable`) is the tool for spectra whose optimum is
proportional.

**The correction to v8.** The exponent is *not* a property of the clock dimension:

| spectrum | exponent | |
|---|---|---|
| (0, 1, φ), d = 3 | −4.03 | one badly approximable ratio |
| (0, 1, φ, φ²), d = 4 | −4.05 | still one ratio: labels ∝ (1, φ, φ²) |
| (0, 1, φ, φ², φ³), d = 5 | −3.99 | still one ratio |
| (0, 1, √2, √3), d = 4 | −3.08 | two independent ratios: simultaneous approximation |
| (0, 1, 2, 3), d = 4 | exact 0 | commensurate: no exponent at all |

So N⁻⁴ is the **one-dimensional** Diophantine rate — the rate whenever the spectrum is generated
by powers of a single quadratic irrational — and in general the rate is set by simultaneous
approximation in $d-2$ dimensions, which continued fractions do not describe. Prop. 1 is
untouched by this: it never used Diophantine structure.

### 5. What tuning the spacing buys (v11)

Prop. 4's $N^{-4}$ law assumes the quantum $\Delta$ is a free design parameter. Real clocks do not
offer a continuum: a cavity mode, an atomic line, a Zeeman splitting are fixed by nature. Two
consequences, both exact:

**Admissibility.** With smallest gap $g$: $\Delta\le g$ suffices for distinct labels, $\Delta\le
R/(d-1)$ is necessary, and $N\ge1+R/\Delta$ levels are required. A quantum above the smallest gap
cannot in general resolve the levels *at all*, however many levels are available.

**What tuning buys.** A pinned quantum is forced to the rounding assignment, so it is confined to
the universal $N^{-2}$ variance and a horizon growing only like $N$. The tuned clock reaches $N^{-4}$
and beats it by a factor growing like $N^2$.

| $N$ | universal bound | pinned quantum | tuned | gain |
|---|---|---|---|---|
| 30 | 7.8e-04 | 4.9e-05 | 5.9e-07 | 83 |
| 300 | 7.3e-06 | 1.1e-06 | 3.9e-11 | 2.8e+04 |
| 3000 | 7.3e-08 | 1.6e-08 | 2.5e-15 | 6.5e+06 |

Measured slopes: universal −2.03, pinned −1.93, tuned −4.13.

The pinned curve is not smooth — at $N=90$ it dips by an order of magnitude because the levels
happen to land almost on integers. That is the same Diophantine luck that governs the whole series,
and it is why the gain is asserted as a *trend* ($N^2$) rather than as monotonicity.

Code: `src/qtime/clock_scales.py`, including a `PHYSICAL_CLOCKS` table (optical, telecom,
21 cm, caesium hyperfine, Zeeman 1 T, 6 GHz cavity, 5 GHz transmon) so the test can be run against
real hardware rather than invented numbers.

## Where it stands

`V*_N` is monotone non-increasing in N and decreases in **plateaus**: a larger clock only
helps when it crosses a convergent, which is why the curve is a staircase and why the
right-hand panel of `figures/fig5_clock_scaling.png` is flat exactly at the convergents.

Computed: three-level spectra in O(N) (`scaling_curve`, N = 10⁶ in 0.6 s), the proportional
family to M = 10⁶ in 0.2 s (`proportional_curve_stable`), exact arithmetic anywhere precision
matters, and d ≥ 4 at N = 10⁵ in 0.1 s. The open end is no longer cost — it is the theorems
listed below.

## Open questions — the v10 candidates

1. **Simultaneous approximation, properly.** $(1,\sqrt2,\sqrt3)$ was measured inside the
   proportional family only. Whether another label family beats it is open, and it is a
   lattice-reduction question (LLL/BKZ). The right statement to aim for: the decay exponent in
   $d$ dimensions is $-(d-2)\cdot 2 - 2$ for "generic independent" ratios, and $-4$ for spectra
   in the rank-one (real-quadratic-field) case. Both ends are now measured; the theorem is not.
2. **Prove the three-level optimiser.** The ±2 window around $\mathrm{nint}(rq)$ reproduces
   enumeration for every spectrum tested up to N = 120, and the nearest integer alone is not
   always optimal (r = φ, q = 4). No proof.
3. **Prove the general optimiser finds the optimum** for d ≥ 4. Currently validated only where
   enumeration is affordable.
4. **Exact algebraic numbers.** Exact rational arithmetic removes round-off but not the
   truncation of an irrational *r*; a genuinely high-precision claim needs $\mathbb{Q}(\sqrt5)$
   rather than a long decimal. Then the constants at p = 10⁹ become checkable.
5. **Non-uniform weights in closed form.** With weights the problem stays a weighted
   least-squares fit but no longer collapses to a two-vector determinant; a determinant
   (Gram-matrix) formula for general d and weights may well exist.

## How to pick this up

```bash
python -m pytest tests/test_clock_design.py -q     # 128 checks, ~5 s
# then, in order: v1, v4, v5, v7, v8, v9
```

## Honest caveats carried forward

- v7's optimum is global only for: fixed finite target spectrum, fixed weights, equally
  spaced clock, distinct ordered labels, free common offset.
- The fidelity bound is sufficient, not necessary; a clock can do better at particular
  times because phase errors wrap.
- None of this says anything about whether time is fundamentally circular. It says how
  accurately a finite equally spaced spectrum can act as a clock, which the manuscript
  needs regardless of the conjecture's truth.
