# Adversarial audit of the new material (v8–v11), 2026-10-18

*Working note. Not paper content.*

The referee report in `REVIEW.md` audited the original draft. It has nothing to say about the
material added since — Sec. II's staircase, Sec. V.C's finite-clock design, the v10 drift result or
v11's admissibility lemma. A referee would audit exactly that, and it is where the errors would be,
because it is the part nobody has read. This note is that audit: what I checked, what was wrong,
what I changed.

Three of the findings were real defects. One of them was a false statement in the manuscript that
had been sitting there since v9.

---

## 1. The admissibility lemma was false as written (fixed)

**Lemma `lem:admissible`** claims: if the clock quantum satisfies `Δ ≤ g` (g = smallest gap), then
rounding every level to the nearest integer gives distinct labels.

I proved it in one line — "consecutive labels differ by ≥ g/Δ ≥ 1 in the rounding, hence strictly
increase" — and that line is wrong. It assumes a monotone rounding convention. `numpy.rint`, which
the code used everywhere, rounds ties to *even*, and round-to-even is not monotone:

| spectrum | Δ | round-to-even | round-half-up |
|---|---|---|---|
| (0, 1.5, 2.5) | 1.0 = g | **(0, 2, 2)** — collision | (0, 2, 3) — distinct |
| (0, 1.5, 2.5, 3.5) | 1.0 = g | **(0, 2, 2, 4)** — collision | (0, 2, 3, 4) |

So the lemma failed at exactly its boundary, on the lattice spectra a referee would try first. A
random scan of 20 000 spectra produced zero failures, which is why no test caught it: the tie set has
measure zero, and hand-constructed counterexamples are the only way to see it.

**Fix.** The lemma now states the convention it needs, `q_j = ⌊(E_j−E_1)/Δ + ½⌋`, and the proof is
correct for it:

> if `u_j = (E_j−E_1)/Δ` then `u_{j+1} − u_j ≥ g/Δ ≥ 1`, and since `q_j ≤ u_j + ½`,
> `u_{j+1} + ½ ≥ q_j + 1`, so `q_{j+1} ≥ q_j + 1`.

The code gained `round_half_up` and every label construction now uses it. Three regression tests pin
the counterexample. The paper also says, in the proof, *why* the convention matters — that is the
honest version, and it costs one sentence.

## 2. "Solved to N = 10⁶ in under a second" was false (fixed)

The scope remark in Sec. V.C claimed the design problem is solvable to `N = 10⁶` in under a second.
Measured, `three_level_optimum` is an `O(N)` scan:

| N | time |
|---|---|
| 10³ | 0.038 s |
| 10⁴ | 0.365 s |
| 10⁵ | 3.84 s |
| 10⁶ | 37.6 s |

Off by a factor of forty. The claim is now the measurement, `0.4 s at N = 10⁴`.

I tried to make the boast true instead of deleting it, by exploiting the fact that the optimum sits at
Fibonacci convergents — for φ the labels at `N = 10⁶` are `(0, 4181, 6765)`. That gives the right
answer for φ in 0.4 ms, but it is **not** a general algorithm:

- for √2 at `N = 11` the optimum is `(0, 7, 10)`, whose ratio `10/7` is the **mediant** of the
  convergents `7/5` and `3/2`;
- extending the candidate set to convergents *and* semiconvergents still misses the optimum for π at
  small `N`.

So the linear scan stays, and the note records why: the optimum can sit at an intermediate fraction,
which is a sharpening of the Proposition 4 picture rather than a defect in it. (The Diophantine
structure that makes the spectrum easy is the same structure that makes the search skip.)

## 3. The golden-ratio constant was "confirmed numerically", not proved (upgraded to a theorem)

The remark claimed the constant `φ⁴/(135 V_E) = φ⁴/60 = 0.1142350327708281` was "confirmed to
3×10⁻⁸ in exact rational arithmetic". Both halves of that are unsatisfying, and a referee would say
so:

- **Floating point fails** past `p ≈ 10³`. The covariance form subtracts two numbers of order `V_E`,
  so at `p = 4181` the computed `p⁴V*` is wrong by 1e-3 — hopeless, when the claim is 1e-8.
- **Rational arithmetic fails differently, and worse.** `φ` is irrational, so the exact computation
  needs `φ` truncated to a rational convergent; the truncation error then *dominates*, and it grows
  with `p`. The existing test documents this: relative error 1e-6 at `p = 4181` and 1e-4 at
  `p = 75025`. The method cannot reach 1e-8 at large labels at all.

**Fix.** Work where `φ` lives: `ℚ(√5)`. The residual is `p − qφ = (p − q/2) − (q/2)√5`, its square is
again in `ℚ(√5)`, and the closed-form denominator is an integer, so `V*` is an exact pair of
`Fraction`s. Nothing is truncated and nothing cancels.

Better, the limit is now *derived* rather than observed. For `q = F_n`, `p = F_{n+1}`:

- `F_{n+1} − φF_n = (−1)^n φ^{−n}` (so the residual squared is `φ^{−2n}`);
- the closed-form denominator `(p+q)² + (2q−p)² + (2p−q)²` equals `F_{n+2}² + F_{n−2}² + L_n²`
  (checked symbolically for `n = 3…29`), which is `∼ (12/5)φ^{2n}`;
- `p⁴ = F_{n+1}⁴ ∼ φ^{4n+4}/25`.

Therefore `p⁴V* = (φ⁴/60)[1 + O(φ^{−2n})]`. The error falls by exactly `φ²` per Fibonacci step,
which is what the exact computation shows:

| p | \|p⁴V* − φ⁴/60\| |
|---|---|
| 5 | 2.5e-03 |
| 34 | 5.3e-05 |
| 233 | 1.1e-06 |
| 1597 | 2.4e-08 |
| 4181 | 3.5e-09 |
| 514229 | 2.3e-13 |

Function: `clock_design.golden_limit_magnitude` (exact pair in, arbitrary-precision modulus out).
Five tests, including one that asserts the float form *does* fail where the exact form succeeds, so
the reason for the extra machinery cannot be quietly lost later.

---

## What I checked and found sound

- **The universal bound** `V*_N ≤ R²/(4(N−1)²)` for `N = 4…128` over five spectra: zero violations.
  It is also provable in one line (every rounding residual lies in an interval of width `Δ`, and a
  distribution on such an interval has variance at most `Δ²/4`), which is now how the v11 remark
  puts it — so the pinned-clock ceiling is a theorem, not a trend.
- **`three_level_optimum`'s heuristic window** (`±2` labels): the optimum is identical at window 2,
  8 and 32 for `N` up to 3000 on both φ and √3, and sits on consecutive Fibonacci convergents. The
  docstring's "validated heuristic" is the right label; the figure's claim survives.
- **Optimiser vs exhaustive enumeration** for d = 3 and d = 4, `N` up to 60, four spectra: zero
  mismatches.
- **Proposition 1 (winding marginalization)** is exact as stated — the coefficient is the Dirichlet
  kernel, `|D_M(ν)| ≤ 1/(M|sin πν|)`, and the `1/M` falloff is correct for fixed `ν ∉ ℤ`.
- **Manuscript lint:** no undefined refs or cites, environments balanced, brace and `$` parity clean,
  nine figure floats.

## What is still open (referee items I did not touch)

1. **`T_C ≈ ħ/ΔE_C` is still not derived** (REVIEW M2). v10 showed the cosmological version of this
   question fails as a prediction; the scope condition is in `notes/drifting_commensurability.md`.
2. **No standard-QM comparison for the circuit-QED protocol** (REVIEW §8.4). The paper states the
   null hypothesis but does not write the competing prediction side by side.
3. **The cosmological section is still the weakest part.** Every statement in it is either known
   (Hawking temperature, Bisognano–Wichmann) or explicitly heuristic. v11 makes Sec. V.C narrower
   and harder to attack; nothing has yet done that for Sec. IX.
4. **Frame 02's ring spacing** in the author's photon model is still unverified — uniform (perspective)
   vs geometric (real discrete scale invariance). That one is a screenshot measurement, not a
   calculation.
