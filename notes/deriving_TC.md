# Deriving T_C: the clock's own phase noise (v14)

*Working note. Not paper content — but this is the calculation the paper needed.*

## The problem it solves

`REVIEW.md` §M2, the referee's central objection:

> "This *real clock* dephasing is already in the literature... so it **cannot by itself tell the conjecture
> apart from standard QM**. The paper must say what the conjecture predicts *beyond* known clock noise."

And the paper's own admission (intro): without a law for `T_C` it is "an interpretation empirically
equivalent to the PW formalism". `T_C` was a free knob, and a free knob makes any prediction
unfalsifiable. The v10 cosmological attempt to pin it failed. This is the surviving route.

## The calculation

Take the clock to be what the experiment actually uses: a cavity mode of frequency ω, linewidth κ,
holding a coherent state of n̄ photons.

1. **The phase diffuses.** A coherent field of n̄ photons has the Schawlow–Townes linewidth
   `Δν = κ/(4π n̄)`. In angular frequency the HWHM is `κ/(4n̄)`. Pure phase diffusion with coefficient
   `D_φ` gives a field correlation `⟨a†(t)a⟩ ∝ e^{−D_φ t/2}`, a Lorentzian of HWHM `D_φ/2`. Equate:
   **`D_φ = κ/(2n̄)`.**

2. **Convert to a clock reading.** The clock reads `t = θ/ω`, so a phase variance `D_φ t` is a timing
   variance `D_φ t/ω²`. Comparing with `σ_t² = T_C t`:

   **`T_C = κ/(2 n̄ ω²)`**

3. **Dephasing.** Substituting into `Γ_jk = T_C ω_jk²/2`:

   **`Γ_jk = κ ω_jk²/(4 n̄ ω²)`,  `T₂ = 4 n̄ ω²/(κ ω_jk²)`**

Every symbol is a device parameter. Nothing is fitted.

## The numbers (transmon 5 GHz, cavity 6 GHz, κ/2π = 1 MHz)

| n̄ | D_φ [rad²/s] | T_C [s] | Γ [s⁻¹] | T₂ [s] |
|---|---|---|---|---|
| 1 | 3.14e6 | 2.2e-15 | 1.09e6 | 9.2e-7 |
| 10 | 3.14e5 | 2.2e-16 | 1.09e5 | 9.2e-6 |
| 100 | 3.14e4 | 2.2e-17 | 1.09e4 | 9.2e-5 |
| 10³ | 3.14e3 | 2.2e-18 | 1.09e3 | **9.2e-4** |
| 10⁴ | 3.14e2 | 2.2e-19 | 1.09e2 | 9.2e-3 |

Measured transmon T₂ is 0.1–1 ms. **At n̄ ≈ 10³ the prediction sits at 0.92 ms — inside the measured
window.** The framework is no longer unfalsifiable: the number is forced, and it lands where
experiments already live.

## What this changes in the paper

This **flips the verdict on the circuit-QED protocol**, and that has to be handled carefully.

- Old: `T_C ~ ħ/ΔE_C` with `ΔE_C = ħω_s` gave `T₂ ≲ 10 ns` → "the most natural identification is
  already excluded by measured coherence times" (abstract (v), §VI item 6).
- New: the derived `T_C` is 2.2e-18 s at n̄ = 10³, seven orders smaller than the dimensional estimate
  `1/ω_s = 3e-11 s`. The same identification now gives `T₂ ≈ 0.9 ms`.

So the *exclusion* belonged to the heuristic, not to the framework. The paper now says exactly that:
the heuristic is excluded, the derived law is testable, and Table III gives the numbers.

## The honest caveats (all in the paper)

1. **The mechanism is textbook.** Schawlow–Townes phase diffusion, plus Egusquiza–Garay–Raya clock
   noise. No new physics. What is new is only that the free parameter is *gone*.
2. **The amplitude must be maintained.** A merely damped coherent state keeps its phase exactly
   (amplitude damping maps coherent states to coherent states), so the phase noise belongs to whatever
   supplies the gain. In steady state spontaneous photons enter at the loss rate, which is why κ
   appears.
3. **It is not a discriminator by itself.** The `n̄⁻¹` and `κ` scalings are the standard ones; ordinary
   clock-noise theory predicts the same. The value is falsifiability, not uniqueness.
4. **The `χ` identification is still hopeless:** `(χ/ω_s)² = 4×10⁻⁸` smaller, i.e. `Γ ≈ 4×10⁻⁵ s⁻¹`.
   That is unchanged, and it is worth stating plainly rather than hiding.

## What is genuinely new to the programme

Before v14 the particle-physics reach of this paper was zero: every quantitative statement contained a
free parameter chosen to make it work. `Γ = κω_jk²/(4n̄ω²)` contains none. That is the difference
between an interpretation and a theory with a number in it. It may well be the only such number in the
paper, and it should be presented as what it is.

## Code

`src/qtime/clock_noise.py`; tests `tests/test_clock_noise.py` (9 checks, including a check of the
`D_φ = κ/(2n̄)` factor against the textbook linewidth formula rather than against itself).
