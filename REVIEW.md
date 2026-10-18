# Referee-Style Review: *The Hilbert-Space Quadrature Conjecture*

**Manuscript reviewed:** "The Hilbert-Space Quadrature Conjecture: Time as the Wave-Function-Generated Quadrature of Circular Phase", S. M. A. Owens. This is the consolidated REVTeX draft (Secs. I–IX plus Appendix A). The review also covers the add-on derivations supplied separately: optimal clock bandwidth ΔE\*, the inflationary back-reaction spectrum, the non-Abelian extension, the simulation code and notebook, the README, and the earlier "structural audit".

**Date:** 18 October 2026

**Status of claims checked:** every equation flagged below was re-derived by hand. Where possible it was also checked numerically; see `tests/` (24 tests, all passing).

**Section numbers** in this review (Sec. I–IX) refer to the *original* draft. The revised manuscript is numbered differently after the second revision: circuit QED is Sec. VI, thermal time Sec. VII, the local relativistic regime Sec. VIII, the cosmological regime Sec. IX, and the discussion Sec. X. See `CHANGES.md`.

---

## 1. Summary verdict

The earlier audit called the draft "publication-ready" and "mathematically airtight". **It is not, and I would not submit it in its current form.** But there is a real, defensible paper inside it. Here is the honest picture:

| Part | Assessment |
|---|---|
| Core idea: Page–Wootters (PW) + covariant clock POVM + Born-rule "quadrature" reading (Secs. III–IV) | **Correct, but mostly already published.** Theorem 1 is the standard PW result for a covariant clock observable. Prior work includes Smith & Ahmadi, *Quantum* **3**, 160 (2019); Höhn, Smith & Lock, *PRD* **104**, 066001 (2021); and Loveridge & Miyadera, *Found. Phys.* **49**, 549 (2019). None of these is cited. |
| The "squaring the circle" geometric reading (Hopf fibre → Born map → time) | **New as an interpretation.** It is a framing, not a new result, and the paper should present it that way. |
| Winding map Π_ω (Prop. 1), said to remove the commensurability problem | **Wrong as stated.** Periodization brings the commensurability condition back (Poisson summation). A correct replacement exists (Zak transform), and it gives a sharper result than the original (see §3, issue M1). |
| Lindblad "theorem" (Thm. 2) | **The algebra is correct, but the result depends on an assumption.** A finite-resolution clock gives a *fixed* suppression factor. Exponential decay appears only if one *assumes* the clock's time error grows diffusively. The identification T_C ~ ħ/ΔE_C is not derived. |
| Circuit-QED proposal | **Not a test of the conjecture as written.** The predicted rate depends on an unjustified choice of energy gap (ħχ). With the qubit's actual splitting, the same formula predicts T₂ < 10 ns, which existing devices already rule out by about 4 orders of magnitude. |
| Relativistic / curved / cosmological extensions (Sec. IX) | **Mixed.** The Klein–Gordon positive-frequency, proper-time and twin-shell results are correct but known. The section also contains several concrete errors: a wrong Hamiltonian, wrong scalings, a sign error, an internal contradiction about Zitterbewegung, and a wrong modular operator. The Hartle–Hawking "theorem" and the near-horizon "derivation" are heuristics presented as theorems. |
| Add-ons not in the consolidated file (ΔE\*, back-reaction spectrum, non-Abelian) | ΔE\* has an algebra slip and a conceptual flaw. The back-reaction "theorem" is not derived, and its main formula has no observable effect. The non-Abelian extension is speculative and is not ready for a paper. |

**Recommendation:** major revision. The paper should present itself as *a geometric/operational synthesis of relational quantum time with covariant clock POVMs, plus a clarified finite-clock dephasing model*. It should not claim to "prove" or "resolve" foundational problems. The revised manuscript in `paper/` makes exactly this change. Every correction is listed in `CHANGES.md`.

---

## 2. What holds up

1. **Covariant Holevo POVM on the Hardy space** (Prop. "Covariant Circular Clock POVM"). Positivity, completeness and shift covariance are proved correctly. Weak convergence of the Pegg–Barnett measures is fine for vectors of finite support. *(One error in the Szegő-kernel formula next to it; see m1.)*
2. **Theorem 1.** For an exact constraint state with commensurate spectrum, the raw clock marginal is uniform and the conditional state obeys the Schrödinger equation exactly. Verified numerically in `test_exact_constraint_state_gives_uniform_marginal_and_exact_unitary_evolution`.
3. **The three-way split** into marginal M_Ψ, joint J_Ψ and conditional ρ_S is useful bookkeeping and answers the "uniform-marginal" objection to PW cleanly.
4. **Spectral-window wave packet** (Prop. 2). The bound ‖Ĉ Ψ_w‖ = ħ‖w′‖ is correct; for a Gaussian, ħ‖w′‖ = ħ/(2ΔT).
5. **The GKSL algebra** in Thm. 2, *given* the diffusive assumption. Verified three independent ways (direct quadrature, closed form, ODE integration).
6. **Klein–Gordon factorization and positive-frequency selection; proper-time Hamiltonian H ≈ E_P + γ⁻¹H_int; twin-shell visibility** τ_grav = √2 ħ/(ΔE_int|ΔN₁₂|). All correct, and all match the cited literature (Smith & Ahmadi 2020; Zych et al. 2011; Pikovski et al. 2015).
7. **Hawking temperature** T_H = ħκ/(2πk_B c) = ħc³/(8πGMk_B). Correct: 6.17×10⁻⁸ K for one solar mass.

---

## 3. Major issues

### M1. The winding map does not avoid the commensurability obstruction (Sec. II.C, Prop. 1)

The proposition claims that periodization, (Π_ω f)(θ) = Σ_m f((θ+2πm)/ω), carries universal-cover physical states to the circle while avoiding the condition E_k/ħω ∈ ℤ.

There are two problems:

- Π_ω is defined on Schwartz functions, but physical states ⟨t|Ψ⟩ = e^{−iĤ_S t/ħ}ψ₀ do not decay, so the series does not converge.
- Interpreted as a distribution, Poisson summation gives

  Σ_m e^{−iE(θ+2πm)/ħω} = e^{−iEθ/ħω} · 2π Σ_l δ(2πE/ħω − 2πl),

  which is **non-zero only when E/ħω ∈ ℤ**. The map restores exactly the obstruction it was meant to remove.

The "partition of unity" step in the proof also fails. A smooth window supported strictly inside one fundamental domain cannot satisfy Σ_m |w(t − mT₀)|² = 1 everywhere.

**Fix applied.** Replace Π_ω with the full **Zak (Bloch–Floquet) transform** (Zf)(θ,κ) = Σ_m e^{−iκm} f((θ+2πm)/ω). This map is unitary, and periodization is just its κ = 0 fibre. Each system energy E then sits in the fibre κ ≡ −2πE/ħω (mod 2π), where the clock "number" spectrum is ℤ + κ/2π. This is a circle clock with twisted boundary conditions. The physical consequence is clean and testable:

> **A cyclic clock read without its winding number keeps coherence only between commensurate energies.** Averaging over M windings multiplies ρ_jk by (1/M) Σ_m e^{−2πimν}, with ν = (E_j − E_k)/ħω. This factor is 1 for ν ∈ ℤ and falls off as 1/M otherwise (`figures/fig2_winding.png`, `test_winding_marginalization_kills_incommensurate_coherence`).

So the operational clock must be the pair (θ, m), which amounts to the universal cover. The circle alone is enough only within one period, or for commensurate spectra. This is a stronger and more honest statement than the original one, and it fits the "circle that squares to a line" metaphor well.

### M2. "Finite clock energy spread ⇒ Lindblad equation" depends on an added assumption (Sec. V, Thm. 2)

- For a Gaussian timing kernel of **fixed** width σ_t, the exact result (the paper's own Eq. exact_gaussian_fourier) is a **time-independent** factor e^{−ω_jk²σ_t²/2}. Exponential decay appears only after substituting σ_t² = T_C·t, a diffusive clock error that is *assumed*, not derived. `figures/fig1_conditioning.png` (b–c) shows the difference: a sharp-but-unsharp readout gives a flat line, and only the diffusion assumption gives decay.
- **Within an exact PW state, the clock's energy spread equals the system's.** For Ψ = Σ α_k |n_k⟩|E_k⟩, the clock's energy distribution is just |α_k|², so "ΔE_C" is not an independent parameter. Even a clock wave packet Σ_n g(n)|n + n_k⟩ factors out of (θ|Ψ⟩ and causes *no* dephasing. What does cause dephasing is (i) unsharp *readout* (a smeared POVM), or (ii) clock noise or clock–system interaction (Smith & Ahmadi 2019). The revised text separates these cases.
- With the natural choice σ_t = ħ/(2ΔE_C), the fixed effect is unmeasurably small: 1 − 3×10⁻¹⁰ for the circuit-QED numbers.
- This "real clock" dephasing is already in the literature: Egusquiza, Garay & Raya, *PRA* **59**, 3236 (1999), and Gambini–Porto–Pullin. It is ordinary quantum mechanics with a noisy clock, so it **cannot by itself tell the conjecture apart from standard QM**. The paper must say what the conjecture predicts *beyond* known clock noise.

**Fix applied.** Thm. 2 is split into (a) an exact proposition: unsharp readout gives fixed Gaussian suppression; and (b) a clearly labelled model assumption: diffusive clock gives GKSL equation. T_C is kept as a phenomenological parameter.

### M3. The circuit-QED protocol does not test the conjecture (Sec. VII)

- The rate Γ = (E₁ − E₀)²/(2ħΔE_C) is evaluated with E₁ − E₀ = ħχ (the dispersive shift, 1 MHz). But the transmon's energy gap relative to the cavity clock is ħω_s (5 GHz). Using χ quietly assumes a frame rotating at ω_s, which is itself defined by an *external* lab clock. That defeats the purpose of a relational test.
- With ΔE = ħω_s, the paper's own formula gives Γ ≈ 1.3×10⁸–4×10⁹ s⁻¹ for n̄ = 10⁴…10, so T₂ ≈ 0.2–8 ns. Real transmons show T₂ ~ 0.1–1 ms. **Under the natural identification, the model is already falsified** (`figures/fig4_cqed.png`, `estimates.circuit_qed_table`).
- With ΔE = ħχ, Γ ≈ 5–170 s⁻¹, which is 1–3 orders of magnitude below ordinary qubit dephasing. That is hard to measure but not impossible. Measurement-induced dephasing grows with n̄, so the opposite n̄^{−1/2} trend would be a usable signature, *if* the identification were justified.
- Heterodyne detection implements a *smeared* covariant phase POVM (the phase marginal of the Husimi Q-function), not the canonical Holevo POVM.
- *(Added in revision 2.)* The dispersive Hamiltonian ħχa†a|1⟩⟨1| holds only for n̄ ≪ n_crit = Δ²/4g² ≈ Δ/4χ. For the paper's own parameters (Δ/2π ≈ 1 GHz, χ/2π ≈ 1 MHz) this gives n_crit ≈ 250, so the proposed sweep n̄ = 10…10⁴ leaves the model's regime. Sweeping three decades inside it needs χ/2π ≲ 2.5 kHz, which pushes the χ-rate down to 10⁻³–10⁻⁵ s⁻¹ (`estimates.critical_photon_number`, `test_dispersive_limit_restricts_photon_sweep`).

**Fix applied.** The original Sec. VII is rewritten as "A laboratory test: circuit QED" (now Sec. VI, directly after the master equation). It states the null hypothesis, shows both identifications, and lists what a discriminating experiment would need.

### M4. Wrong modular operator, and the "equivalence" claim is too strong (Sec. VI, Prop. "Modular Time Equivalence")

- For the algebra 1_C ⊗ B(H_S) and the vector Ψ, the modular operator is **Δ_Ψ = ρ_C⁻¹ ⊗ ρ_S**, not ρ_C ⊗ ρ_S⁻¹. The stated form gives the flow with the opposite sign from the one used in the proof. This was checked by building Tomita's S = J Δ^{1/2} directly from its definition (`test_modular_operator_is_rhoC_inverse_tensor_rhoS`).
- The "equivalence" holds only when the Schmidt spectrum is exactly Gibbs, |α_k|² ∝ e^{−βE_k}. In that case it is the standard KMS fact underlying Connes–Rovelli, not a new equivalence.

**Fix applied.** Δ corrected; the result is restated as "consistency with thermal time for Gibbs Schmidt spectra".

### M5. Internal contradiction: Hardy clock vs. Zitterbewegung (Sec. IX.A vs. IX.B)

Sec. IX.A argues that the Hardy projection *removes* the negative-frequency sheet. Sec. IX.B then has the same Hardy clock *damp* Zitterbewegung, which needs interference between positive- and negative-energy components. Both cannot be true for the same clock: with a Hardy clock there are no negative-energy components left to interfere.

**Fix applied.** The Zitterbewegung result is now stated only for a clock with a two-sided spectrum (the universal cover). It is noted that the Hardy restriction removes Zitterbewegung outright, which is the Foldy–Wouthuysen/Newton–Wigner choice. Also, the summary and slides write Γ_Zitt ≥ 2m²c²/(ħΔE_C); that has the wrong dimensions. The correct expression is 2m²c⁴/(ħΔE_C).

### M6. Curved-spacetime Hamiltonian and near-horizon scaling (Sec. IX.D, IX.G)

- Ĥ_grav was written as √N (c²P² + M²c⁴) N, which has dimension energy² and is not symmetric. The correct form is **Ĥ_grav = N^{1/2} √(c²P̂_Σ² + M̂²c⁴) N^{1/2}**.
- The text sets r = r_s + ρ (coordinate distance) and then uses N ≈ κρ/c², which holds only for **proper** distance. Because the rate goes as Γ ∝ 1/N(r_C), it scales as **ρ_proper⁻¹**, equivalently (r − r_s)^{−1/2}. The claimed "ρ^{−1/2}" mixes the two definitions (`test_near_horizon_lapse_is_linear_in_proper_distance`).
- In the near-horizon proof, the boost/Killing relation should read K = c Ĥ_Killing/(ħκ), not c³/(ħκ). In the add-on, T_H = ħκ/(2πk_B) is missing a factor of c.
- The "derivation" of T_H puts the answer in by hand. The temperature comes entirely from Bisognano–Wichmann (Minkowski) or Sewell (1982, Schwarzschild/Kruskal); nothing in the quadrature framework produces it. Also, the wedge algebra is type III₁, so "ρ_wedge = Tr_ext|Ψ⟩⟨Ψ|" and K = −ln ρ_wedge are only formal. Finally, the Γ divergence comes from the clock's *redshifted energy spread*, not from "thermal blueshift"; these are unrelated mechanisms.

**Fix applied.** Hamiltonian, scaling and conversion factor corrected. The section is relabelled "Consistency with Bisognano–Wichmann/Sewell", with the type III caveat.

### M7. Minisuperspace section: sign and dimension errors (Sec. IX.E)

- The scale-factor kinetic term is missing the 1/a factor of the ADM Hamiltonian, H ⊃ −p_a²/(2m₀a). In the isolated equation, the potential term has the wrong sign. In Ĥ_cosmo, the prefactor 2m₀/a³ does not match the previous line's 2a³.
- With the signs corrected, Ĥ²_cosmo contains −2a⁶V. For V > 0 this operator is not positive, so the square root is not well defined in general. The φ-clock is exact only for a **massless** scalar (V = 0), which is the known Wheeler–DeWitt/loop-quantum-cosmology model of Ashtekar, Pawłowski & Singh (2006). For V(φ) ≠ 0 the Hamiltonian depends on time and the square-root factorization is only approximate.
- φ ∈ ℝ with an unbounded momentum p_φ. That is the universal-cover clock, not the Hardy-circle clock the paper builds elsewhere.
- "Tomonaga–Schwinger" is the wrong name; that term refers to many-fingered time in field theory.
- σ_H² = (∂_φH)²Δφ² + γ_cosmo/Δφ² is an ansatz with an undefined γ_cosmo. It cannot be called a prediction.

**Fix applied.** Equations corrected; the exact case (V = 0) is presented; the Hubble-variance formula is moved to "open problems".

### M8. The Hartle–Hawking "theorem" is not a theorem (Sec. IX.F)

The state |Ψ_NB⟩ = ∫dθ_E |θ_E) ⊗ ∫_cap^{θ_E} 𝒟g 𝒟φ e^{−I_E} |geom(θ_E)⟩ is formal. "∫_cap^{θ_E}" is undefined. The Euclidean action is not a generator of translations. The Hartle–Hawking wave function depends on boundary data, not on a Euclidean time label, because Euclidean time is integrated over. The clock space is called both L²(S¹_E) and Hardy. The "equator at θ_E = π/2" mixes the S¹ angle with the hemisphere's polar angle.

**Fix applied.** Downgraded to a clearly labelled *heuristic proposal*, with these caveats stated.

### M9. The paper overstates novelty and rigour throughout

The phrasings "we prove", "resolves", "exact", "mathematically airtight", "eliminates the problem of time" and "resolving the historical Klein–Gordon negative-probability problem" do not fit what is shown. (The positive-frequency restriction is standard; Newton–Wigner localization is not Lorentz covariant and spreads faster than light, per Hegerfeldt 1974.) A referee in this area will know Smith–Ahmadi, Höhn–Smith–Lock and Egusquiza et al. Leaving them out will read as either unawareness or overclaiming.

**Fix applied.** The introduction now has an explicit "relation to prior work" paragraph and a "what is new / what is not" list. "Conjecture" is kept for the conjecture itself; "theorem" is reserved for statements that are actually proved.

---

## 4. Minor issues

| # | Location | Issue | Fix |
|---|---|---|---|
| m1 | Eq. (Szegő kernel) | Σ_{n≥0} e^{inx} = πδ(x) + 1/(1 − e^{i(x+i0)}) counts the delta twice; the +i0 prescription already contains it. | Σ_{n≥0} e^{inx} = 1/(1 − e^{i(x+i0)}) = πδ_{2π}(x) + ½ + (i/2) P.V. cot(x/2). Checked numerically. |
| m2 | Sec. II, Thm. 1 | Sign conventions: C_circ = ħωN − H_S is equivalent to a clock entering with the opposite sign. A semi-bounded clock represents only system energies in a window (after shifting H_S). | Remark added: the Hardy restriction is a *selection rule*, not a free lunch. |
| m3 | Thm. 1 | Needs E_k/ħω ∈ ℕ₀ (not merely ℤ); stated, but the abstract implies generality. | Scope made explicit. |
| m4 | Prop. 2 | Unitary conditional evolution holds by construction for |Ψ_w⟩, so the proposition is a consistency check, not a derivation. | Said so. |
| m5 | Sec. VI, Wigner | Fine. Note that "squaring" here is the ordinary Wigner transform. | – |
| m6 | Sec. IX.C | Twin-paradox visibility factor applies to a two-level internal superposition (n, m); state this. | Clarified. |
| m7 | Abstract | Lists "inflationary backreaction", but no such section exists in the manuscript. | Removed from abstract; discussed under open problems. |
| m8 | Audit, slides | Claim "Sections I–IX.H"; the consolidated file has IX.A–G. ΔE\* and back-reaction are missing. | Structure fixed; ΔE\* included as a remark. |
| m9 | Bibliography | Kuchař proceedings editor is **G.** Kunstatter, not S. Missing references: Smith–Ahmadi 2019; Höhn–Smith–Lock 2021; Loveridge–Miyadera 2019; Egusquiza–Garay–Raya 1999; Zak 1967; Bisognano–Wichmann 1975/76; Sewell 1982; Ashtekar–Pawłowski–Singh 2006; Newton–Wigner 1949; Hegerfeldt 1974; Busch–Grabowski–Lahti 1995; Unruh 1976; Hawking 1975. | Added. |
| m10 | Source file | Ends with `end{document}` (missing backslash). The PDF export turned inline math into Unicode fragments, so the .tex shipped inside the PDFs does not compile. | Rebuilt from scratch; compiles. |
| m11 | Throughout | Name: "Bloch–Floquet winding quotient" describes plain periodization. With the Zak transform, "Bloch–Floquet" becomes accurate. | Renamed. |

---

## 5. Add-on material (not in the consolidated manuscript)

**Optimal clock bandwidth ΔE\*.**
- The minimizer ΔE\* = |ΔE_jk|/(2^{1/4}√|ΔN₁₂|) is correct *for the model Γ = a/x + bx*.
- The minimum value is wrong. It should be **Γ_min = 2^{1/4}|ΔE_jk| |ΔN₁₂|^{1/2}/ħ**, not 2^{1/4}|ΔE_jk||ΔN₁₂|/ħ.
- Conceptually, the formula adds an exponential (Markovian) rate to the inverse of a *Gaussian* decay time, and it combines decoherence of two different quantities (system energy coherence vs. clock centre-of-mass coherence). Treat it as a heuristic scaling only.
- For scale: a 1 m vertical separation on Earth gives ΔN ≈ 1.1×10⁻¹⁶, so ΔE\* ≈ 8×10⁷ × the system gap. Gravity is irrelevant for lab clocks.
- Included in the revision as a labelled remark.

**Back-reaction power spectrum.**
- ⟨δφ²⟩ = ∫₀^{aH} d³k/(2π)³ |u_k|² diverges logarithmically in the infrared; H²/4π² is the variance *per e-fold*, not the total.
- The "proof" says only that convolution "introduces" the exponential.
- The factor exp[−H²(k/aH)²/(8π²Δφ²)] equals the k-independent constant exp[−H²/(8π²Δφ²)] at horizon crossing and tends to 1 on super-horizon scales. So it produces no tilt and no running, contrary to what is claimed.
- Not included. Listed as an open problem.

**Non-Abelian extension.**
- The group "phase states" |g) = Σ √d_R D^{(R)}_{mn}(g)|R,m,n⟩ are fine (Peter–Weyl). But the constraint η_a Ĵ^a + Ĥ is not gauge invariant unless η_a is itself covariant.
- Γ_YM = g²C₂(R)/(2ħΔE_YM) is not derived: g² is dimensionless, so the formula has the wrong dimensions unless an energy scale is supplied.
- Not included.

---

## 6. Code review

| File | Problem | Fix |
|---|---|---|
| `quadrature_time_simulation.py` / notebook | It is labelled "exact bipartite conditioning", but **no bipartite state is ever built**. The code convolves hand-written ideal qubit states with a classical kernel. The parameter called "energy spread Δn" is really the *readout resolution*. For non-integer `omega_ratio` (allowed by the notebook slider: 0.5, 1.5, …), e^{−iωθ} is not periodic on the circle, so the result has a discontinuity at 2π. That is the commensurability issue showing up in the code itself. | `src/qtime/clocks.py` builds |Ψ⟩ ∈ ℂ^N ⊗ ℂ^d explicitly, constructs sharp and smeared covariant POVMs, and computes conditional states by real partial traces. Commensurability is enforced, with an error raised otherwise. |
| Same | The Lindblad "analytical envelope" is a constant line, which contradicts the paper's claim of exponential decay. The code was right; the claim was wrong (see M2). | `figures/fig1_conditioning.png` (c) compares the two honestly. |
| `simulate_rel_clock_interferometer` | The physics is correct. But the internal spectrum is equally spaced, so V(t) **fully revives** every T_rev = 2π/(ω₀|Δ|) ≈ 1.29. The plotted window t ∈ [0, 15] contains about 11 revivals. The Gaussian "decay envelope" overlaid on it holds only before the first revival, so this is reversible dephasing, not decoherence. | `interferometer.py` adds the revival period and suggested time step; `figures/fig3_interferometer.png`. |
| Notebook JSON | The JSON shipped in the PDF lost its `**` operators in the export, so it would not run. | Regenerated: `notebooks/quadrature_time.ipynb`. |

---

## 7. About the earlier audit and the chat assistance

The earlier "comprehensive structural audit" rated the draft "Publication-Ready", "mathematically airtight", and "structurally complete, mathematically rigorous, and physically profound". Those ratings were not justified. They missed every item in §3. The first chat's claim that the synthesis is "distinctly yours" is fair for the *framing*, but it should have pointed out that PW with covariant POVMs had already been worked out in 2019–2021. I would treat the chat-generated "drop-in theorems" (Hartle–Hawking, near-horizon, back-reaction, non-Abelian) with particular caution. They have the form of results without the content.

## 8. Path to submission

1. Use the revised manuscript (`paper/quadrature_time.pdf`) as the base. It keeps your framing and conjecture, corrects the errors, and states novelty accurately.
2. Choose venue to match the scope. The conceptual synthesis plus the winding-number result and the clarified finite-clock model fits *Foundations of Physics*, *Quantum*, or *Classical and Quantum Gravity* (relational-dynamics papers). It does not fit PRD/PRL as a "resolution of the problem of time".
3. The most promising original direction is to develop **M1** (winding-number marginalization / Zak fibres) into a proper result about cyclic clocks. It is new, correct, simple, and directly testable with a two-level system timed by a cyclic clock without a cycle counter.
4. Before claiming an experimental signature, write the standard-QM prediction for the identical setup and show where the conjecture departs from it.
5. Get an expert in relational quantum dynamics to read it before submission.
