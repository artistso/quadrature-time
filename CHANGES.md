# Change Log

## v11: which clocks can be relational clocks, and what tuning buys (2026-10-18)

The previous entries treated the clock quantum $\Delta=\hbar\omega$ as a free design parameter.
It is not, for any real clock. This entry states the scope condition that follows, and it
sharpens Prop. 4 rather than adding to it.

### Manuscript (Sec. V.C)

| Addition | Where |
|---|---|
| **Lemma (admissibility).** With smallest gap $g$: $\Delta\le g$ suffices for distinct labels; $\Delta\le R/(d-1)$ is necessary; $N\ge1+R/\Delta$ levels are needed. | `lem:admissible` |
| **Remark (what tuning buys).** A *pinned* quantum is forced to the rounding assignment, so it is confined to the universal $N^{-2}$ and a horizon growing like $N$. The $N^{-4}$ law requires the freedom to sit at a convergent, and then beats the pinned clock by a factor growing like $N^2$. Stated as the scope condition both propositions should be read with. | `rem:tuning` |
| New Fig. 9 (`fig:tuning`): the two scalings and the $N^2$ gain. The wobble in the pinned curve is Diophantine luck, not noise, and the caption says so. | `scripts/make_figures.py::fig_tuning` |

### Code (v11)

| Addition | Where |
|---|---|
| New module `src/qtime/clock_scales.py`: `design_rule`, `admissible` (sufficient test first, then an offset scan between the sufficient and necessary bounds), `levels_needed`, `mean_variance_at_quantum` (offset-averaged, since a clock's zero is not physical), `horizon_at_quantum`, `tuning_gain`, and a table of `PHYSICAL_CLOCKS` with reference quanta. | new |
| 15 new tests, suite 127 -> **142**. | `tests/test_clock_design.py` |

Two things I got wrong on the way, recorded because they are instructive:

1. My first example spectrum was $(0,0.5,1.0)$ eV, which is *commensurate* -- every quantum fits it
   exactly and the horizon is infinite. The phenomenon only appears for irrational ratios; the
   module's examples use $(0,1,\varphi)$.
2. The gain from tuning is not monotone in $N$: at $N=90$ it dips by an order of magnitude because
   the pinned levels happen to land almost exactly on integers. That is the same arithmetic that
   governs the whole series, showing up as scatter, and it is now asserted as a *trend* ($N^2$)
   rather than as monotonicity.

## The spiral staircase, and a scope condition from an expanding universe (2026-10-18)

### Manuscript

| Addition | Where |
|---|---|
| New opening to Sec. II: **the clock as a spiral staircase**. Walking along it is unrolled time (the universal cover); from above it is the circle, with the lap unrecoverable, which *is* the winding number; from the side it is a sinusoid, the oscillation a clock represents. Line, circle and oscillation are three projections of one helix, and each projection is a choice of what to forget -- the line forgets nothing and costs a ground state (Sec. III), the circle forgets the lap and costs coherence (Prop. 2), the oscillation keeps neither. | `quadrature_time.tex`, start of Sec. II |
| New Fig. 3 (`fig:staircase`): the helix, its two projections, and the top view showing four instants at one face point (ties directly to Fig. 1). | `scripts/make_figures.py::fig_staircase` |
| New Fig. 4 (`fig:staircaseconical`): the vanishing point. A narrowing staircase is a conical helix; its top view is a logarithmic spiral, self-similar under a constant factor per lap, with a turn rate that diverges at the axis. What looks like a walker disappearing into a centre is a fixed point of the projection. Caption claims are asserted in the code against a numerical integration, not against themselves. | `.py::fig_staircase_conical` |

Figure numbers shift again from here on; all cross-references are by label.

### Code

| Addition | Where |
|---|---|
| **v10, first step: commensurability drift.** `scale_factor` (de Sitter default, matter+Lambda inverted from the exact integral), `closed_mode_frequency_s3` ($\omega_k=c\sqrt{k^2-1}/a$ on a round $S^3$), `winding_envelope` ($|D_M|$), `commensurability_ratio`, `drift_rate`, `drift_over_efolds`, `first_unprotected_efolds`. | `clock_design.py` |
| Six new tests (suite 127 -> **133**). | `tests/test_clock_design.py` |
| **The result is negative, and is written up as such:** since $\omega\propto1/a$ gives $\dot\nu=+\nu H$, a coherence protected at early times is destroyed in **under a hundredth of an e-fold**, independently of mode number, $M$ and coupling. So the drift is a sharp cutoff, not a slow sweep -- a scope condition (commensurability can be used within an era, never across one), not a prediction. Replaces the "worth chasing" note of the previous commit. | `notes/drifting_commensurability.md` |

## Author sketch -> Figure 1 (2026-10-18)

An author sketch ("2 AM no longer exists") turned out to be a literal illustration of
winding marginalization, so it is now **Fig. 1** (`fig:cover`) with a paragraph in Sec. II.D.

| Change | Where |
|---|---|
| New paragraph in Sec. II.D: the clock face as a covering map, why a reading needs its winding number, and why a *global* relabeling is a rotation of the face (inert) while a relabeling applied *after a seam* is a cut of the cover that deletes a label (DST). | `quadrature_time.tex`, after Prop. 2 |
| New `\label{sec:zak}` on the Zak-transform subsection, referenced from that paragraph. | Sec. II.C |
| New Fig. 1: (a) one reading, four instants; (b) relabeling function $g(t)$, constant for a global shift, a step for DST. The caption's numbers are asserted in the plotting code (4 occurrences of 2 AM in 96 h; 23 distinct labels on the seam day; no extra jumps under a global shift), so the figure cannot drift from its caption. | `scripts/make_figures.py::fig_clock_cover`, `figures/fig_cover.*` |
| Sketches archived with a reading of each; the second (time-dilation vortex) is recorded as *not yet physics*. | `archive/sketches/`, `archive/README.md` |

Figure numbers shift by one from here on (the conditioning figure becomes Fig. 2, and so on);
all cross-references are by label, so nothing else needed changing.

## v9 of the finite-clock series, and a new subsection in the manuscript (2026-10-18)

Closes the two open items from v8 and folds the results into the paper.

### Manuscript (Sec. V.C, new)

| Addition | Where |
|---|---|
| New subsection **"How many levels does a clock need?"**: Lemma (energy-mismatch fidelity bound), Proposition (universal bound, Eq. 25), Proposition (three-level closed form and the N^-4 law, Eqs. 26-27), plus a scope-and-numerics remark and a new figure. | `paper/quadrature_time.tex`, `sec:levels`; figure `fig:clockscaling` |
| Abstract gains clause (iv); "What is new here" gains item (d) (the old (d) becomes (e)); plan-of-the-paper and discussion updated; appendix "Numerical model" now lists the new checks; new reference Hurwitz (1891); new `lemma` theorem environment. | same file |
| Copy of `fig5_clock_scaling.pdf` into `paper/figures/` for the manuscript build. | `paper/figures/` |

The subsection keeps the paper's tone: the bound is stated as unconditional but loose, the
closed form as exact but restricted (three levels, uniform weights), and the extremal claim about
the golden ratio as a *numerical finding* with Hurwitz supplying only the exponent-side bound,
because the constant also carries r^4/V_E(r).

**`paper/quadrature_time.pdf` is stale**: there is no TeX in this environment, so the PDF has not
been regenerated. Run `cd paper && pdflatex quadrature_time && pdflatex quadrature_time`.

### Code (v9)

| Addition | Where |
|---|---|
| **Exact rational arithmetic.** `gap_ratio_fraction` (exact decimal value, so 2.5 is exactly 5/2), `exact_variance`, `three_level_closed_form_exact`. Commensurate spectra now give variance exactly zero instead of 1e-34; the covariance form `V_E - C²/V_q` is shown to be unusable past p ~ 2000, which is why exact arithmetic is required rather than optional. | `clock_design.py` |
| **Numerically stable residual form** `mismatch_variance_stable`, used for all optimiser scoring and for the vectorised proportional scan. | same |
| **General optimiser for d >= 4** `search_optimum`: grid starts in the (spacing, offset) design plane, fit-then-round fixed points, coordinate and window descent, spectrum-derived seed families. Validated against exhaustive enumeration on 300 randomised instances (d = 3 to N = 120, d = 4 to N = 36), zero mismatches; d = 4, N = 10^5 in 0.1 s. Documented as a heuristic with a known large-N coverage limit. | same |
| **Proportional family** `proportional_curve` (exact, moderate M) and `proportional_curve_stable` (vectorised, M = 10^6 in 0.2 s), plus `decay_exponent`. | same |
| **Correction to v8.** The exponent is not universal: powers of one quadratic irrational give N^-4 in d = 3, 4, 5; (1, sqrt2, sqrt3) gives N^-3.1. | `notebooks/finite_clock/*_v9.ipynb` |
| New notebook | `notebooks/finite_clock/Finite_Relational_Clock_Exact_Arithmetic_v9.ipynb` (self-verifying, 17 checks) |
| New figure | `figures/fig6_clock_dimension.{png,pdf}`, added to `scripts/make_figures.py` |
| Tests | 103 checks in `tests/test_clock_design.py`; suite total 24 -> 127 |

## v8 of the finite-clock series (2026-10-18) — new results, no change to the manuscript

Closes the question v7 left open: how the optimal mismatch variance V*_N decays with the
number of clock levels. Nothing in `paper/` is touched; this is the start of the next
paper, not a revision of the current one.

| Addition | Where |
|---|---|
| **Proposition 1 (universal bound).** Rounding each shifted energy to its nearest label gives \|δ_j\| ≤ Δ/2, so Popoviciu gives Var_p(δ) ≤ Δ²/4; the range R must fit in N−1 steps, so Δ ≥ R/(N−1). Hence V*_N ≤ R²/(4(N−1)²) and T(ε) ≥ 2(N−1)√ε/R. Spectrum-free, unconditional, no number theory. Verified: 0 violations for N = 4…128 over five spectra, rounding construction ≤ 0.222 of the bound. | `clock_design.universal_bound`, `universal_horizon`; notebook v8 §2 |
| **Proposition 2 (three-level closed form).** For E = (0,1,r), uniform weights: V*(q,p) = (p−rq)²/[(p+q)²+(2q−p)²+(2p−q)²]. The numerator is the rational-approximation error of the gap ratio, the denominator is 27V_E up to O(q). Verified against the density-matrix computation to 1e-9. | `clock_design.three_level_closed_form`; notebook v8 §3 |
| **The N⁻⁴ law and its constant.** Substituting \|r−p/q\| = κ/q² along convergents gives V*_N ≃ κ²r⁴/(27V_E)·N⁻⁴, so the guaranteed horizon grows like N², not N. At Fibonacci labels the golden-ratio constant φ⁴/(135V_E) = 0.1142350327 is matched to 3e-8 by p = 4181; measured exponents at convergent labels in 10³–10⁶ are −4.000 for φ, 1+√2, √2 and −4.06 for e−1. | notebook v8 §4–5; `figures/fig5_clock_scaling` |
| **The Diophantine reformulation.** Hurwitz's theorem (κ ≤ 1/√5, equality at φ) makes the golden ratio the extremal, hardest three-level spectrum. Stated as a numerical result, not a proof: the constant also carries r⁴/V_E(r). | notebook v8 §4, §8 |
| **O(N) exact optimiser.** For three levels, scan the middle label and a ±2 window around nint(rq); reproduces exhaustive enumeration for every spectrum tested up to N = 120 and reaches N = 10⁶ in 0.6 s (enumeration would need ~10¹¹ assignments). The nearest integer alone is *not* always optimal (r = φ, q = 4), which is why the window is needed. Validated, not proved. | `clock_design.three_level_optimum`, `scaling_curve` |
| New module | `src/qtime/clock_design.py` (~470 lines), exported from `qtime` |
| New tests | `tests/test_clock_design.py`, 72 checks; suite total 24 → **96** |
| New notebook | `notebooks/finite_clock/Finite_Relational_Clock_Asymptotics_v8.ipynb` (self-verifying: its last section re-runs 17 checks) |
| New figure | `figures/fig5_clock_scaling.{png,pdf}`, regenerated by `scripts/make_figures.py` |
| `notebooks/finite_clock/README.md` | Rewritten around the three results, with five concrete v9 candidates (d ≥ 4; proving the optimiser; the full extremal statement; exact integer arithmetic beyond the double-precision floor at p ≈ 2000; folding Props. 1–2 into the manuscript) |

Two caveats recorded in the module and the notebook rather than glossed:
the closed form needs uniform weights (with weights the problem stays a weighted
least-squares fit but no longer collapses to a two-vector determinant), and the ±2 window
in the optimiser is validated rather than proved.

## Repository maintenance (2026-10-18) — no physics content changed

A GitHub web upload of the whole working directory had flattened the repository: the
`paper/`, `src/qtime/`, `tests/`, `scripts/`, `figures/` and `notebooks/` directories
were dissolved and their contents written to the repository root, together with a
partial dump of `.git/objects` (38-character hex filenames) and the chat-generated PDFs
and transcripts from earlier sessions.

| Change | Detail |
|---|---|
| Layout restored | `paper/`, `src/qtime/`, `tests/`, `scripts/`, `figures/`, `notebooks/` recreated from the revision-2 tree recorded in `archive/bundles/quadrature-time.bundle` (commit `569be92`). All 29 files were confirmed byte-identical to the flattened copies first; `pytest` is 24/24 again and `scripts/make_figures.py` regenerates all four figures. |
| Duplicates removed | Root-level copies of `clocks.py`, `estimates.py`, `interferometer.py`, `modular.py`, `make_figures.py`, `test_clocks.py`, `test_physics.py`, `__init__.py`, `quadrature_time.tex|pdf|ipynb`, `fig1…fig4.*` (all byte-identical to their proper-path twins). |
| Stray git objects removed | The 37 hex-named files were loose git objects (verified by zlib header); every one of their contents is already in history. |
| Chat artifacts archived | 37 PDFs, `.txt` transcripts and images moved to `archive/chat-notes/`, with an inventory and a "check `REVIEW.md` before citing" warning in `archive/README.md`. The bundle moved to `archive/bundles/`. |
| Finite-clock series gathered | The four `Finite_*Relational_Clock_*` notebooks moved to `notebooks/finite_clock/`. |
| `README.md` | Repository-layout table updated; new paragraph pointing at the finite-clock series. |
| `notebooks/finite_clock/README.md` | **New.** Documents the v1 → v4 → v5 → v7 chain, the reduction of the design problem to a weighted least-squares fit, the universal `V*_N ≤ R²/(4(N−1)²)` bound (verified, 0 violations for N = 4…128 over six spectra), and the open asymptotic question. |
| `…_Fidelity_Bound_v5.ipynb` | Title cell said "v4" (copy-paste from the notebook it supersedes) and the closing section was numbered 7 twice; both fixed. Content unchanged. |

`notebooks/finite_clock` note: v6 is referenced by v7 ("This strengthens v6's grid
search") but was never committed; the chain in this repository is v1, v4, v5, v7.

## Revision 2 (response to the NotebookLM critique; see `NOTEBOOKLM_RESPONSE.md`)

| Change | Where |
|---|---|
| Circuit-QED section moved directly after the master equation and retitled "A laboratory test: circuit QED" | old rev-1 Sec. VII → **Sec. VI** |
| Thermal-time section now stands alone; Wigner–Weyl subsection moved into the quadrature section | rev-1 Sec. VI → **Sec. VII**; Wigner → **Sec. IV.A** |
| Extensions split by regime: **Sec. VIII** "Local relativistic regime: from flat spacetime to the horizon" (KG → Dirac → proper time → static spacetimes + ΔE\* → near-horizon) and **Sec. IX** "Cosmological and boundary regime" (minisuperspace → no-boundary conjecture → open directions) | rev-1 Sec. VIII |
| Near-horizon subsection moved from after Hartle–Hawking to the end of the local-gravity ladder | VIII.E |
| New "Plan of the paper" paragraph and a stress-test introduction to Sec. VIII | Sec. I, VIII |
| New interpretive paragraphs tying "squaring the circle" to the formal objects: after Def. 2 (the Born map is the squaring step), after Thm. 1 (flat squared measure = no preferred moment), in II.D (winding bookkeeping, not interference, straightens the circle), in IV.A (circle → ordering parameter, system → phase-space arena), with cautions | Secs. II.D, IV, IV.A |
| **New physics caveat:** dispersive limit n_crit ≈ Δ/4χ ≈ 250 for the paper's parameters, so the proposed n̄ sweep to 10⁴ leaves the regime; item 5 plus an explicit protocol paragraph with backgrounds and what a null result would mean; shading added to Fig. 3 | Sec. VI |
| Hubble-variance and back-reaction directions now in a dedicated "Open directions in this regime" subsection (still open, not claimed) | IX.C |
| Code: `estimates.critical_photon_number`, `estimates.chi_for_dispersive_sweep`, new test (24 tests total); figure title no longer hard-codes a section number | `src/`, `tests/`, `scripts/` |

## Revision 1 (original rebuild)

Section headings in this part use the **revision-1** numbering (before the reordering above).

#### Consolidated draft → revision 1

"Old" = the consolidated REVTeX draft from the chat PDFs (Secs. I–IX.G + Appendix A) plus its add-ons. "New" = `paper/quadrature_time.tex`. Issue codes (M1…, m1…) refer to `REVIEW.md`.

### Framing
| Old | New | Why |
|---|---|---|
| "We formulate and establish…", "we prove", "resolves four foundational tensions", "mathematically self-consistent foundation" | "We develop a geometric and operational reading…; contributions are interpretive and clarifying" | M9: the core results are the PW mechanism with covariant POVMs (Smith–Ahmadi 2019; Höhn–Smith–Lock 2021; Loveridge–Miyadera 2019) |
| No discussion of prior work | New paragraphs: *Relation to prior work*, *What is new here*, *What is not claimed* | M9 |
| Abstract lists "inflationary backreaction" | Removed (there is no such section); moved to open problems | m7 |

### Section II: Kinematics
| Old | New | Why |
|---|---|---|
| Prop. 1 "Bloch–Floquet winding quotient homomorphism" (periodization Π_ω) said to avoid commensurability | **Prop. 1 (Zak transform):** the line is a direct integral of twisted circles; periodization is the κ = 0 fibre. **Cor. 1:** each energy lives in fibre κ_E ≡ −2πE/ħω | M1: periodization brings back E/ħω ∈ ℤ (Poisson summation) |
| Partition-of-unity "isometric lift" | Removed (not needed; the Zak transform is unitary) | M1 (step was false) |
| — | **Prop. 2 (winding marginalization):** a cyclic clock read without its winding number keeps only commensurate coherences, factor D_M(ν) | New result; replaces the incorrect claim |
| — | Remark 1 on sign conventions: the Hardy restriction is a selection rule | m2 |
| Universal cover presented as the fix, with no caveat | Notes that P_t is unbounded, so this ideal clock has no ground state | honesty |

### Section III: POVM
| Old | New |
|---|---|
| (θ\|θ′) = πδ + 1/(1 − e^{i(x+i0)}) | (θ\|θ′) = 1/(1 − e^{i(x+i0)}) = πδ + ½ + (i/2) P.V. cot(x/2)  (m1) |
| Weak convergence "for all vectors" | "on finitely supported vectors" |
| — | Remark: twisted-fibre version of the POVM |

### Section IV: Conjecture and tri-fold structure
- Theorem 1 kept; attribution added (PW; Smith–Ahmadi; Höhn et al.); scope n_k ∈ ℕ₀ stated (m3).
- Prop. on spectral windows kept; Gaussian value ħ/2ΔT given; noted that unitarity holds by construction (m4).
- Conjecture now refers to winding number counting.

### Section V: Finite clocks (M2)
| Old | New |
|---|---|
| Thm. 2 "Finite-Clock Double-Commutator Master Equation" (finite ΔE_C ⇒ exact Lindblad) | **Prop. 5:** unsharp readout gives a *static* factor e^{−ν²σ²/2}. **Assumption 1:** diffusive clock σ_t² = σ₀² + T_C t. **Prop. 6:** GKSL equation under that assumption |
| T_C ~ ħ/ΔE_C stated as the rate | Labelled a dimensional heuristic, not derived; T_C is phenomenological |
| — | Remark: in an exact PW state the clock's energy spread is not an independent parameter |
| — | Attribution: Egusquiza–Garay–Raya 1999, Milburn, Gambini–Porto–Pullin |

### Section VI: Modular flow (M4)
- Δ_Ψ = ρ_C ⊗ ρ_S⁻¹ → **Δ_Ψ = ρ_C⁻¹ ⊗ ρ_S** (verified by building the Tomita operator numerically).
- "Modular Time Equivalence" → "Modular flow for Gibbs Schmidt spectra"; restricted to the Gibbs case.

### Section VII: Experiment (M3)
- "Experimental Protocol" → "Experimental considerations".
- Added Table I comparing ΔE = ħχ with ΔE = ħω_s. The latter predicts T₂ < 10 ns and is excluded.
- Added the null-model requirement, the size of the static effect, and the heterodyne-smearing caveat.

### Section VIII (old IX): Extensions
| Subsection | Change |
|---|---|
| Klein–Gordon | Uses a line clock with H_C = −P_t ≥ 0 (avoids commensurability). Notes this is the standard Newton–Wigner restriction and cites Hegerfeldt. "Resolves the KG negative-probability problem" removed. |
| Dirac | Contradiction resolved (M5): a Hardy clock removes Zitterbewegung; damping only for a two-sided clock under Assumption 1. Uses m²c⁴ (not m²c²). |
| Proper time | Kept. Two-level scope stated (m6). Revivals noted. |
| Static spacetimes | Ĥ_grav corrected to N^{1/2}√(c²P² + M²c⁴)N^{1/2} (M6). Twin-shell result attributed to Zych/Pikovski. Horizon rate marked model-dependent. Near-horizon scaling fixed: **∝ ρ_proper⁻¹ = (r − r_s)^{−1/2}**. "Thermal blueshift" wording removed. |
| ΔE\* (add-on) | Added as Remark 4. Γ_min corrected to 2^{1/4}\|ΔE\|\|ΔN\|^{1/2}/ħ. Caveats and lab numbers added. |
| Minisuperspace | WDW equation corrected: 1/a factor, sign of the V term, prefactor (M7). Exact case V = 0, p = 1 (Θ = −m₀⁻¹∂²_{ln a}), attributed to Ashtekar–Pawłowski–Singh 2006. "Tomonaga–Schwinger" removed. Hubble-variance ansatz moved to open problems. |
| Hartle–Hawking | "Theorem" → **Conjecture 2 (heuristic)**, with the reasons it is not yet a result (M8). |
| Near horizon | "Derivation" → "consistency with Bisognano–Wichmann/Sewell". K = cĤ/(ħκ) (not c³). Type III₁ caveat. Sewell 1982 cited. |

### Not carried over
- **Back-reaction power spectrum (IX.F add-on).** IR-divergent variance mislabelled; no derivation; the exponential factor is k-independent at horizon crossing (no tilt). Listed as an open problem.
- **Non-Abelian clocks (Deliverable 2).** Gauge covariance of the constraint not shown; Γ_YM has the wrong dimensions. Listed as an open problem.
- **Full Python listing in the appendix.** Replaced by a description of the tested package in `src/qtime`.

### Bibliography
- Fixed: Kuchař proceedings editor **G.** Kunstatter; Pauli Handbuch details.
- Added: Smith–Ahmadi 2019; Höhn–Smith–Lock 2021; Loveridge–Miyadera 2019; Busch–Grabowski–Lahti 1995; Egusquiza–Garay–Raya 1999; Zak 1967; Newton–Wigner 1949; Hegerfeldt 1974; Foldy–Wouthuysen 1950; Ashtekar–Pawłowski–Singh 2006; Bisognano–Wichmann 1975, 1976; Sewell 1982; Hawking 1975; Unruh 1976.

### Code
- `quadrature_time_simulation.py` → `src/qtime/clocks.py`. It now builds a real bipartite state, real POVMs and real partial traces. The original convolved hand-written states and mislabelled the readout resolution as "energy spread".
- `simulate_rel_clock_interferometer` → `src/qtime/interferometer.py`. Adds revival period and time-step helper; the plot shows revivals instead of implying monotone decay.
- New: `src/qtime/modular.py` (brute-force Tomita operator), `src/qtime/estimates.py` (all numbers quoted), `tests/` (23 tests).
- Notebook rebuilt (`notebooks/quadrature_time.ipynb`); the PDF-exported JSON was broken.
