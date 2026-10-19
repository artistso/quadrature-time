# Change Log

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
