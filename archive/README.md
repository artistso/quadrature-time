# archive/

Reference material, kept for provenance but **not** part of the paper or the package.
Nothing here is imported, tested, or cited by `src/`, `tests/` or `paper/`.

## `chat-notes/`

Drafts, critiques and transcripts produced in chat sessions alongside the manuscript:

- **Derivation drafts** — optimal clock bandwidth ΔE\* (superscripted "Derive the exact
  optimal clock energy bandwidth ΔE…"), the inflationary back-reaction power spectrum,
  the Hartle–Hawking no-boundary proposal, the near-horizon limit, the stationary
  curved-spacetime extension, the non-Abelian clock proposal, the quantum-optics /
  circuit-QED design note, and a second Hilbert-Space Quadrature Conjecture extension.
  Each exists as a matched `.pdf` + `.txt` pair.
- **Critiques and audits** — "Perform a comprehensive structural review…", "Has anybody
  proposed this yet–", "How Wave Functions Generate Time.pdf" (the NotebookLM audio
  overview transcript), `notebooklm.txt` (its verbatim transcript), and the earlier
  "Revised Rigorous Architecture" notes.
- **Assembled drafts and exports** — "Provide the complete, final consolidated REVTeX
  4…", "Merge our 3-page paper and this Revised Rigorous Architecture…", "I say let's do
  all three", "Mr. Owens, we may proceed rigorously…", "Write a complete, formal 3-page
  foundational physics…", "Write out the formal mathematical proposition…", "Quadrature
  Time Simulation Notebook", "Write a Python simulation modeling a relativistic clock…".
- **Assets** — `Notes_261018_034424.jpg` (handwritten notes) and
  `The_Hilbert-Space_Quadrature_Conjecture.png` (summary graphic).

**Status:** superseded. Every claim in these files was triaged in `REVIEW.md`
(issues M1–M9, m1–m11) and the surviving, corrected versions live in `paper/`.
Three of them were found to be wrong and are *not* in the paper: the winding-periodization
proposition (replaced by the Zak transform), the Lindblad "theorem" (split into a
proposition plus an assumption), and the back-reaction "theorem" (moved to open
problems). The ΔE\* minimiser survives as Remark 4 but with a corrected minimum value.
Do not cite anything in this folder without checking `REVIEW.md` first.

## `bundles/`

`quadrature-time.bundle` — a complete git bundle of the repository at the revision-2
release (commit `569be92`, "Revised manuscript, referee review, and tested qtime
package"). It is a snapshot, not a source of truth: `main` has moved on. Kept only as a
recoverable copy of the pre-reorganisation tree.

## Note on the 2026-10-18 reorganisation

A GitHub web upload had flattened this repository: `paper/`, `src/qtime/`, `tests/`,
`scripts/`, `figures/` and `notebooks/` were dissolved, their contents dumped at the
repository root next to a partial copy of `.git/objects` (38-character hex filenames) and
these chat artifacts. File contents were byte-identical to the bundle, so nothing was
lost; the layout was restored from it and the strays were removed. See `CHANGES.md`,
"Repository maintenance".

## `sketches/`

Hand drawings made by the author while the finite-clock series was being written. They are raw
material, not claims; nothing here is cited by the paper.

**`2026-10-18_dst_seam.jpg` — "2 AM no longer exists."** Daylight saving time drawn as an
operator on clock labels: FALL and SPRING arrive at the work day, LATE and EARLY branch off it,
the labels 2 AM and 3 AM sit on either side of a line, and the pair
$\langle t_a e^{i\theta} \mid t_b e^{i\theta_2}\rangle$ is written underneath. Despite the
joke framing this is a correct and sharp statement about the manuscript's core topic, and it is
now **Figure 1** (`fig:cover`) with the surrounding paragraph in Sec. II.D:

- the reading map $t\mapsto\theta/2\pi$ is a covering of $\mathbb{R}$ onto $S^1$, so "2 AM"
  names four distinct instants in four days and needs the winding number $m$;
- a *global* relabeling is a rotation of the face — physically inert, no $m$ lost, undetectable;
- a relabeling applied *after a seam* (DST) is not a rotation but a cut of the cover: the seam
  hour is never displayed that day, so the clock shows 23 distinct hours, not 24.

The arithmetic is checked in `make_figures.py` (4 occurrences of 2 AM in 96 h; 23 labels on the
seam day; zero extra jumps under a global shift), so the caption cannot drift from the figure.

**`2026-10-18_time_dilation_vortex.jpg` — "slowest time / fastest spin".** A qualitative picture of
a vortex carrying a time-dilation spectrum: $\gamma$, $N$, $\Delta N$, entropy/action, a
$Δ \cdot 2\pi^2$ factor along $S_N$, a thermal $\pm\cos$ component, "fastest confit",
"longest time of vortex", and the whole thing repeated inside a circle labelled *Universe* — the
same structure at a larger scale. **Status: not used, and not yet physics.** It is a picture of
structure, in the spirit of §IX.C's open directions, and it would need to be turned into a
concrete claim (which sector? which clock constraint? what is the spectrum $S_N$?) before it could
appear anywhere. Kept because the composition is right even where the notation is not yet defined.
