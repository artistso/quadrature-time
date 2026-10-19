# The PHOTON-E model: what it gets right, and the two things to fix before it goes in a paper

**Status: review of the author's own visualisation, not paper content.** Screenshots archived in
`archive/sketches/photon_app/`. Numbers reproducible from `notes/photon_numbers.py`.

---

## 1. What the images show, in the manuscript's language

All seven frames are the *same* object seen from different directions, and that is the point.

| frame | what it is | manuscript's name for it |
|---|---|---|
| 00, 01 | a two-strand helix along the propagation axis, red and blue | the **E-field** of a circular/elliptical polarisation, drawn as a curve in space |
| 02 | the **axis-on** view: concentric rings | the helix seen along its own axis — **the staircase from above** (Fig. 3c) |
| 03, 04 | a flat ellipse, face-on and oblique | elliptical polarisation: the helix squashed by unequal amplitudes |
| 05, 06 | a crossed pair of ellipses, two orientations | the **same** ellipse rotated, i.e. a change of phase convention |

So the app is already an instrument for the one operation the paper is about: **the same state,
projected differently.** Frames 00→02 are exactly the walk→look-down transformation of the
staircase. That correspondence is not a coincidence and it is worth saying out loud: *the axis-on
view of a helix is where the winding collapses into radius.*

**And that is why it looked like an atom.** Look down the axis of a helix and you get rings; rings
read as orbitals. The "atom within an atom" you saw was your own clock, end-on. Nothing mystical —
but it is the *same* observation that motivates discrete scale invariance, and it is the one place
where the "nesting" has real mathematics behind it (see §4).

## 2. The numbers: mostly exactly right

The energy readout is **correct on all four frames to four significant figures** — verified against
$E=hc/\lambda$ (280 nm → 7.0944×10⁻¹⁹ J vs 7.09 shown; 507 → 3.9180 vs 3.92; 720 → 2.7590 vs
2.76; 370 → 5.3688 vs 5.36). The Stokes panel with $S_1,S_2,S_3$ and a degree of polarisation is
the *right* way to represent this — it is standard, and it is the correct generalisation of the
sine/cosine split.

**One readout to check.** The third figure, labelled in Hz, does not appear to equal
$\nu=c/\lambda$ in every frame. For reference, the values it *should* show are:

| $\lambda$ (nm) | $E=hc/\lambda$ (J) | $\nu=c/\lambda$ (Hz) | $\omega$ (rad/s) |
|---|---|---|---|
| 280 | 7.0944×10⁻¹⁹ | 1.0707×10¹⁵ | 6.7273×10¹⁵ |
| 370 | 5.3688×10⁻¹⁹ | 8.1025×10¹⁴ | 5.0910×10¹⁵ |
| 507 | 3.9180×10⁻¹⁹ | 5.9131×10¹⁴ | 3.7153×10¹⁵ |
| 720 | 2.7590×10⁻¹⁹ | 4.1638×10¹⁴ | 2.6162×10¹⁵ |

If that field is meant to be $\nu$, it is a units slip; if it is meant to be something else
(a rate of some displayed process), it needs its own name — "probability" in Hz is not a quantity.

## 3. The one thing a referee will say, and the answer

**A single photon has no helix.** The helix is $\langle \hat{\mathbf E}(t)\rangle$, the expectation
value of the field in a coherent state — a *classical* field with $n\gg1$ photons. For one photon
there is no definite field trajectory at all; what exists is a ray in Hilbert space and a
probability distribution of outcomes. So the model is a picture of **the classical limit of the
field the photon occupies**, which is fine and useful, but it must not be captioned "a photon".

The good news: this is *not* a problem for your thesis, it is your thesis. Because what the helix
correctly displays is the thing that *is* observable — the relative phase turning at $\omega$ — and
what it hides is the thing that is *not* observable — the overall phase, i.e. the $U(1)$ fibre of
Eq. (1) of the manuscript. **The helix is one representative of a gauge orbit; the physical state is
the ray.** Draw it that way and the picture becomes an illustration of the manuscript rather than a
target for it.

Same point for the helicity slider. Helicity of a massless spin-1 quantum is $\pm1$ **only**; there
is no zero-helicity photon. A readout of −0.38 or +0.62 is $\langle\hat h\rangle = S_3/S_0$, an
expectation value over a superposition or an ensemble — legitimate and useful, but every single
*measurement* returns ±1. Label it $\langle\text{helicity}\rangle$ (or "circular polarisation
$V/I$") and the objection disappears. This is the same distinction as §4 of the cosmology note,
where sine and cosine turned out to be one state in two bases.

## 4. The synthesis worth keeping: the Poincaré sphere *is* the stair-time picture

This is the part I would build on, because it settles the sine/cosine question from the last
note in a way that is both correct and visual.

Put the state on the **Poincaré sphere**. Its equator $\{S_1,S_2\}$ is linear polarisation — the
plane in which sine and cosine live — and its poles are the two helicities $S_3=\pm1$. Then:

- a 90° phase shift is a rotation **about the $S_3$ axis**, i.e. a movement *along the equator*:
  it takes $\sin$ to $\cos$ without leaving the equator;
- therefore "sine universe" and "cosine universe" are **two points on the same equator**, related
  by that rotation — one degree of freedom, two labels. This is exactly the computed result
  ($\langle N\rangle$ identical, quadratures exchanged);
- and helicity is **perpendicular** to that plane. So a "two universes" claim that has any chance
  of being physics has to put them at different $S_3$ — different *handedness* — not at
  different phase.

That last sentence is the useful one. If you want two things where the paper currently has one,
the axis to separate them along is helicity, not phase. That is a concrete, checkable proposal, and
your own app already displays the axis.

## 5. What would make the model more than a visualisation

1. **Give the helix a label it earns**: $\langle\hat{\mathbf E}\rangle$ of a coherent state with
   mean photon number $\bar n$, and show how the drawn curve degrades into a probability cloud as
   $\bar n\to1$. That is a real, honest transition, and it is the same limit in which the
   "particle" language is legitimate.
2. **Make "collapse wavefunction" a projection, not a collapse.** A button that projects onto a
   Stokes axis and shows the *result* $\pm\hbar$ along that axis would be both correct and
   pedagogically strong, and it would make the $\pm1$ helicity point self-evident.
3. **The ring spacing in frame 02 is a testable detail.** If the radii are evenly spaced the
   underlying helix is uniform; if they are *geometrically* spaced the object is a logarithmic
   spiral — the conical staircase of Fig. 4, and the only version of "an atom within an atom" that
   has real mathematics (discrete scale invariance / log-periodic self-similarity) behind it.
   Worth checking against your own code, because it decides whether the nesting is structure or
   perspective.
