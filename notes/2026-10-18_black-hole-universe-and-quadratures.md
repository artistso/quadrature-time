# The original picture, sorted

**Status: working note, not paper content.** Nothing here is a claim of the manuscript. It exists
so that the ideas from the 2026-10-18 conversation are recorded with their *status* attached:
what is already known, what is checkable today, and what still needs a definition before it can be
called physics. Every number below is reproducible from `notes/cosmology_numbers.py`.

---

## 1. The pieces, stated separately

The original proposal is not one claim; it is at least six, and they have very different status.

| # | The claim | Status |
|---|---|---|
| a | A clock is the representation of an oscillation: a particle is also a wave, and its phase is attached to the Hilbert-space lattice. | **Already the manuscript.** A massive particle carries an internal clock of frequency $\omega = mc^2/\hbar$; that is exactly the proper-time clock of §VIII.C, and the De Broglie/Compton clock of Zych–Pikovski–Brukner. Well founded. |
| b | The clouds that met at the Big Bang were so massive that we are *inside* a black hole, surrounded by larger ones. | **Established literature.** Pathria (1972); Smolin's black-hole selection (1994–97); Popławski's torsion bounce (2010); Gaztañaga's black-hole universe (2022), which fits the expansion history without dark energy. Not new — and *must not be presented as new*. |
| c | Two universes in a phase space of two intersecting geodesics, one sine and one cosine — the natural frequency of the universe from the Big Bang. | **Needs a decision.** See §4. As literally stated it is probably a *redundancy*, not two universes. |
| d | Expansion from *within* a collapsed object (supernova-like), not soup ejected into empty space; matter trapped, accumulating charge and energy, unable to escape the growing expansion. | **Partly standard, partly well-posed.** Bounce cosmology is standard. The "trapped radiation" part has a precise and beautiful cosmological counterpart (§3.2). |
| e | A 3D model of a photon, scaled up, looks like an atom; scaled again, an atom within an atom, trailing out very long. | **Metaphor, not yet physics.** The nearest real structure is discrete scale invariance / self-similar cosmological models (Oldershaw). Worth keeping as a *hypothesis generator*, not a claim. |
| f | A clock looks flat; seeing its actual length changes its shape; you cannot see both at once. | **Already formalised — in this manuscript.** And it now has a number (§3.3). |

---

## 2. The one thing that is exact, and it is a gift

Your intuition that a sufficiently massive region *is* a black hole is not numerology. At critical
density, the Schwarzschild radius of the mass inside a sphere equals the radius of that sphere —
**exactly**, for every radius:

$$M(R) = \rho_c \tfrac{4}{3}\pi R^3,\qquad
r_s = \frac{2GM}{c^2} = \frac{2G}{c^2}\cdot\frac{3H_0^2}{8\pi G}\cdot\frac{4}{3}\pi R^3
= \frac{H_0^2R^3}{c^2},$$

and at $R = R_H = c/H_0$ this gives $r_s = R_H$ identically. Numerically, with $H_0 = 70$ km/s/Mpc:

| $R$ | $M$ inside | $r_s$ | $r_s/R$ |
|---|---|---|---|
| $0.5\,R_H = 6.98$ Gly | $1.11\times10^{52}$ kg | 1.75 Gly | 0.250000000000 |
| $R_H = 13.97$ Gly | $8.90\times10^{52}$ kg | 13.97 Gly | **1.000000000000** |
| $2\,R_H = 27.9$ Gly | $7.12\times10^{53}$ kg | 111.7 Gly | 4.000000000000 |

That is the Pathria observation, and it is why "we live inside a black hole" keeps being taken
seriously: the Hubble sphere of a critical-density universe *is* its own Schwarzschild radius.

**But the honest caveat, which the popular versions always skip.** The observable universe is not
the Hubble sphere: it extends to about **46.5 Gly**, 3.3 $R_H$, because of what expansion did to the
light on its way to us. Using the standard mass $M_{\rm obs}\approx1.5\times10^{53}$ kg, the
Schwarzschild radius is **23.6 Gly** — comparable to, but only 0.51 of, the comoving observable
radius. So the correct statement is:

> the mass of the observable universe is within a factor of two of the mass that would make it a
> black hole, and the critical-density Hubble sphere is exactly at its own Schwarzschild radius.

That is a real coincidence, of order unity, and it is the one your picture rests on. It is *not*
the statement "we are provably inside a black hole", and a referee will make the distinction
whether or not the popular press does.

---

## 3. Three things that are checkable today

### 3.1 The photon's long climb out of the Sun, and the universe's

You said photons take a very long time to get from the core to the surface. That is right, and it
is pure random walk: $t \simeq R^2/(\lambda c)$.

| mean free path $\lambda$ | time to escape $R_\odot$ |
|---|---|
| 1 mm | 5.1 × 10⁴ yr |
| 1 cm | 5.1 × 10³ yr |
| 10 cm | 5.1 × 10² yr |

(light crosses the Sun in 2.3 s, so this is a factor of order $10^{8}$–$10^{10}$ — the literature
quotes 10⁴–10⁵ yr.) **The cosmological counterpart is exact and is usually missed:** the universe
was opaque from the Big Bang until **recombination at $t\approx3.8\times10^{5}$ yr, $z\approx1100$**.
The CMB is the universe's photosphere. The ratio of the two diffusion eras is about **3.8** — the
same order. Your "trapped until the surface lets go" picture is not an analogy for the early
universe; it *is* the early universe, and it is already measured.

### 3.2 The trade-off you described, in numbers

"The clock looks flat, but its real length is different, and you can't see both at once." In the
manuscript's language that is the circle ($S^1$, the flat face, one reading) versus its universal
cover ($\mathbb R$, the unrolled length, the winding number) — Fig. 1 and Proposition 2. But it is
*also* a quantitative trade-off that falls out of the finite-clock work:

**a clock with $N$ levels reads to $1/N$ of a cycle, and holds a target infidelity for a time
$T \propto N^2$.**

| $N$ | reading resolution $1/N$ | horizon $T$ ($\varepsilon=10^{-2}$) | $T\,(1/N)^2$ |
|---|---|---|---|
| 10 | 0.100 | 19 | 0.19 |
| 100 | 0.010 | 2344 | 0.23 |
| 1000 | 0.001 | 2.88 × 10⁵ | 0.29 |
| 3000 | 0.00033 | 2.00 × 10⁶ | 0.22 |

So **a ten-fold sharper reading costs a hundred-fold of horizon.** Sharpening the face leaves you
less of the length. That is your sentence, with an exponent on it, and it is already in the
manuscript as the $N^{-4}$ variance law.

### 3.3 Sine and cosine are a quarter period apart, not two universes

Computed directly for a single mode: apply a 90° phase shift and the quadratures exchange,
$\langle x\rangle \leftrightarrow \langle p\rangle$, while

$$\langle N\rangle = 1.250000000000 \quad\text{before and after.}$$

Same energy, same photon number, same entropy. And since $\cos\omega t = \sin(\omega t + \pi/2)$,
a "cosine universe" is the "sine universe" a quarter period later. In *your own manuscript* the
overall phase is the unobservable $U(1)$ fibre — so sine and cosine are **the same physical state
at two different times**, not two universes.

This matters because it is the difference between a redundancy and a claim. If the two "universes"
are the two quadratures of one mode, a referee will correctly say that you have described one
degree of freedom in two bases.

---

## 4. The decision that unblocks everything

"Two universes, in a phase space of two intersecting geodesics, one sine one cosine" has three
possible meanings, and they lead to completely different physics:

1. **Two quadratures of one oscillation.** Then it is a redundancy (above). Nothing to test; the
   statement is a change of basis.
2. **Two modes / two fields whose worldlines cross.** Then "intersecting geodesics" is a genuine
   event — a collision — and the crossing has consequences: interference, or a Bogoliubov mixing of
   the two, which is observable.
3. **Two observers or sectors in one curved spacetime.** Then it is real, and in principle
   measurable (this is the twin-shell interference of §VIII.D).

Pick one. On (1) the idea folds into the existing $U(1)$ picture; on (2) or (3) it becomes a model
you can write down and falsify.

---

## 5. What I would actually do next

**(i) Make the black-hole universe a constraint, not a picture.** The interior of a Schwarzschild
black hole is Kantowski–Sachs, not FRW. Impose your clock constraint on *that* minisuperspace and
ask what "trapped radiation building up" is. This is the same exercise as §IX.A (which already has
one exact case, Prop. 10) with a different metric — a weekend of algebra, and it either closes or
it doesn't.

**(ii) The one genuinely new-looking consequence.** In a closed, expanding universe the natural
frequency is set by the geometry: $\omega_k = c\sqrt{k^2+1}/a(t)$. But the whole commensurability
apparatus (Props. 1–2, Fig. 2) is built on a *fixed* clock quantum $\hbar\omega$. If $\omega$
falls as $a$ grows, then **commensurability drifts**: levels that were locked at early times
unlock later, and the winding-marginalization factor $D_M(\nu)$ with $\nu = \Delta E/\hbar\omega(t)$
sweeps through its zeros. That is a calculable, model-dependent, potentially observable statement
that no part of the current paper makes — and it comes from joining *your* two ideas (the
oscillation is the clock; the clock is in an expanding closed geometry). This is the thing I would
chase.

**(iii) Retire "atom within atom" or replace it.** The real structure behind self-similar scaling
is discrete scale invariance (log-periodic behaviour), which is testable; "an atom within an atom"
is not, and putting it in the paper would cost you credibility on the parts that are solid.

---

## 6. Cautions, in the order a referee would raise them

1. **Novelty.** (b) is 1972–2022 literature. Cite it or you lose the argument before it starts.
2. **"Mass velocity of the growing expansion."** Superluminal recession is not a speed violation —
   space expands, nothing moves through it locally. So "nothing can escape" has to be a *horizon*
   statement (a trapped surface, or the Hubble radius crossing), not a speed comparison.
3. **Two universes.** See §4 — as stated it is most likely a basis choice.
4. **The vortex drawing.** The composition is right; the notation is not yet defined. It needs
   $\Delta$, $S_N$ and the constraint named before it can be anything but a sketch.
