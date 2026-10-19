# Commensurability drifts in an expanding universe

**Status: v10, first step. A negative result, and a small clean one.** Code:
`src/qtime/clock_design.py` (v10 block), tests: `tests/test_clock_design.py` (6 checks),
numbers reproducible from the same module.

## The idea

Everything in the manuscript assumes a fixed clock quantum $\hbar\omega$. But if the clock is the
oscillation of the universe itself, the geometry supplies the frequency, and in an expanding
universe that frequency is not fixed. On a round $S^3$ of curvature radius $a$ the scalar
eigenvalues are $(k^2-1)/a^2$, $k\ge2$, so

$$\omega_k(a) = \frac{c\sqrt{k^2-1}}{a}.$$

Because $\omega\propto a^{-1}$, the commensurability ratio $\nu=\Delta E/\hbar\omega$ that controls
the winding-marginalization factor $D_M(\nu)$ of Proposition 2 is **time dependent**:

$$\dot\nu = -\nu\,\frac{d\ln\omega}{dt} = +\nu H(t), \qquad\text{i.e.}\qquad \nu(t) = \nu_0\,a(t).$$

Commensurability, the condition that makes a cyclic clock a faithful timekeeper, is therefore not
a property of the spectrum alone in an expanding universe. It is a *state* that decays.

## The result, and it is negative

The coherence factor is

$$|D_M(\nu)| = \left|\frac{\sin\pi M\nu}{M\sin\pi\nu}\right|,$$

which is 1 at integer $\nu$ and zero at half-integers. Since $\nu$ grows like $a$, a coherence that
is exactly protected at early times is destroyed after an amount of expansion that does **not**
depend on the mode number, on $M$, or on the strength of the coupling:

| $\nu_0$ | e-folds until $\lvert\nu-\text{integer}\rvert>0.05$ | e-folds until $\lvert D_{200}\rvert<0.5$ |
|---|---|---|
| 1 | 0.049 | 0.003 |
| 2 | 0.025 | 0.002 |
| 3 | 0.017 | 0.001 |
| 10 | 0.005 | < 0.001 |

**Under a hundredth of an e-fold.** The reason is structural, not numerical: the drift is
*proportional* to $\nu$, so higher harmonics (which are the ones that carry more coherence weight
at late times) die faster, and the fraction of a winding period travelled is
$\ln a$ regardless of everything else. A universe that expands at all cannot protect an
incommensurate coherence by the mechanism of Proposition 2.

**So this is not the observable I hoped it would be.** The interesting question it was meant to
answer — "does the drift leave a signature?" — is answered in the negative at leading order: the
drift is so fast that the effect is a *sharp cutoff*, not a slow sweep. What survives is the
statement of when commensurability can be used at all: only within an era, never across one. That
is worth having, because it is exactly the kind of scope condition the manuscript insists on
elsewhere (cf. Assumption 1 and the dispersive limit $n_{\rm crit}$ of Sec. VI).

## What would make it non-trivial

1. **Squeezing.** In de Sitter, modes are squeezed, not merely redshifted; the two quadratures
   acquire different variances. Since the sine/cosine split is exactly a quadrature split (§5 of
   the cosmology note), a squeezing calculation is the right place for a two-component picture to
   earn its keep. That is derivable and does not need the winding machinery.
2. **A clock with a gap, not a mode.** If the clock quantum comes from a *gap* (a two-level atom,
   $\hbar\omega$ fixed by atomic physics) rather than from the geometry, then $\nu$ does **not**
   drift, and the manuscript's results stand unchanged. So the drift result is really a condition:
   it says which clocks can be relational clocks over cosmological times, and which cannot.
3. **Bounces.** Across a bounce $\omega$ diverges and the drift rate is momentarily undefined;
   a non-singular bounce (Popławski, Ashtekar–Pawłowski–Singh) is the only setting where "the
   frequency of the universe" can be tracked through the event at all.

## Honest summary

I proposed this as the thing worth chasing. It is worth having as a **scope condition**, and it is
not worth having as a prediction: the leading effect is that commensurability cannot survive
expansion, and no slowly-varying signature comes out of it. Recorded as a negative result rather
than quietly dropped.
