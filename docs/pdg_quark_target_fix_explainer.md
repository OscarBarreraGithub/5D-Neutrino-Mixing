# PDG Quark-Target Fix — A Physics Explainer

This note explains, in physics language rather than code language, a recent
change to the quark-sector acceptance gate of the parameter scan. It is
written for the project owner and for future collaborators who pick up the
quark side of the analysis cold. The intent is pedagogical: read it once,
top to bottom, and you should know what was wrong, what is now right, and
what is still on the to-do list.

## 1. The collaborator's complaint, in physical terms

The scan we run accepts a candidate $(r,\ \mathrm{overall\ scale},\ \Lambda_{\mathrm{IR}},\ c_Q,\ c_u,\ c_d)$
in two stages. First we ask: does this point reproduce the Standard Model
quark spectrum and the CKM matrix? Only if it passes that "spectrum gate"
do we then evaluate the $\Delta F = 2$ Wilson coefficients
($\varepsilon_K$, $\Delta M_K$, $\Delta M_{B_d}$, $\Delta M_{B_s}$,
$\Delta M_{D^0}$) to see whether the KK-gluon contribution sits inside the
experimental budget.

The complaint was about that first gate, and it was a fair one. The old
gate compared the model's quark masses to a hard-coded target table whose
provenance was not traceable to PDG, declared to live "at $\mu = 3\,\mathrm{TeV}$,"
and applied a single uniform $10\%$ log-residual tolerance to every quark.
That is wrong on three counts.

The first count is the most physically alarming. One accepted point gave
$m_s \approx 0.065\,\mathrm{GeV}$ at the gate's reference scale, whereas
the PDG $\overline{\mathrm{MS}}$ value for the strange quark, evolved to
the same scale, is around $0.055\,\mathrm{GeV}$. That is an 18% deviation
in $m_s$. PDG knows $m_s$ to better than 1%. Why does this matter for
$\Delta F = 2$? Because the strange-quark mass enters $\varepsilon_K$
through the chiral enhancement factor

$$ r_\chi \;=\; \left( \frac{m_K}{m_s + m_d} \right)^{\!2} \;\approx\; 25, $$

so a fractional error in $m_s$ propagates approximately linearly into the
predicted $\varepsilon_K$. An 18% error in $m_s$ is, to first
approximation, an 18% error in the very observable whose tightness drives
most of the new-physics constraint. The acceptance gate was, in other
words, allowing into the survivor set points that misrepresented
$\varepsilon_K$ by more than the experimental tightness of $\varepsilon_K$
itself. The downstream conclusions about which $r$ and which localization
geometries pass the FCNC budget would be biased by this leak.

The second count is the use of "$\mu = 3\,\mathrm{TeV}$" as the reference
scale for the mass targets. The 3 TeV scale really does belong somewhere in
this code — it is the reference at which we anchor $\alpha_s$ for running
the $\Delta F = 2$ Wilson coefficients down to hadronic scales. But that
is a *Wilson-coefficient* concern. There is no PDG table of quark masses
at 3 TeV; PDG quotes light-quark masses at $\mu = 2\,\mathrm{GeV}$, and
the heavy quarks at their own self-consistent scales $m_q(m_q)$. Pretending
that the hard-coded targets were "at 3 TeV" elided the running step that
should have been done explicitly.

The third count is that one tolerance does not fit all six quarks. PDG
knows $m_b$ to roughly 0.2%; it knows $m_u$ to roughly 17%. Demanding 10%
agreement on $m_u$ is meaningless (the input is fuzzier than the gate),
while demanding 10% agreement on $m_b$ is comically loose (the gate is
much fuzzier than the input). A single threshold cannot simultaneously
reflect the actual experimental knowledge across the spectrum.

## 2. What "PDG $\overline{\mathrm{MS}}$ at a consistent scale" means

A quark mass is not a physical observable in the way the proton mass is.
It is a renormalized parameter of the QCD Lagrangian, and its numerical
value depends on what scheme you renormalize in and at what scale you
quote the result. PDG quotes $\overline{\mathrm{MS}}$ ("MS-bar") masses,
which is the convention universally used in flavor phenomenology. The
$\overline{\mathrm{MS}}$ mass runs with scale according to the QCD
renormalization group:

$$
\mu^2\, \frac{d\, m_q}{d\,\mu^2} \;=\; -\, \gamma_m\!\bigl(\alpha_s(\mu)\bigr)\, m_q ,
$$

where $\gamma_m$ is the mass anomalous dimension, known through four loops.
The qualitative behavior is that $m_q(\mu)$ shrinks logarithmically as
$\mu$ increases. For example, $m_b(m_b) \approx 4.18\,\mathrm{GeV}$, but
$m_b$ evaluated at the top-mass scale is about $2.74\,\mathrm{GeV}$ —
the same physical quark, the same Lagrangian parameter, just expressed at
a different scale.

The other subtlety is *flavor-threshold matching*. The number of active
quark flavors $n_f$ in the running depends on which quarks are
"integrated in" at a given scale. Below $m_b$ the running uses $n_f = 4$;
above it uses $n_f = 5$. Crossing the threshold requires a finite matching
relation that ties $m_q(m_b^-,\ n_f=4)$ to $m_q(m_b^+,\ n_f=5)$, and
similarly for $m_c$. The repo now uses Chetyrkin–Kniehl–Steinhauser-style
matching at the heavy-quark thresholds, evaluated to two or three loops.

"Consistent scale" simply means: pick *one* renormalization scale for the
whole acceptance comparison, evolve every PDG input to that scale before
comparing it to anything the model produces, and evolve the model's
outputs to the same scale. Mixing scales — comparing the model's
$m_b$-at-3-TeV to PDG's $m_b(m_b)$ without running, for instance — is the
mistake the old gate made implicitly.

## 3. Why $\mu_{\mathrm{common}} = m_t(m_t) \approx 163.5\,\mathrm{GeV}$

Three candidate scales are reasonable on their face: the low scale
$\mu = 2\,\mathrm{GeV}$ where PDG quotes the light quarks; the
electroweak scale set by the top quark, $m_t(m_t)$; and the KK scale
$\mu = 3\,\mathrm{TeV}$ used elsewhere in this code.

The low scale is rejected because it is *below* both heavy-quark
thresholds. To compare the bottom quark there we would have to integrate
$m_b$ all the way down through two thresholds, which adds avoidable
matching error and forces a non-perturbative-leaning $\alpha_s$ into the
mass-running calculation. There is no physics gain.

The 3 TeV scale is rejected for a different reason. The model's quark
masses are generated by overlap integrals with an electroweak-scale Higgs;
they are physically defined "at the EW scale" in the sense that this is
the scale of the matching between the 5D theory and the SM effective
theory. Quoting them at 3 TeV requires extrapolating QCD running an extra
three decades for no reason. Worse, it invites confusion with the Wilson
coefficient anchor (see §6 below).

The top-mass scale $m_t(m_t)$ threads the needle. It sits comfortably
above all the QCD thresholds, so the running is fully perturbative; it is
a natural "edge of the SM" scale where electroweak matching is conducted
in standard-model phenomenology; and it is far enough from 3 TeV that
nobody will accidentally conflate the mass-target scale with the
Wilson-coefficient anchor. The repo therefore uses
$\mu_{\mathrm{common}} = m_t(m_t) = 163.5\,\mathrm{GeV}$ for *mass targets*,
identified in code as `qcd.constants.M_TOP_MS`.

## 4. Why per-quark 2σ tolerances, with numbers

The right way to ask whether a model point matches PDG is to ask whether
its predicted $m_q$ falls inside the PDG uncertainty band, *quark by
quark*. PDG's relative uncertainties are nearly invariant under
one-loop multiplicative running, so the relative-2σ tolerances are quoted
once and reused at any common scale. The numbers in use, all expressed at
$\mu_{\mathrm{common}}$:

| quark | central at $\mu_{\mathrm{common}}$ | $2\sigma$ relative |
|-------|------------------------------------|--------------------|
| u     | $\sim 0.0012$ GeV                  | $0.45$             |
| d     | $\sim 0.0028$ GeV                  | $0.18$             |
| s     | $\sim 0.055$ GeV                   | $0.017$            |
| c     | $\sim 0.626$ GeV                   | $0.0072$           |
| b     | $\sim 2.74$ GeV                    | $0.0033$           |
| t     | $162.5$ GeV → $163.5$              | $0.0086$           |

In addition there is a numerical floor of $0.3\%$ per quark to absorb
optimizer round-off; nothing is gated tighter than this even where PDG
itself is.

The contrast with the old global $10\%$ is dramatic for the quarks where
PDG is sharp. The strange tolerance has tightened by a factor of six; the
bottom tolerance by a factor of thirty. Conversely, $m_u$ and $m_d$ have
*loosened*, which is correct: PDG is the entity with the noisy
measurement, and demanding sub-PDG agreement there was meaningless.

## 5. CKM: a single scale-free table

For the CKM matrix the situation is simpler. Between the electroweak scale
and the hadronic scales relevant to $\Delta F = 2$ (any
$\mu_{\mathrm{had}} \gtrsim m_W$), the CKM elements are scale-free at
one-loop QCD. This is the standard PDG convention: there is no QCD
running that mixes flavor labels in the unitary CKM matrix. We therefore
freeze CKM at $m_Z$, use it directly at $\mu_{\mathrm{common}}$ for the
gate, and use the same numbers when feeding $\Delta F = 2$ matching at
$\mu_{\mathrm{had}}$. Per-element 2σ tolerances follow PDG 2024 directly:
$|V_{us}|$ at $0.4\%$, $|V_{cb}|$ at $4\%$, $|V_{ub}|$ at $8\%$, and the
Jarlskog invariant $J$ at $\sim 6\%$.

## 6. Two scales, kept orthogonal on purpose

This is the single most important conceptual point of the change.

There are now two distinct reference scales that live in the code, and
they answer two different physics questions:

- $\mu_{\mathrm{common}} = 163.5\,\mathrm{GeV}$ is the **mass-target
  scale**. It exists only to compare quark masses between PDG and the
  model's prediction.
- $\mu_{\mathrm{ref}} = 3\,\mathrm{TeV}$ is the **Wilson-coefficient
  reference scale**. It is the scale at which the KK-gluon-induced
  $\Delta F = 2$ Wilson coefficients are matched and at which $\alpha_s$
  is anchored for the subsequent RG flow down to $\mu_{\mathrm{had}}$.

These are not redundant; they are not aliases; and the change explicitly
preserves the second while overhauling the first. A registry-level test
in the code asserts that the WC reference scale is still 3 TeV after the
change — a regression guard against a future agent accidentally
"cleaning up" the apparent duplication. Conflating them would be subtle
to spot and disastrous for the FCNC predictions, because the KK-gluon
amplitude carries an $\alpha_s$ at the matching scale and the
short-distance Wilson coefficients run between $\mu_{\mathrm{ref}}$ and
$\mu_{\mathrm{had}}$ over almost three decades.

## 7. Plan B and the audit gate

The new tolerances are tighter than the old ones — substantially tighter
for $s$, $c$, and $b$. There is a real possibility that the production
re-scan will return zero accepted points. The natural reflex is to
loosen tolerances. We are explicitly forbidding the silent version of
this reflex.

Instead, before any tolerance is relaxed beyond 2σ, the code requires a
per-quark *median residual audit*. For each quark $Q$, the median
log-residual across the production scan is computed and committed to
disk. If $Q$'s median is well below 2σ, the rejection of a marginal
point is a *fit-precision* issue and modestly relaxing the floor is
defensible. If $Q$'s median sits *at* or *above* 2σ, that is *physics
tension*: the model genuinely cannot produce $Q$ at the PDG value with
the localization machinery as currently constrained, and any relaxation
must be accompanied by a written paragraph in
`docs/quark_scan_assumptions_compact.tex` saying so. A null result, in
other words, would itself be an interesting physics result, and we want
to be in a position to report it as one rather than to hide it by
loosening the gate.

## 8. What is implemented and what is deferred

Implemented now: PDG 2024 inputs, the mass-running module with threshold
matching, the 2σ-per-quark and per-CKM-element tolerances wired into the
acceptance classifier, the orthogonal-scale guard in the modern-input
registry, the assumptions-document update, and tests for round-trip
consistency, threshold continuity, charm-decoupling, and a synthetic
regression fixture that pins the $m_s = 0.065\,\mathrm{GeV}$ failure
mode.

Deferred to compute jobs: the Phase-0 calibration scan, which would
re-run the optimizer with the new objective but with the acceptance gate
*disabled*, and would produce empirical 95th-percentile log-residuals
that could replace the 0.3% floor with a measured one. Until that job is
run, the production tolerance vector is the deterministic fallback —
PDG 2σ per quark, with a flat $0.3\%$ floor. The full production re-scan
under the new gate is also deferred to the user.

## 9. Worked example: the strange quark

PDG 2024 quotes
$m_s(\mu=2\,\mathrm{GeV},\ n_f=4) = 93.5 \pm 0.8\,\mathrm{MeV}$.
The walk to the common scale goes:

1. Start at $\mu = 2\,\mathrm{GeV}$ in $n_f = 4$ with the PDG central
   value, $93.5\,\mathrm{MeV}$.
2. Run upward in $n_f = 4$ to the bottom threshold,
   $\mu = m_b \approx 4.18\,\mathrm{GeV}$. The mass shrinks because
   $\gamma_m$ is positive and $\alpha_s$ is largish here.
3. Apply finite matching at $m_b$ to switch $n_f: 4 \to 5$.
4. Run upward in $n_f = 5$ from $m_b$ to $163.5\,\mathrm{GeV}$, where
   the running slows because $\alpha_s$ is smaller.
5. Result: $m_s(163.5\,\mathrm{GeV}) \approx 0.055\,\mathrm{GeV}$.

The 2σ tolerance is the PDG relative uncertainty,
$2 \times 0.8/93.5 \approx 1.7\%$, applied to the central at the common
scale, giving an acceptance window of roughly
$0.054 \le m_s(163.5\,\mathrm{GeV}) \le 0.056\,\mathrm{GeV}$.

The offending accepted point had $m_s \approx 0.065\,\mathrm{GeV}$ at
the gate's reference scale — about $18\%$ above the central. That is
roughly $11\sigma$ in PDG units. With the new gate it is decisively
rejected, exactly as it should be, and a synthetic test fixture pins
this behavior so the protection cannot regress silently.

## 10. Open questions worth a second look

Three things are worth thinking about once Phase-0 and the full re-scan
have run.

The first is geometric. The bottom quark's tolerance is now $0.3\%$,
which is genuinely tight. If the new accepted set is heavily depopulated
of points whose limiting failure mode is $m_b$, that is informative:
the geometry has trouble producing $m_b$ at PDG precision, and one
might ask whether a particular range of $c_b$ (or of the
top-bottom localization split) is being squeezed out. That would be a
physics result, not just a fit-quality complaint.

The second is whether the new acceptance set still spans the
phenomenologically interesting range of the warp-ratio parameter $r$ and
the overall scale. If the accepted manifold collapses to a thin sliver,
the FCNC conclusions become statistically fragile and the scan grid may
need refinement near that sliver.

The third is the interaction between the spectrum gate and the
$\Delta F = 2$ gate. With sharper $m_s$, $m_c$, $m_b$ in particular, the
KK-gluon contributions to $\varepsilon_K$ and $\Delta M_{B_{d,s}}$ are
now evaluated against more reliable Standard Model inputs, so the
implied bound on the KK-gluon mass and on the Yukawa-overlap structure
should be quantitatively stronger. It will be worth comparing the
post-fix bound on $\Lambda_{\mathrm{IR}}$ to the pre-fix bound — not to
"check" anything, but because that comparison is itself a measurement of
how much the leakage in the old gate was diluting the FCNC constraint.
