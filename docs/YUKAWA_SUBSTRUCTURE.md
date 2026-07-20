# Yukawa substructure: what we asked, what we found, what it means

**Status:** explainer, written 2026-07-20. Summarizes the July 2026 Yukawa
perturbation study (the "F7" study), the independent per-element anatomy
result, the limitations identified by the 2026-07-15 audit, and the full
study the results motivate. Companion decisions live in
`docs/OPEN_QUESTIONS.md` (item 8).

---

## 1. The question

The original ask was: do the 5D Yukawa matrices of surviving scan points have
substructure? Concretely, four parts:

1. Perturb each Yukawa entry and find how far it can move before the point
   stops fitting, including which observable fails first and how badly.
2. Add random sub-percent Gaussian noise to the complete Yukawa matrices and
   measure survival.
3. Treat the local response as a point-dependent linear map on Yukawa space.
4. Study the geometry of those maps: do the allowed points lie on a
   lower-dimensional submanifold, and if so, does that manifold correspond to
   a recognizable flavor symmetry?

The physics stake: a point that survives epsilon_K by an accidental phase
cancellation and a point that survives because of a structural mechanism (an
approximate symmetry, an alignment) can sit at the same place in an exclusion
plot. They are physically very different models. Local response geometry is
how you tell them apart.

## 2. Background: why Yukawa structure controls epsilon_K in RS

In the warped setup, 4D Yukawa couplings are products of anarchic O(1)
5D Yukawas and geometric overlap factors, y_ij = Y5_ij f_Qi f_qj, so the
observed mass and CKM hierarchies come from the f factors while the 5D
entries stay unstructured. The same overlaps produce flavor-violating
KK-gluon couplings. For epsilon_K the dangerous object is the imaginary part
of the LR Wilson coefficient C4 (chirally and RG enhanced), roughly
Im C4 ~ (g_s*^2 / M_KK^2) (f_Q1 f_d1)(f_Q2 f_d2) x O(1) phases.
Whether a given point passes therefore depends on the detailed magnitudes
and phases of the down-sector Yukawa entries: exactly the substructure this
study probes. Conventions and lanes: `docs/MODEL_CONVENTIONS.md`.

## 3. What was actually done (the F7 study)

Code: `scripts/yukawa_perturbation_study.py` (per-class RNG seeding made
process-stable on 2026-07-20 via `scripts/reproducible_seeds.py`). All runs
at M_KK = 3 TeV under the legacy perturbative coupling convention.

Four base-point classes were selected:

| Class | Construction | Why it survives epsilon_K |
|---|---|---|
| phase-tuned anarchy | largest-abs(C4) anarchic survivor found in a finite search (`scripts/instrument_epsK_phase.py`) | accidental phase cancellation in Im C4 |
| typical anarchy | magnitude-suppressed anarchic survivor (`scripts/anarchic_bauer_s1.py` machinery) | small abs(C4), no tuning |
| rank-one / U(2) toy | structured Yukawas with a U(2)-like texture (`scripts/rankone_u2_lane.py`) | partial directional protection |
| real-Y_d toy ("Nelson-Barr") | Y_d constrained real | Im C4 = 0 by construction |

Three experiments were run on each:

1. **Multiplicative real noise** on Y_d entries at increasing sigma;
   survival fraction of the epsilon_K gate recorded.
2. **Finite-difference gradient** of Im C4 with respect to the 18 real
   components of Y_d; from it a local "tuning radius" (how far you can move
   along the worst direction before the gate fails).
3. **PCA/SVD of normalized gradient directions** over small ensembles per
   class, to count how many independent dangerous directions exist.

## 4. Results

| Class | Noise survival | Local tuning radius | Reading |
|---|---|---|---|
| phase-tuned anarchy | below 50 percent survival already at sigma about 5e-3; median failure reaches about 36x the bound at large noise | 7.9e-4 | extremely fragile; the cancellation is razor thin |
| typical anarchy | robust at percent level; about 55 percent at 10 percent noise | 0.15 | protected by small magnitudes, not by tuning |
| rank-one / U(2) toy | 100 percent through 3 percent noise, 85 percent at 10, 63 percent at 30 | 0.024 | magnitude suppression plus partial directional protection |
| real-Y_d toy | 100 percent under all real multiplicative noise | infinite in real directions | exact by construction: real noise preserves real Y_d |

The cleanest geometric result is the gradient ensemble of the real-Y_d
class: the gradients of Im C4 occupy exactly the nine imaginary directions
of Y_d, with a sharp singular-value cliff after component nine. That is a
direct numerical detection of a nine-real-dimensional protected flat: the
reality condition is visible in the response geometry without being put in
by hand. The U(2) toy shows a softer, mostly imaginary-sensitive spectrum
with no sharp cliff, consistent with partial rather than exact protection.

Headline: the original intuition is supported. A point that survives by
phase cancellation is destroyed by sub-percent Yukawa noise, while a point
protected by a structural condition is insensitive to perturbations that
respect the condition. Fragility under noise is a practical symmetry
detector.

## 5. An independent clue: the down sector moves collectively

A separate utility, `scripts/yukawa_per_element_anatomy.py`, asks a related
question of the old 83,961-point optimizer sample
(`data/accepted_points_with_yukawas.csv`, gitignored but preserved locally
with provenance): do accepted points suppress the six off-diagonal Y_d
entries with one collective dial, or tune each entry independently?

- Y_u off-diagonal within-point log spread: median 0.534 dex
  (95th percentile 0.630 dex).
- Y_d off-diagonal within-point log spread: median 0.148 dex
  (95th percentile 0.339 dex).

The down-sector off-diagonals of accepted points move together about 3.6x
more tightly than the up-sector ones. That is what you would expect if one
approximate low-dimensional suppression spurion controls the down sector
(the sector epsilon_K cares about), rather than six unrelated tunings. It
is a hypothesis generator, not a result: the sample is a pre-audit
existence/envelope fit, and the statement is basis-dependent.

## 6. What these results do NOT establish (audit findings)

The 2026-07-15 audit (section 8) reviewed the study and classified it as a
useful proof of concept that does not yet answer the original question. The
substantive limitations, distilled:

1. **Frozen background.** The bulk profiles f_Q, f_u, f_d were held fixed
   while Yukawas moved. The study measures the frozen local response, not
   the distance to the refitted solution set. "Does it still fit?" requires
   re-profiling the nuisance and bulk parameters after each perturbation.
2. **Perturbation class too narrow.** Real multiplicative noise on Y_d only:
   it cannot explore generic CP phases and barely moves near-zero entries.
   The real-Y_d immunity is therefore exact by construction for this noise
   class (a consistency check, not a discovery), and the study is not a
   demonstration of a complete Nelson-Barr model.
3. **One observable.** Only Im C4 was differentiated. Masses, CKM, Delta m_K,
   B and D mixing, EW, and perturbativity gates are absent from the gradient;
   the pass count used only the epsilon_K ratio.
4. **PCA is not manifold dimension.** The reported "effective dimension" is
   the span of gradient covectors across an ensemble. A single scalar
   constraint generically has a 17-dimensional pointwise kernel in the
   18-dimensional real Y_d space; calling the tuned point "isolated" is too
   strong. The interesting object is how that kernel intersects the
   mass/CKM tangent space and the other active constraint normals.
5. **Selection and convention caveats.** The phase-tuned point is a
   best-of-budget selection; results are at M_KK = 3 TeV under the legacy
   coupling; step sizes and Euclidean norms are coordinate-dependent; quark
   rephasing/weak-basis redundancy was not quotiented. (The seed
   reproducibility defect flagged by the audit was fixed on 2026-07-20.)

None of the numerical radii above are release results. They are internal
diagnostics under a stated convention.

## 7. The study that answers the question

Specified in full in section 9 of
`docs/REPOSITORY_AUDIT_AND_RESEARCH_ROADMAP_2026-07-15.md`. In brief:

- Work in the 36-real-dimensional space (Re Y_u, Im Y_u, Re Y_d, Im Y_d)
  plus nuisance/profile parameters, with the rephasing redundancy fixed or
  quotiented.
- Define a whitened residual vector over ALL fit targets and constraints
  (masses, CKM including the CP phase, epsilon_K, Delta m_K, B_d, B_s, D,
  EW, perturbativity), each normalized by a named uncertainty.
- Report two sensitivities per point: the frozen response (mechanism
  diagnostic) and the profiled response after re-minimizing nuisances (the
  actual "does it still fit" answer).
- Bracket-and-bisect per-entry failure distances in both additive and
  fractional directions, both signs, with the failing observable and its
  pull recorded; random noise with real and complex Gaussian ensembles,
  quoting survival distributions over many independently accepted points
  rather than one best point.
- Build the Jacobian of the residual map, whiten, eliminate nuisance
  directions (Schur complement), and SVD: stiff directions, soft/tangent
  directions, effective codimension, and the observable composition of each
  normal. At constraint boundaries use tangent cones (the boundary is a
  stratified set, not one smooth manifold).
- Beyond linear response: directional curvature along soft directions,
  predictor-corrector continuation, principal angles between local soft
  subspaces across points, and viable tube volume as a naturalness measure.
  The information-geometry language (stiff/sloppy "hyperribbons") is
  established in Quinn et al., arXiv:2111.07176; the new content would be
  the gauge-invariant application to warped flavor and the recovery of
  protecting spurions from data.

## 8. Why this matters for the paper

The proposed paper direction (audit section 10, directions A, B, C) turns
this machinery into the central result: infer the approximate protecting
generators (CP, U(2), FPR/MFV spurions) from Jacobian kernels of surviving
points instead of assuming them; measure how fast controlled
symmetry-breaking spurions regenerate the LR kaon operator in the exact FPR
limit; and produce a naturalness phase diagram that separates decoupling,
magnitude suppression, and phase cancellation. The four point classes above
already show, at proof-of-concept level, that these mechanisms have
quantitatively different response signatures. The down-sector collectivity
of section 5, if it survives basis fixing and profiling on a corrected
ensemble, is the data-driven bridge from the scan to an analytic flavor
ansatz.

## 9. Where everything lives

| Item | Path |
|---|---|
| Perturbation study driver (parts A/B/C) | `scripts/yukawa_perturbation_study.py` |
| Stable per-class seeding | `scripts/reproducible_seeds.py` |
| Phase-tuned point construction | `scripts/instrument_epsK_phase.py` |
| Anarchic baseline machinery | `scripts/anarchic_bauer_s1.py`, `scripts/run_rs_anarchy.py` |
| Rank-one / U(2) toy | `scripts/rankone_u2_lane.py` |
| Per-element anatomy of accepted points | `scripts/yukawa_per_element_anatomy.py` |
| Accepted-point sample (local, gitignored) | `data/accepted_points_with_yukawas.csv` + provenance JSON |
| Earlier single-entry stability probe | `scripts/perturb_yukawa_stability.py` |
| Audit review of the study | `docs/REPOSITORY_AUDIT_AND_RESEARCH_ROADMAP_2026-07-15.md`, sections 8-9 |
