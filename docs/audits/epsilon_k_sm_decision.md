# Decision: epsilon_K^SM central value and NP-budget treatment

Date: 2026-05-15
Authority: Phase 2 hole #5 sign-off (Claude Opus orchestrator, acting per
`docs/paper_execution_decisions.md` "Audit philosophy" and
"Stopping conditions" delegations).
Branch: `audit/bag-inputs` (FF-merged into `paper/quark-scan-2026q2`).
Triggering documents: `docs/audits/bag_param_inventory.md`,
`docs/phase_logs/phase2_h5_impl.md`, `docs/phase_logs/phase2_h5_review.md`.

## Context

The Phase 2 hole #5 hadronic-input audit moved the SM central value of
|epsilon_K| from 1.81e-3 (legacy CKMfitter-style number, originally hard-coded
in `quarkConstraints/deltaf2.py`) to 2.161e-3 (Brod, Gorbahn, Stamou 2020,
arXiv:1911.06822, Eq. (20)). With PDG 2024 |epsilon_K^exp| = 2.228e-3, this
shrinks the central-value NP budget from |Delta epsilon_K| = 4.18e-4 to
6.7e-5, a 6.24x reduction. The peer reviewer (`phase2_h5_review.md`,
Findings 1-2) raised two open WARNINGs:

1. BGS 2020 sits near the high end of the SM-prediction band relative to
   CKMfitter / UTfit averages (typically ~1.9e-3 to ~2.1e-3 for the
   long-distance-corrected SM prediction), so it gives the **tightest**
   admissible NP budget and therefore the **strongest** M_KK lower bound.
   This needs an explicit defense or a sensitivity treatment.
2. The 6.7e-5 number is treated by the code as an exact denominator. BGS
   2020's published total uncertainty is sigma_BGS = 0.18e-3, far larger
   than the 6.7e-5 central gap. PDG experimental uncertainty alone is
   1.1e-5, about 16% of the central gap. Neither is currently propagated
   into the M_KK bound.

## Decision

**Option (a): Accept BGS 2020 as the central value, and add a 1-sigma
sensitivity band derived from the combined SM-theory and experimental
uncertainty.** The paper's headline epsilon_K-driven M_KK lower bound
becomes a **band quote**, not a single number.

### Why (a) and not (b) or (c)

- **(b) Switch to UTfit ~ 2.0e-3 as the central.** Rejected. The
  orchestrator's audit philosophy (`docs/paper_execution_decisions.md`,
  "Audit philosophy") is to **prefer narrowing claims over silently
  changing numbers, unless the existing convention is contradicted by
  literature consensus.** BGS 2020 is the most recent fully NNLO
  perturbative-QCD SM prediction with explicit long-distance and
  kappa_epsilon corrections; it is not contradicted by CKMfitter / UTfit,
  it simply sits at the high end of a band of mutually-consistent
  predictions. Replacing it with UTfit would be a silent convention switch
  with no scientific basis other than wanting a looser bound. That is the
  wrong way to be conservative.
- **(c) Band only, no central.** Rejected as poor presentation: a paper
  that wants to be read needs a quotable number, and a band quote with an
  explicit central is the standard way to provide one in
  flavor-constraint literature (e.g., 1303.5877, 1911.06822, FLAG 2024
  itself).
- **(a) Central + band.** This is the standard treatment of an
  asymmetric, theory-dominated systematic. It preserves the most recent
  modern perturbative input as the central value, but it forces the paper
  to honestly disclose that the headline number depends on the SM
  prediction in a non-trivial way.

## Concrete prescription for the methodology note

The methodology note shall report epsilon_K-driven bounds as a band:

```
M_KK^min (50%, g_s* = 3, epsilon_K only) = X^{+a}_{-b} TeV
```

with the band constructed from the following sources, combined in
quadrature unless a more careful profile likelihood is performed later:

| Source                                          | 1-sigma value          |
|---|---|
| SM theory (BGS 2020, sigma_BGS)                 | 0.18e-3 on epsilon_K^SM |
| Experimental (PDG 2024, sigma_exp)              | 0.011e-3 on epsilon_K^exp |
| SM-choice sensitivity (UTfit/CKMfitter span)    | ~0.15e-3 on epsilon_K^SM |

The **central** quoted bound uses |exp - SM(BGS 2020)| = 6.7e-5.

The **upper-band edge** (looser bound; less NP room) uses the **smaller**
of |exp - SM| evaluated at SM = (BGS central - 1 sigma combined) and the
absolute experimental error floor. For BGS central - 1 sigma combined =
2.161e-3 - sqrt(0.18^2 + 0.011^2 + 0.15^2)e-3 = 2.161e-3 - 0.235e-3 =
1.926e-3, giving an alternative budget |2.228 - 1.926|e-3 = 3.0e-4 (~ 4.5x
larger than central; close to the legacy 4.18e-4 number).

The **lower-band edge** (tighter bound; less NP room) uses SM = (BGS
central + 1 sigma combined) = 2.396e-3, giving |2.228 - 2.396|e-3 = 1.7e-4
absolute value (but the sign flips, meaning the central SM is allowed to
exceed the experimental value within uncertainty; the prescription in this
case is to use the experimental error alone, ~1.1e-5, treating the
constraint as effectively saturated by the experiment).

The methodology note shall:

1. State the central NP budget is 6.7e-5 (BGS 2020 minus PDG 2024).
2. State the band budget runs from ~1e-5 (lower edge) to ~3e-4 (upper
   edge).
3. Quote the corresponding M_KK band, using the leading
   M_KK ~ 1 / sqrt(|Delta epsilon_K|) scaling: band half-width on M_KK is
   roughly 0.5 * d log|Delta epsilon_K| / d log(M_KK) ~ 0.5 of the
   fractional budget uncertainty. A budget range from 1e-5 to 3e-4 maps to
   an M_KK range of roughly sqrt(3e-4 / 6.7e-5) ~ 2.1x on the upper edge
   and sqrt(1e-5 / 6.7e-5) ~ 0.39x on the lower edge of the central bound.
4. Add a single-sentence physics-decision note citing this file.

## Implications for scan re-runs

Per `docs/paper_execution_decisions.md`, the re-run policy triggers when an
audit changes a Wilson coefficient by >10% at M_KK = 3 TeV. The
INVALIDATION_GATE was tripped by hole #5 already; the central-budget shift
is 6.24x, vastly above threshold.

However: the scan re-run is **deferred** until hole #6 (Wilson-RG audit)
also signs off. Re-running RUNA twice (once after #5, once after #6) is
wasteful if hole #6 changes Wilson normalizations or thresholds, since
that would require a third re-run. The orchestrator's stopping conditions
allow this: scan re-runs feed into the final manifest, and the final
manifest is unsigned until both audits pass.

When the scan **is** re-run, it must be run at all three budget edges
(central, +1 sigma, -1 sigma) so the M_KK band can be quoted at
fixed percentile thresholds. The driver is `scripts/run_rs_anarchy.py`;
add a CLI flag `--epsilon-k-budget {central,low,high}` or accept a
numeric override.

## Follow-up requirements before paper claim freeze

1. Edit `docs/quark_scan_methodology_note.tex` to quote the M_KK
   epsilon_K-driven bound as a band, citing this decision file.
2. Confirm hole #6 (Wilson-RG audit) before re-running RUNA.
3. Re-run RUNA at three budget edges; record outputs in
   `docs/artifact_manifest.md`.
4. Add a sensitivity-band figure to the figure set if a single-band quote
   is not enough.

## Sign-off references

- `docs/phase_logs/phase2_h5_impl.md` (implementation log, dc9c498 +
  82a96f0 + 695f35e).
- `docs/phase_logs/phase2_h5_review.md` (peer review,
  APPROVE-WITH-CONCERNS).
- `docs/phase_logs/phase2_h5_signoff.md` (this sign-off,
  PASS-WITH-FOLLOWUP).
- `docs/paper_execution_decisions.md` (orchestrator delegation
  authority).
