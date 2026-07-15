# Note to collaborators: full-repository audit and fixes (July 2026)

## Short version

We commissioned a full, independent physics-and-code audit of the whole repository (every module,
script, notebook, derivation, and the experimental anchors). It found a set of real errors, most of
which we had not caught because they were hiding in convention mismatches, stale certification docs,
and self-referential validation that could not see its own bugs. We have now fixed or explicitly
dispositioned every finding (7 criticals, 36 majors, and roughly 50 minor notes), each fix produced by
one agent and independently re-derived and re-checked by a second before it was accepted. The test
suite is green. A per-finding ledger with the disposition and rationale for each item lives in
`docs/audits/full_repo_audit_2026-07/`.

## The headline number: what actually changed for the epsilon_K floor

The single most important correction is to how we quote the flavor (epsilon_K) production floor. The old
core code gated new physics on the bare central gap between the measured and Standard-Model epsilon_K
(6.7e-5), with no uncertainty propagated. That is a 0.44-sigma cut: absurdly tight, and it was also
inconsistent with the catalog path, which used a different (4.5x looser) number for the same observable.
So the same physics was being scored two different ways depending on which entry point you used, and
neither matched the one-sigma band our own decision note (`docs/audits/epsilon_k_sm_decision.md`)
mandated.

We unified both paths onto a single, sign-aware, one-sigma band policy. The consequence is that the
epsilon_K-driven floor is now honestly a band rather than a single number: roughly 3.3 to 3.6 TeV on the
loose (NP-raises-epsilon_K) edge, about 6.3 to 7 TeV central, and about 4.8 to 5.4 TeV on the other edge.
The old "7 TeV" was the central value of a mislabeled, over-tight gate, not a hard wall.

Separately, we found and hardened a coupling-convention issue (M-13): the production scan uses the
correct physical KK mass on the mass axis, but it uses the perturbative 4D gluon coupling rather than the
volume-enhanced RS coupling (a factor of about sqrt(2 pi k r_c) ~ 8.5). If we adopt the physical RS
coupling, the epsilon_K floor rises by roughly that factor, to order 59 TeV. We did NOT silently change
this number: production still runs the legacy perturbative convention, the convention is now typed and
labeled so the two cannot be confused, and the choice of whether to re-quote at the physical coupling
(and re-run the scan) is flagged as an open decision for the group. Please look at this one; it is a real
publication-level choice.

## Why these were subtle (and where we simply blundered)

Most of the serious errors share a few root causes, and they are worth internalizing so we do not repeat
them:

1. Convention splits. The same quantity was computed in two code paths that were never cross-checked
   against each other: epsilon_K (core vs catalog, 4.5x apart), the mu-to-e-gamma dipole (the 4e-8
   Perez-Randall prefactor was calibrated at Lambda_IR but the catalog path fed it the physical KK mass, a
   factor 2.45^4 ~ 36 in the branching ratio), and three different values of the gluon coupling g_s* living
   under one shared label. Once a quantity has two homes, they drift. We have collapsed these onto shared,
   asserted, and typed policies so a divergence now fails loudly.

2. Factor-of-two normalizations. A cluster of genuine bugs were plain factors of two or sign errors buried
   in normalization conventions: the top FCNC dipole widths (x2 too large), the Higgs coupling profile (an
   unfixed twin of a Zbb bug we had already fixed elsewhere, wrong sign for UV-localized states), the
   mu-to-3e interference coefficient (2 sqrt 2 e instead of the Kuno-Okada 8e, and the wrong sign), the
   tau three-body and radiative rates (missing the tau leptonic branching fraction, a factor 5.65), the
   Z-prime C9/C10 matching (x2 and a sign), and the B-to-K* kernel. Each is small in isolation but each
   moved a real constraint.

3. Stale certification and circular validation. Our paper-facing review documents still certified the
   pre-fix (swapped, doubled) matrix elements as "confirmed," even though the code had already been fixed;
   and the "modern" verifiers compared fields written by the same builder, so they were structurally
   incapable of catching any of the above. We corrected the review docs to match the code and added honest
   banners stating exactly what the self-consistency checkers do and do not verify.

4. The FPR "Lane C" reproduction was broken end to end. This is the closest thing to a straight blunder:
   the Wilson-coefficient renormalization-group step ran in the inverted direction with an untransposed
   anomalous-dimension matrix, a Fierz sign was frozen wrong, the matrix elements were a factor of 4 off,
   the volume factor was missing, the default path structurally zeroed the dominant operator, and the
   seed-to-profile map inverted the localization hierarchy. Its own validation was circular, so nothing
   flagged it. We fixed every convention error we could pin without the paper, de-circularized the tests
   against an independent closed-form oracle, and quarantined the lane (it feeds no headline number) with
   an explicit note that the exact Table-I calibration is still pending the arXiv source, which is not in
   the repo.

## What else was corrected

The existence (oblique S,T) floor was documented as 18 to 20 TeV but the shipped anchors actually solve to
15.96 TeV; we reconciled the docs to the code and flagged the anchor citation as needing to be pinned
before publication. The QCD running had a wrong 3-loop decoupling constant and a spurious top-threshold
matching (m_t reconciled to 162.5 GeV). The KK Bessel solver silently returned wrong tower masses beyond
the fifth root. Collider SSM benchmarks were being applied as hard vetoes to states whose couplings are
suppressed, so we downgraded them to advisory. The anarchic "Lane A" reproduction walls were about a
factor of two too high from an uncanceled sqrt(2) in the Yukawa bridge plus a missing Wolfenstein A. And
several reproducibility problems (a cited results file that is not in the tree, un-bannered stale figures
and notebooks asserting the retracted 25 to 30 TeV Zbb floor) are now labeled honestly.

## How it was fixed, and how much to trust it

Every code finding went through the same loop: one agent researched and implemented the fix, a second,
independent agent adversarially re-derived and re-executed the physics (not just re-ran the tests) before
the fix was accepted, and only then was it committed. Findings that we deliberately did NOT change are
recorded with the reason (either they need a dedicated physics pass, or they are convention choices, or
the audit itself vindicated the original code on a second look). The two items that most deserve your eyes
are the epsilon_K coupling-convention re-quote (7 vs ~59 TeV) and the EW anchor citation, both flagged as
open decisions rather than quietly resolved.

Everything is on one branch for review before it touches main.
