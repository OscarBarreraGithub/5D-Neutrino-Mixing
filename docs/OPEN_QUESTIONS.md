# Open questions and pending physics decisions

**Status:** live document, created 2026-07-20. This is the single place where
unresolved physics decisions are collected. Items marked **DECISION (PI)**
block downstream work and need input from Lisa; items marked OPEN are known
gaps that agents or students can work on once prioritized.

---

## 1. DECISION (PI): KK-gluon coupling normalization (audit P0-2 / M-13)

**This is the main blocker for any new production scan and for quoting a
defensible Lane-B floor.**

The production scan evaluates every Delta F = 2 constraint with the legacy
perturbative KK-gluon coupling,

- `coupling_policy_id = perturbative_4d_legacy`: g_s* = g_s(M_KK), about 1.0.

The physically motivated RS value is volume-enhanced:

- `coupling_policy_id = rs_volume_sqrt2L_physical`: g_s* = g_s(M_KK) x sqrt(2L),
  about 8.5 g_s, with L = log(k/Lambda_IR) about 37.

Both policies are implemented in `quarkConstraints/couplings.py`; the scan
currently hardcodes the legacy one. Since the NP Wilson coefficients scale as
g_s*^2 / M_KK^2, adopting the physical coupling moves the epsilon_K wall by
roughly sqrt(2L), i.e. the verified ~6.3-7 TeV Lane-B wall would become
roughly ~59 TeV by a fixed-axis rescale (not an exact floor; the exact number
requires a rerun). Details and provenance: the M-13 warning block in
`docs/FLOOR_SUMMARY.md` and audit item P0-2 in
`docs/REPOSITORY_AUDIT_AND_RESEARCH_ROADMAP_2026-07-15.md`.

**Question for Lisa, concretely:**

> For the paper, do we (a) keep the explicitly labeled legacy perturbative
> coupling g_s* ~ g_s as a conservative, convention-tagged choice, (b) adopt
> the volume-enhanced physical RS coupling g_s* = g_s sqrt(2L) and rerun the
> production scan, or (c) match the coupling convention of a specific
> reference paper (Csaki-Falkowski-Weiler 0804.1954, Bauer et al. 0912.1625,
> or FPR 0710.1869) so our floors are directly comparable to its figures?
> The choice rescales the epsilon_K floor by a factor of about 8.5, so it is
> a physics decision, not a code decision.

Recommendation prepared by the audit: pick the convention first (Gate 1),
then rerun (Gate 5). No scan compute should be spent before this is fixed.

## 2. DECISION (PI): default down-sector alignment for Lane C (residual C-2)

The exact FPR (0710.1869) implementation needs a default choice for the
RH-down alignment model (how V5KM acts, and what carries the residual Q4/Q5
LR content). As of the 2026-07-20 fix the default kaon observable fails
closed on nonzero Q4/Q5, so this choice is now explicit rather than silent.
Related: Bauer's S2 (RH-down U(3) degeneracy) is a more targeted epsilon_K
cure than the LH V5KM rotation, because epsilon_K is dominated by the LR
(C4) operator. Which alignment limit is the paper's headline model?

## 3. OPEN: exact FPR Table-I benchmark calibration (Gate 3)

Reproduce arXiv:0710.1869 Table I / Eq. (3) inputs exactly in Lane C and
verify the ~2 TeV statement under the paper's own conventions. Blocked in
part by items 1 and 2.

## 4. OPEN: Lane-A corrected median rerun (P0-3)

The code fixes of July 2026 moved the exploratory Lane-A epsilon_K median
from about 30 TeV to about 5.7 TeV under the stated legacy convention
(commit 2c9989b). Headline documents have been reconciled to say the
corrected pre-rerun value; a full corrected Lane-A production rerun (with
the coupling decision of item 1 applied) is still needed before quoting any
new median as a result.

## 5. OPEN: RS-EW coupling machinery for the catalog (G1)

98 of 103 catalog constraints carry a documented NP-side proxy because
`ParameterPoint` does not expose the full RS EW sector (KK W/Z profiles,
Z-fermion coupling shifts, lepton profiles). Building the shared RS-EW
machinery would upgrade most of the beauty/kaon/charm/top rare-decay and
Z-pole families at once. See `.orchestration/NEEDS_HUMAN_PHYSICS.md`
(categories G1, G2, G3) for the full structured list.

## 6. OPEN: physics questions from the review pass

From `reports/physics_reviews/open_questions.md`:

- How should the lepton-sector c parameters be chosen? Can it be done the
  same way as the quark sector (fit to masses plus mixing)?
- The perturbativity gate bounds |Y| < 4 from above. Is there a principled
  lower bound on how small 5D Yukawa entries may be?

## 7. OPEN: benchmark_quark_mfv gates fail honestly (CI stays red at that step)

With lint green (2026-07-20), CI now actually executes the benchmarks, and
`scripts/benchmark_quark_mfv.py` exits 1: its frozen MFV benchmark point
fails the `paper h_RS < 0.30` gate and the Delta F = 2 gate (K-system ratio
about 4.96, dominant C4_LR). This is NOT a July regression: the identical
failure reproduces at the pre-audit-fix commit (verified 2026-07-20 in a
clean worktree at `a393c8b^`). The benchmark's expectations are stale
relative to the audited epsilon_K budget and matching. Decision needed:
recalibrate the benchmark point/gates against a stated convention, or retire
this gate in favor of the Gate-4 validation matrix. Deliberately NOT patched
by loosening tolerances.

## 8. OPEN: Yukawa substructure program (the profiled study)

The completed noise/gradient study is a proof of concept, not the fit
boundary. The full program (36-dimensional perturbation space, profiled
refits, whitened Jacobians, tangent cones) is specified in
`docs/YUKAWA_SUBSTRUCTURE.md` and section 9 of the July audit. It is the
core of the proposed paper direction and should be scoped with Lisa.
