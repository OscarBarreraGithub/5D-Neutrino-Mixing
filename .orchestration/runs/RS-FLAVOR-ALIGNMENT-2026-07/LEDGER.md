# RS-FLAVOR-ALIGNMENT-2026-07 — research ledger

**Goal.** Anarchic RS forces M_KK high because epsilon_K's KK-gluon LR operator Q_4
has an O(1) CP phase. Low-M_KK survival requires a phase/off-diagonal cancellation.
We want to (1) characterize that cancellation and turn it into a *prediction*, and
(2) find a better-motivated-than-anarchy Yukawa structure ("alignment") that
suppresses epsilon_K and predicts something distinctive — ideally a novel step
beyond the known cures (RH-down U(3) "S2", LH V_5KM "FPR").

## The crux (physics)
C_4 ~ (g_s*^2 / M_KK^2) (G_L)_12 (G_R)_12,
  (G_L)_12 ~ f_Q1 f_Q2 e^{i phi_L},  (G_R)_12 ~ f_d1 f_d2 e^{i phi_R}.
epsilon_K ∝ Im(C_4) ∝ sin(phi_L + phi_R).  Delta m_K ∝ Re(C_4) ∝ cos(phi_L+phi_R).
Survival at low M_KK  <=>  phi_L + phi_R ≈ 0 or pi (Im small, |C_4| still large).
Custodial is blind to Q_4 (color-octet, EW-singlet).

## Two questions
- **Q1 (cancellation → prediction).** If Im(C_4^K)≈0 is forced, |C_4^K| stays large,
  so Delta m_K should saturate; and the CP violation rotated out of down-12 must
  reappear (up-sector/D0, EDMs, B phases, K->pi nu nu, eps'/eps). Which is the
  sharpest smoking gun?
- **Q2 (submanifold / better draw).** Codimension & symmetry of the survival set;
  fine-tuning vs natural region; a few-parameter motivated draw that suppresses
  epsilon_K AND predicts something. The "big novel step."

## Loop protocol (dual sign-off gate)
Each hypothesis Hn runs on its own git branch `align/Hn-<slug>`:
1. IDEA — Opus thinker (+ Codex independent) propose; orchestrator (Claude main) picks.
2. IMPLEMENT — Codex worker writes the test/generator on the branch; Claude reviews.
3. TEST — run the forward MC / analysis; store figures+numbers under this run dir.
4. CROSS-REVIEW — Codex and Claude each independently assess (real effect? bug?
   prediction falsifiable?). Both must APPROVE or the finding is PARKED, not merged.
5. VERDICT — CONFIRMED / KILLED / PARKED, recorded below. Merge only CONFIRMED tooling.

**No result is "real" until Codex and Claude have independently signed off.**

## Perspectives in flight (this session)
- Explore agent — code map of the epsilon_K/C_4/phase computation (API to build tests).
- Opus thinker — physics theory + hypothesis menu (Claude perspective).
- Codex worker — independent physics theory + hypothesis menu (cross-review).

## FINDINGS
### F1 — eps_K survival is MAGNITUDE-dominated, not phase-cancellation (PLAUSIBLE, needs Codex sign-off)
Data: existing `scan_outputs/anarchic_reproduction/anarchic_bauer_S1.parquet` (990k draws,
112k mass+CKM-consistent). Among 474 eps_K survivors at M_KK<=3 TeV:
- **57% MAGNITUDE-suppressed** (small |C_4|: ratio_dm_K < 0.1x typical) — kills eps_K AND Delta m_K together.
- **~7% PHASE-aligned** (large |C_4|, sin(Phi)~0: ratio_dm_K >= failer-typical) — the Delta m_K-SATURATING tail; these are the "fine-tuned off-diagonals."
- corr(log ratio_eps_K, log ratio_dm_K) = **+0.82** (magnitude-dominated; phase-cancellation alone would give ANTI-correlation).
- survivor eps_K/dm_K ratio (~tan Phi) is 27x smaller than failers (phase-skew real but subdominant).
**Reframing:** the naive "you tune the phase to cancel eps_K" is the MINORITY path (~7%). Flat
anarchy mostly survives by drawing accidentally small flavor-violating MAGNITUDES (small rotations).
This overturns the working premise and Opus-H1's "all survivors are Delta m_K-loud." Caveat: uses
ratio_dm_K as a |C_4| proxy; the instrumented run (true Phi) will confirm/refute.

### F2 — Nelson-Barr / CP-only lane works numerically (VALIDATED, Codex-B spec + Claude build)
NB draw (Y_d real anarchic, CP via up-sector rank-one spurion; eta_leak = down leakage).
At M_KK=3 TeV, eta_leak=0: |sin Phi_12| median 1.2e-16 (down phase exactly cancels),
|C_4| RETAINED at anarchic scale (log10 -15.9 vs flat -15.5 => Delta m_K large),
**100% eps_K pass** (414/414) vs flat ~1%, CP displaced to up sector (|sin PhiD_12|~0.64),
CKM alive (|J|/J_pdg~1.4). This is the defensible-new corner (D): eps_K-safe WITH anarchic
down magnitude (unlike S2/FPR/Redi-Weiler which suppress magnitude). OPEN decisive test
(running, job bhozq792x): does I_d (EDM CP invariant) stay small for NB, or does up-sector
CP leak back via Y_uY_u^dag? That is the make-or-break for the "EDM-quiet" claim.

### NOVELTY BOUNDARY (Codex-D) — what we can honestly stake
- Do NOT claim broad "CP-only alignment fixes eps_K+EDM" (Redi-Weiler 1106.6357, Minimal CP
  Violation) or "warped Nelson-Barr" (Girmohanta 2203.09002) as new.
- DEFENSIBLE new lane: CP-only kaon alignment KEEPING ANARCHIC DOWN MAGNITUDE (Delta m_K stays
  large) + the (Delta m_K, nEDM) discriminator + the F1 survival taxonomy. Differentiator MUST
  be anarchic-down-magnitude + kaon/EDM phenomenology. Radiative stability needs a spurion proof.
- Rank-one/U(2) lane (Codex-E2): structural magnitude suppression, floor ~2-3 TeV; WARNING the
  c-fit must enforce c_Q1=c_Q2, c_d1=c_d2 or U(2) silently breaks. Decisive metric: M50_all
  (min M_KK with >=50% pass-all) across {flat | RH-align | rank-one+U(2) | +CP-sequester}.

## Hypotheses (reframed by F1: magnitude-alignment [known: S2/FPR] vs phase-alignment [novel])
| ID | Claim (one line) | Branch | Status |
|----|------------------|--------|--------|
| H1 | The ~7% phase-aligned tail is real with the TRUE Phi and is Delta m_K-saturating (+ EDM-loud) | align/instrument-phase | RUN IN FLIGHT (500k instrumented draws) |
| H2 | **BIG NOVEL STEP:** IR-brane spontaneous CP / phase-only alignment keeps magnitudes anarchic (masses+CKM still predicted) but controls phases => eps_K-safe + Delta m_K-loud + EDM-QUIET. Discriminates from S2 (mag kill) and flat-tail (EDM-loud) in the (Delta m_K, d_n) plane | tbd | DESIGN |
| H3 | The dominant (magnitude) survival mode IS the S2/FPR alignment locus in disguise (small V_dL/V_dR 12) — quantify the survival submanifold = codim structure | tbd | DESIGN |
| H4 | Benchmark S2 (RH-down U(3)) & FPR (V_5KM) with the instrumented generator: confirm they act on the |C_4| MAGNITUDE axis | tbd | DESIGN |
| H5 | One-parameter Y_d = r Y_u + delta partial alignment ties eps_K suppression to CKM ratios (distinctive vs S2/FPR) | tbd | DESIGN |

**Dependency/GAP:** no neutron/quark EDM estimator in the code yet. H1/H2's sharpest smoking gun
(EDM) needs a chromo-dipole estimator built from the KK structure. TASK before H1/H2 can close.

## THE STEER (correct deep-research report, 2026-07-02) — 3-part structural change
The literature verdict: the fix is NOT prettier alignment of two i.i.d.-anarchic Yukawas.
It is (1) **replace the anarchic ENSEMBLE** with sequentially-generated LOW-RANK spurions
(Y_f ~ y3 n3n3d + eps y2 n2n2d + eps' y1 n1n1d; Greljo-Thomsen) under an accidental
**intermediate** light-family symmetry (U(2)^5; Isidori/Barbieri) — the sweet spot BETWEEN
anarchy and MFV; (2) **factorize/localize** flavor breaking away from the deep-IR
Higgs/compositeness sector (flavor branes, multi-scale; Vecchi); (3) **sequester CP** so
kaon CPV is spurion-suppressed not phase-cancelled (CP-odd bulk scalar; Cheung-Fitzpatrick-
Randall / Nelson-Barr). Each piece exists separately; NO calculable 5D construction combines
all three -> that gap is our target. Keeps the RS virtue (geometric mass hierarchy), fixes
eps_K structurally. Refs: FPR 5D-MFV, Santiago min-flavor-protection, Csaki-Perez-Surujon-
Weiler flavor-shining, A4, Bauer-Malm-Neubert extra-SU(3), Greljo-Thomsen rank-one,
Vecchi flavor-branes, Cheung-Fitzpatrick-Randall CP-sequester, Isidori U(2)^5.

**Convergence:** our F1 + Opus + Codex independently reached the CP-sequester axis
(Nelson-Barr); the report ADDS the ensemble axis (rank-one) + intermediate-symmetry framing.

## CODEX FLEET (5 agents, disjoint, ANALYSIS/SPEC only — running 2026-07-02)
- A (bd5q7zltp): EDM estimator spec (the missing observable / sharpest smoking gun).
- B (btxrpy9j1): Nelson-Barr / single-CP-spurion draw spec (CP-sequester axis).
- C (bm064ggg0): tuning-volume analysis, explains F1 (P_pass(R_K) arcsin law).
- D (btfflykqs): novelty / literature boundary audit (what is genuinely new).
- E2 (bmhkuall2): the rank-one / U(2) "wrong ensemble" axis (report's deepest point).
- (old cx_E on the WRONG eigenquestions report was NEVER launched — user corrected it.)
Synthesis target: build EDM estimator + rank-one & NB draw modes into instrument_epsK_phase.py,
then the decisive floor-vs-mechanism comparison {anarchy | rank-one | NB | S2 | FPR} in
(eps_K floor, Delta m_K, d_n) space, and the great report.

## Log
- 2026-07-02: run opened; 3 independent perspectives launched (map + Opus theory + Codex theory).
- 2026-07-02: **F1** found from existing parquet (magnitude-dominated survival). Built
  `scripts/instrument_epsK_phase.py` (retains complex G_L/G_R, true Phi) on branch
  `align/instrument-phase`; smoke-tested; launched 500k-draw instrumented run (M_KK 2,3 TeV).
  Codex theory job still running.
