## Full-repository physics + code audit: fixes

Addresses every finding in the July 2026 full-repo audit (`docs/audits/full_repo_audit_2026-07/`, imported
from the Fable review): **7 criticals (C-1..C-7), 36 majors (M-1..M-36), plus ~60 minor notes** (Section 4 +
Section 8.4). Each code finding went through the same loop, run serially: **Codex researched and implemented
the fix, then an independent Claude agent adversarially re-derived and re-executed the physics** (not just
re-ran tests) before it was accepted. Nothing was silently dropped: findings not changed are recorded with a
rationale.

Full test suite: **1806 passed, 1 skipped** (the skip is a pre-existing environmental fixture). Test count
rose from 1767 as regression/oracle tests were added throughout.

### Disposition ledgers (read these first)
- `docs/audits/full_repo_audit_2026-07/FIX_LEDGER.md` — every numbered finding, its disposition, the fix
  commit, and the audit verdict.
- `docs/audits/full_repo_audit_2026-07/MINOR_FINDINGS_LEDGER.md` — every Section-4/8.4 minor (8 FIXED, 50 DOC,
  2 WONT-FIX, 1 ALREADY-CORRECT, 2 ESCALATE-then-resolved).
- `docs/audits/full_repo_audit_2026-07/COLLABORATOR_MESSAGE.md` — a plain-language summary for the group.

### Two OPEN DECISIONS flagged for you (not resolved here, on purpose)
1. **epsilon_K physical-coupling re-quote (7 vs ~59 TeV).** Production uses the physical KK mass but the
   *legacy perturbative* gluon coupling (g_s* ~ 1), not the RS volume-enhanced coupling (g_s* ~ 8.5 g_s).
   Adopting the physical coupling raises the epsilon_K floor by ~sqrt(2 pi k r_c) ~ 8.5, to roughly 59 TeV
   (needs a scan rerun for the exact number). M-13 HARDENS the convention (typed, asserted, three distinct
   coupling-policy ids) and keeps production numbers bit-identical; whether to re-quote is your call.
2. **EW001 oblique anchor citation.** The anchors are a plausible recent U=0 fit but the exact "PDG 2025
   Table 10.8" attribution is unpinned; flagged "pin before publication."

### Notable physics changes to headline numbers
- **epsilon_K production floor is now a BAND** (M-1/M-2/M-4): the old core gate was the bare 0.44-sigma
  central gap and was 4.5x tighter than the catalog path. Unified onto one sign-aware 1-sigma band policy
  (per `epsilon_k_sm_decision.md`): raise-edge ~3.3-3.6 TeV, central ~6.3-7 TeV, lower-edge ~4.8-5.4 TeV.
- **EW existence floor reconciled to 15.96 TeV** (C-3): docs said 18-20 TeV; the shipped anchors solve to
  15.96 TeV. Docs corrected to the code.
- **Lane C (FPR) convention-corrected but QUARANTINED** (C-1,C-2,C-5,C-6,M-8,M-13,M-27,M-28): RG direction +
  transpose, Fierz sign, matrix elements /4, volume factor, xi_g, seed->profile orientation, LR/C4 capability;
  de-circularized against an independent BBL oracle. Lane C feeds no headline number; exact FPR Table-I
  calibration remains pending arXiv:0710.1869 (not in-repo).
- Mechanical factor-2/sign fixes: top FCNC dipole x2 (M-9/M-31), mu->3e interference (M-29), tau leptonic-BR
  (M-30), B->K* kernel (M-32), R(D) proxy (M-33), Higgs B-profile B1-twin sign (M-10), Z' C9/C10 (RD-01).
- QCD (M-19/M-20/M-21): CKS 3-loop d3, equal-nf matching no-op, m_t(m_t)=162.5 GeV, sigma_t. KK Bessel deep
  tower roots (C-4). Lane A anarchic walls ~2x (M-11/M-12). mu->e gamma convention unified (M-14).
- Budget/veto logic (M-5/M-15/M-16/M-17/M-18): B022 SM point passes, gauge-equivalent reported seed,
  single overall_scale, hard-partial vetoes, fail-closed on malformed input.
- Docs/reproducibility (M-22/M-23/M-24/M-25/M-34/M-35/M-36): pre-B3 certification corrected, stale
  parquet/Zbb-floor references and notebooks bannered, honesty banners on self-consistency checkers,
  website joint-veto floor.

### Scope note
This branch is based on `align/instrument-phase`, so it also brings that session's research work to main,
**including the Yukawa perturbativity work** (`|Y*|<3` enforcement and the F7 Yukawa-perturbation study), as
requested, so this is one consolidated PR. The audit-fix commits are cleanly separated from the research
commits in the history.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
