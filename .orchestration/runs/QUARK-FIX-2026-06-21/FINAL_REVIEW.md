# QUARK-FIX — orchestrator final review (2026-06-22)

End-to-end verification by the orchestrator after dual-signoff (Opus-level; Codex
backup still owed before push).

## Change set (exactly the 6 approved items, nothing extraneous)
- B2 epsilon_K SVD->PDG rephasing: quarkConstraints/fit.py (+124)
- B3 dF=2 O4/O5 unswap + 1/(2m_M), all 6 sites: deltaf2.py, modern/phenomenology.py
- B1 CGHNP Zbb retranslation: rs_ew_couplings.py
- M2 EW001 dT physical-M_KK convention: oblique_stu.py
- M1 T010 R_b -> T011-style NP-shift budget: T010.py
- M5 tag substring fix: T001.py, T002.py, run_full_catalog_scan.py, build_scan_explorer.py
- Tests: 13 test files updated with literature-anchored ABSOLUTE pins (fail on old code).

## Verification
- Focused suite on fixed modules: 95 passed (epsilon_k, deltaf2, zbb, T010, EW001, K001, B003, quark_fit).
- Full suite at impl-review: 1767 passed / 1 skip (independent Opus reviewer, re-derived numbers).
- Scan provenance: both scans logged "git HEAD f94a517 + uncommitted QUARK-FIX working tree" => fixes ARE live in the runs.

## Physics result (minimal model, 100k quark-only, 99,979 evaluated)
Per-constraint own M_KK floor (sharp 100%->0% walls):
  S,T,U (EW001, PROXY)     20 TeV  <- leading (M2 made dT ~6x stronger)
  epsilon_K (K001, rigor)   7 TeV  <- leading RIGOROUS / flavor constraint
  Z->bb (T010, rigor)       5 TeV  <- was the ~25-30 TeV (108 TeV @1sigma) ARTIFACT
  Delta m_s (B003)          2 TeV
  D0, Delta m_d, Z->bb Ab  <=1 TeV (never bind)
Combined inclusive floor: 20 TeV (proxy-driven, S,T,U).
Strict (rigorous-only) floor: 7 TeV (epsilon_K).
Dropping Z->bb leaves the combined floor at 20 TeV => Z->bb is NOT the driver. CONFIRMED.
Collider cut m1>=5.5 TeV never binds (far below 20 TeV). 
mu->e gamma (L001): lepton-sector, evaluated only on the ~5% lepton-viable (perturbative
  seesaw) draws; veto 45%->16% over 1->20 TeV; reported as a separate panel, NOT ANDed
  into the quark/EW floor.

## Consistency with the audit's predicted post-fix picture: MATCHES.
- B1+M1: minimal Zbb artifact collapses (108 TeV -> 5 TeV). [audit: "minimal floor is an artifact"]
- B2+B3: epsilon_K becomes a real ~7 TeV floor. [audit: corrected epsilon_K strengthens, becomes dominant rigorous]
- M2: EW001 ~6x stronger -> 20 TeV inclusive (proxy). [audit M2: EW001 dT x6 underestimated]

## Deliverables
- scripts/build_constraint_matrix.py (JSONL -> per-draw pass/fail parquet, evaluated flag)
- notebooks/constraint_explorer.ipynb (interactive toggle, precompute-then-AND, no recompute)
- Robust matrix: scan_outputs/wq_quarkonly_20260622T090807/constraint_matrix.parquet
- Lepton matrix: scan_outputs/fix100k_minimal_20260622T080053/constraint_matrix.parquet

## NOT done (deferred per user "100k for now" scope)
- Custodial re-run (minimal only requested).
- STATE_OF_PROJECT.md / methodology-note floor restatement (still cite old artifact numbers).
- Per-item commits + Codex pre-push backup review. NOTHING pushed.
