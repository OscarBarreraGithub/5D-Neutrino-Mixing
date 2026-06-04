1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.
4. Phase formulas are standard: β and β_s definitions match UT conventions; recompute gives 2β=45.09992°, sin2β=0.7083388935, φ_s=-0.0379451028 rad, within PDG-scale targets.
5. Wiring/SM limit: B002=0.7083388934693238 and B004=-0.037945102831894784 come from `repo_default_ckm_phases()`, not YAML central literals; zero NP reproduces SM.
6. NP/M12: phase shift remains `arg(1 + M12_NP/M12_SM)` using complex QCD-running DeltaF=2 M12; no DeltaF2/M12 files touched.
7. Flags/isolation/numerics: CKM-phase NEEDS-HUMAN flag removed from B002/B004 runtime diagnostics/notes; other historical PARTIAL sidecar metadata kept; scalar result fields are real finite floats.
8. Tests: focused `22 passed`; full `pytest tests/ -q` => `1696 passed, 1 skipped in 786.95s`.
9. W4-OK