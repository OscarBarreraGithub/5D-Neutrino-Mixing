# Slice 3 prompt — RS-EW / Zbb / custodial / oblique

You are a skeptical theoretical-physics code reviewer auditing a Randall-Sundrum (warped extra dimension) quark-flavor codebase at /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. READ-ONLY on the codebase — do not modify any repo code. You MAY run python to check numbers.

Your slice: RS electroweak sector — Z→bb̄, Z FCNCs, custodial protection, oblique parameters. Files to audit in depth:
- quarkConstraints/rs_ew_couplings.py (the big one: z_delta_g_L/R matrices, minimal RS Zbb shift, Casagrande-style m_b²/(2Λ_IR²) fermion-KK admixture, custodial P_LR protection branch around lines 1532-1583)
- quarkConstraints/oblique_stu.py (S, T predictions)
- flavor_catalog_constraints/primary/top_higgs_ew/T010.py (Z→bb̄ minimal), T011.py, T014.py (Z FCNC widths), EW001.py (S,T fit)
- tests/test_rs_ew_phase6a_zbb_fermion_mixing.py, tests/test_rs_ew_custodial_pr1.py, tests/constraints/.../test_T010/T011/T014/EW001
Context docs: docs/STATE_OF_PROJECT.md, .orchestration/PHASE2_PROGRAM_LEDGER.md, .orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md, docs/quark_scan_methodology_note.tex.

Physics to verify:
1. Z-coupling shift from gauge-KK mixing: standard RS result (Casagrande-Goertz-Haisch-Neubert-Pfoh 0807.4937; Agashe et al.) is δg ∝ g_SM · (m_Z²/M_KK²) · L · [f²(c) overlap] with L = log(1/ε) ≈ 37 for IR-localized fermions, plus a non-enhanced universal piece. Check: powers of m_Z and the scale in the denominator — Λ_IR vs physical M_KK (a mix-up = factor 2.45² ≈ 6 in the constraint, directly moving the headline 25-30 TeV minimal floor); presence of the L enhancement; sign/direction of δg_Lb vs the known RS prediction for R_b; EW quantum numbers (b_L: T³ = −½, Q = −⅓; g_L^SM,b = (−½ + s_w²/3)·g/c_w).
2. Fermion-KK-mixing piece m_b²/(2Λ_IR²): formula source and whether bottom-only application is consistent.
3. Custodial P_LR (Agashe-Contino-Da Rold-Pomarol hep-ph/0605341): protection zeroes the gauge-KK δg_Lb at tree level for T³_L = T³_R embeddings, NOT δg_Rb. Implementation "zeroes protected down-left diagonal entries" — all down-type or only b? Is 'pr1_minimal_offdiag' (keep minimal off-diagonals for T014) conservative and apples-to-apples for the minimal-vs-custodial comparison?
4. T010/T011 gating: observable compared (R_b/A_b/δg room), anchor values, one- vs two-sided.
5. EW001/oblique: minimal-RS T ∝ (v²/M_KK²)·L and S without custodial protection — formulas, Peskin-Takeuchi normalization, custodial branch behavior, S-T fit values/correlation. Docs call this a proxy; check order-of-magnitude and scaling anyway.
6. T014 Z FCNC width: Γ(Z→qq̄′) factors (N_c, chirality sum, phase space) and the experimental room.

The headline claim under test: minimal floor 25-30 TeV driven by T010/Zbb; custodial strict floor 2-3 TeV. Scrutinize anything that could shift these.

DISTINGUISH documented approximations/proxies from ACTUAL MISTAKES.

OUTPUT PROTOCOL: Write your FULL structured report (per finding: file:line, code expression, correct result + why, numerical impact on floors, severity BLOCKER/MAJOR/MINOR, confidence; plus "verified correct" list) to
/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-3-rsew.md
using the Write tool. Then return ONLY a summary ≤30 lines: finding count by severity, one line per BLOCKER/MAJOR, overall verdict.
