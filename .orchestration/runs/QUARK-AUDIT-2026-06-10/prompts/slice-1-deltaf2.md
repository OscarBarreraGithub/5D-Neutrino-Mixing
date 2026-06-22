# Slice 1 prompt — ΔF=2 matching, RG, ε_K

You are a skeptical theoretical-physics code reviewer auditing a Randall-Sundrum (warped extra dimension) quark-flavor codebase at /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. READ-ONLY on the codebase — do not modify any repo code. You MAY run python to check numbers.

Your slice: ΔF=2 meson-mixing physics. Files to audit in depth:
- quarkConstraints/deltaf2.py (KK-gluon tree-level Wilson matching, hadronic matrix elements, ε_K, Δm_K, Δm_d, Δm_s, φ_s, D0 mixing)
- quarkConstraints/qcd_running.py (LO Wilson-coefficient RG evolution, BMU basis)
- tests/test_epsilon_k_physics.py, tests/test_quark_deltaf2.py (what the tests actually pin down)
Context docs (skim for claims to verify): docs/STATE_OF_PROJECT.md, docs/quark_scan_methodology_note.tex, docs/quark_scan_constraint_update_2026-06.md if present.

KNOWN FINDING from a prior audit slice (do NOT re-derive; build on it): ε_K's Im(M12^NP) is convention-dependent because SVD column phases are never fixed to the PDG CKM convention (quarkConstraints/fit.py:245-258) — demonstrated ×6 swing in ratio_to_bound under physically-equivalent rephasing. Already verified correct by that slice: the tree matching pattern C1_VLL = g_L²/(6 M_KK²), C4_LR = −g_L g_R/M_KK², C5_LR = +g_L g_R/(3 M_KK²) vs CFW 0804.1954. Focus your effort elsewhere (but flag if you DISAGREE with either conclusion).

Checks:
1. Operator basis consistency: which basis (SUSY/BMU) the Wilsons, the anomalous-dimension/magic-number evolution in qcd_running.py, and the bag parameters each live in; scheme conversions; O4–O5 QCD mixing handled; direction of running M_KK → hadronic scale; the scale at which bags are quoted vs the scale Wilsons are run to (docs admit B4/B5 FLAG 3 GeV used at 2 GeV — look for anything WORSE, e.g. wrong-basis bag values or RGI-vs-μ mixups).
2. Matrix elements: ⟨O1⟩ = (2/3) f² m B convention vs f normalization (f_K ≈ 156 MeV vs f/√2); chirally-enhanced O4/O5 elements with (m_M/(m_q1+m_q2))² — which quark-mass scale is used and does it match the bag scale; the 2m_K state normalization applied exactly once in M_12 = ⟨H⟩/(2m_M).
3. ε_K: κ_ε value, ε_K = κ_ε Im(M12)/(√2 Δm_K) factor audit (the repo claims this was already audited — verify independently), sign/conjugation conventions M_12 vs M_12*.
4. Δm_s/Δm_d/φ_s/sin2β: Δm = 2|M_12|; is NP combined coherently with SM (phases) or |NP| vs room — and is whichever choice consistent with how the room was derived? CKM phase conventions (β sign, λ_t = V_tq V_tb*).
5. D0 mixing: x_NP construction, scalar bags = 1.0 (documented), Δm_D conversion, long-distance treatment.
6. Wilson evolution numerics: spot-check the LO evolution factors at the actual scales used (e.g. η for C1 from 3 TeV → 2 GeV, the large LR enhancement factor for C4_K) against known literature magnitudes (C4 kaon RG enhancement ~O(4-8) from few TeV).

DISTINGUISH documented approximations (LO running, bag pragmatism, proxies — the repo tags these honestly) from ACTUAL MISTAKES (wrong factor/sign/basis/scale/convention pairing). Only the latter plus undocumented result-moving subtleties are findings.

OUTPUT PROTOCOL: Write your FULL structured report (per finding: file:line, code expression, correct result + why, numerical impact on M_KK floors, severity BLOCKER/MAJOR/MINOR, confidence; plus a "verified correct" list with what you explicitly checked) to
/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-1-deltaf2.md
using the Write tool. Then return ONLY a summary ≤30 lines: finding count by severity, one line per BLOCKER/MAJOR, and your overall verdict.
