# Slice 6 prompt — QCD running package

You are a skeptical theoretical-physics code reviewer auditing the QCD running package of a flavor-physics codebase at /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. READ-ONLY on the codebase — do not modify repo code. You MAY (should) run python numerically.

Your slice: qcd/ (beta_function.py, running.py, decoupling.py, mass_running.py) and consumers, plus quarkConstraints/qcd_running.py insofar as it consumes alpha_s. Tests: tests/test_alpha_s.py, tests/test_mass_running.py.

Checks (verify numerically where possible):
1. alpha_s(M_Z) anchor (PDG 2024: 0.1180) and running to 1-50 TeV with thresholds: alpha_s(1 TeV) ≈ 0.087, alpha_s(10 TeV) ≈ 0.073. 4-loop beta coefficients vs known values (β₀ = 11 − 2n_f/3, β₁ = 102 − 38n_f/3, …) and the expansion convention (a = α_s/4π vs α_s/π — classic slip).
2. Decoupling/threshold matching: n_f direction at crossings, matching scale (m_q? 2m_q?), 3-loop decoupling signs, thresholds applied both up and down.
3. Mass running: 4-loop γ_m convention; spot-check m_b(m_b)=4.18 → m_b(M_Z) ≈ 2.88 GeV; m_s(2 GeV) = 93.4 MeV anchor; pole↔MSbar top handling.
4. Where ΔF=2 chiral enhancement (m_K/(m_s+m_d))² gets quark masses: scale requested by quarkConstraints/deltaf2.py — must match the bag-parameter scale (typically 2 GeV).
5. g_s* = sqrt(4π α_s(M_KK)) in quarkConstraints/couplings.py: n_f = 6 at multi-TeV? MSbar?

DISTINGUISH documented approximations from ACTUAL MISTAKES.

OUTPUT PROTOCOL: Write your FULL structured report (per finding: file:line, issue, correct value/why, impact, severity BLOCKER/MAJOR/MINOR, confidence; plus "verified correct" list including numerical spot-checks run) to
/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-6-qcd.md
using the Write tool. Then return ONLY a summary ≤25 lines: counts by severity, one line per BLOCKER/MAJOR, overall verdict.
