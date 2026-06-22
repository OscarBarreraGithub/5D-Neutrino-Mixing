# Slice 5 prompt — scan harness & conventions

You are a skeptical code+physics reviewer auditing the production scan harness of a Randall-Sundrum quark-flavor codebase at /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. READ-ONLY on the codebase — do not modify repo code. You MAY run python.

Your slice: the pipeline behind the headline results (minimal RS physical-M_KK floor 25-30 TeV; custodial strict 2-3 TeV, inclusive 7 TeV). Audit in depth:
- scripts/run_full_catalog_scan.py (parameter draws, M_KK/Λ_IR conversion DEFAULT_XI_KK=2.4487, quark-only allowlist, strict vs inclusive classification, tag_result rigorous/proxy/partial/stub, seeding)
- scripts/wq_quarkonly_1m_plan.py (grid r ∈ {0.05,0.1,0.25,0.5,1.0}, M_KK ∈ {1..50} TeV, 20k draws/cell, seeds)
- scripts/build_wq_quarkonly_comparison.py (minimal/custodial row pairing)
- flavor_catalog/website/scripts/build_scan_explorer.py (floor extraction, FLOOR_THRESHOLD=0.5)
- scan_outputs/wq_quarkonly_1M_20128400/analysis/analysis_report.md + custodial twin (sanity-check quoted numbers)
- tests/test_full_catalog_scan_harness.py
Context: docs/STATE_OF_PROJECT.md, .orchestration/PHASE2_PROGRAM_LEDGER.md, docs/quark_scan_constraint_update_2026-06.md if present.

Checks:
1. Scale-convention leaks: scan axis is PHYSICAL M_KK; Λ_IR = M_KK/ξ_KK feeds the model. Trace every scale entry point per constraint family (ΔF=2 1/M², Zbb m_Z²/scale², dipoles, collider proxies, oblique) — each must use its intended convention exactly once. One constraint internally using Λ_IR where it means physical M_KK shifts its floor ×2.45. README warns the LFV default is ξ=1 — confirm no quark constraint inherited that.
2. 878,707 evaluated of 1,000,000 rows: find the drop reason and whether drops bias survival (do failures cluster at low M_KK? which denominator do veto fractions use?).
3. Strict vs inclusive: confirm classification code matches docs; 'partial'/'stub' can never veto; tag_result string-matching brittleness; exceptions → stub → silently not vetoing — are exceptions counted/reported anywhere?
4. Determinism/pairing: do minimal and custodial runs draw IDENTICAL Yukawa matrices per (r, M_KK, seed) — verify the draw path, not just seed labels.
5. Floor extraction: grid {1,2,3,5,7,10,15,20,30,50} TeV with veto-fraction ≤ 0.5 threshold — how were "25-30" and "2-3 TeV" derived; does prose anywhere claim "excluded below X" where the math only supports "≥50% of draws vetoed"? Threshold sensitivity (0.5 vs 0.05).
6. Yukawa draw measure: anarchic distribution (uniform in what? range?), perturbativity cuts consistent across both runs.
7. Spot-check analysis-report numbers against comparison artifacts if cheap.

DISTINGUISH documented conventions from ACTUAL MISTAKES.

OUTPUT PROTOCOL: Write your FULL structured report (per finding: file:line, issue, why it matters, impact on headline floors, severity BLOCKER/MAJOR/MINOR, confidence; plus "verified correct" list) to
/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-5-harness.md
using the Write tool. Then return ONLY a summary ≤30 lines: counts by severity, one line per BLOCKER/MAJOR, overall verdict.
