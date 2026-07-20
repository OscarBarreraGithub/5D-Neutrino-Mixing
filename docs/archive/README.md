# docs/archive: superseded and historical documentation

Nothing in this directory describes the current state of the project.
It is kept because the project deliberately preserves its history, including
what went wrong and why. For current documentation start at
`docs/README.md`.

## Layout

| Path | What it is | Why it is here |
|---|---|---|
| `superseded/quark_scan_constraint_update_2026-06.md` | The "Z->bb is dominant, 25-30 TeV floor" note | Central claim was the B1 bug; retracted with banner |
| `superseded/artifact_manifest.md` | Seal of the May 2026 `quarkscan-paper-rc1` Delta F = 2 artifacts | Pre-audit, LEGACY banner; rc1 floors are not current |
| `superseded/custodial_problem_statement.md` | June 2026 custodial problem statement | Framing partly retracted after the B1 fix |
| `superseded/paper_execution_roadmap.md` | May 2026 paper-branch roadmap (v2) | Superseded by the July 2026 audit roadmap |
| `superseded/paper_execution_decisions.md` | May 2026 orchestrator decision log | Historical process record |
| `superseded/quark_scan_consolidation_report.tex` | Phase-2 consolidation report source | Pre-audit numbers; never rebuilt post-audit |
| `superseded/RS_vs_SM_flavor_reassessment_plan.md` | Plan for the RS-vs-SM flavor note | Plan completed; the note itself is `docs/RS_vs_SM_flavor_note.tex` |
| `phase_logs/` | 48 impl/review/signoff logs from the May 2026 constraint-catalog build | Historical process records of the dual-signoff gate |
| `overleaf/RS_Flavor_Constraints.pdf` | The pre-pivot lepton-sector Overleaf paper (April 2026) | The project pivoted to the quark sector in May 2026 |

## Related historical trees

- `.orchestration/` holds the multi-agent build ledgers and run logs;
  its May 2026 consolidation cluster is under
  `.orchestration/archive/2026-05_consolidation/`.
- `results/figures/quark_pre_audit_constants/` holds figure exports produced
  with pre-audit constants; they are superseded by `results/figures/quark/`.
- The physics review sources formerly at repo-root `review_local/` now live at
  `reports/physics_reviews/` (they are current, not archived; historical
  documents that mention `review_local/` refer to that directory).
