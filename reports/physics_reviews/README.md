# Physics reviews (per-observable derivation checks)

This directory (formerly the repo-root `review_local/`) holds the internal,
equation-by-equation physics reviews written while auditing each constraint
family. Each `.tex` file states the observable, the formulas the repo
implements, the literature source and convention for every ingredient, and any
discrepancies found. LaTeX sources are tracked; compiled PDFs are not
(build locally with `pdflatex <file>.tex`).

These are working reviews, not polished deliverables. The polished
collaborator-facing material lives in `reports/collaborator_2026-06/`.

| File | Covers |
|---|---|
| `physics_summary.tex` | Cross-observable summary of the review pass |
| `constraint_formulas.tex` | Master formula sheet for the scan constraints |
| `constraint_formulas_review.md` | Review notes targeting the formula sheet |
| `deltaF2_framework_review.tex` | The general Delta F = 2 matching/RG framework |
| `epsilon_k_review.tex` | epsilon_K (kaon CP violation) |
| `delta_m_d_review.tex` / `delta_m_s_review.tex` | B_d and B_s mass differences |
| `d0_mixing_review.tex` | D0 mixing (Gedalia et al. funnel) |
| `z_to_bb_review.tex` | Z -> bb (including the B1 bug post-mortem context) |
| `stu_review.tex` | Oblique S, T, U |
| `custodial_review.tex` | Custodial (SU(2)_R / P_LR) protection |
| `collider_kk_review.tex` | Direct KK collider bounds |
| `mu_e_gamma_review.tex` | mu -> e gamma (lepton sector) |
| `rs_flavor_anarchy_review.tex` | Anarchic RS flavor baseline (lane A) |
| `open_questions.md` | Physics questions collected during the reviews |

Vintage: most reviews are from 2026-06-21; several were refreshed on
2026-07-15 alongside the independent repository audit
(`docs/REPOSITORY_AUDIT_AND_RESEARCH_ROADMAP_2026-07-15.md`).
