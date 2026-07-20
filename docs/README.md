# Documentation guide (start here)

This page is the map of all project documentation: what is current, what
each document is for, which paper and convention each piece of physics comes
from, and where the corresponding code lives. If you are a collaborator
seeing the repository for the first time, read this page, then the three
documents in the first table below, in order.

The project: constraints on Randall-Sundrum warped extra dimensions from
flavor and electroweak physics. Originally a lepton-sector scan (hence the
repository name); since May 2026 the active program is quark-sector, built
around audited Delta F = 2 constraints, RS electroweak observables, and a
103-constraint flavor/collider catalog with a production scan harness.

**Precedence rule:** where any two documents disagree, `FLOOR_SUMMARY.md`
and `reports/collaborator_2026-06/CONTENT.md` win, then `STATE_OF_PROJECT.md`.
Anything under `docs/archive/` is historical and never wins.

---

## 1. Current state (read these first)

| Doc | What it tells you |
|---|---|
| [`FLOOR_SUMMARY.md`](FLOOR_SUMMARY.md) | The one-page canonical KK-scale floor table, with every number tagged by lane and convention |
| [`MODEL_CONVENTIONS.md`](MODEL_CONVENTIONS.md) | The three flavor-model lanes (A anarchic, B production, C FPR ideal) and the conventions that must never be blended |
| [`STATE_OF_PROJECT.md`](STATE_OF_PROJECT.md) | Component-by-component status: what is rigorous, what is approximate, what was actually scanned |
| [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) | Open latent issues in the code and their blast radius |
| [`OPEN_QUESTIONS.md`](OPEN_QUESTIONS.md) | Pending physics decisions (including the ones waiting on the PI) and open work items |
| [`REPOSITORY_AUDIT_AND_RESEARCH_ROADMAP_2026-07-15.md`](REPOSITORY_AUDIT_AND_RESEARCH_ROADMAP_2026-07-15.md) | The July 2026 independent audit: findings (P0-1..P0-7), sequenced execution plan (Gates 0-6), and proposed research directions |
| [`YUKAWA_SUBSTRUCTURE.md`](YUKAWA_SUBSTRUCTURE.md) | The Yukawa perturbation/substructure study: question, results, limitations, and the full program |

## 2. What was done, from which paper, and where it lives

Full per-constraint literature provenance (paper, figure, experimental
inputs, how to map our scan data) is in
[`references/REFERENCES.md`](../references/REFERENCES.md); PDFs of every
source paper are in `references/papers/`, BibTeX in `references/refs.bib`.
The table below is the coarse map.

| Physics | Primary literature | Code | Convention / provenance notes | Review doc |
|---|---|---|---|---|
| Warp geometry, f-factors, KK spectra | Randall-Sundrum; Agashe et al. hep-ph/0412089 | `warpConfig/`, `solvers/` | c = M5/k, alpha = abs(c + 1/2); `derivations/conventions.tex` | `derivations/` |
| Quark mass/CKM fit, MFV spurions | Fitzpatrick-Perez-Randall 0710.1869 | `quarkConstraints/` (fit, model, couplings) | `MODEL_CONVENTIONS.md` (C_Q = r C_u + C_d, bulk-mass map) | `reports/physics_reviews/constraint_formulas.tex` |
| KK-gluon Delta F = 2 matching, epsilon_K | Csaki-Falkowski-Weiler 0804.1954; Bauer et al. 0912.1625; Blanke et al. 0809.1073 | `quarkConstraints/` Delta F = 2 core + `qcd/` running | coupling policy ids in `couplings.py`; the pending normalization decision is `OPEN_QUESTIONS.md` item 1 | `reports/physics_reviews/epsilon_k_review.tex`, `deltaF2_framework_review.tex` |
| Lane A anarchic baseline | Bauer-Casagrande-Haisch-Neubert 0912.1625 (S1-S4) | `scripts/anarchic_bauer_s1.py`, `scripts/run_rs_anarchy.py` | lane A only; never quote as "our model" | `reports/physics_reviews/rs_flavor_anarchy_review.tex` |
| Lane C exact FPR / paper mode | FPR 0710.1869 | `quarkConstraints/paper_0710_1869/` (quarantined, frozen contracts) | `quarkConstraints/PAPER_0710_1869.md`; Q1-only vs Wilson surface split | audit module report 04 |
| B_d, B_s mixing | Blanke et al. 0809.1073 (Table 3, amplitude conventions) | `flavor_catalog_constraints/primary/beauty/` (B001-B004) | SM M12 phase from (V_tq* V_tb)^2, fixed 2026-07-20 | `reports/physics_reviews/delta_m_d_review.tex`, `delta_m_s_review.tex` |
| D0 mixing | Gedalia-Grossman-Nir-Perez 0906.1879 | charm constraints (C001/C002) | | `reports/physics_reviews/d0_mixing_review.tex` |
| Z -> bb | Casagrande et al. 0807.4937 (Fig. 8) | Zbb adapter (T010/T011) | post-B1-fix floor ~5 TeV; the old 25-30 TeV claim was a bug, see archive | `reports/physics_reviews/z_to_bb_review.tex` |
| Oblique S, T, U | Casagrande et al. 0807.4937 (Fig. 4); Archer 1201.1561 | `quarkConstraints/oblique_stu.py` (EW001) | lane-independent existence floor, code-verified 15.96 TeV | `reports/physics_reviews/stu_review.tex` |
| Custodial protection | Agashe et al. custodial SU(2)_R x P_LR literature | `custodial_rs_su2r` model | `CUSTODIAL_PROVENANCE.md` (equation-level ledger), `CUSTODIAL_LITERATURE_CATALOG.md`, `CUSTODIAL_EXPLAINER.tex`, `CUSTODIAL_PAPER_MAP.tex` | `reports/physics_reviews/custodial_review.tex` |
| Direct KK collider bounds | CMS 2603.23454; ATLAS 2512.17856, 2102.13405; Lillie-Randall-Wang hep-ph/0701166 | collider constraints (CR001-CR014) | mass-edge recasts, proxy-flagged | `reports/physics_reviews/collider_kk_review.tex` |
| mu -> e gamma (lepton) | Agashe-Blechman-Petriello hep-ph/0606021; Perez-Randall 0805.4652 | `flavorConstraints/` (L001) | scan default is the MEG II 2025 limit 1.5e-13 | `reports/physics_reviews/mu_e_gamma_review.tex` |
| Lepton masses, seesaw, PMNS | Perez-Randall 0805.4652 | `neutrinos/`, `yukawa/`, `diagonalization/`, `scanParams/` | follow-up scope, pre-pivot machinery | root `README.md` package list |
| 103-constraint catalog | per-constraint anchors cite their sources | `flavor_catalog_constraints/` (95 primary + 8 secondary), website in `flavor_catalog/website/` | proxy status per constraint in `.orchestration/NEEDS_HUMAN_PHYSICS.md` | `flavor_catalog/CATALOG_METHODOLOGY.tex` |
| QCD running | standard 4-loop + thresholds; audited July 2026 | `qcd/` | m_t(m_t) = 162.5 GeV MS-bar (audit M-20) | audit module report 03 |
| RS-vs-SM flavor comparison | | `RS_vs_SM_flavor_note.tex` (+ PDF) | verdict and framing in the note itself | |

## 3. Reports and reviews

- [`../reports/collaborator_2026-06/`](../reports/collaborator_2026-06/):
  the polished June 2026 collaborator report (`report.pdf`, `sendout.pdf`,
  `notes.pdf`, custodial 2x2 comparison); its physics source of truth is
  `CONTENT.md`. A self-contained cut-down variant is under `local/`.
- [`../reports/physics_reviews/`](../reports/physics_reviews/): internal
  per-observable derivation reviews (formerly repo-root `review_local/`);
  see its README for the file map.
- [`audits/`](audits/): audit ledgers. The July 2026 full-repo audit
  deliverables are under `audits/full_repo_audit_2026-07/`
  (compendium, fix ledger, minor-findings ledger); older May-vintage
  ledgers (cfw_*, epsilon_k_sm_decision, ...) are historical context.

## 4. Paper sources in this repo

- `quark_scan_methodology_note.tex` (+ PDF): legacy Delta F = 2-only
  methodology note, pre-audit, carries a SUPERSEDED banner. Kept as a paper
  skeleton; its floors are not current.
- `quark_scan_assumptions_compact.tex` (+ PDF): compact assumptions
  walkthrough.
- `CUSTODIAL_EXPLAINER.tex`, `CUSTODIAL_PAPER_MAP.tex`,
  `RS_vs_SM_flavor_note.tex` (+ PDFs): June 2026 custodial and RS-vs-SM
  deliverables.
- `../flavor_catalog/catalog_master.tex`, `CATALOG_METHODOLOGY.tex`: the
  catalog document set.
- `../derivations/`: step-by-step LaTeX derivations (conventions, Bessel/KK
  modes, RS EW gauge KK coupling).

## 5. History and archive

- [`archive/`](archive/): superseded docs (each with a banner), the May 2026
  phase logs, and the pre-pivot lepton-sector paper. See `archive/README.md`
  for what is there and why.
- `../.orchestration/`: ledgers and run trails of the multi-agent build
  (constraint rebuild, audits, reviews). Live ledgers at top level;
  May 2026 consolidation cluster under `.orchestration/archive/`.
- `pdg_quark_target_fix_explainer.md`: May 2026 physics explainer of the
  PDG quark-target gate fix (pedagogical, kept).
- `kk_collider_problem_statement.md`: problem statement for making collider
  bounds rigorous (open work, see `OPEN_QUESTIONS.md` item 5/G3).

## 6. Results and artifacts

- `../results/figures/quark/`: current post-audit figure exports.
- `../results/figures/quark_pre_audit_constants/`: SUPERSEDED (see its README).
- `../results/paper_0710_1869/`: frozen Lane-C verifier artifacts.
- `../artifacts/`: LEGACY rc1 manifest seal (pre-audit) and collaborator
  point-set CSVs.
- `../scan_outputs/`, `../data/`, `../logs/`: local, gitignored heavy trees;
  provenance JSONs accompany the data files that matter.
