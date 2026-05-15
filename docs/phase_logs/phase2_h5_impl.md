# Phase 2 Hole #5 Implementation Log: hadronic bag-input audit

Date: 2026-05-15
Branch: `audit/bag-inputs`
Base branch: `paper/quark-scan-2026q2`

## Inventory summary

Audited `quarkConstraints/deltaf2.py` scalar hadronic and Delta m/epsilon_K inputs for K, B_d, B_s, and D systems. The full table is in `docs/audits/bag_param_inventory.md`.

Counts over 36 physical input rows:

- Green: 28
- Yellow: 4
- Red: 4

Yellow rows were documented only: `B_1_BD`, `B_4_BD`, `B_1_BS`, and `B_4_D` drift by 5-10% relative to the selected contemporary source.

Red rows were updated: `B_1_K`, `B_4_K`, `B_5_K`, and `EPSILON_K_SM`.

## Source citations summary

Primary sources used:

- FLAG Review 2024, arXiv:2411.04268: kaon, B, and D decay constants; light-quark masses; `B_K`; kaon BSM bag-parameter averages.
- PDG 2024, Review of Particle Physics, Phys. Rev. D 110, 030001: meson masses, experimental mass splittings, quark masses, and measured `epsilon_K`.
- HPQCD 2019, arXiv:1907.01025: B_d/B_s MSbar-NDR bag parameters at `mu=m_b` and B-mixing SM reference values.
- ETM 2015, arXiv:1505.06639: comparable D0 bag parameters in MSbar at 3 GeV.
- Buras, Guadagnoli, Isidori 2010, arXiv:1002.3612: `kappa_epsilon = 0.94` provenance.
- Brod, Gorbahn, Stamou 2020, arXiv:1911.06822: SM `epsilon_K = 2.161e-3` central value.
- Buras, Misiak, Urban 2001, hep-ph/0102316: operator-basis background only; no Wilson/RG convention changes were made in this hole.

## Code changes

Commit `dc9c498` (`physics(deltaf2): update kaon inputs to FLAG 2024 and BGS 2020`):

- Replaced `B_1_K = 0.717` with `0.5503`, the FLAG 2024 `B_K^MSbar(2 GeV)` value.
- Replaced `B_4_K = 0.78` with `0.903`, the FLAG 2024 Nf=2+1 BSM kaon average.
- Replaced `B_5_K = 0.57` with `0.691`, the FLAG 2024 Nf=2+1 BSM kaon average.
- Replaced `EPSILON_K_SM = 1.81e-3` with `2.161e-3`, Brod-Gorbahn-Stamou 2020.
- Added `test_audited_deltaf2_hadronic_constants_match_selected_sources` with 0.1% relative checks.

No Wilson-coefficient running, RG threshold, operator-basis, CFW comparison, or scan-output files were changed.

## Methodology-note update

Commit `82a96f0` (`docs(paper): document hadronic input provenance`):

- Added `docs/audits/bag_param_inventory.md` with the full provenance and drift table.
- Added appendix subsection `Hadronic input provenance` to `docs/quark_scan_methodology_note.tex`.
- Updated the methodology note's epsilon_K budget sentence from `4.18e-4` to `6.7e-5` to match the audited constants.
- Rebuilt `docs/quark_scan_methodology_note.pdf`; `pdfinfo` reports `Pages: 16`.

## Verification

- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q tests/test_quark_deltaf2.py` -> `5 passed in 5.10s`.
- `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` run twice from `docs/` -> exit 0 both times; no unresolved references.
- `rm -f *.aux *.log *.out` run from `docs/` after the PDF rebuild.

## Invalidation gate

INVALIDATION_GATE_TRIPPED: `B_1_K`, `B_4_K`, `B_5_K`, and `EPSILON_K_SM` change epsilon_K amplitudes or ratios by more than 10% at `M_KK = 3 TeV`; affected final-claim scan is RUNA, because RUNA's final p50/p95/figure claims depend on pre-audit `epsilon_k_ratio` values.

Amplitude probe details with `Im C_i = 1 / (3 TeV)^2`:

- `B_1_K`: one-operator `O1_VLL` epsilon_K estimate shifts by 23.25% relative to the pre-audit amplitude.
- `B_4_K`: one-operator `O4_LR` epsilon_K estimate shifts by 15.77% relative to the pre-audit amplitude.
- `B_5_K`: one-operator `O5_LR` epsilon_K estimate shifts by 21.23% relative to the pre-audit amplitude.
- `EPSILON_K_SM`: central-value NP budget changes from `4.18e-4` to `6.70e-5`, increasing `epsilon_k_ratio` by a factor of `6.238805970149273`.

The exact Wilson/RG convention and matrix-element prefactor audit remains hole #6. The duplicate modern backend constants were noted as out of scope and not touched.

Hole #5 ready for peer review.
