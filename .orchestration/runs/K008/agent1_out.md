Plan executed: extend K006 dilepton machinery, add K008 adapter/constraint, add focused tests.

Physics: ISU/Buras-Mescia-Smith CPV decomposition with direct CP `A_eff = Im(lambda_t Y_t + Y_NP)/(Y_t 1e-4)`. Sources: `K008.yaml`, https://arxiv.org/abs/hep-ph/0404127, https://arxiv.org/abs/hep-ex/0309072.

Reused: `quarkConstraints/rare_kaon_dilepton.py` Y-function RS proxy and existing adapter. Built: append-only `K_L -> pi0 e+e-` result path.

Budget: HARD `2.8e-10`, from `K008.yaml` PDG/KTeV 90% CL limit.

SM validation: direct CP `4.8363e-12` within YAML `3e-12..6e-12`; constructive total `3.8006e-11` within `8e-12..4.5e-11`.

NEEDS-HUMAN-PHYSICS: ChPT interference sign/CPC total treatment and full RS EW/Z/Z'/electron matching remain flagged; HARD verdict uses direct-CP short-distance term.

Files changed: `quarkConstraints/rare_kaon_dilepton.py`, `flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton.py`, `flavor_catalog_constraints/primary/kaon/K008.py`, `tests/constraints/primary/kaon/test_K008.py`.

Tests: `test_K008.py` 10 passed; `python -m pytest tests/constraints/ -q` 285 passed.