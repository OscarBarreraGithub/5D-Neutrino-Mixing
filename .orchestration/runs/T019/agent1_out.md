Plan: studied examples/scaffold/T019 YAML, added T019, tested targeted and full suite.

Physics: pure-NP `BR(h -> e tau) = m_h(|Y_eτ|² + |Y_τe|²)/(8π Γ_h)`, SM set to `0.0`.

Reused: existing `physics_adapters.higgs_lfv` and `quarkConstraints.higgs_lfv`; no new adapter or T018 changes.

Budget: `2.0e-3` from `ATLAS2023:T019:etau_limit` in `T019.yaml` (`95% CL`, arXiv:2302.05225); PDG ATLAS entry matches.

Validation: for `Y_eτ=2e-4+1e-4i`, `Y_τe=0.5e-4i`, independent core/manual BR is `6.428398173640323e-05`, ratio `0.03214199086820161`.

Gap: flagged `NEEDS-HUMAN-PHYSICS` because full RS off-diagonal charged-lepton Higgs-Yukawa matching is not on `ParameterPoint`.

Files changed: `flavor_catalog_constraints/primary/top_higgs_ew/T019.py`, `tests/constraints/primary/top_higgs_ew/test_T019.py`.

Tests: T019 focused `12 passed`; full `python -m pytest tests/constraints/ -q` -> `582 passed`.