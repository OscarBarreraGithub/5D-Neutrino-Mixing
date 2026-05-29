Plan executed: scaffold/YAML study, adapter proxy, CR009 constraint, focused tests, full suite.

Physics: CR009 compares `Lambda_RS` proxy to high-mass DY `llqq` contact-scale lower limits from `CR009.yaml`.

Reused/built: reused `collider_resonance` lower-bound comparison; appended only a small DY contact-scale wrapper.

Budget: active HARD limit is `36.4 TeV`, `PDG2025:CR009:CMS_RR_constructive_range_endpoint`, PDG 2025 quoting CMS arXiv:2103.02708.

Validation: no SM value applies; independent core check at `Lambda_proxy=40.0 TeV` gives ratio `36.4/40.0 = 0.91`, pass.

NEEDS-HUMAN-PHYSICS: full helicity/operator matching, couplings, interference, EFT validity, and binned DY likelihood remain proxy-only.

Files changed: `flavor_catalog_constraints/physics_adapters/collider_resonance.py`, `flavor_catalog_constraints/primary/collider_rs/CR009.py`, `tests/constraints/primary/collider_rs/test_CR009.py`.

Tests: `test_CR009.py` `8 passed`; full `python -m pytest tests/constraints/ -q` `898 passed`.

Worktree note: unrelated untracked orchestration/other constraint files are present and were left untouched.