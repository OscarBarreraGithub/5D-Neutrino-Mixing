Plan complete: studied scaffold/K001/K004/K019/K009/K010, implemented K020, added tests, ran full suite.

Physics: pure-NP HARD bound for `K+ -> pi+ mu+ e-`; reused rare-kaon `s->d ll`/K019 LFV proxy, built full-q2 `K->pi` form-factor three-body adapter.

Budget: `1.3e-11` from `K020.yaml` row `PDG2025:K020:Kplus_piplus_mup_em_limit`.

Validation: SM LFV rate `0.0`; sample independent recomputation matches `BR = 1.323883001e-12`, ratio `0.101837`.

NEEDS-HUMAN-PHYSICS: off-diagonal e-mu lepton coupling is not on `ParameterPoint`; requires explicit `lepton_mass_basis_couplings` proxy.

Files: `flavor_catalog_constraints/physics_adapters/rare_kaon_lfv_semileptonic.py`, `flavor_catalog_constraints/secondary/kaon/K020.py`, `tests/constraints/secondary/kaon/test_K020.py`.

Tests: `test_K020.py` 11 passed; `python -m pytest tests/constraints/ -q` 1002 passed.