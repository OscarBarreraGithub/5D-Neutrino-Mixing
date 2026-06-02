Plan: inspected scaffold/K001/K004/K019 YAML, verified K019 registration, added secondary tests, ran full constraints suite.

Physics: `K_L -> e mu` is treated as SM `0.0`, pure-NP HARD upper bound using `s -> d l l` rare-kaon dilepton `kappa_mu` normalization with explicit LFV `e-mu` proxy.

Reused/built: reused `rare_kaon_dilepton`; K019 uses append-only `rare_kaon_lfv_dilepton` adapter for unequal-lepton LFV two-body rate.

Budget: `BR < 4.7e-12` from `flavor_catalog/processes/secondary/kaon/K019.yaml` PDG 2025 / BNL E871 90% CL.

Validation: independent test recomputes nonzero LFV BR from underlying core Wilsons; safe point `1.7882686599234883e-12`, ratio `0.38048`; excluded point ratio `38.048`.

NEEDS-HUMAN-PHYSICS: full off-diagonal charged-lepton neutral-current RS matching is missing; explicit `lepton_mass_basis_couplings` proxy is flagged in diagnostics.

Files: added `tests/constraints/secondary/kaon/test_K019.py`; K019 implementation/adapter are present at requested paths.

Tests: `tests/constraints/secondary/kaon/test_K019.py` -> `11 passed`; `tests/constraints/` -> `959 passed`.