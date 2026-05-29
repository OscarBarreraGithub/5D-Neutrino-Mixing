Plan executed: inspect examples/scaffold -> add K010 adapter/constraint -> add tests -> run full suite.

Physics: `BR(K_S -> pi0 e+e-) = 5.2e-9 |a_S^eff|^2`, with `|a_S|` derived from K010 full-rate YAML. Sources: K010 PDG/NA48 YAML data from https://pdgprod.lbl.gov/pdgprod/pdgLiveJson/api/listings/S012.10 and a_S formula from https://arxiv.org/abs/hep-ph/0308008.

Reused K008 rare-kaon-dilepton y7V/y7A proxy; built only a new K_S a_S adapter. RS matching is flagged `NEEDS-HUMAN-PHYSICS`.

Budget: HARD, direction-aware `+2.9e-9 / -2.4e-9`; max budget `2.9e-9`, from `K010.yaml extrapolated_total_rate`.

Validation: zero-NP a_S baseline `5.80e-9`; `|a_S| = 1.056117709`. Safe point ratio `0.00296`; excluded test point ratio `50.61`.

Files changed: [adapter](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton_ks.py), [constraint](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/kaon/K010.py), [test](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/kaon/test_K010.py).

Tests: `python -m pytest tests/constraints/primary/kaon/test_K010.py -q` -> 10 passed. `python -m pytest tests/constraints/ -q` -> 407 passed.