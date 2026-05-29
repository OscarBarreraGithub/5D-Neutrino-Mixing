Plan done: read scaffold/examples/YAML, implemented K012, added tests, ran full suite.

Physics: K012 reuses K006 `s -> d mu+mu-` short-distance Y-function machinery and rescales by `tau_KS/tau_KL = 0.0017501954652071933`. Total rate remains LD-dominated.

Reused vs built: reused `rare_kaon_dilepton`; added append-only adapter `rare_kaon_dilepton_kshort_mumu.py`.

Budget: HARD limit `2.1e-10` from `K012.yaml` `PDG2025:K012:headline_limit`.

SM validation: K012 SD proxy gives `1.4368757234288928e-12`; YAML has only approximate total SM context `5e-12`, not an SD anchor.

NEEDS-HUMAN-PHYSICS: full K_S/K_L projection, LD two-photon treatment, EW KK/Z/Z' and muon-sector matching remain proxy gaps.

Files changed: `flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton_kshort_mumu.py`, `flavor_catalog_constraints/primary/kaon/K012.py`, `tests/constraints/primary/kaon/test_K012.py`.

Tests: `python -m pytest tests/constraints/ -q` → `440 passed in 6.66s`.