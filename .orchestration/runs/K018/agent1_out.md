Plan executed: load K018 anchors, build core extraction, wrap via adapter, implement constraint, add tests, run suite.

Physics: K_l3 uses `|V_us| = (|V_us| f_+(0)) / f_+(0)` from K018.yaml: PDG2025 `0.21656(35)` and FLAG2024 `f_+(0)=0.9698(17)`.

Built: `quarkConstraints/ckm_extraction.py` and `physics_adapters/ckm_extraction.py`. Reused: scaffold `load_anchor` via the EW002/EW003 list-anchor pattern.

Budget: SOFT, `0.00053` from `K018.yaml` `PDG2025:K018:Vus_from_Kl3.uncertainty_total`.

Validation: extracted `|V_us| = 0.223303773974`, pull vs YAML derived `|V_us|=0.22330` is `0.00712070569 sigma`.

NEEDS-HUMAN-PHYSICS: flagged for missing RS charged-current W/W'/KK, right-handed-current, lepton/G_F, and radiative/isospin matching inputs.

Files changed: `quarkConstraints/ckm_extraction.py`, `flavor_catalog_constraints/physics_adapters/ckm_extraction.py`, `flavor_catalog_constraints/primary/kaon/K018.py`, `tests/constraints/primary/kaon/test_K018.py`.

Tests: `python -m pytest tests/constraints/ -q` -> `482 passed in 13.53s`.