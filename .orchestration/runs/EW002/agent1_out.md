Plan executed: inspected scaffold/examples, implemented EW002, added focused tests, ran full constraint suite.

Physics: data-level CKM first-row sum vs SM unitarity target, `Delta_CKM = sum - 1`.

Sources: EW002.yaml PDG `first_row_sum`, `Vud_superallowed`, `Vus_kaon_average`.

Machinery: built new `physics_adapters/ckm_unitarity.py`; no RS physics core existed.

Budget: `0.0007` from `EW002.yaml` `combined_quoted_elsewhere`.

Validation: sum `0.9983`, `Delta_CKM = -0.0017`, pull `2.428571`.

Severity/result: `SOFT`; current data fails the 1-sigma budget as a non-vetoing SM tension.

Gap: `NEEDS-HUMAN-PHYSICS` flagged for missing grounded RS charged-current/G_F matching; EW002 lacks standalone `|Vub|` block, so central uses the quoted PDG sum.

Files changed: `flavor_catalog_constraints/primary/top_higgs_ew/EW002.py`, `flavor_catalog_constraints/physics_adapters/ckm_unitarity.py`, `tests/constraints/primary/top_higgs_ew/test_EW002.py`.

Tests: `python -m pytest tests/constraints/ -q` -> `125 passed`.