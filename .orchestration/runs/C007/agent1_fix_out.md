Option b for q2: relabeled C007 as a smooth full-q2 proxy; no LHCb window/acceptance mask is implied, and NEEDS-HUMAN-PHYSICS is retained.

Form factors: added YAML/snapshot anchors and `load_anchor` wiring for f_+(0), m_pole, a, b; exposed 8.96% f_+(0) normalization uncertainty in diagnostics. Source: https://doi.org/10.1140/epjc/s10052-017-4888-4

Module: `quarkConstraints/rare_charm_semileptonic.py` now emits explicit q2-treatment/window-applied diagnostics.

Tests: `test_C007.py -q` -> 11 passed; `tests/constraints/ -q` -> 352 passed.