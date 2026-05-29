Touched `quarkConstraints/rare_b_nunu.py`, `flavor_catalog_constraints/primary/beauty/B022.py`, and `tests/constraints/primary/beauty/test_B022.py`.

Fixed charged-B formula: `BR = BR_LD + (BR_SM_total - BR_LD) * R_K`, with `BR_LD=6.09e-7`; SM validates at `5.58e-6`.

Corrected benchmark: `R_K=4.002115` now gives `BR=2.050351497441e-5` instead of naive `2.233180316982e-5`.

B022 value anchors now route through `load_anchor`; test forces missing/mismatched `load_anchor` and gets `AnchorError`. `NEEDS-HUMAN-PHYSICS` Z-like proxy flag retained.

Tests: `test_B022.py -q` -> `11 passed`; `tests/constraints/ -q` -> `160 passed`.