Plan: implemented L010 append-only with a tau-e adapter, constraint, and tests.

Physics: pure-NP `BR(tau -> 3e)` with SM fixed to `0.0`; dipole + vector contact + constructive dipole-contact interference envelope.

Reused/built: reused `quarkConstraints.lfv_three_body`; built only `lfv_three_body_taue.py` to pin `initial=tau`, `final=e` and reject mu-e/tau-mu spurion aliases.

Budget: `2.7e-8` from `L010.yaml` `primary_current_limit` / Belle 2010 PDG-selected limit. Belle II is retained as prospect context, not the current veto budget.

Validation: tau->e dipole conversion factor checked as `0.0104975811678193`; proxy numerics independently recomputed against the shared core after tau-e normalization.

NEEDS-HUMAN-PHYSICS: full RS lepton neutral-current, EW KK/Z/Z', box, and tau-dipole matching are still proxy inputs and flagged in diagnostics.

Files changed: [lfv_three_body_taue.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/lfv_three_body_taue.py), [L010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charged_lepton/L010.py), [test_L010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/charged_lepton/test_L010.py).

Pytest: `test_L010.py` 16 passed; charged-lepton family 52 passed; `tests/constraints/` excluding unrelated `test_B025.py` 606 passed. Full `tests/constraints/` is blocked by an unrelated B025 collection-time `AnchorError` on missing `joint_fit_context.rd_value`, left untouched per isolation rules.