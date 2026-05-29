Changed L001 unevaluated handling: missing/invalid lepton dipole input now returns no prediction, `evaluated=False`, an `unevaluated_reason`, and the explicit `NOT EVALUATED — ...` note.

Evaluated proxy path still carries `NEEDS-HUMAN-PHYSICS`; real BR-vs-`1.5e-13` verdict unchanged.

Tests added for unevaluated missing/invalid cases.

Pytest: L001 `10 passed`; requested `tests/constraints/` `126 passed in 3.36s`.