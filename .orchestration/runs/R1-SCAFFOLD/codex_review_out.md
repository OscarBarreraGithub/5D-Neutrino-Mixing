1. BLOCKER `flavor_catalog_constraints/base.py:158`: `ConstraintResult` only rejects complex. Probe: complex raised `TypeError`; `NaN`, `Inf`, string numeric, bad `severity`, and non-bool `passes` were accepted. Fix: enforce `Severity`, `bool`, `float`, and `math.isfinite`.
2. BLOCKER `flavor_catalog_constraints/anchors.py:173,260`: `load_anchor` cannot validate `value_id`, expected `block_key`, units, or CL. Probe: missing block raised `AnchorError`, but no `value_id` + wrong units/CL accepted as `block_key='wrong_existing_block'`. Fix: add typed `value_id/confidence_level` fields and expected-* validation.
3. BLOCKER `flavor_catalog_constraints/base.py:104`, `flavor_catalog_constraints/point_builder.py:70`: `ParameterPoint` is shallow-frozen only. Probe: `p.raw = ...` raised `FrozenInstanceError`; `p.extras[...] = ...` succeeded. Fix: wrap extras in immutable mapping.
4. SHOULD-FIX `flavor_catalog_constraints/registry.py:181`: `reset_for_tests()` leaves imported modules cached; probe after reset rediscovery returned `0` constraints with `0` import failures. Fix: unload leaf modules or make reset private/safer.
5. SHOULD-FIX `flavor_catalog_constraints/TEMPLATE.py:79`: template calls bare `load_anchor`, so copied constraints won’t enforce units/CL/value-id unless authors add bespoke checks. Fix after anchor API hardening.

Probe pass points: registry smoke found `103` constraints, `0` import failures; exception isolation worked (`K001` failed with `RuntimeError`, `K002` continued); discovery order probe was stable 3/3. `build_from_quark_couplings` probe was deterministic for a dummy `M_KK`.

Tests: `python -m pytest tests/constraints/ -q` → `1054 passed in 34.01s`.

SCAFFOLD-NEEDS-FIXES.