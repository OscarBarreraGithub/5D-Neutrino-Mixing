# DUAL-REVIEW of the scaffold-hardening fix (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Do NOT modify code — findings + verdict. A prior codex review found 3 BLOCKERS in the scaffold; a fix was applied. Verify the fix is correct, complete, and BACKWARD-COMPATIBLE (did not weaken anything or break the 103 constraints).

Changed files: flavor_catalog_constraints/base.py, anchors.py, registry.py, TEMPLATE.py, tests/constraints/test_contract.py.

VERIFY (probe yourself):
1. **NaN/Inf + type guard** (base.py ~110/171): `ConstraintResult` numeric fields now reject NaN AND Inf (not just complex); `passes` must be bool; `severity` must be `Severity`. Probe: construct results with NaN ratio, Inf budget, non-bool passes, bad severity → each must raise. Confirm complex still rejected. Confirm a NORMAL valid result still constructs.
2. **Anchor validation** (anchors.py ~181/271): new OPTIONAL `expected_value_id`/`expected_block_key`/`expected_units`/`expected_confidence_level` params raise `AnchorError` on mismatch WHEN PROVIDED, and are NO-OP when omitted (so the 103 existing calls are unaffected). Probe a wrong block_key with the param set (raises) and without (unchanged old behavior).
3. **Immutable extras** (base.py/point_builder): mutating `ParameterPoint.extras`/`raw` now raises. Probe. Confirm reads still work and constraints still evaluate.
4. **reset_for_tests** (registry.py ~168): after `reset_for_tests(); discover()` the registry repopulates (103, not 0). Probe.
5. **No regression / no weakening**: confirm the changes didn't silently relax any existing guard; run `python -m pytest tests/constraints/ -q` (expect ≥1061 passed) and report counts. Confirm the new tests in test_contract.py actually exercise the failure paths (not trivial).
6. **Backward-compat**: spot-check 2-3 existing constraints still load anchors + evaluate unchanged.

OUTPUT (<=16 lines): numbered findings (BLOCKER/SHOULD-FIX/NIT) with file:line + the ACTUAL probe results, pytest counts. End with: SCAFFOLD-FIX-OK or SCAFFOLD-FIX-NEEDS-FIXES.
