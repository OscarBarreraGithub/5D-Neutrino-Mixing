# SCAFFOLD HARDENING (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. First a SHORT plan, then implement. The constraint-framework scaffold was retro-reviewed: codex found 3 BLOCKERS + 2 SHOULD-FIX; Opus concurred they are real hardening gaps. Fix them to make the framework bulletproof for a 100M-point scan. A codex reviewer and an Opus reviewer will check your work (both must approve).

ABSOLUTE CONSTRAINTS (correctness, not optional):
- BACKWARD-COMPATIBLE: `python -m pytest tests/constraints/ -q` MUST stay green (currently 1054 passed). Do NOT force-refactor the 103 existing constraints. New validation must be OPTIONAL/additive so existing `load_anchor` calls keep working.
- Keep the numeric-real contract; keep auto-discovery + per-constraint exception isolation in `evaluate_all` intact.

FIXES:
1. **base.py `ConstraintResult.__post_init__`** — numeric fields must be real AND FINITE floats: keep complex rejection, ADD `math.isfinite` rejection of NaN/Inf (a NaN `ratio` silently bypasses the `ratio<=1` veto). Also enforce `passes` is `bool` and `severity` is a `Severity`. Raise `TypeError`/`ValueError` with a clear message. Add/extend a unit test proving NaN/Inf/non-bool/bad-severity now raise.
2. **anchors.py `load_anchor`** — add OPTIONAL, backward-compatible keyword params (e.g. `expected_value_id`, `expected_block_key`, `expected_units`, `expected_confidence_level`) that, WHEN PROVIDED, raise `AnchorError` on mismatch. Model `value_id`/CL fields if the loader doesn't already surface them. When NOT provided, behavior is unchanged (so the 103 constraints keep passing). Add tests probing each mismatch raises.
3. **point_builder.py `ParameterPoint`** — make `extras` (and `raw`) immutable to consumers (e.g. expose as `MappingProxyType`) so a constraint cannot mutate a shared point during a long scan loop. Reads must be unaffected; existing constraints already only read, so the suite must stay green. Add a test that mutation attempts raise.
4. **registry.py `reset_for_tests()`** — currently leaves modules cached so a post-reset `discover()` returns 0 constraints. Make reset correct (re-discovery repopulates) OR make it safe/clearly-scoped. Add a test that reset→discover repopulates the registry.
5. **TEMPLATE.py** — update the anchor-loading example to demonstrate the new optional validation params, so future constraints inherit defense-in-depth.

OUTPUT (<=16 lines): the short plan (1-3 lines), then changes (file:line), a PROBE showing each new guard now rejects bad input (NaN ratio, wrong block_key, extras mutation, reset→discover), and pytest counts (must be ≥1054 passed). End with: FIX-DONE.
