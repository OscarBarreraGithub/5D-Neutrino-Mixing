# B007 FIX PASS (agent1) — you already implemented B007; the code reviewer found ONE blocker. Fix it, keep everything else intact, re-run tests. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Constraint: B007 — B0→ee / Bs→ee (secondary/beauty). Your files:
- flavor_catalog_constraints/secondary/beauty/B007.py
- flavor_catalog_constraints/physics_adapters/rare_b_electronic.py
- tests/constraints/secondary/beauty/test_B007.py

## BLOCKER to fix (agent3, code review)
In `_limit_anchor_from_value_row(...)` (around B007.py:257), the mismatched-anchor probe does NOT fail loudly: when `load_anchor` returns an anchor whose `block_key` differs from the requested `block_key`, the code accepts it silently. This violates the loud-fail anchor contract.

**Fix:** assign the returned anchor, verify `anchor.block_key == block_key` (the requested block), and raise `AnchorError` on mismatch (mirror how other constraints validate the resolved anchor — e.g. value_id/units mismatch already raise). Then ADD a regression test near tests/constraints/secondary/beauty/test_B007.py:208 that monkeypatches `load_anchor` to return a wrong `block_key` and asserts `AnchorError` is raised.

## NIT (optional, low priority — only if trivial)
Virtual-anchor loading temporarily monkeypatches `anchor_scaffold.load_pdg_block` (B007.py:241). Works and restores in finally; leave as-is unless a clean scaffold helper is obvious. Do NOT over-engineer.

## Rules
- Do NOT change the physics (agent2 PHYSICS-OK: SM Bs→ee=8.540e-14, Bd→ee=2.450e-15 must stay). Do NOT modify the reused rare_b_dilepton muon functions.
- Keep numeric ConstraintResult fields real floats; complex Wilsons in diagnostics.
- Re-run: `python -m pytest tests/constraints/secondary/beauty/test_B007.py -q` AND `python -m pytest tests/constraints/ -q`. Report counts.

OUTPUT (≤12 lines): what you changed (file:line), the new regression test, pytest counts. End with: FIX-DONE.
