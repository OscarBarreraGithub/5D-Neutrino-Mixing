You are agent1 fixing B003 (Delta m_s). Read the code review: .orchestration/runs/B003/agent3_out.md (physics .orchestration/runs/B003/agent2_out.md was PHYSICS-OK). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.
FIX:
1. (SHOULD-FIX) Top-level ConstraintResult.sm_prediction and .experimental currently hold half-M12 values, not the OBSERVABLE Delta m_s. Set them to the Delta m_s values (sm_delta_m_gev / experimental_delta_m_gev); keep the M12 / half-M12 quantities in diagnostics. predicted should remain the NP quantity consistent with the budget; document units in notes.
2. (SHOULD-FIX) Validate the anchor before budget construction: require units == "ps^-1" and observable matches "Delta m_s"; raise AnchorError (loud) on mismatch. Add a probe/test.
3. (NIT) Fix the misleading error message for the missing auxiliary anchor (it says pdg_or_equivalent; make it auxiliary-specific).
Keep changes ISOLATED to B003.py and its test (and B003.yaml only if needed). Run `python -m pytest tests/constraints/primary/beauty/test_B003.py -q` and `python -m pytest tests/constraints/ -q` until green.
OUTPUT (<=12 lines): changes, the new top-level field values vs diagnostics, pytest counts, git diff --stat. No full files.
