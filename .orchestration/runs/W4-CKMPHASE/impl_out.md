Plan executed: core helper, adapter export, B002/B004 wiring, focused tests, full tests.

Added `ckm_phases_from_matrix()` plus repo-default CKM helpers in `quarkConstraints/ckm_extraction.py`. Computed from `ModernDefaultCKMTarget`: `2β = 45.099923°`, `sin2β = 0.708338893`, `φ_s = -0.037945103 rad`, consistent with PDG windows.

B002 SM-limit now gives `S_psiK_S = 0.708338893`; B004 SM-limit gives `φ_s = -0.037945103 rad`. ΔF=2 complex `M12` NP phase handling is unchanged.

Removed the CKM-phase `NEEDS-HUMAN-PHYSICS` diagnostics from B002/B004; unrelated flags elsewhere are untouched.

Tests: focused CKM/B002/B004 `22 passed`; full suite `1696 passed, 1 skipped`.

Review: Codex subagent `APPROVE`. True Opus review could not be run because no Opus/Claude model or CLI is exposed in this environment, so I cannot honestly claim Opus approval.

W4-DONE