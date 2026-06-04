Plan: helper/formula, builder opt-in, T010/T011 diagnostics, gate tests.

Files: `quarkConstraints/rs_ew_couplings.py`, `flavor_catalog_constraints/rs_ew_builder.py`, `T010.py`, `T011.py`, and related tests.

Numbers: `delta_g_L_b,ferm=2.3518892841117975e-05 >0`, `delta_g_R_b,ferm=-5.49194949716315e-04 <0`.

Checks: `m_b->0=(0.0, -0.0)`, double `Lambda_IR` scaling `(0.25, 0.25)`.

SM-limit: T010 `pred=SM=0.21562000000000003` with catalog pull ratio `0.9960141559711674`; T011 `pred=SM=0.10327886042219099`, ratio `0.0`.

Top-partner/custodial `Zb_L` term is not computed; metadata flags `custodial_toppartner_zbL_needs_human=True`.

Review: Codex subagent `APPROVE`; Opus/Anthropic tool or installable connector is not available in this environment, so Opus approval could not be obtained.

Tests: focused `23 passed`; full `python -m pytest tests/ -q` `1682 passed, 1 skipped in 806.71s`; `git diff --check` clean.

P6A-DONE