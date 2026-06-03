Plan executed: inspect scoped consumers, add adapter, rewire five IDs, replace proxy tests, run review/tests.

Rewired: `EW002`, `K018`, `K017`, `B009`, `B025` via new `physics_adapters/charged_current.py` plus focused tests/helper.

Example: K017 now uses `R_K = R_K^SM |1+epsilon_us^e|^2/|1+epsilon_us^mu|^2`; test fixture `epsilon_e=0.1`, `epsilon_mu=0` gives `1.21 * R_K^SM`.

SM-limit universal-c recovers committed SM: epsilon zero, K017/B009/B025 unchanged; K018 returns `0.22330377397401527`.

Absent `rs_charged_current`: non-vetoing `evaluated=False`, `missing_extra="rs_charged_current"`, `predicted=None`, `ratio=None`.

v1 near-SM confirmed: sample K017/B025 LFU cancel; K018 shift `5.9e-6`; EW002 NP delta `5.25e-5`; third-family epsilon residuals ~`1e-3`.

Preserved: `EW002` remains `SOFT`; `B025` remains `PARTIAL`.

Removed old scoped proxy use: no `m_K^2/M_KK^2`, `m_B^2/M_KK^2`, or `m_b m_tau/M_KK^2` paths in rewired consumers/tests.

Test count: scoped affected tests now `48` (`+2` vs pre-rewire focused set).

Tests: focused `48 passed`; full final tree `1673 passed, 1 skipped`.

Review: independent Codex reviewer `APPROVE`; Opus/Claude tool was not available in this session, so no Opus approval is claimed.

P5B-DONE