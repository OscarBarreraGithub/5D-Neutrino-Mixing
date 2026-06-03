Plan: rewire only T010/T011/T012 to 3a `rs_ew_couplings`, add explicit old-point degradation, update focused tests.

Files: [T010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T010.py:486), [T011.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T011.py:663), [T012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T012.py:351).

Zbb path: `delta_g_L=+2.580074997922e-03`, `delta_g_R=-1.636688349408e-03`; shifted `gL=-4.202532583354e-01`, `gR=+7.552997831726e-02`.

Zcc path: `delta_g_L=-4.269875187085e-05`, `delta_g_R=+1.969187886049e-03`; shifted `gL=+3.456239679148e-01`, `gR=-1.523641454473e-01`.

SM-limit: max `delta_g=0`; `predicted==sm_prediction` for T010/T011/T012; ratios `0.996014155971`, `0`, `0.089724223634`.

Absent old-style point: T010/T011/T012 all `passes=True`, `predicted=None`, `ratio=None`, `diagnostics.evaluated=False`, `missing_extra=rs_ew_couplings`.

Tests: focused `27 passed`; full `python -m pytest tests/ -q` -> `1646 passed, 1 skipped in 808.06s`.

Review: Codex local review APPROVE; Opus approval not claimable here because no Opus tool is available and the independent subagent review attempt hit 429.

P3B-DONE