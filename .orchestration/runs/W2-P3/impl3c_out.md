Plan: inspect 3c/3b, rewire T014 only, replace proxy tests with rigorous/SM/absent coverage, run suite/review.

Rewired: [T014.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:70), evaluate path at [T014.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:435), adapter call at [T014.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/secondary/top_higgs_ew/T014.py:462).

Rigorous bs: `z_delta_g_L_d[1,2]=(-2.0387380225097076e-04-2.2036569174673724e-08j)`, `z_delta_g_R_d[1,2]=(-1.9683644415245964e-04+0j)` -> `BR=1.32030759087837e-07`.

SM-limit: `BR=0.0`, pass. Absent old-style point: `evaluated=False`, `missing_extra=rs_ew_couplings`, non-vetoing.

Tests: T014 stayed `11 -> 11` items; proxy numerics replaced by rigorous builder offdiag test; zero proxy test replaced by SM-limit builder test; missing path now old-style absent; proxy KK scaling replaced by direct-coupling no-rescale.

Pytest: focused `11 passed`; full `1646 passed, 1 skipped in 804.62s`.

Review: Codex sub-agent APPROVE. Opus approval was not obtainable because no Opus/Anthropic connector or model is exposed/installable in this session, so I am not claiming Opus approval.

P3C-DONE