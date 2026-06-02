1. NIT: ΔF=2 M12 amplitude/running checks are N/A for B008; this is a ΔB=1 leptonic BR. Relevant rate uses `|P|^2+|S|^2` with `C10-C10'` and `m_l=m_tau`: `quarkConstraints/rare_b_dilepton.py:468`, `:485`, tau override `flavor_catalog_constraints/physics_adapters/rare_b_tauonic.py:57`.
2. NIT: No `*_with_running` path is expected/available for this C10-only leptonic proxy; C10 has no QCD running effect in this v1 rate, and C9 is diagnostic-only for leptonic BR: `rare_b_dilepton.py:559`, `:623`.
3. NIT: Budget is defensible. Code uses PDG one-mode 95% CL limits as HARD budgets and also enforces total BR <= limit: Bs `6.8e-3`, Bd `2.1e-3`; SM fractions are tiny, `1.14e-4` and `1.06e-5` of limits: `B008.py:396`, `:466`.
4. NIT: Anchors match yaml/snapshots: SM Bs `7.73(0.49)e-7`, Bd `2.22(0.19)e-8`; formula gives Bs `7.737985e-7`, Bd `2.190973e-8`: `B008.yaml:83`, `:104`, `:168`.
5. NIT: Severity/units/notes are appropriate: HARD, branching-fraction units, and RS C9/C10 proxy is marked NEEDS-HUMAN-PHYSICS: `B008.py:68`, `:84`, `:483`, `:641`.

PHYSICS-OK