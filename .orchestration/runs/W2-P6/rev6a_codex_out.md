1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: modified tracked files are isolated to builder/couplings + T010/T011 + tests; two untracked `.orchestration` notes are present.
4. Formula checked against [Casagrande arXiv:0807.4937](https://arxiv.org/abs/0807.4937): signs and row/column `B_d/B_Q` assignment match the requested minimal non-custodial ZMA convention.
5. Recomputed sample: `dgL=+2.3518892841117975e-05`, `dgR=-5.49194949716315e-04`; signs are `dgL>0`, `dgR<0`.
6. Scale convention confirmed: `Lambda_IR=spectrum.lambda_ir_gev=1224.2351231258115 GeV`, not physical first gauge mass `3000.000000000015 GeV`.
7. `m_b=0`: fermion piece `0.0/-0.0`, mixed matrices equal gauge-only; doubling `Lambda_IR` gives `0.25` for both shifts.
8. No double-count: total equals gauge plus fermion; e.g. `L22 0.00259701058 -> 0.00262052948`, `R22 -0.00166042276 -> -0.00220961771`.
9. Top-partner/custodial term is deferred: flag true, no `delta_g_L_b_top_partner` silently added.
10. T010/T011 read gauge+fermion `z_delta_g`; SM-limit recovers `T010=0.21562000000000003`, `T011=0.10327886042219099`, real float fields.
11. `pytest tests/ -q`: `1682 passed, 1 skipped in 802.13s`.
P6A-OK