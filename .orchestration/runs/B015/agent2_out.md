1. NIT - Δm/CP amplitude check is N/A: B015 is a ΔB=1 inclusive BR, not an M12 observable; verdict uses a positive rate shape with `|C_i|^2` and `Re(C7 C9*)` at `quarkConstraints/rare_b_dilepton.py:1319-1335`, then a BR pull at `flavor_catalog_constraints/primary/beauty/B015.py:491-495`.

2. NIT - QCD running is appropriate for this ΔB=1 proxy: no ΔF=2 `*_with_running`/2 GeV path should be used. C7/C8 runs to `mu_b=4.8 GeV` via `compute_bsgamma_wilsons` at `rare_b_dilepton.py:1431-1435` and `bsgamma.py:346-357`; for `MKK=3 TeV`, `U77=0.512`, `U78=0.119`, so `C7_low=0.631*C7_matching` if `C8=C7`.

3. NIT - Budget is uncertainty-aware and defensible under the prompt’s proxy rule: low-q2 exp `1.58±0.37e-6`, SM μμ `1.62±0.09e-6`, residual `0.04e-6`, combined sigma `0.381e-6`, proxy `0.486e-6`, total `0.9068e-6`; see `B015.py:408-445`.

4. NIT - Anchors match YAML/snapshots: HFLAV total `5.84±0.69e-6`, low-q2 `1.58±0.37e-6`, SM μμ `1.62±0.09e-6`; see `B015.yaml:95-195`, snapshots `hflav...txt:11-16`, `huber...txt:20-28`, loader `B015.py:337-377`.

5. NIT - SM validation is correct: normalized inclusive SM prediction is exactly `1.62e-6` using C7/C9/C10 shape; see `rare_b_dilepton.py:1450-1478`, `B015.py:509-514`, `test_B015.py:244-275`.

6. NIT - Severity, units, and limitations are consistent: HARD, BR dimensionless with `10^-6` scaling, q2 in `GeV^2`, masses/scales in GeV, and `NEEDS-HUMAN-PHYSICS` is explicit; see `B015.py:1-38`, `B015.py:502-504`, `B015.py:677-684`.

PHYSICS-OK