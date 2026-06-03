Added `tests/test_rs_ew_phase4a.py:144` LFV toy; asserts `z_delta_g_L_e[e,mu]=1.212164466423143e-05+0j`.
LFV Wilsons asserted: `s_to_d e_mu C9=2.3656625045856276e-04+6.666935783109433e-05j`, `C10=-2.3656625045856276e-04-6.666935783109433e-05j`.
Scaling check: `rate(6000)/rate(3000)=0.0625`, matching `(3000/6000)^4=0.0625`.
Added neutral contact write-lock assertion at `tests/test_rs_ew_phase4a.py:328`.
Fixed Phase-4a docstrings at `quarkConstraints/rs_ew_couplings.py:50` and `:409`.
Pytest: `tests/test_rs_ew_phase4a.py -q` -> `7 passed`; `python -m pytest tests/ -q` -> `1653 passed, 1 skipped`.
P4A-FIX-DONE