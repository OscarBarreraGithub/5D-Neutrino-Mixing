Changed `rare_b_dilepton.py` to compute `A_DeltaGamma` from the actual complex `P,S` amplitudes and use it in the time-integrated BR; B005 diagnostics/tests updated, RS proxy flag retained.

SM BR(B_s→μμ): `3.6483515801441554e-09`, with `A_DeltaGamma=1`.

Range now: `A_DeltaGamma=-1..+1`; for `y_s=0.0645`, time factor `0.9394081728511038..1.068947087119188`.

Pytest: `test_B005.py -q` = `11 passed`.

Full `tests/constraints/ -q` is not green due unrelated existing failures: `B022` import `NameError: load_pdg_block`; ignoring that test shows `146 passed, 2 failed` with unrelated `B011` plus contract discovery of `B022`.