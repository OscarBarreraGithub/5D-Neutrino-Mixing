Scoped edits: `C002.py` and `test_C002.py`; existing deltaf2 D0 complex-`M12` helper used unchanged.

Changes: budget now comes from HFLAV `phiM=Arg(M12)` snapshot interval; q/p and `phi_D` route through `load_anchor`; raw NP argument diagnostic relabeled; NEEDS-HUMAN-PHYSICS caveat kept.

New HARD budget: `phiM_95=(-1.48, 1.35) deg`, dimensionless `0.02582800048591857`, so `|Im M12^NP| <= 8.4741669594298825e-17 GeV`.

Cross-check safe: `M12=-9.852835695838563e-17+0j`, `pred=0`, `ratio=0`, pass.

Cross-check finite: `M12=-8.488938395120150e-17-7.679675331923575e-17j`, `pred=0.023406508174104161`, `ratio=0.906245459723659`, pass.

Cross-check excluded: `M12=-3.395575358048060e-14-3.071870132769430e-14j`, `pred=9.3626032696416637`, `ratio=362.49818388946358`, fail.

Pytest: `test_C002.py -q` = 8 passed; `tests/constraints/ -q` = 91 passed.