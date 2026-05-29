Changed `test_T003.py`: removed hand-coded BR formula; cross-check now recomputes dipoles from couplings and calls `quarkConstraints.top_fcnc.photon_dipole_branching_fraction`.

Cross-check: `Constraint.evaluate=6.0684210348288108e-06`, core `BR=6.0684210348288108e-06`, `diff=0`, `ratio=0.40188218773700735`.

Pytest: `test_T003.py -q` -> `9 passed`; `tests/constraints/ -q` -> `352 passed`.