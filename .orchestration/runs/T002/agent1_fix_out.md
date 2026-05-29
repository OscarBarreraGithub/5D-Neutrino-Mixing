Changed `test_T002.py`: removed the local manual width helper; cross-check now uses `quarkConstraints.top_fcnc.compute_top_z_fcnc_proxy` + `z_vector_branching_fraction`.
Cross-check: core BR = constraint predicted = `4.817352363171973e-07`; ratio = `7.769923166406408e-03`.
Proxy vectors: `X_L=9.239087104177778e-04+2.309771776044444e-04j`, `X_R=0+2.771726131253333e-04j`.
`python -m pytest tests/constraints/primary/top_higgs_ew/test_T002.py -q`: 9 passed.
`python -m pytest tests/constraints/ -q`: 297 passed.