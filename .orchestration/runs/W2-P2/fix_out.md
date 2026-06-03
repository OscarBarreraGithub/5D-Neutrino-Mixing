Root fix: [rs_ew_spectrum.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_ew_spectrum.py:140) now scans ordered gauge-NN cross-product brackets; build uses it at line 289 and rejects duplicate/regressing roots at line 256.

Corrected roots x1-x6: 2.450509663814, 5.567547286728, 8.701957144582, 11.840260538785, 14.980018189298, 18.120466878758.

Corrected Q=1024 overlaps: a(0.2)=21.860547338064457, a(0.65)=-1.457139933026502.

New tests assert strict increasing/unique tower and independent bracketed overlap recomputation: [test_rs_ew_spectrum.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/test_rs_ew_spectrum.py:99).

Pytest: `tests/test_rs_ew_spectrum.py -q` = 7 passed; `tests/ -q` = 1645 passed, 1 skipped.

PHASE2-FIX-DONE