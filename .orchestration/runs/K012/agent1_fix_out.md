Implemented K012 Im projection.

Formula: `BR_SD(K_S->mumu)_{l=0} = (tau_KS/tau_KL) kappa_mu [Im(lambda_t Y_t + Y_NP)/lambda^5 - Im(lambda_c) P_c(Y)/lambda]^2`, i.e. `Im[-lambda_c Y_c + lambda_t C10]^2`.

SM SD: `1.864194959306828e-13`; YAML `5e-12` remains total-SM context, not SD validation.

Real-only NP test: `left_sd=0.2` real gives `1.864194959306828e-13`, passes, ratio `8.877119e-04`.

RS proxy still flags `NEEDS-HUMAN-PHYSICS`.

Pytest: `tests/constraints/primary/kaon/test_K012.py -q` -> `12 passed`; `tests/constraints/ -q` -> `464 passed`.

Source checked: https://arxiv.org/abs/2104.06427.