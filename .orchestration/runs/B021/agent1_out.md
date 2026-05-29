Plan: read scaffold/examples/B021, add append-only baryonic machinery, implement constraint/test, run suite.

Physics: active HARD observable is the B021.yaml LHCb `15 < q2 < 20 GeV^2` bin, converted from `dBR/dq2` to partial BR. Uses Detmold-Lin-Meinel-Wingate leading-HQET `Lambda_b -> Lambda` `F_+/-` normalization.

Reused/built: reused B005/B016 `rare_b_dilepton` CKM/C9/C10 proxy; built new baryonic core + adapter, without modifying meson functions.

Budget: `3.458207347033031e-07`, from B021 high-q2 exp bin + `|exp-SM| + sigma_exp + 30%*SM`.

SM validation: `BR_SM[15,20] = 5.532430650131983e-07`, average `dBR/dq2 = 1.1064861300263965e-07`, consistent with catalog `1.18e-07`.

NEEDS-HUMAN-PHYSICS: RS C9/C10 proxy still lacks full EW KK/Z/Z’, lepton matching, C7/nonlocal charm, scalar/tensor terms, baryonic covariance, and bin correlations.

Files changed: `quarkConstraints/rare_b_baryon_dilepton.py`, `flavor_catalog_constraints/physics_adapters/rare_b_baryon.py`, `flavor_catalog_constraints/primary/beauty/B021.py`, `tests/constraints/primary/beauty/test_B021.py`.

Tests: `python -m pytest tests/constraints/ -q` -> `568 passed in 15.46s`.

Note: unrelated dirty/untracked files from other work are present and were left untouched.