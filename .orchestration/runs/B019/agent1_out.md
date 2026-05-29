Plan: B019 implemented as the B018-style LFU ratio for beauty K* mode, with adapter/core, constraint, and tests.

Physics/sources: `B019.yaml` HFLAV Dec. 2025 `R_K* central-q2 = 1.028 +/- 0.074`; BIP 2016 SM/QED `R_K* = 1.00 +/- 0.01_QED`.

Reused/built: reused `rare_b_dilepton` C9/C10 Wilson matching; added append-only `quarkConstraints/rare_b_kstar_dilepton.py` plus `physics_adapters/rare_b_kstar_dilepton.py`.

Budget: HARD `|1.028 - 1.00| + 0.074 + 0.01 = 0.112`, sourced from `B019.yaml`.

SM validation: proxy SM `R_K* = 1.000`; K* central-q2 SM bin normalization `3.676601736716339e-08`; SM pass ratio `0.25`.

NEEDS-HUMAN-PHYSICS: RS C9/C10 and K* form-factor normalization are documented proxies; full EW/lepton/QED/C7/nonlocal-charm/angular covariance matching is missing.

Files changed: `B019.py`, `test_B019.py`, `rare_b_kstar_dilepton.py`, `physics_adapters/rare_b_kstar_dilepton.py`.

Tests: B019 targeted `11 passed`; full `python -m pytest tests/constraints/ -q` `463 passed`.

Note: unrelated dirty/untracked files are present in the worktree and were left untouched.