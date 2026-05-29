Plan done: studied K001/scaffold/K004/RS couplings; added rare-kaon core + adapter; added K004 + tests; ran suite.

SM: Buras/BGS formula with `kappa_+=5.173e-11*(lambda/0.225)^8`, `P_c=0.404`, `X_t=1.481`, repo CKM target. Sources: K004 YAML, [arXiv:1503.02693](https://arxiv.org/abs/1503.02693), [arXiv:2105.02868](https://arxiv.org/abs/2105.02868).

RS proxy: `X_NP=(Delta_sd^L+Delta_sd^R) Delta_nu/(g_SM^2 M_KK^2)`, using mass-basis `left_down/right_down` overlaps and a Z-like `g/(2 c_W)` coupling.

SM-limit validation: `BR_SM = 8.498092727794516e-11`.

Budget: NA62 `9.6^{+1.9}_{-1.8}e-11` plus BV SM `8.60±0.42e-11`; combined `+1.9458674e-11 / -1.8483506e-11`.

Created: `quarkConstraints/rare_kaon_snd.py`, `physics_adapters/rare_kaon.py`, `primary/kaon/K004.py`, `test_K004.py`.

Cross-checks: SM point ratio `0.596`; safe NP `left=1e-5` ratio `0.606` passes; large NP `left=5e-2` ratio `98.52` fails.

Tests: `python -m pytest tests/constraints/ -q` -> `65 passed`; K004 file -> `9 passed`.

`git diff --stat` tracked-only currently shows pre-existing `flavor_catalog/processes/beauty/B001.yaml | 17 +`; new K004 files are untracked.

NEEDS-HUMAN-PHYSICS: full RS electroweak KK/Z/Z-prime tower and neutrino-coupling matching are not available on `ParameterPoint`; v1 uses the documented Z-like proxy.