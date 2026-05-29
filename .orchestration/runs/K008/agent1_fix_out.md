Implemented K008 fix with append-only additions to `quarkConstraints/rare_kaon_dilepton.py`, plus adapter/K008/test routing.

y7 structure now: direct CP uses `2.67 * (V^2 + A^2)` with `y7V=0.73`, `y7A=-0.68`; interference uses vector-only `8.91 * |a_S| * |V|`. Source checked: https://arxiv.org/abs/hep-ph/0404127

RS y7V/y7A remains a proxy and is flagged `NEEDS-HUMAN-PHYSICS`; no single axial `Y_NP` rescale is used for K008.

Diagnostics now state semileptonic QCD running is not applied, LO factor `1.0`, effect `0.0`; CPC diagnostic uses YAML `4.7e-13`.

Tests: `test_K008.py` 11 passed; `tests/constraints/ -q` 298 passed.