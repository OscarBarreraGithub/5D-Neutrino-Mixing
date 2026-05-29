1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation OK: `T019.py` imports only scaffold/base/registry plus `physics_adapters.higgs_lfv`; no T018 import, and no diff in existing T018/adapter/core files.

Numerical cross-check: `Y_eτ=2e-4+1e-4j`, `Y_τe=5e-5j` gives constraint BR `6.428398173640321e-05`, core BR `6.428398173640321e-05`, manual BR `6.428398173640323e-05`; ratio `0.0321419908682016`, pass `True`.

Safe/fail probe: `Y_eτ=2e-4` gives BR `4.897827179916435e-05`, ratio `0.024489135899582175`, pass `True`; `Y_eτ=1e-2` gives BR `0.12244567949791091`, ratio `61.22283974895545`, pass `False`.

Contract/determinism OK: numeric result fields are real floats, complex Yukawas stay in diagnostics, missing input returns non-crashing unevaluated result, repeated evaluate returns identical result.

Anchor probe OK: YAML limit `0.002` via `pdg_or_equivalent.values[2]`; mismatched `load_anchor` block raises `AnchorError`.

Registry import smoke OK: `fcc.get("T019")` resolves `top_higgs_ew`, HARD, `BR(h -> e tau)`.

Pytest: `python -m pytest tests/constraints/ -q` -> `626 passed in 14.42s`.

CODE-OK