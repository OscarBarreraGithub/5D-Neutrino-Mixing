1. No BLOCKER / SHOULD-FIX / NIT findings for CR004.

Evidence:
- Isolation: CR004 imports only the collider adapter; adapter diff is append-only, with no existing function bodies changed.
- Anchor: YAML values loaded via `load_anchor`; missing value id and mismatched block both raise `AnchorError`.
- Numerics: core cross-check agrees exactly.
  - Safe: `m=1.8 TeV`, limit `1.57 TeV`, ratio `0.872222222222`, pass.
  - Excluded: `m=1.2 TeV`, limit `1.57 TeV`, ratio `1.30833333333`, fail.
- Contract/determinism: numeric result fields are real `float`; repeated evaluates equal; point extras unchanged; missing couplings return non-crashing `passes=True`, `predicted=None`, `ratio=None`.
- Registry smoke: `CR004 collider_rs HARD m(B pair -> tW/bZ/bH)`.
- Pytest: `tests/constraints/primary/collider_rs/test_CR004.py -q`: `8 passed`; `python -m pytest tests/constraints/ -q`: `678 passed in 14.84s`.

CODE-OK