1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation OK: K008 imports physics only via `physics_adapters.rare_kaon_dilepton`; AST check found no changed/removed pre-existing symbols in core/adapter, only appended π0ee symbols.

Anchor OK: `load_anchor` path used; probes raised `AnchorError` for missing observable and unit mismatch.

Numerical cross-check: budget `2.8e-10`; SM direct `4.836298329689e-12` within YAML `3e-12..6e-12`; SM constructive total `3.800575583687e-11` within YAML `8e-12..45e-12`.

Safe point `left=1e-5j`: constraint `4.898960071204e-12`, independent manual `4.898960071204e-12`, ratio `0.01749628597`, passes `True`.

Excluded point `left=5e-2j`: constraint `5.358778045830e-09`, independent manual `5.358778045830e-09`, ratio `19.1384930208`, passes `False`.

Contract/determinism OK: real `ConstraintResult` numeric fields, complex amplitudes in diagnostics, missing-couplings path non-crashing, repeated evaluate equal/no coupling mutation.

Tests: `tests/constraints/primary/kaon/test_K008.py` `10 passed`; full `python -m pytest tests/constraints/ -q` `296 passed`.

Registry smoke: `K008 kaon HARD BR(K_L -> pi0 e+ e-)_direct CP`.

CODE-OK