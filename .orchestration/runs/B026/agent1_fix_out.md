Changed B026 to raise `AnchorError` when `load_anchor.block_key` mismatches the expected candidate; added monkeypatch regression for both anchor loads.

Pytest B026: `13 passed`.
Pytest `tests/constraints/ -q`: `680 passed`.