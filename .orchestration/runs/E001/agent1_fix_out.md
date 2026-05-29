Updated E001 to require `experimental.block_key == "canonical_limit"` and raise `AnchorError` otherwise.

Added superseded ACME monkeypatch coverage and switched expected EDM values to `quarkConstraints.edm` core helpers. NP `NEEDS-HUMAN-PHYSICS` proxy unchanged.

Pytest: E001 focused `12 passed`; full `tests/constraints/ -q` `510 passed in 14.85s`.