Option (a): added YAML SM-SD theory anchor `0.79 +/- 0.12e-9` from Gorbahn-Haisch 2006: https://arxiv.org/abs/hep-ph/0605203

Changed K006 to load that anchor via `load_anchor`, validate the repo SM formula against it, and keep the SD-bound verdict `HARD`.

Updated K006 tests to derive SM validation from the YAML anchor, not a hardcoded `8.2098e-10`.

Tests: K006 focused `11 passed`; `tests/constraints/ -q` `244 passed`.