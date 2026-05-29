Added shifted-Zcc independent core cross-check in `test_T012.py`: manual `dgL/dgR`, direct `quarkConstraints.zpole` recompute, compare `R_c/A_c/ratio`.

Shifted point: `R_c=0.18385094638346`, `A_c=0.707756438361834`, `ratio=3.9169821278201`.

Pytest: `test_T012.py -q` -> 11 passed; `tests/constraints/ -q` -> 297 passed.