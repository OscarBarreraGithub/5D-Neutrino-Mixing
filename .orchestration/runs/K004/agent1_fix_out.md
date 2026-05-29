1. BLOCKER: `K004.py` now routes experimental, SM, and validation SM block selection through scaffold `load_anchor`.
2. SHOULD-FIX: K004 tests now derive budget and SM-anchor expectations from YAML/constraint anchor; retained SM-limit numeric regression.
3. Optional nits: added `(1+Delta_EM)=0.997`; corrected SM-anchor citation wording to Buras-Venturini 2022, with BGS 2021 as validation.
4. SM-limit BR: `8.472598449611133e-11`.
5. `python -m pytest tests/constraints/ -q`: `65 passed in 2.54s`.
6. `git diff --stat`: empty because these K004 files are untracked in this checkout; B001 untouched.