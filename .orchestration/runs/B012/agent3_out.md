1. BLOCKER [B012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B012.py:146): `_load_b012_anchor` accepts a valid-looking `Anchor` with `block_key="wrong_block"`; my probe returned `NO_RAISE`. Missing-anchor probe did raise `AnchorError`. Fix: validate neutral/charged `block_key` against the expected candidate tuples and add the monkeypatch mismatch test near [test_B012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B012.py:147).

Isolation OK: B012 imports physics only via `flavor_catalog_constraints.physics_adapters.bsgamma`; adapter diff is append-only 118 insertions; `quarkConstraints/bsgamma.py` untouched.

Contract OK: missing-extra path gives `passes=True, predicted=None, ratio=None`; evaluated numeric fields are `float`, complex C7 values stay in diagnostics; deterministic probe equal and no coupling mutation.

Numerics: SM/normalization `4.163e-05`, budget `9.2e-07`, validates against `4.3e-5` within 5%.

Independent core cross-check: constraint `4.113330056975298e-05`, direct `quarkConstraints.bsgamma` `4.113330056975298e-05`, abs diff `0`, ratio_to_SM `0.988068714142517`.

Safe/excluded: right=0.05 predicts `4.207862477952566e-05`, ratio `0.487635629919197`, pass; right=0.08 predicts `4.27784794355857e-05`, ratio `1.248347212593149`, fail.

Pytest: `tests/constraints/primary/beauty/test_B012.py -q` = 11 passed; `python -m pytest tests/constraints/ -q` = 193 passed.

Registry smoke: 19 constraints, B012 registered, family `beauty`, import_failures `{}`.

CODE-NEEDS-FIXES