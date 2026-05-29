1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation OK: `L023.py` imports only the `neutrino_trident` adapter, no other constraint; no pre-existing relevant adapter/base/registry diffs found.

Numerical cross-checks: CCFR anchor `0.82 +/- 0.28`, 95% CL `z=1.9599639845`, budget `0.5487899157`, allowed band `[0.2712100843, 1.3687899157]`.
SM: predicted `1.0`, manual `1.0`, ratio `0.3279943652`, pass.
Shift `(dv,da)=(0.2,-0.1)`: predicted `1.1317562937`, manual `1.1317562937`, ratio `0.5680794868`, pass.
Safe `dv=0.1`: predicted `1.0839410435`, manual `1.0839410435`, ratio `0.4809509723`, pass.
Excluded `dv=1.0`: predicted `2.0306646356`, manual `2.0306646356`, ratio `2.2060621033`, fail.

Contract probes OK: missing couplings returns non-crashing unevaluated result; numeric fields are real/non-complex; repeated evaluate is deterministic and did not mutate input. Anchor probes fail loudly for mismatched block and missing CCFR.

Tests: `python -m pytest tests/constraints/primary/charged_lepton/test_L023.py -q` -> `12 passed`; `python -m pytest tests/constraints/ -q` -> `910 passed`. Registry smoke imports `L023 charged_lepton PRIMARY HARD`.

CODE-OK