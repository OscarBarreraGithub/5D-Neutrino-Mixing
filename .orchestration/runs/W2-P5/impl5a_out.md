Plan executed: add charged W diagonalization, build charged-current extra, add 5a gates, run full suite.

Added [rs_charged_current.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_charged_current.py) with frozen `RSChargedCurrentCouplings`; charged diagonalization lives in [rs_ew_spectrum.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_ew_spectrum.py) as `RSEWSpectrum.charged_w_diagonalization`.

Extended [rs_ew_builder.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/rs_ew_builder.py) with opt-in `include_charged_current=True`; added `rs_charged_current` to `KNOWN_EXTRA_KEYS`. No constraint rewires.

Asserted numbers: root `2.4505096638137482`, `m_w=77.50672575578955`, `m_wprime=3071.5336376392315`, `eta_W=-1.0`.

Sample `delta_g_W_ud_L[0,1]=-8.169915813726022e-07+3.4880867617547286e-06j`; `epsilon[0,1,0]=2.6427715369627687e-05-1.2332357024682375e-05j`.

`delta_G_F/G_F=-1.8334382815196566e-05`; universal max `delta_g_W_ud_L=0.0`, universal max `epsilon=3.3881317890172014e-21`.

Tests: `tests/test_rs_charged_current_phase5a.py` -> `7 passed`; RS-EW slice -> `28 passed`; full `python -m pytest tests/ -q` -> `1671 passed, 1 skipped`.

Review: Codex re-review APPROVE after fixes. Opus is not available as a real tool/model in this session, so no Opus approval is claimed.

P5A-DONE.